# -*- coding: utf-8 -*-
"""
compare_scene.py
================================================================================
Headless comparison scene for benchmarking two unilateral-Lagrangian-constraint
implementations on the same two-tube CTR scene:

    mode='cpulc'   ContactPointsUnilateralConstraint reading directly from BCM
                   (mu = 0.2 -- DIFFERENT from the original ctr_two_tubes.py
                    which uses mu = 0; we standardise to mu = 0.2 in BOTH
                    modes here for a clean comparison).

    mode='feeder'  Stock UnilateralLagrangianConstraint<Vec3d> fed each step
                   by a ContactFeeder bridging from BCM data.

Mode dispatch is via the environment variable CTR_COMPARE_MODE.  Output path
is via CTR_COMPARE_OUT.  Defaults: 'cpulc' and './comparison_<mode>.npz'.

Geometry, SSIM, BCM, MOs, solvers, dt, RequiredPlugins -- byte-identical between
modes.  Only the constraint block differs.

Usage:
    runSofa --start -g batch compare_scene.py
    (with CTR_COMPARE_MODE / CTR_COMPARE_OUT set in the environment)

The original ctr_two_tubes.py / gui.py / live_monitor.py / init_monitoring.py
files are NOT touched.  Their helpers are imported.
"""

import math
import os
import sys

import Sofa
import Sofa.Core

# ---- Make sure the original module is importable -----------------------------
# If the existing scene files live next to this one, the directory is already
# on sys.path when runSofa loads us; otherwise the user can prepend it via
# the CTR_COMPARE_SOURCE env var.
_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)

_extra = os.environ.get('CTR_COMPARE_SOURCE')
if _extra and _extra not in sys.path:
    sys.path.insert(0, _extra)

# ---- Imports from the existing scene module ---------------------------------
# These are import-safe (top-level only defines constants and functions; the
# actual scene is built inside ctr_two_tubes.createScene which we do NOT call).
from ctr_two_tubes import (
    T1_PARAMS, T2_PARAMS, MAX_K,
    ALARM_DISTANCE, ALGORITHM, DEFAULT_NORMAL,
    compute_concentric_offset, add_cosserat_tube, add_tube_visual,
    CTRController,
)

from compare_controllers import (
    Schedule, ScriptedActuator, ComparisonRecorder,
)


# =============================================================================
#  Top-of-file tunables.  Change these to alter the comparison schedule.
# =============================================================================

DT                       = 1e-3               # fixed throughout the run
MU                       = 0.2                   # both modes

INIT_STEPS               = 5_000
HOLD_STEPS               = 200

TRANS_INCREMENT_M        = 1e-3                  # 1 mm
TRANS_N_INCREMENTS       = 60                     # -> 30 mm total
TRANS_DELTA_PER_STEP     = 5e-6                  # 1 um/step  (1000 steps/incr)

ROT_INCREMENT_RAD        = math.radians(1.0)     # 1 deg
ROT_N_INCREMENTS         = 0                     # -> 30 deg total
ROT_DELTA_PER_STEP       = math.radians(0.05)    # 0.01 deg/step (100 steps/incr)


# =============================================================================
#  createScene
# =============================================================================

def createScene(root_node):
    # --------------------------------------------------------------------
    #  Mode dispatch
    # --------------------------------------------------------------------
    mode = os.environ.get('CTR_COMPARE_MODE', 'cpulc').strip().lower()
    if mode not in ('cpulc', 'feeder'):
        raise RuntimeError(f"CTR_COMPARE_MODE must be 'cpulc' or 'feeder', "
                           f"got {mode!r}")

    out_path = os.environ.get(
        'CTR_COMPARE_OUT', f'./comparison_{mode}.npz'
    )
    print(f"[compare_scene] mode = {mode!r}", flush=True)
    print(f"[compare_scene] out  = {out_path!r}", flush=True)

    # --------------------------------------------------------------------
    #  Required plugins (identical to ctr_two_tubes.py)
    # --------------------------------------------------------------------
    root_node.addObject('RequiredPlugin', pluginName=[
        'Cosserat',
        'Sofa.Component.AnimationLoop',
        'Sofa.Component.Constraint.Lagrangian.Correction',
        'Sofa.Component.Constraint.Lagrangian.Model',
        'Sofa.Component.Constraint.Lagrangian.Solver',
        'Sofa.Component.Constraint.Projective',
        'Sofa.Component.LinearSolver.Direct',
        'Sofa.Component.Mapping.Linear',
        'Sofa.Component.Mapping.NonLinear',
        'Sofa.Component.Mass',
        'Sofa.Component.ODESolver.Backward',
        'Sofa.Component.Setting',
        'Sofa.Component.StateContainer',
        'Sofa.Component.Topology.Container.Constant',
        'Sofa.Component.Visual',
        'Sofa.GL.Component.Rendering3D',
        'Sofa.GUI.Component',
    ])

    root_node.gravity = [0., 0., 0.]
    root_node.dt      = DT

    root_node.addObject('DefaultVisualManagerLoop')
    root_node.addObject('FreeMotionAnimationLoop')
    root_node.addObject('BackgroundSetting', color=[1.0, 1.0, 1.0, 0])

    constraint_solver = root_node.addObject(
        'BlockGaussSeidelConstraintSolver',
        name='ConstraintSolver',
        tolerance=1e-5,
        maxIterations=500,
    )

    root_node.addObject('Camera',
                        position=[0.5, -0.3, 0.3],
                        lookAt=[0.08, 0.0, 0.0])
    root_node.addObject('VisualStyle',
                        displayFlags='showVisualModels hideBehaviorModels '
                                     'hideCollisionModels '
                                     'hideBoundingCollisionModels '
                                     'hideForceFields '
                                     'hideInteractionForceFields '
                                     'hideWireframe '
                                     'hideMechanicalMappings')

    # --------------------------------------------------------------------
    #  Tubes (identical to ctr_two_tubes.py)
    # --------------------------------------------------------------------
    x_t3 = compute_concentric_offset(T1_PARAMS, T2_PARAMS)

    t1_base_mo, _, _, t1_frame_node, _ = add_cosserat_tube(
        root_node, T1_PARAMS,
        init_strategy='natural',
    )
    t2_base_mo, _, _, t2_frame_node, _ = add_cosserat_tube(
        root_node, T2_PARAMS,
        x_offset=x_t3,
        init_strategy='conform_to_outer',
        outer_params=T1_PARAMS,
    )

    add_tube_visual(t1_frame_node, T1_PARAMS, color=T1_PARAMS['color'])
    add_tube_visual(t2_frame_node, T2_PARAMS, color=T2_PARAMS['color'])

    # --------------------------------------------------------------------
    #  Contact pipeline -- SSIM + BCM identical between modes
    # --------------------------------------------------------------------
    intersection_node = root_node.addChild('IntersectionNode')

    t1_MO = t1_frame_node.FramesMO
    t2_MO = t2_frame_node.FramesMO

    ssim = intersection_node.addObject(
        'SphereSweptIntersectionMethod',
        name='ssim',
        beam1Frames=t1_MO.getLinkPath() + '.position',
        beam2Frames=t2_MO.getLinkPath() + '.position',
        beam1Velocities=t1_MO.getLinkPath() + '.velocity',
        beam2Velocities=t2_MO.getLinkPath() + '.velocity',
        radius1=T1_PARAMS['rex'],
        radius2=T2_PARAMS['rex'],
        innerRadius1=T1_PARAMS['rin'],
        innerRadius2=T2_PARAMS['rin'],
        algorithmType=ALGORITHM,
        contactConfiguration='nested',
        defaultNormal=DEFAULT_NORMAL,
        broadPhaseMarginFactor=1.5,
    )

    contact_output = t1_frame_node.addChild('contactOutput')
    t2_frame_node.addChild(contact_output)

    contactMO = contact_output.addObject(
        'MechanicalObject', template='Vec3d',
        name='contactMO_gap',
        position=[[0., 0., 0.]] * 2 * MAX_K,
    )

    bcm = contact_output.addObject(
        'BeamContactMapping',
        name='bcm',
        input1=t1_MO.getLinkPath(),
        input2=t2_MO.getLinkPath(),
        output=contactMO.getLinkPath(),
        ssim=ssim.getLinkPath(),
        mappingMode='contactPoints',
    )

    # --------------------------------------------------------------------
    #  Constraint block -- THIS is where the modes differ
    # --------------------------------------------------------------------
    if mode == 'cpulc':
        contact_output.addObject(
            'ContactPointsUnilateralConstraint',
            name='cpuc',
            mu=MU,                                          # = 0.2 (was 0 in original)
            contactTriads=bcm.getLinkPath() + '.contactTriads',
            gapSign=bcm.getLinkPath() + '.gapSign',
        )
        print(f"[compare_scene] using ContactPointsUnilateralConstraint, "
              f"mu={MU}", flush=True)
    else:  # mode == 'feeder'
        contact_output.addObject(
            'UnilateralLagrangianConstraint',
            template='Vec3d',
            name='ulc',
            object1=contactMO.getLinkPath(),
            object2=contactMO.getLinkPath(),
        )
        contact_output.addObject(
            'ContactFeeder',
            name='feeder',
            distances     = bcm.getLinkPath() + '.distances',
            constraint    = '@ulc',
            alarmDistance = ALARM_DISTANCE,
            mu            = MU,
            contactTriads = bcm.getLinkPath() + '.contactTriads',
            gapSign       = bcm.getLinkPath() + '.gapSign',
        )
        print(f"[compare_scene] using UnilateralLagrangianConstraint + "
              f"ContactFeeder, mu={MU}", flush=True)

    # --------------------------------------------------------------------
    #  Schedule, ScriptedActuator, ComparisonRecorder
    # --------------------------------------------------------------------
    trans_steps_per_incr = int(round(TRANS_INCREMENT_M / TRANS_DELTA_PER_STEP))
    rot_steps_per_incr   = int(round(ROT_INCREMENT_RAD / ROT_DELTA_PER_STEP))

    schedule = Schedule(
        init_steps           = INIT_STEPS,
        hold_steps           = HOLD_STEPS,
        trans_n_increments   = TRANS_N_INCREMENTS,
        trans_steps_per_incr = trans_steps_per_incr,
        rot_n_increments     = ROT_N_INCREMENTS,
        rot_steps_per_incr   = rot_steps_per_incr,
    )
    print(f"[compare_scene] schedule: total_steps={schedule.total_steps}, "
          f"total_snapshots={schedule.total_snapshots}", flush=True)
    print(f"[compare_scene]   init={schedule.init_steps}, "
          f"trans cycle={schedule.trans_cycle_len} "
          f"({trans_steps_per_incr} act + {HOLD_STEPS} hold) x {TRANS_N_INCREMENTS}, "
          f"rot cycle={schedule.rot_cycle_len} "
          f"({rot_steps_per_incr} act + {HOLD_STEPS} hold) x {ROT_N_INCREMENTS}",
          flush=True)

    actuator = root_node.addObject(ScriptedActuator(
        name='ScriptedActuator',
        root_node=root_node,
        t1_base_mo=t1_base_mo,
        t2_base_mo=t2_base_mo,
        t2_x_offset=x_t3,
        schedule=schedule,
        trans_delta_per_step=TRANS_DELTA_PER_STEP,
        rot_delta_per_step=ROT_DELTA_PER_STEP,
        set_pose_fn=CTRController._set_pose,
    ))

    # The recorder needs Tube_3 frame curvilinear abscissae.  Read them
    # from the DiscreteCosseratMapping object (same trick as
    # InitializationMonitor in ctr_two_tubes.createScene).
    t2_curv_abs_frames = list(
        t2_frame_node.cosseratMapping.curv_abs_output.value
    )

    intersection_node.addObject(ComparisonRecorder(
        name='ComparisonRecorder',
        root_node=root_node,
        actuator=actuator,
        schedule=schedule,
        constraint_solver=constraint_solver,
        t2_MO=t2_MO,
        bcm=bcm,
        contact_mo=contactMO,
        t2_frame_curv_abs=t2_curv_abs_frames,
        out_path=out_path,
        mode_label=mode,
    ))

    # Auto-start in batch mode -- there's no GUI to click Initialize.
    root_node.animate = True

    return root_node
