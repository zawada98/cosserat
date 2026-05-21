# -*- coding: utf-8 -*-
"""
SphereSweptIntersectionMethod – demonstration scene
====================================================

Beam 1 – COMPLETELY FIXED  (acts as a rigid obstacle along +X)
Beam 2 – PARALLEL to Beam 1, spring-driven by distributed frame targets.

Beam 2 starts GAP_Z mm directly above Beam 1.  During initialization the
distributed control frames descend to the nested-contact clearance so the
whole beam settles onto the inner wall.  In control phase those targets impose
the no-slip rolling orbit and axial spin homogeneously along the length.

Live data collection + 3D plots
---------------------------------
An SSIMPlotter controller records at every step:
  • simulation time  t
  • curvilinear parameter on Beam 1  s1*
  • curvilinear parameter on Beam 2  s2*
  • signed gap  δ  (< 0 ⇒ interpenetration)

At the end of the simulation it generates two 3D plots:
  (1)  δ  vs  (s1 , t)   saved as  delta_vs_s1_time.png
  (2)  δ  vs  (s2 , t)   saved as  delta_vs_s2_time.png

HOW TO STOP THE SIMULATION AND SAVE THE PLOTS
----------------------------------------------
  1. Press  [Q]  inside the runSofa viewer window.
  2. Press  [P]  to save plots without stopping.
  3. Let the simulation run until MAX_STEPS is reached.

Run with:
    runSofa concentric_straight_beams_rolling_spring.py
"""


import matplotlib
matplotlib.use("Agg")          # non-interactive backend – safe inside runSofa

import math
import Sofa
import Sofa.Core


import os

from live_monitor   import LiveContactMonitor
from gui2 import CTRGuiBridgeStraightOuter
try:
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment
    _HAS_OPENPYXL = True
except ImportError:
    _HAS_OPENPYXL = False

SCENE_DIR = os.path.dirname(os.path.abspath(__file__))


# ──────────────────────────────────────────────────────────────────────────────
#  Scene parameters
# ──────────────────────────────────────────────────────────────────────────────
BEAM_LENGTH      = 120.0    # [mm]  total length of each beam
NB_SECTIONS      = 5        # number of Cosserat sections per beam
NB_FRAMES        = 10       # number of output Rigid3d frames per beam
RADIUS1_EX       = 6.0      # [mm]  cross-section radius
RADIUS1_IN       = 5.0      # [mm]  cross-section radius
RADIUS2_EX       = 2.0      # [mm]  cross-section radius
RADIUS2_IN       = 1.0      # [mm]  cross-section radius
YOUNG_MODULUS    = 3.0e6    # [Pa]
POISSON_RATIO    = 0.49
STIFFNESS        = 1e8  # control-point spring stiffness for Beam 2 base
DT               = 1e-3   # [s]   time step
MAX_STEPS        = 500      # stop automatically after this many steps
MAX_K = NB_FRAMES
ALARM_DISTANCE = 3.0 * RADIUS2_EX

# Keep damping small on this scene.  A large stiffness-proportional Rayleigh
# term damps the high-stiffness control spring itself, which makes the proximal
# base lag far behind the moving control point.
RAYLEIGH_STIFFNESS = 1e-3
RAYLEIGH_MASS      = 0

class JiggleRecorder(Sofa.Core.Controller):
    def __init__(self, t2_frames_mo, out_path, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mo = t2_frames_mo
        self.out_path = out_path
        self.rows = []

    def onAnimateEndEvent(self, event):
        t = self.getContext().time.value
        tip = list(self.mo.position.value[-1])[:3]  # last frame center
        mid = list(self.mo.position.value[len(self.mo.position.value)//2])[:3]
        self.rows.append((t, *tip, *mid))

    def onSimulationInitDoneEvent(self, event):
        pass

    def __del__(self):
        try:
            with open(self.out_path, 'w') as f:
                f.write("t,tip_x,tip_y,tip_z,mid_x,mid_y,mid_z\n")
                for r in self.rows:
                    f.write(",".join(str(v) for v in r) + "\n")
        except Exception:
            pass

# ──────────────────────────────────────────────────────────────────────────────
#  Geometry
# ──────────────────────────────────────────────────────────────────────────────
#
#  Beam 1 – along global +X, centred at Y=0, Z=0  (fully fixed obstacle)
#    base at [0, 0, 0], quaternion [0,0,0,1]
#
#  Beam 2 – cantilever, clamped at base, free end falls under gravity (−Z).
#    base at [0, GAP_Y, GAP_Z].
#
#    GAP_Z > 0  : Beam 2 starts above Beam 1.
#    GAP_Y ≠ 0  : lateral eccentricity. Contact normal will lie in the (Y,Z)
#                 plane, tilted by atan2(GAP_Y, effective vertical gap at
#                 contact). GAP_Y = 0 reproduces the original axisymmetric
#                 (purely vertical) configuration.
#
# Initial vertical gap between the two parallel beam axes
GAP_Z = 0.0    # [mm]   must be > 2*RADIUS to start without interpenetration
GAP_Y = 0.0

# In the nested configuration, Beam 2 rests against the inner wall of Beam 1
# when the centerline offset reaches RADIUS1_IN - RADIUS2_EX.  During the
# initialization phase the distributed frame-control targets are lowered to
# this Z so the whole inner beam settles onto the wall.
INNER_BEAM_SETTLE_Z = GAP_Z - (RADIUS1_IN - RADIUS2_EX)
INNER_BEAM_SETTLE_STEP = 0.05  # [mm/step]

# Rolling geometry in the Y-Z cross-section.  The center of Beam 2 moves on a
# circle of radius ROLLING_ORBIT_RADIUS inside the inner wall of Beam 1.  The
# no-slip rolling spin is dtheta = -(orbit_radius / inner_radius) * dphi.
ROLLING_ORBIT_RADIUS = RADIUS1_IN - RADIUS2_EX
ROLLING_SPIN_RATIO = ROLLING_ORBIT_RADIUS / RADIUS2_EX

class ContactDebugRecorder(Sofa.Core.Controller):
    def __init__(self, t1_frames_mo, t2_frames_mo, t2_control_mo, t2_base_mo,
                 bcm, out_path, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.t1_frames_mo = t1_frames_mo
        self.t2_frames_mo = t2_frames_mo
        self.t2_control_mo = t2_control_mo
        self.t2_base_mo = t2_base_mo
        self.bcm = bcm
        self.out_path = out_path
        self._file = None

    def onSimulationInitDoneEvent(self, event):
        self._file = open(self.out_path, 'w')
        self._file.write(
            "t,"
            "control_x,control_y,control_z,"
            "base_x,base_y,base_z,"
            "control_minus_base_x,control_minus_base_y,control_minus_base_z,"
            "control_base_distance,"
            "min_centerline_distance,same_index_min_centerline_distance,"
            "estimated_nested_clearance,estimated_same_index_nested_clearance,"
            "min_bcm_gap,max_bcm_gap,num_bcm_gaps\n"
        )
        self._file.flush()

    def onAnimateEndEvent(self, event):
        if self._file is None:
            return

        t = float(self.getContext().time.value)
        control = list(self.t2_control_mo.position.value[0])[:3]
        base = list(self.t2_base_mo.position.value[0])[:3]
        diff = [control[i] - base[i] for i in range(3)]
        diff_norm = math.sqrt(sum(v * v for v in diff))

        f1 = [list(p)[:3] for p in self.t1_frames_mo.position.value]
        f2 = [list(p)[:3] for p in self.t2_frames_mo.position.value]

        min_center_dist = float('inf')
        for p1 in f1:
            for p2 in f2:
                d = math.sqrt(sum((p2[i] - p1[i]) ** 2 for i in range(3)))
                min_center_dist = min(min_center_dist, d)

        same_min_center_dist = float('inf')
        for p1, p2 in zip(f1, f2):
            d = math.sqrt(sum((p2[i] - p1[i]) ** 2 for i in range(3)))
            same_min_center_dist = min(same_min_center_dist, d)

        estimated_clearance = RADIUS1_IN - RADIUS2_EX - min_center_dist
        estimated_same_clearance = RADIUS1_IN - RADIUS2_EX - same_min_center_dist

        bcm_gaps = []
        try:
            for d in self.bcm.distances.value:
                gap = float(d[0])
                if abs(gap) < 1e8:
                    bcm_gaps.append(gap)
        except Exception:
            pass

        min_bcm_gap = min(bcm_gaps) if bcm_gaps else float('nan')
        max_bcm_gap = max(bcm_gaps) if bcm_gaps else float('nan')

        row = [
            t,
            *control,
            *base,
            *diff,
            diff_norm,
            min_center_dist,
            same_min_center_dist,
            estimated_clearance,
            estimated_same_clearance,
            min_bcm_gap,
            max_bcm_gap,
            len(bcm_gaps),
        ]
        self._file.write(",".join(str(v) for v in row) + "\n")
        self._file.flush()

    def __del__(self):
        try:
            if self._file is not None:
                self._file.close()
        except Exception:
            pass


class CTRController(Sofa.Core.Controller):
    """
    GUI controller for the rolling spring-driven scene.

    This keeps the spring-driven architecture from
    concentric_straight_beams_driven_spring.py, but the spring target is
    distributed over all mapped inner-beam frames instead of applied only at
    the proximal base.  The GUI's Z-rotation target is interpreted as the
    rolling/orbit angle phi, then every target frame is constrained to:

        y = d cos(phi0 + phi)
        z = d sin(phi0 + phi)
        rx = -(d / r_inner) phi

    where d = RADIUS1_IN - RADIUS2_EX.
    """

    DT_RAMP_PER_STEP = 1.02

    def __init__(self,
                 root_node,
                 t1_base_mo,
                 t2_base_mo,
                 t2_control_mo,
                 gui_bridge,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.root_node     = root_node
        self.t1_base_mo    = t1_base_mo
        self.t2_base_mo    = t2_base_mo
        self.t2_control_mo = t2_control_mo
        self.gui           = gui_bridge

        # Outer: pin to construction-time pose forever.
        self._t1_rest_pose = read_rigid_pose(t1_base_mo)

        # Inner: filled on init->control transition (NOT at construction).
        self._t2_rest_poses = None
        self._rolling_phi0 = None

        self._t2_tx = 0.0
        self._rolling_phi = 0.0

        try:
            dt0 = float(self.root_node.dt.value)
        except Exception:
            dt0 = 1e-4
        self._dt_current = dt0
        self._dt_target  = dt0

        self._step       = 0
        self._prev_phase = 'waiting'

        self._init_start_wall = None
        self._INIT_TIMEOUT_S = 5.0

        # ------------------------------------------------------------------
    def onAnimateBeginEvent(self, event):
        self._step += 1
        snap  = self.gui.snapshot()
        phase = snap['phase']

        # 1) Outer tube: pin every step regardless of phase.
        write_rigid_pose(self.t1_base_mo, self._t1_rest_pose)

        if phase == 'initializing':
            import time as _time
            if self._init_start_wall is None:  # first entry
                self._init_start_wall = _time.monotonic()
            elapsed = _time.monotonic() - self._init_start_wall
            if elapsed >= self._INIT_TIMEOUT_S:
                print(f"[CTRController] init timeout ({self._INIT_TIMEOUT_S:.1f} s "
                      f"wall clock) — forcing control phase at step {self._step}")
                self.gui.signal_init_complete()  # flips phase + enables GUI widgets
                # phase will read 'control' on the next snapshot(); nothing more
                # needs to happen this step.

        # 2) Detect init -> control transition and snapshot the inner base
        if self._prev_phase != 'control' and phase == 'control':
            self._t2_rest_poses = read_rigid_poses(self.t2_control_mo)
            x, y, z, qx, qy, qz, qw = self._t2_rest_poses[0]
            self._rolling_phi0 = math.atan2(z, y)
            print(f"[CTRController] init -> control at step {self._step}; "
                  f"captured inner-tube distributed control pose: "
                  f"pos=({x:+.4f}, {y:+.4f}, {z:+.4f}) "
                  f"quat=({qx:+.4f}, {qy:+.4f}, {qz:+.4f}, {qw:+.4f}); "
                  f"rolling phi0={math.degrees(self._rolling_phi0):+.2f} deg")
        self._prev_phase = phase

        # 3) Control phase: dt ramp + DOF advancement + write control point.
        if phase == 'control':
            dt_req = self.gui.consume_dt_request()
            if dt_req is not None:
                self._dt_target = float(dt_req)
                print(f"[CTRController] dt target -> {self._dt_target:.6g} s "
                      f"at step {self._step}; ramping from "
                      f"{self._dt_current:.6g} s "
                      f"at <= {(self.DT_RAMP_PER_STEP - 1) * 100:.1f}% per step")

            if abs(self._dt_current - self._dt_target) > 1e-15:
                ratio = self._dt_target / self._dt_current
                if ratio > self.DT_RAMP_PER_STEP:
                    self._dt_current *= self.DT_RAMP_PER_STEP
                elif ratio < 1.0 / self.DT_RAMP_PER_STEP:
                    self._dt_current /= self.DT_RAMP_PER_STEP
                else:
                    self._dt_current = self._dt_target
                try:
                    self.root_node.dt = self._dt_current
                except Exception as e:
                    print(f"[CTRController] failed to write dt: {e!r}")

            max_t = float(snap['translation_step_m'])
            max_r = float(snap['rotation_step_rad'])

            self._t2_tx = self._step_toward(self._t2_tx,
                                            float(snap['t2_tx_target_m']),  max_t)
            # Reinterpret the GUI's Z-rotation slider as the rolling orbit
            # angle phi.  The controller owns Y, Z, and axial spin.
            self._rolling_phi = self._step_toward(
                self._rolling_phi,
                float(snap['t2_rz_target_rad']),
                max_r,
            )

            poses = self._compose_inner_poses()
            self._write_control_point_poses(poses)
        else:
            self._settle_inner_control_z()

        # 'waiting' / 'initializing': do not drive the control point from GUI
        # targets.  The inner base remains spring-driven from its current
        # control-point rest pose until control phase begins.

    def _compose_inner_poses(self):
        orbit_angle = self._rolling_phi0 + self._rolling_phi
        new_y = ROLLING_ORBIT_RADIUS * math.cos(orbit_angle)
        new_z = ROLLING_ORBIT_RADIUS * math.sin(orbit_angle)

        rolling_spin = -ROLLING_SPIN_RATIO * self._rolling_phi
        q_twist = self._quat_x(rolling_spin)

        poses = []
        for rest_pose in self._t2_rest_poses:
            new_x = rest_pose[0] + self._t2_tx * 1000.0
            q_rest = tuple(rest_pose[3:7])
            q_out = self._quat_normalize(self._quat_mul(q_rest, q_twist))
            poses.append([
                new_x, new_y, new_z,
                q_out[0], q_out[1], q_out[2], q_out[3],
            ])
        return poses

    # ---- (unchanged: _step_toward, _quat_x/z/mul/normalize) --
    def _settle_inner_control_z(self):
        poses = read_rigid_poses(self.t2_control_mo)
        for pose in poses:
            pose[2] = self._step_toward(float(pose[2]),
                                        INNER_BEAM_SETTLE_Z,
                                        INNER_BEAM_SETTLE_STEP)
        self._write_control_point_poses(poses)

    def _write_control_point_poses(self, poses):
        """SOFA-thread write to spring targets, not to simulated frames."""
        write_rigid_poses(self.t2_control_mo, poses)

    @staticmethod
    def _step_toward(current, target, max_step):
        delta = target - current
        if delta >  max_step: return current + max_step
        if delta < -max_step: return current - max_step
        return target

    @staticmethod
    def _quat_x(theta):
        h = 0.5 * theta
        return (math.sin(h), 0.0, 0.0, math.cos(h))

    @staticmethod
    def _quat_z(theta):
        h = 0.5 * theta
        return (0.0, 0.0, math.sin(h), math.cos(h))

    @staticmethod
    def _quat_mul(a, b):
        ax, ay, az, aw = a
        bx, by, bz, bw = b
        return (
            aw*bx + ax*bw + ay*bz - az*by,
            aw*by - ax*bz + ay*bw + az*bx,
            aw*bz + ax*by - ay*bx + az*bw,
            aw*bw - ax*bx - ay*by - az*bz,
        )

    @staticmethod
    def _quat_normalize(q):
        n2 = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]
        if n2 < 1e-24:
            return (0.0, 0.0, 0.0, 1.0)
        inv = 1.0 / math.sqrt(n2)
        return (q[0]*inv, q[1]*inv, q[2]*inv, q[3]*inv)


def extract_contact_points(parent_node, contactMO, MAX_K):
    even_indices = list(range(0, 2 * MAX_K, 2))
    odd_indices  = list(range(1, 2 * MAX_K, 2))

    node_A = parent_node.addChild('Pc_A')
    contactMO_A = node_A.addObject(
        'MechanicalObject', template='Vec3d',
        name='contactMO_A',
        position=[[0., 0., 0.]] * MAX_K)
    node_A.addObject(
        'SubsetMapping', template='Vec3d,Vec3d',
        input=contactMO.getLinkPath(),
        output='@contactMO_A',
        indices=even_indices,
        handleTopologyChange=False)

    node_B = parent_node.addChild('Pc_B')
    contactMO_B = node_B.addObject(
        'MechanicalObject', template='Vec3d',
        name='contactMO_B',
        position=[[0., 0., 0.]] * MAX_K)
    node_B.addObject(
        'SubsetMapping', template='Vec3d,Vec3d',
        input=contactMO.getLinkPath(),
        output='@contactMO_B',
        indices=odd_indices,
        handleTopologyChange=False)

    return contactMO_A, contactMO_B
# ──────────────────────────────────────────────────────────────────────────────
#  Cosserat beam builder
# ──────────────────────────────────────────────────────────────────────────────

def _make_section_params(nb_sections, length):
    sec_len = length / nb_sections
    return (
        [[0., 0., 0.]] * nb_sections,
        [sec_len] * nb_sections,
        [i * sec_len for i in range(nb_sections + 1)],
    )


def _make_frame_params(nb_frames, length):
    fl       = length / nb_frames
    frames   = [[i * fl, 0., 0.,  0., 0., 0., 1.] for i in range(nb_frames + 1)]
    curv_abs = [i * fl for i in range(nb_frames + 1)]
    edge_indices = [[i, i + 1] for i in range(nb_frames)]
    return frames, curv_abs, edge_indices

def read_rigid_pose(mo):
    """Return the first Rigid3d pose as a plain 7-value Python list."""
    pose = list(mo.position.value[0])
    if len(pose) != 7:
        raise ValueError(f"Expected one Rigid3d pose with 7 values, got {pose!r}")
    return pose

def read_rigid_poses(mo):
    """Return all Rigid3d poses as plain 7-value Python lists."""
    poses = [list(pose) for pose in mo.position.value]
    for pose in poses:
        if len(pose) != 7:
            raise ValueError(f"Expected Rigid3d poses with 7 values, got {pose!r}")
    return poses

def write_rigid_pose(mo, pose):
    """
    Robocop-style Rigid3d update: mutate SOFA data through writeable().

    The GUI thread never calls this directly; CTRController calls it from
    SOFA's animation callback, and the RestShapeSpringsForceField reads the
    resulting controlPointMO pose as its external rest shape.
    """
    rigid_pose = list(pose)
    if len(rigid_pose) != 7:
        raise ValueError(f"Expected a Rigid3d pose with 7 values, got {pose!r}")
    with mo.position.writeable() as position:
        position[0] = rigid_pose

def write_rigid_poses(mo, poses):
    """Write all Rigid3d poses through SOFA's writeable data API."""
    rigid_poses = [list(pose) for pose in poses]
    for pose in rigid_poses:
        if len(pose) != 7:
            raise ValueError(f"Expected Rigid3d poses with 7 values, got {pose!r}")
    with mo.position.writeable() as position:
        for i, pose in enumerate(rigid_poses):
            position[i] = pose

def add_base_control_point(parent_node, name, base_pos, base_quat):
    control_node = parent_node.addChild(name + '_base_control')
    control_mo = control_node.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='controlPointMO',
        position=[list(base_pos) + list(base_quat)],
        showObject=True,
        showObjectScale=3.0)
    return control_mo

def add_frame_control_points(parent_node, name, frames, base_pos, base_quat):
    control_node = parent_node.addChild(name + '_frame_controls')
    bx, by, bz = base_pos
    qx, qy, qz, qw = base_quat
    control_positions = [
        [bx + frame[0], by + frame[1], bz + frame[2], qx, qy, qz, qw]
        for frame in frames
    ]
    control_mo = control_node.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='controlFramesMO',
        position=control_positions,
        showObject=True,
        showObjectScale=2.0)
    return control_mo

def add_cosserat_beam(parent_node, name, base_pos, base_quat,
                      nb_sections, nb_frames, length, radius, radius_in,
                      young_modulus, poisson_ratio, stiffness,
                      beam_number, fully_fixed=False, base_control_mo=None,
                      base_control_active_directions=None):
    """
    Build a Cosserat beam and attach it to *parent_node*.

    Parameters
    ----------
    fully_fixed : bool
        • False  (default) – beam is clamped at its base via a
          RestShapeSpringsForceField and carries a UniformMass so gravity
          acts on it (the tip is free to deform / fall).
        • True  – ALL degrees of freedom are frozen with
          FixedProjectiveConstraint on both the rigid-base MO and the
          cosserat-coordinate MO.  No mass is added.  The beam acts as a
          perfectly rigid, immovable obstacle.

    Returns
    -------
    frame_node : Sofa.Core.Node
        The child node that holds the mapped Rigid3d frames (FramesMO).
    """
    bx, by, bz     = base_pos
    qx, qy, qz, qw = base_quat

    sec_pos, sec_len, curv_in = _make_section_params(nb_sections, length)
    frames, curv_out, edge_indices = _make_frame_params(nb_frames, length)

    beam_node       = parent_node.addChild(name)
    solver_node = beam_node.addChild('solverNode')
    solver_node.addObject('EulerImplicitSolver',
                          rayleighStiffness=RAYLEIGH_STIFFNESS,
                          rayleighMass=RAYLEIGH_MASS)
    solver_node.addObject('SparseLDLSolver', name='solver',
                          template='CompressedRowSparseMatrixd')
    #solver_node.addObject('BlockGaussSeidelConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    solver_node.addObject('GenericConstraintCorrection', linearSolver = solver_node.solver.getLinkPath())

    rigid_base_node = solver_node.addChild('rigidBase')

    base_mo = rigid_base_node.addObject(
        'MechanicalObject', template='Rigid3d', name='RigidBaseMO',
        position=[bx, by, bz, qx, qy, qz, qw],
        showObject=True, showObjectScale=3.0)

    if fully_fixed:
        # ── Fix the rigid base completely ──────────────────────────────────
        rigid_base_node.addObject(
            'FixedProjectiveConstraint', name='fixBase',
            indices=[0])
    else:
        rigid_base_node.addObject('UniformMass', name='baseMass', totalMass=0.001)

    if (not fully_fixed) and base_control_mo is not None:
        if base_control_active_directions is None:
            base_control_active_directions = [1, 1, 1, 1, 1, 1, 1]

        rigid_base_node.addObject(
            'RestShapeSpringsForceField',
            name='base_control_spring',
            stiffness=stiffness,
            angularStiffness=stiffness,
            external_rest_shape=base_control_mo.getLinkPath(),
            external_points=[0],
            mstate='@RigidBaseMO',
            points=[0],
            template='Rigid3d',
            activeDirections=base_control_active_directions)

    coord_node = solver_node.addChild('cosseratCoordinate')
    coord_mo   = coord_node.addObject(
        'MechanicalObject', template='Vec3d',
        name='cosserat_state', position=sec_pos)
    coord_node.addObject(
        'BeamHookeLawForceField',
        crossSectionShape='circular',
        length=sec_len, radius=radius,
        innerRadius = radius_in,
        youngModulus=young_modulus,
        poissonRatio=poisson_ratio)

    if fully_fixed:
        # ── Fix all cosserat strains (zero = straight reference config) ────
        all_section_ids = list(range(nb_sections))
        coord_node.addObject(
            'FixedProjectiveConstraint', name='fixStrains',
            indices=all_section_ids)

    frame_node = solver_node.addChild('mappedFrames')

    frames_mo = frame_node.addObject(
        'MechanicalObject', template='Rigid3d', name='FramesMO',
        position=frames, showObject=True, showObjectScale=2.0)

    if not fully_fixed:
        # Only the free beam needs mass; the fixed beam has none.
        frame_node.addObject('UniformMass', totalMass=0.1)

    frame_node.addObject(
        'DiscreteCosseratMapping', name='cosseratMapping',
        curv_abs_input=curv_in, curv_abs_output=curv_out,
        input1='@../cosseratCoordinate/cosserat_state',
        input2='@../rigidBase/RigidBaseMO',
        output='@FramesMO', debug=False)

    edges_positions = [[x, y, z] for x, y, z, *_ in frames]

    return frame_node, frames, base_mo

def add_visual_model(framesNode, frames, rex, color=(0.85, 0.15, 0.15, 1.0)):
    """
    Build a tube-shaped visual model around the mapped Rigid3d frames of a
    Cosserat beam.

    Parameters
    ----------
    framesNode : Sofa.Core.Node
        The node holding the FramesMO (output Rigid3d frames of the beam).
    frames : list
        The list of frame positions+quaternions (used here only for its length).
    rex : float
        Outer radius of the cross-section [mm] – defines the tube radius.
    color : sequence of 4 floats, optional
        RGBA color of the tube surface, each component in [0, 1].
        Default is red (0.85, 0.15, 0.15, 1.0).
        Accepts any iterable of 4 floats (tuple, list, numpy array …).
    """
    N = len(frames)
    n_sides = 60
    TWO_PI = 2.0 * math.pi

    def _ring_positions(r):
        return [
            [0.0,
                r * math.cos(TWO_PI * k / n_sides),
             r * math.sin(TWO_PI * k / n_sides),
             ]
            for k in range(n_sides)
        ]

    def _tube_quads(n_frames, n_s):
        quads = []
        for i in range(n_frames - 1):
            for k in range(n_s):
                k1 = (k + 1) % n_s
                v0, v1 = i * n_s + k, i * n_s + k1
                v2, v3 = v0 + n_s, v1 + n_s
                quads.append([v0, v1, v2, v3])
        return quads

    rigid_idx = [i for i in range(N) for _ in range(n_sides)]
    quads = _tube_quads(N, n_sides)
    safe_name = "Tube"


    for r, color_list, suffix in [
        (rex, list(color), 'outer')
    ]:
        ring_pos = _ring_positions(r) * N
        color_str = " ".join(str(v) for v in color_list)

        vis = framesNode.addChild(f'visual_{suffix}_{safe_name}')

        # Topology (rest positions, connectivity)
        vis.addObject('MeshTopology', name='topo',
                      position=ring_pos, quads=quads)

        vis.addObject('MechanicalObject', name='visMO',
                      template='Vec3d', position=ring_pos)
        vis.addObject('RigidMapping',
                      input='@../FramesMO',
                      output='@visMO',
                      rigidIndexPerPoint=rigid_idx,
                      globalToLocalCoords=False)

        # OglModel in child node, driven by IdentityMapping from visMO
        ogl = vis.addChild('ogl')
        ogl.addObject('OglModel', name='oglModel',
                      src='@../topo',
                      color=color_str)
        ogl.addObject('IdentityMapping',
                      input='@../visMO',
                      output='@oglModel')


# ── Output directory (same folder as the scene file) ─────────────────────────
SCENE_DIR = os.path.dirname(os.path.abspath(__file__))


def createScene(root_node: Sofa.Core.Node):
    """
    Scene layout
    ============

    Beam 1 (fixed obstacle)
    -----------------------
    • Lies along the global +X axis at Z = 0.
    • Base at [0, 0, 0], no rotation.
    • ALL DOFs frozen with FixedProjectiveConstraint – it will never move.
    • No mass.

    Beam 2 (spring-driven beam settling under gravity)
    -----------------------------------------
    • Lies along the global +X axis – PARALLEL to Beam 1.
    • Base clamped at [0, 0, GAP_Z] (directly above Beam 1's base), same
      orientation quaternion [0,0,0,1].
    • Pulled at its base by a RestShapeSpringsForceField whose external rest
      shape is a GUI-driven Rigid3d control point.
    • Carries a UniformMass; gravity (−Z) bends it while the base-control
      spring lowers the proximal end to the nested contact clearance.
    """

    root_node.gravity = [0., 0., -9810.]   # mm/s² gravity
    root_node.dt      = DT

    root_node.addObject('RequiredPlugin', pluginName=[
        'Sofa.Component.Constraint.Projective',
        'Sofa.Component.LinearSolver.Direct',
        'Sofa.Component.Mass',
        'Sofa.Component.ODESolver.Backward',
        'Sofa.Component.SolidMechanics.Spring',
        'Sofa.Component.Visual',
        'Sofa.GL.Component.Rendering3D',
        'Cosserat',
        'Sofa.Component.StateContainer',
        'Sofa.Component.Setting',
        'Sofa.Component.Constraint.Lagrangian.Correction',
        'Sofa.Component.Mapping.Linear',
        'Sofa.Component.Mapping.NonLinear',
        'Sofa.Component.Topology.Container.Constant',
        "Sofa.Component.Collision.Detection.Intersection",
        "Sofa.Component.Collision.Response.Contact",
        "Sofa.Component.Constraint.Lagrangian.Correction",
        "Sofa.Component.Constraint.Lagrangian.Solver",
        "Sofa.GUI.Component", "Sofa.Component.Collision.Geometry", "Sofa.Component.LinearSolver.Direct",
        "Sofa.Component.Mapping.Linear", "Sofa.Component.MechanicalLoad",
        "Sofa.Component.StateContainer", "Sofa.Component.ODESolver.Backward",
        "Sofa.Component.AnimationLoop",
        "Sofa.Component.Collision.Detection.Algorithm",
        "Sofa.Component.Topology.Container.Dynamic",
        'Sofa.Component.Constraint.Lagrangian.Model'
    ])

    root_node.addObject('DefaultVisualManagerLoop')
    root_node.addObject('FreeMotionAnimationLoop')
    root_node.addObject('BackgroundSetting', color=[0.05, 0.05, 0.12, 1.0])

    # Camera positioned to see both beams from a 3/4-angle view
    root_node.addObject('Camera',
                        position=[-40, -80, 120],
                        lookAt=[50, 0, 0])

    root_node.addObject('VisualStyle',
                        displayFlags='showVisualModels showBehaviorModels '
                                     'showCollisionModels '
                                     'hideBoundingCollisionModels '
                                     'hideForceFields '
                                     'hideInteractionForceFields '
                                     'hideWireframe '
                                     'showMechanicalMappings')

    root_node.addObject('BlockGaussSeidelConstraintSolver',tolerance=1e-5, maxIterations=500)

    # ── Beam 1 – horizontal along +X, COMPLETELY FIXED ───────────────────────
    #
    #   Axis:  X ∈ [0, L]   Y = 0   Z = 0
    #   This beam never deforms or translates; it is a rigid obstacle.
    #
    beam1_framesNode, beam1_frames, base_mo1 = add_cosserat_beam(
        root_node, 'Beam1_Fixed',
        base_pos   = [0., 0., 0.],
        base_quat  = [0., 0., 0., 1.],
        nb_sections = NB_SECTIONS,
        nb_frames   = NB_FRAMES,
        length      = BEAM_LENGTH,
        radius      = RADIUS1_EX,
        radius_in   = RADIUS1_IN,
        young_modulus  = YOUNG_MODULUS,
        poisson_ratio  = POISSON_RATIO,
        stiffness      = STIFFNESS,
        beam_number    = 1,
        fully_fixed    = True,
    )


    add_visual_model(beam1_framesNode, beam1_frames,
                     RADIUS1_IN, color=(0.55, 0.20, 0.75, 0.35))

    add_visual_model(beam1_framesNode, beam1_frames,
                     RADIUS1_EX, color=(0.55, 0.20, 0.75, 0.35))

    # ── Beam 2 – along +X (parallel to Beam 1), spring-driven base descends ─
    #
    #   Same orientation as Beam 1 (quaternion [0,0,0,1]).
    #   Base at [0, 0, GAP_Z] so both beams share the same X axis range and
    #   Beam 2 starts GAP_Z mm directly above Beam 1.
    #
    #   Gravity acts in -Z while the initialization controller lowers the
    #   distributed frame-control targets so the full beam settles onto Beam 1.
    #   In control phase, every mapped frame follows the same rolling orbit
    #   and no-slip axial spin through a RestShapeSpringsForceField.  The
    #   controller does not overwrite the simulated frame DOFs directly.
    #
    beam2_framesNode, beam2_frames, base_mo2 = add_cosserat_beam(
        root_node, 'Beam2_Cantilever',
        base_pos   = [0., GAP_Y, GAP_Z],
        base_quat  = [0., 0., 0., 1.],
        nb_sections = NB_SECTIONS,
        nb_frames   = NB_FRAMES,
        length      = BEAM_LENGTH,
        radius      = RADIUS2_EX,
        radius_in   = RADIUS2_IN,
        young_modulus  = YOUNG_MODULUS,
        poisson_ratio  = POISSON_RATIO,
        stiffness      = STIFFNESS,
        beam_number    = 2,
        fully_fixed    = False,
        base_control_mo=None,
    )
    beam2_control_mo = add_frame_control_points(
        root_node,
        'Beam2',
        frames=beam2_frames,
        base_pos=[0., GAP_Y, GAP_Z],
        base_quat=[0., 0., 0., 1.])
    beam2_framesNode.addObject(
        'RestShapeSpringsForceField',
        name='rolling_frame_control_spring',
        stiffness=STIFFNESS,
        angularStiffness=STIFFNESS,
        external_rest_shape=beam2_control_mo.getLinkPath(),
        external_points=list(range(NB_FRAMES + 1)),
        mstate='@FramesMO',
        points=list(range(NB_FRAMES + 1)),
        template='Rigid3d',
        activeDirections=[1, 1, 1, 1, 1, 1, 1])
    add_visual_model(beam2_framesNode,beam2_frames,
                     RADIUS2_EX, color=(0.95, 0.85, 0.15, 1.0))
    add_visual_model(beam2_framesNode, beam2_frames,
                     RADIUS2_IN, color=(0.95, 0.85, 0.15, 1.0))

    intersection_node = root_node.addChild('IntersectionNode')

    # ---- GUI bridge ---------------------------------------------------------
    gui_bridge = CTRGuiBridgeStraightOuter(
        root_node=root_node,
        max_tx_m=0.04,
        # The rolling controller owns Y, Z, and axial spin.  The reused GUI's
        # Z-rotation slider is interpreted as the rolling/orbit angle phi.
        max_ty_m=0.0,
        max_tz_m=0.0,
        max_rx_deg=0.0,
        max_rz_deg=720.0,
        init_dt=float(root_node.dt.value),
        default_trans_step_um=50.0,
        default_rot_step_deg=5.0,
        default_control_dt=1e-3,
        dt_min=1e-6, dt_max=1e-1,
        rot_step_max_deg=20.0,
    )

    beam1_MO = beam1_framesNode.FramesMO
    beam2_MO = beam2_framesNode.FramesMO

    ssim=intersection_node.addObject(
        'SphereSweptIntersectionMethod',
        name='ssim',
        beam1Frames=beam1_MO.getLinkPath() + '.position',
        beam2Frames=beam2_MO.getLinkPath() + '.position',
        beam1Velocities = beam1_MO.getLinkPath() + '.velocity',
        beam2Velocities=beam2_MO.getLinkPath() + '.velocity',
        radius1=RADIUS1_EX,
        radius2=RADIUS2_EX,
        innerRadius1 = RADIUS1_IN,
        innerRadius2 = RADIUS2_IN,
        contactConfiguration = "nested",
        defaultNormal = "0 0 -1",
    )

    contact_output = beam1_framesNode.addChild('contactOutput')
    beam2_framesNode.addChild(contact_output)

    contactMO = contact_output.addObject(
        'MechanicalObject', template='Vec3d',
        name='contactMO_gap',
        position=[[0., 0., 0.]] * 2*MAX_K,
    rest_position = [[0., 0., 0.]] * 2 * MAX_K)


    bcm = contact_output.addObject(
        'BeamContactMapping',
        name='bcm',
        input1=beam1_MO.getLinkPath(),
        input2=beam2_MO.getLinkPath(),
        output=contactMO.getLinkPath(),
        ssim = ssim.getLinkPath(),
        mappingMode='contactPoints'
    )

    cpuc = contact_output.addObject(
        'ContactPointsUnilateralConstraint2',
        name='cpuc',
        mu = 0,
        contactTriads = bcm.getLinkPath() + '.contactTriads',
        gapSign = bcm.getLinkPath() + '.gapSign',
        activationTolerance = 1e-4,
        )

    root_node.addObject(CTRController(
        name='CTRController',
        root_node=root_node,
        t1_base_mo=base_mo1,
        t2_base_mo=base_mo2,
        t2_control_mo=beam2_control_mo,
        gui_bridge=gui_bridge,
    ))

    t2_curv_abs_frames = list(
        beam2_framesNode.cosseratMapping.curv_abs_output.value
    )

    intersection_node.addObject(LiveContactMonitor(
        name='LiveMonitor',
        t2_MO=beam2_framesNode.FramesMO,
        bcm=bcm,
        contact_mo=contactMO,
        t2_frame_curv_abs=t2_curv_abs_frames,
        bridge=gui_bridge,
        every_n_steps=20,
        contact_constraint=cpuc,
        force_unit_scale=1e-3,  # scene units are kg-mm-s; kg*mm/s^2 = 1e-3 N
    ))

    intersection_node.addObject(JiggleRecorder(
        name='JiggleRec',
        t2_frames_mo=beam2_framesNode.FramesMO,
        out_path=os.path.join(SCENE_DIR, 'jiggle.csv'),
    ))

    intersection_node.addObject(ContactDebugRecorder(
        name='ContactDebugRecorder',
        t1_frames_mo=beam1_framesNode.FramesMO,
        t2_frames_mo=beam2_framesNode.FramesMO,
        t2_control_mo=beam2_control_mo,
        t2_base_mo=base_mo2,
        bcm=bcm,
        out_path=os.path.join(SCENE_DIR, 'contact_debug.csv'),
    ))
    return root_node
