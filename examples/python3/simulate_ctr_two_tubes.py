# -*- coding: utf-8 -*-
"""
simulate_ctr_two_tubes.py
=========================
SOFA / Cosserat Plugin — Concentric Tube Robot (CTR) with two pre-curved tubes.

Tubes (all SI units)
--------------------
  Tube_1  (outer) : L=17 cm, R_crv=12 cm, arc=60 deg,    OD=3 mm,   ID=2.7 mm
  Tube_2  (inner) : L=33 cm, R_crv=6 cm,  arc=143.2 deg, OD=0.8 mm, ID=0.54 mm

NODE TOPOLOGY (single-parent, mirrors PrecurvedTube exactly)
------------------------------------------------------------
  tube_node  (EulerImplicit + SparseLDL — the solver scope)
    ├── <n>_rigid_base   (Rigid3d base DOF + proximal BC)
    ├── <n>_coss_state   (Vec3d strain DOFs + BeamHooke)
    └── <n>_frames       (child of tube_node ONLY — single parent)
          ├── FramesMO   (Rigid3d output frames)
          ├── DiscreteCosseratMapping
          │     input1 = "@../<n>_coss_state/cosserat_state"
          │     input2 = "@../<n>_rigid_base/cosserat_base_mo"
          └── <n>_visu
                ├── MeshTopology (ring_pos, quads — rest shape)
                ├── visMO        (Vec3d — mechanical target of RigidMapping)
                ├── RigidMapping (FramesMO → visMO, rigidIndexPerPoint)
                └── ogl
                      ├── OglModel       (src=../topo)
                      └── IdentityMapping (visMO → OglModel)

WHY SINGLE-PARENT?
  SOFA's VisualUpdateVisitor traverses a *tree*.  When frame_node has two
  parents (rigid_base AND coss_state), the visitor reaches it via rigid_base,
  processes its visual children, then arrives again via coss_state — and skips
  it (already-visited flag).  The visual update never fires a second time, so
  OglModel positions freeze at t=0.  Making frame_node a child of tube_node
  only (single parent) eliminates the double-visit and the visual chain updates
  every step, exactly as in PrecurvedTube.

  DiscreteCosseratMapping still reaches both inputs via relative paths:
    "@../<n>_coss_state/cosserat_state"   (strain DOFs)
    "@../<n>_rigid_base/cosserat_base_mo" (base frame)
  Both are siblings of frame_node under tube_node, so the paths resolve
  correctly and force back-propagation (applyJT) works as expected.
"""

import math
import Sofa
import Sofa.Core


# =============================================================================
#  TUBE PHYSICAL PARAMETERS
# =============================================================================

T1_PARAMS = {
    'name':          'Tube_1',
    'tube_number':   1,
    'str_length':    0.17,
    'crv_radius':    0.12,
    'crv_angle_deg': 60.0,
    'rex':           15e-4,
    'rin':           13.5e-4,
    'E':             6e10,
    'v':             0.33,
    'density':       6450,
    'nb_sections':   10,
    'nb_frames':     20,
    'color':         [0.15, 0.50, 1.00, 1.0],   # blue
}

T2_PARAMS = {
    'name':          'Tube_3',
    'tube_number':   3,
    'str_length':    0.33,
    'crv_radius':    0.06,
    'crv_angle_deg': 143.239,
    'rex':           4e-4,
    'rin':           2.7e-4,
    'E':             6e10,
    'v':             0.33,
    'density':       6450,
    'nb_sections':   20,
    'nb_frames':     40,
    'color':         [1.00, 0.45, 0.10, 1.0],   # orange
}

N_CIRCLE = 10   # points per cross-sectional ring


# =============================================================================
#  GEOMETRY HELPERS
# =============================================================================

def compute_tube_geometry(p):
    """
    Build Cosserat discretization for a pre-curved tube.

    Layout: [BASE] --- straight (L_str) ---+--- curved arc (L_crv) --- [TIP]

    init_states  = all zeros   → tube starts straight along X at t=0
    rest_states  = kappa on curved sections → BeamHooke drives it to natural shape
    """
    L     = p['str_length']
    R     = p['crv_radius']
    theta = math.radians(p['crv_angle_deg'])
    ns    = p['nb_sections']
    nf    = p['nb_frames']

    L_crv = R * theta
    L_str = max(0.0, L - L_crv)
    kappa = 1.0 / R

    if L_str < 1e-12:
        n_str, n_crv = 0, ns
    else:
        frac_str = L_str / L
        n_str    = max(1, round(frac_str * ns))
        n_crv    = max(1, ns - n_str)
        n_str    = ns - n_crv

    ls = L_str / n_str if n_str > 0 else 0.0
    lc = L_crv / n_crv

    section_lengths = [ls] * n_str + [lc] * n_crv
    rest_states     = [[0., 0., 0.]] * n_str + [[0., 0., kappa]] * n_crv
    init_states     = [[0., 0., 0.]] * ns

    sec_curv_abs = [0.0]
    s = 0.0
    for sl in section_lengths:
        s += sl
        sec_curv_abs.append(round(s, 10))

    lf           = L / nf
    frm_curv_abs = [round(i * lf, 10) for i in range(nf + 1)]
    # All frames along X at t=0 (identity quaternion)
    frame_positions = [[s, 0., 0., 0., 0., 0., 1.] for s in frm_curv_abs]

    # One edge per consecutive frame pair: (0,1), (1,2), ..., (nf-1, nf)
    # nf+1 frames → nf edges
    edge_indices = [[i, i + 1] for i in range(nf)]


    return (section_lengths, rest_states, init_states,
            sec_curv_abs, frame_positions, frm_curv_abs,
            edge_indices)


def tube_mass(p):
    ri, re = p['rin'], p['rex']
    return p['density'] * math.pi * (re**2 - ri**2) * p['str_length']


def build_tube_quads(n_frames, N):
    """Quad faces for a cylindrical surface (n_frames rings × N points)."""
    quads = []
    for i in range(n_frames - 1):
        for j in range(N):
            a = i * N + j
            b = i * N + (j + 1) % N
            c = (i + 1) * N + (j + 1) % N
            d = (i + 1) * N + j
            quads.append([a, b, c, d])
    return quads


def build_ring_positions(r, n_sides, n_frames):
    """
    Rest positions for the cylindrical surface: n_frames copies of a ring of
    n_sides points in the local YZ plane.  RigidMapping applies each frame's
    Rigid3d transform on top of these every timestep.
    """
    TWO_PI = 2.0 * math.pi
    ring = [[ 0.0,
             r * math.cos(TWO_PI * k / n_sides),
             r * math.sin(TWO_PI * k / n_sides)] for k in range(n_sides)]
    return ring * n_frames


def add_collision_model(parent_node, positions, edge_indices, tube_number,tube_params):

    # Create collision node
    collision_node = parent_node.addChild('CollisionNode')

    # Add topology and collision components
    collision_node.addObject(
        'EdgeSetTopologyContainer',
        name="collisEdgeSet",
        position=positions,
        edges=edge_indices
    )
    collision_node.addObject(
        'EdgeSetTopologyModifier',
        name="colliseEdgeModifier"
    )
    collision_node.addObject('MechanicalObject',
                             name="CollisionDOFs",
                             template='Vec3d')

    collision_node.addObject(
        'LineCollisionModel',
        bothSide="1",
        group=str(tube_number),
        proximity=tube_params["rex"]
    )
    collision_node.addObject(
        'PointCollisionModel',
        bothSide="1",
        group=str(tube_number),
    )
    collision_node.addObject('RigidMapping',
                             input='@../FramesMO',
                             output='@CollisionDOFs')

# =============================================================================
#  ADD ONE COSSERAT TUBE
# =============================================================================

def add_cosserat_tube(root_node, p):
    """
    Build the full Cosserat beam hierarchy for one pre-curved tube.

    Critical design: frame_node is a child of tube_node ONLY (single parent).
    See module docstring for the rationale.
    """


    name = p['name']

    (section_lengths, rest_states, init_states,
     sec_curv_abs, frame_positions, frm_curv_abs, edge_indices) = compute_tube_geometry(p)

    re, ri = p['rex'], p['rin']
    nf     = p['nb_frames']
    mass   = tube_mass(p)

    # ── Rigid base (sibling of coss_state and frame_node under tube_node) ───
    tube_node = root_node.addChild(name)
    tube_node.addObject('UncoupledConstraintCorrection')
    rigid_base = tube_node.addChild(name + '_rigid_base')
    base_mo = rigid_base.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='cosserat_base_mo',
        position=[[0., 0., 0., 0., 0., 0., 1.]],
        showObject=True,
        showObjectScale=0.001,
    )
    rigid_base.addObject('PartialFixedProjectiveConstraint',
                         name='proximal_bc',
                         fixedDirections=[0, 1, 1, 0, 1, 1],
                         indices=[0])

    I = math.pi / 4.0 * (p["rex"] ** 4 - p["rin"] ** 4)  # second moment of area
    L_avg = sum(section_lengths) / len(section_lengths)
    compliance = L_avg / (p['E'] * I)



    # ── Cosserat strain state (sibling of rigid_base and frame_node) ─────────
    coss_state = tube_node.addChild(name + '_coss_state')
    coss_mo = coss_state.addObject(
        'MechanicalObject',
        template='Vec3d',
        name='cosserat_state',
        position=init_states,
        rest_position=rest_states,
    )
    coss_state.addObject(
        'BeamHookeLawForceField',
        name='beam_force',
        crossSectionShape='circular',
        length=section_lengths,
        radius=re,
        innerRadius=ri,
        youngModulus=p['E'],
        poissonRatio=p['v'],
        template='Vec3d',
    )

    # ── Output frames — SINGLE parent (tube_node) ────────────────────────────
    # Do NOT addChild from rigid_base or coss_state here.
    # DiscreteCosseratMapping reaches both inputs via relative sibling paths.
    frame_node = tube_node.addChild(name + '_frames')

    frames_mo = frame_node.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='FramesMO',
        position=frame_positions,
        showObject=True,
        showObjectScale=0.001
    )
    frame_node.addObject('UniformMass', name='mass', totalMass=mass)
    frame_node.addObject(
        'DiscreteCosseratMapping',
        name='cosseratMapping',
        curv_abs_input=sec_curv_abs,
        curv_abs_output=frm_curv_abs,
        # Relative paths — both targets are siblings of frame_node under tube_node
        input1='@../' + name + '_coss_state/cosserat_state',
        input2='@../' + name + '_rigid_base/cosserat_base_mo',
        output='@FramesMO',
        debug=False,
        radius=re,
    )

    positions_3d = [[x, y, z] for x, y, z, *_ in frame_positions]
    add_collision_model(parent_node=frame_node,
                        tube_number=p["tube_number"],
                        positions=positions_3d,
                        edge_indices= edge_indices,
                        tube_params = p)

    n_frames_total = nf + 1
    ring_pos  = build_ring_positions(re, N_CIRCLE, n_frames_total)
    rigid_idx = [f for f in range(n_frames_total) for _ in range(N_CIRCLE)]
    quads     = build_tube_quads(n_frames_total, N_CIRCLE)
    color_str = " ".join(str(v) for v in p['color'])

    visu_node = frame_node.addChild(name + '_visu')
    visu_node.addObject('MeshTopology', name='topo',
                        position=ring_pos, quads=quads)
    visu_node.addObject('MechanicalObject', name='visMO',
                        template='Vec3d', position=ring_pos)
    visu_node.addObject('RigidMapping',
                        input='@../FramesMO',
                        output='@visMO',
                        rigidIndexPerPoint=rigid_idx,
                        globalToLocalCoords=False)

    ogl_node = visu_node.addChild('ogl')
    ogl_node.addObject('OglModel', name='oglModel',
                       src='@../topo', color=color_str)
    ogl_node.addObject('IdentityMapping',
                       input='@../visMO',
                       output='@oglModel')

    return base_mo, coss_mo, tube_node


# =============================================================================
#  CONTROLLER
# =============================================================================

class CTRController(Sofa.Core.Controller):
    """
    Two-phase actuation for the CTR.

    Phase 1 (ROTATE_STEPS steps) : rotate  Tube_2 around X
    Phase 2 (TRANS_STEPS  steps) : translate Tube_1 along +X
    """

    ROTATE_STEPS = 200
    TRANS_STEPS  = 200
    ROT_RATE_DEG = 1.0      # deg / step
    TRANS_RATE_M = 5e-4     # m   / step
    TIME_OUT = 200

    def __init__(self, t1_base_mo, t2_base_mo, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.t1_base_mo    = t1_base_mo
        self.t2_base_mo    = t2_base_mo
        self._step         = 0
        self._total_angle  = 0.0
        self._total_tx     = 0.0
        self._rot_rate_rad = math.radians(self.ROT_RATE_DEG)

    def _set_rotation_x(self, mo, angle_rad):
        half = angle_rad * 0.5
        s, c = math.sin(half), math.cos(half)
        with mo.position.writeable() as pos:
            p = list(pos[0])
            p[3], p[4], p[5], p[6] = s, 0.0, 0.0, c
            pos[0] = p

    def _set_translation_x(self, mo, tx):
        with mo.position.writeable() as pos:
            p = list(pos[0])
            p[0] = tx
            pos[0] = p

    def onAnimateBeginEvent(self, event):
        self._step += 1
        n = self._step

        # if 1 <= n <= self.ROTATE_STEPS:
        #     self._total_angle += self._rot_rate_rad
        #     self._set_rotation_x(self.t2_base_mo, self._total_angle)
        #     if n % 50 == 0:
        #         print(f"[CTR | step {n:4d}] Phase 1 – "
        #               f"Tube_2 rotation: {math.degrees(self._total_angle):6.1f} deg")

        # self._set_translation_x(self.t1_base_mo, 0)
        # self._set_translation_x(self.t2_base_mo, 0)
        #
        # if self.TIME_OUT < n <= self.TRANS_STEPS + self.TIME_OUT:
        #     self._total_tx += self.TRANS_RATE_M
        #     self._set_translation_x(self.t1_base_mo, self._total_tx)
        #     phase_step = n - self.ROTATE_STEPS
        #     if phase_step % 50 == 0:
        #         print(f"[CTR | step {n:4d}] Phase 2 – "
        #               f"Tube_1 insertion: {self._total_tx * 1e3:6.1f} mm")
        #
        # elif n == self.ROTATE_STEPS + self.TRANS_STEPS + 1:
        #     print("[CTR] Actuation complete — free simulation continues.")


# =============================================================================
#  CREATE SCENE
# =============================================================================

def createScene(root_node):

    root_node.addObject('RequiredPlugin', pluginName=[
        'Cosserat',
        'Sofa.Component.AnimationLoop',
        'Sofa.Component.Constraint.Projective',
        'Sofa.Component.LinearSolver.Direct',
        'Sofa.Component.Mapping.Linear',        # IdentityMapping
        'Sofa.Component.Mapping.NonLinear',     # RigidMapping
        'Sofa.Component.Mass',
        'Sofa.Component.ODESolver.Backward',
        'Sofa.Component.StateContainer',
        'Sofa.Component.Topology.Container.Constant',  # MeshTopology
        'Sofa.Component.Visual',
        'Sofa.GL.Component.Rendering3D',        # OglModel
        'Sofa.Component.Setting',
        'Sofa.Component.Collision.Detection.Algorithm',  # CollisionPipeline, BruteForceBroadPhase, BVHNarrowPhase
        'Sofa.Component.Collision.Detection.Intersection',  # LocalMinDistance
        'Sofa.Component.Collision.Response.Contact',  # RuleBasedContactManager, FrictionContactConstraint
        'Sofa.Component.Collision.Geometry',  # LineCollisionModel, PointCollisionModel
        'Sofa.Component.Topology.Container.Dynamic',  # EdgeSetTopologyContainer, EdgeSetTopologyModifier
        'Sofa.Component.Constraint.Lagrangian.Solver',  # GenericConstraintSolver
        'Sofa.Component.Constraint.Lagrangian.Correction',
        'MultiThreading',
    ])

    root_node.gravity = [0., 0., 0.]
    root_node.dt      = 0.001

    root_node.addObject('DefaultVisualManagerLoop')
    root_node.addObject('FreeMotionAnimationLoop')
    root_node.addObject('BackgroundSetting', color=[1.0,1.0,1.0,0])

    root_node.addObject('Camera',
                        position=[2.5, -1.2, 1.0],
                        lookAt=[1.2, 0.0, 0.0])

    root_node.addObject('VisualStyle',
                        displayFlags='showVisualModels hideBehaviorModels '
                                     'hideCollisionModels '
                                     'hideBoundingCollisionModels '
                                     'hideForceFields '
                                     'hideInteractionForceFields '
                                     'hideWireframe '
                                     'hideMechanicalMappings')

    root_node.addObject('CollisionPipeline')
    root_node.addObject('BruteForceBroadPhase')
    root_node.addObject('BVHNarrowPhase')
    root_node.addObject('LocalMinDistance',
                        alarmDistance=4e-4, contactDistance=1e-4)

    root_node.addObject('RuleBasedContactManager',
                        response='FrictionContactConstraint',
                        responseParams='mu=0.2')
    root_node.addObject('BlockGaussSeidelConstraintSolver', name='ConstraintSolver', tolerance=1e-5,
                                 maxIterations=5e2)
    root_node.addObject('EulerImplicitSolver',
                        name='odesolver',
                        rayleighStiffness=0.2,
                        rayleighMass=0.1,
                        firstOrder=False)

    root_node.addObject('SparseLDLSolver',
                        name='Solver',
                        template='CompressedRowSparseMatrixd')


    t1_base_mo, _, _ = add_cosserat_tube(root_node, T1_PARAMS)
    t2_base_mo, _, _ = add_cosserat_tube(root_node, T2_PARAMS)

    root_node.addObject(CTRController(
        name='CTRController',
        t1_base_mo=t1_base_mo,
        t2_base_mo=t2_base_mo,
    ))

    for p in (T1_PARAMS, T2_PARAMS):
        L_crv = p['crv_radius'] * math.radians(p['crv_angle_deg'])
        L_str = p['str_length'] - L_crv
        print(f"[CTR] {p['name']:8s}  "
              f"L={p['str_length']*1e3:.0f} mm  "
              f"L_str={L_str*1e3:.1f} mm  "
              f"L_crv={L_crv*1e3:.1f} mm  "
              f"kappa={1/p['crv_radius']:.2f} rad/m  "
              f"mass={tube_mass(p)*1e3:.2f} g")

    return root_node