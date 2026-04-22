# -*- coding: utf-8 -*-
"""
simulate_ctr_two_tubes.py
=========================
SOFA / Cosserat Plugin — Concentric Tube Robot (CTR) with two pre-curved tubes.

Tubes (all SI units)
--------------------
  Tube_1  (outer) : L=17 cm, R_crv=12 cm, arc=60 deg,    OD=3 mm,   ID=2.7 mm
  Tube_3  (inner) : L=33 cm, R_crv=6 cm,  arc=143.2 deg, OD=0.8 mm, ID=0.54 mm

INITIAL CONFIGURATION — CONCENTRIC PLACEMENT
---------------------------------------------
  The global frame sits at Tube_1's rigid base (world origin).
  Tube_1 takes its natural configuration freely from t=0 (no change).

  Tube_3 is the longer tube and is placed behind Tube_1 by X mm along -X so
  that no part of Tube_3 protrudes beyond Tube_1's tip at t=0.  This avoids
  the ill-posed question of what shape the protruding section should adopt.

  Offset derivation (both tubes straight along +X at t=0):
    Tube_1 base @ x=0,    tip @ x = +L1
    Tube_3 base @ x=-X,   tip @ x = -X + L3

  Tip-coincidence condition  ->  X = L3 - L1  (> 0 because L3 > L1)

  With the parameters below:
    X = 0.33 - 0.17 = 0.16 m = 160 mm
    Tube_3 rigid base at x = -0.16 m

  compute_concentric_offset() encodes this computation and is called once
  in createScene().  Only Tube_3 receives the offset; Tube_1 is untouched.

NODE TOPOLOGY (single-parent, mirrors PrecurvedTube exactly)
------------------------------------------------------------
  tube_node  (EulerImplicit + SparseLDL -- the solver scope)
    +-- <n>_rigid_base   (Rigid3d base DOF + proximal BC)
    +-- <n>_coss_state   (Vec3d strain DOFs + BeamHooke)
    +-- <n>_frames       (child of SolverNode ONLY -- single parent)
          +-- FramesMO   (Rigid3d output frames)
          +-- UniformMass
          +-- DiscreteCosseratMapping
                input1 = "@../<n>_coss_state/cosserat_state"
                input2 = "@../<n>_rigid_base/cosserat_base_mo"
          +-- <n>_visu
                +-- MeshTopology (ring_pos, quads -- rest shape)
                +-- visMO        (Vec3d -- mechanical target of RigidMapping)
                +-- RigidMapping (FramesMO -> visMO, rigidIndexPerPoint)
                +-- ogl
                      +-- OglModel       (src=../topo)
                      +-- IdentityMapping (visMO -> OglModel)

WHY SINGLE-PARENT?
  SOFA's VisualUpdateVisitor traverses a *tree*.  When frame_node has two
  parents (rigid_base AND coss_state), the visitor reaches it via rigid_base,
  processes its visual children, then arrives again via coss_state -- and skips
  it (already-visited flag).  The visual update never fires a second time, so
  OglModel positions freeze at t=0.  Making frame_node a child of SolverNode
  only (single parent) eliminates the double-visit and the visual chain updates
  every step, exactly as in PrecurvedTube.

WHY UncoupledConstraintCorrection (not GenericConstraintCorrection)?
  BeamHookeLawForceField acts on Vec3d strain DOFs and does NOT implement
  addKToMatrix().  SparseLDLSolver therefore assembles a degenerate stiffness
  matrix, causing GenericConstraintCorrection to silently produce wrong
  compliance values and corrupt the constraint solve.
  UncoupledConstraintCorrection uses a diagonal compliance approximation
  (L_avg / (E*I) per section) and bypasses the degenerate assembly entirely.

CUSTOM CONTACT PIPELINE (replaces classic SOFA collision pipeline entirely)
--------------------------------------------------------------------------
  The entire classic pipeline (CollisionPipeline, BruteForceBroadPhase,
  BVHNarrowPhase, LocalMinDistance, RuleBasedContactManager,
  LineCollisionModel, PointCollisionModel) is removed.  Contact is handled by:

    SphereSweptIntersectionMethod (SSIM)  [DataEngine]
      For CTR internal contact: radius1 = rin_1 (inner wall of Tube_1),
      radius2 = rex_3 (outer wall of Tube_3).
      gap = rin_1 - d_centreline - rex_3
         > 0  ->  clearance   (no contact)
         = 0  ->  touching
         < 0  ->  penetrating

      NOTE: SSIM must be configured for internal-contact mode so that it
      computes gap = r1 - d - r2 (not the external formula d - r1 - r2).
      The radius arguments passed here are the physically correct surfaces;
      the sign convention in SphereSweptIntersectionMethod.cpp must match.

    BeamContactMapping (BCM)  [Multi2Mapping<Rigid3d,Rigid3d,Vec3d>]
      mode = 'gap'  (single output MO, avoids SceneCheckMapping conflict).
      input1 = Tube_1/SolverNode/Tube_1_frames/FramesMO
      input2 = Tube_3/SolverNode/Tube_3_frames/FramesMO
      out[k] = delta[k] = Pc_B[k] - Pc_A[k]

    UnilateralLagrangianConstraint<Vec3d>  (ULC)
      object1 = contactMO_ref  (fixed zero, MAX_K DOFs)
      object2 = contactMO_gap  (BCM output)
      Violation: dfree = delta_n_free[k]; enforces delta_n >= 0.

    ContactFeeder  [BaseObject + AnimateBeginEvent]
      Reads SSIM distances -> clear() + addContact() on ULC each step.
      Activates pair k when distances[k][0] < ALARM_DISTANCE.

  Scene graph
  -----------
  root
  +-- FreeMotionAnimationLoop
  +-- BlockGaussSeidelConstraintSolver
  +-- Tube_1/
  |   +-- SolverNode/
  |       +-- EulerImplicitSolver + SparseLDLSolver + UncoupledConstraintCorrection
  |       +-- Tube_1_rigid_base/  cosserat_base_mo + PartialFixedProjectiveConstraint
  |       +-- Tube_1_coss_state/  cosserat_state + BeamHookeLawForceField
  |       +-- Tube_1_frames/      FramesMO + UniformMass + DiscreteCosseratMapping
  |           +-- Tube_1_visu/    MeshTopology + visMO + RigidMapping + ogl/
  +-- Tube_3/
  |   +-- SolverNode/
  |       +-- EulerImplicitSolver + SparseLDLSolver + UncoupledConstraintCorrection
  |       +-- Tube_3_rigid_base/  cosserat_base_mo + PartialFixedProjectiveConstraint
  |       +-- Tube_3_coss_state/  cosserat_state + BeamHookeLawForceField
  |       +-- Tube_3_frames/      FramesMO + UniformMass + DiscreteCosseratMapping
  |           +-- Tube_3_visu/    MeshTopology + visMO + RigidMapping + ogl/
  +-- contact_node/
      +-- SphereSweptIntersectionMethod  (ssim)
      +-- contactMO_ref  (Vec3d, MAX_K zero DOFs -- ULC object1 zero reference)
      +-- contactMO_gap  (Vec3d, K DOFs    -- BCM sole output, delta[k] = Pc_B-Pc_A)
      +-- BeamContactMapping  (bcm)  mappingMode='gap'
      +-- UnilateralLagrangianConstraint  (ulc)
      +-- ContactFeeder  (feeder)
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
    'str_length':    0.17,        # total arc length [m]
    'crv_radius':    0.12,        # curvature radius [m]
    'crv_angle_deg': 60.0,
    'rex':           15e-4,       # outer radius [m]
    'rin':           13.5e-4,     # inner (lumen) radius [m]
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
    'str_length':    0.33,        # total arc length [m]
    'crv_radius':    0.06,        # curvature radius [m]
    'crv_angle_deg': 143.239,
    'rex':           4e-4,        # outer radius [m]
    'rin':           2.7e-4,      # inner (lumen) radius [m]
    'E':             6e10,
    'v':             0.33,
    'density':       6450,
    'nb_sections':   20,
    'nb_frames':     40,
    'color':         [1.00, 0.45, 0.10, 1.0],   # orange
}

N_CIRCLE = 10   # points per cross-sectional ring for visual model


# =============================================================================
#  CONTACT PIPELINE PARAMETERS
# =============================================================================

# "ALGO_1": segment-to-segment closest pair (NB_SEC_1 x NB_SEC_3 pairs max)
# "ALGO_2": node-to-segment Newton-Raphson  ((NB_FRM_1+1) x NB_SEC_3 pairs max)
ALGORITHM = "ALGO_1"

# Pre-allocated size of contactMO_ref (zero-reference MO for ULC).
# Must be >= maximum number of contact pairs SSIM can ever produce.
# ALGO_1 upper bound: nb_sections_1 x nb_sections_3
# ALGO_2 upper bound: (nb_frames_1 + 1) x nb_sections_3
MAX_K = max(
    T1_PARAMS['nb_sections'] * T2_PARAMS['nb_sections'],          # ALGO_1: 10x20 = 200
    (T1_PARAMS['nb_frames'] + 1) * T2_PARAMS['nb_sections'],      # ALGO_2: 21x20 = 420
)   # -> 420

# CTR internal contact geometry:
#   The relevant surfaces are the inner wall of Tube_1 (radius = rin_1)
#   and the outer wall of Tube_3 (radius = rex_3).
#   Physical wall-to-wall clearance when coaxial:
#     gap_wall = rin_1 - rex_3 = 13.5e-4 - 4e-4 = 9.5e-4 m
#
# ContactFeeder activates a constraint for pair k when gap[k] < ALARM_DISTANCE.
# Setting ALARM_DISTANCE = gap_wall + margin catches near-contact before
# penetration: feeder activates as soon as any eccentricity consumes the gap.
_GAP_WALL      = T1_PARAMS['rin'] - T2_PARAMS['rex']   # 9.5e-4 m
ALARM_DISTANCE = _GAP_WALL + 0.5e-4                    # 1.0e-3 m  (0.5 mm margin)


# =============================================================================
#  GEOMETRY HELPERS
# =============================================================================

def compute_concentric_offset(p_outer, p_inner):
    """
    Compute the signed x-offset [m] to apply to the inner tube's rigid base
    so that the inner tube is fully housed inside the outer tube at t=0,
    with their tips coinciding along the X-axis.

    Geometry (both tubes straight along +X at t=0)
    -----------------------------------------------
      Outer (Tube_1):  base @ x = 0          tip @ x = +L_outer
      Inner (Tube_3):  base @ x = -X         tip @ x = -X + L_inner

      Tip-coincidence:  -X + L_inner = L_outer
                         X = L_inner - L_outer   (> 0 for a valid CTR)

    With the parameters in this file:
      L_outer = 0.17 m,  L_inner = 0.33 m
      X       = 0.16 m = 160 mm
      -> inner rigid base placed at x = -0.16 m

    Note: 'str_length' in the param dict is the TOTAL arc length
    (straight + curved), which is exactly what is needed here.

    Parameters
    ----------
    p_outer : dict   Param dict of the shorter (outer) tube -- T1_PARAMS.
    p_inner : dict   Param dict of the longer  (inner) tube -- T2_PARAMS.

    Returns
    -------
    float  Signed x-offset for the inner tube's rigid base [m].
           Always negative (inner base is behind outer base).

    Raises
    ------
    ValueError  If the inner tube is not strictly longer than the outer tube.
    """
    L_outer = p_outer['str_length']
    L_inner = p_inner['str_length']

    if L_inner <= L_outer:
        raise ValueError(
            f"Inner tube '{p_inner['name']}' (L = {L_inner * 1e3:.1f} mm) must be "
            f"strictly longer than outer tube '{p_outer['name']}' "
            f"(L = {L_outer * 1e3:.1f} mm) for a valid CTR concentric placement.  "
            f"Got X = {(L_inner - L_outer) * 1e3:.1f} mm <= 0."
        )

    X = L_inner - L_outer      # 0.16 m for this scene
    print(
        f"[CTR] Concentric offset: X = {X * 1e3:.1f} mm  "
        f"-- '{p_inner['name']}' rigid base at x = {-X * 1e3:.1f} mm  "
        f"(tips coincide at x = +{L_outer * 1e3:.1f} mm)"
    )
    return -X   # negative: inner base sits behind the outer base


def compute_tube_geometry(p, x_offset=0.0):
    """
    Build Cosserat discretization for a pre-curved tube.

    Layout: [BASE] --- straight (L_str) ---+--- curved arc (L_crv) --- [TIP]

    init_states  = all zeros   -> tube starts straight along X at t=0
    rest_states  = kappa on curved sections -> BeamHooke drives it to natural shape

    Parameters
    ----------
    p        : dict   Tube parameter dictionary.
    x_offset : float  World-frame x-translation of the rigid base [m].
                      0.0 for Tube_1 (default -- untouched).
                      compute_concentric_offset(T1_PARAMS, T2_PARAMS) for Tube_3.

    Notes on x_offset
    ------------------
    sec_curv_abs and frm_curv_abs are intrinsic arc-length coordinates measured
    from the tube's own base.  They are consumed by DiscreteCosseratMapping and
    must always start at 0 -- they are NEVER offset.

    Only frame_positions (the world-frame Rigid3d initial positions fed to
    FramesMO) and the rigid base position in cosserat_base_mo are shifted by
    x_offset.  After the first DiscreteCosseratMapping::apply() call these
    positions are overwritten by the mapping; the offset ensures geometric
    consistency during SOFA's init() phase.
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

    # Intrinsic arc-length coordinates -- NOT offset, always start at 0.
    sec_curv_abs = [0.0]
    s = 0.0
    for sl in section_lengths:
        s += sl
        sec_curv_abs.append(round(s, 10))

    lf           = L / nf
    frm_curv_abs = [round(i * lf, 10) for i in range(nf + 1)]

    # World-frame initial positions: x_offset shifts Tube_3's base behind Tube_1.
    # frm_curv_abs[i] is the intrinsic arc length from the tube's own base;
    # world x-coordinate of frame i = x_offset + frm_curv_abs[i].
    frame_positions = [[x_offset + s, 0., 0., 0., 0., 0., 1.]
                       for s in frm_curv_abs]

    # One edge per consecutive frame pair: (0,1), (1,2), ..., (nf-1, nf)
    edge_indices = [[i, i + 1] for i in range(nf)]

    return (section_lengths, rest_states, init_states,
            sec_curv_abs, frame_positions, frm_curv_abs,
            edge_indices)


def tube_mass(p):
    ri, re = p['rin'], p['rex']
    return p['density'] * math.pi * (re**2 - ri**2) * p['str_length']


def build_tube_quads(n_frames, N):
    """Quad faces for a cylindrical surface (n_frames rings x N points)."""
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
    ring = [[0.0,
             r * math.cos(TWO_PI * k / n_sides),
             r * math.sin(TWO_PI * k / n_sides)] for k in range(n_sides)]
    return ring * n_frames


# =============================================================================
#  ADD ONE COSSERAT TUBE  (no classic collision model)
# =============================================================================

def add_cosserat_tube(root_node, p, x_offset=0.0):
    """
    Build the full Cosserat beam hierarchy for one pre-curved tube.

    The classic CollisionNode (LineCollisionModel / PointCollisionModel /
    EdgeSetTopologyContainer) is NOT added here; contact is handled entirely
    by the SSIM custom pipeline in contact_node.

    Critical design: frame_node is a child of SolverNode ONLY (single parent).
    See module docstring for the rationale.

    Parameters
    ----------
    root_node : Sofa node  Scene root.
    p         : dict       Tube parameter dictionary.
    x_offset  : float      World-frame x-position of the rigid base [m].
                           0.0 (default) for Tube_1 -- no change.
                           Pass compute_concentric_offset(T1_PARAMS, T2_PARAMS)
                           for Tube_3 to place it concentrically behind Tube_1.

    Returns
    -------
    base_mo    : MechanicalObject<Rigid3d>   Rigid base DOF (for controller).
    coss_mo    : MechanicalObject<Vec3d>     Cosserat strain DOFs.
    tube_node  : Sofa.Core.Node              Top-level tube node.
    frame_node : Sofa.Core.Node              Node holding FramesMO -- passed to
                                             SSIM/BCM for contact wiring.
    """
    name = p['name']

    # edge_indices unused (no CollisionNode in this pipeline).
    (section_lengths, rest_states, init_states,
     sec_curv_abs, frame_positions, frm_curv_abs,
     _edge_indices) = compute_tube_geometry(p, x_offset)

    re, ri = p['rex'], p['rin']
    nf     = p['nb_frames']
    mass   = tube_mass(p)

    # ---- Per-section compliance for UncoupledConstraintCorrection -----------
    # BeamHookeLawForceField does NOT implement addKToMatrix(), so
    # GenericConstraintCorrection cannot assemble a valid compliance matrix.
    # UncoupledConstraintCorrection with an explicit diagonal compliance value
    # bypasses the degenerate SparseLDLSolver assembly entirely.
    I_sec      = math.pi / 4.0 * (re ** 4 - ri ** 4)   # second moment of area [m^4]
    L_avg      = sum(section_lengths) / len(section_lengths)
    compliance = L_avg / (p['E'] * I_sec)               # [1 / (Pa * m^3)]

    # ---- Solver scope --------------------------------------------------------
    tube_node   = root_node.addChild(name)
    solver_node = tube_node.addChild('SolverNode')
    solver_node.addObject('EulerImplicitSolver',
                          name='odesolver',
                          rayleighStiffness=0.2,
                          rayleighMass=0.1,
                          firstOrder=False)
    solver_node.addObject('SparseLDLSolver',
                          name='Solver',
                          template='CompressedRowSparseMatrixd')

    solver_node.addObject('GenericConstraintCorrection')

    # ---- Rigid base ----------------------------------------------------------
    rigid_base = solver_node.addChild(name + '_rigid_base')
    base_mo = rigid_base.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='cosserat_base_mo',
        # x_offset applied so Tube_3's kinematic root starts at x = -160 mm.
        # Tube_1 receives x_offset = 0.0 -> position unchanged.
        position=[[x_offset, 0., 0., 0., 0., 0., 1.]],
        showObject=True,
        showObjectScale=0.001,
    )
    # Allow translation along X (insertion) and rotation around X (axial twist).
    # Fix Y/Z translation and Y/Z rotation.
    rigid_base.addObject('PartialFixedProjectiveConstraint',
                         name='proximal_bc',
                         fixedDirections=[0, 1, 1, 0, 1, 1],
                         indices=[0])

    # ---- Cosserat strain state -----------------------------------------------
    coss_state = solver_node.addChild(name + '_coss_state')
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

    # ---- Output frames -- SINGLE parent (SolverNode) ------------------------
    # Do NOT addChild from rigid_base or coss_state here.
    # DiscreteCosseratMapping reaches both inputs via relative sibling paths.
    frame_node = solver_node.addChild(name + '_frames')

    frame_node.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='FramesMO',
        # frame_positions already carry x_offset (set in compute_tube_geometry).
        position=frame_positions,
        showObject=True,
        showObjectScale=0.001,
    )
    frame_node.addObject('UniformMass', name='mass', totalMass=mass)
    frame_node.addObject(
        'DiscreteCosseratMapping',
        name='cosseratMapping',
        # curv_abs_input / curv_abs_output are intrinsic arc lengths measured
        # from the tube's own base -- they are NEVER offset.
        curv_abs_input=sec_curv_abs,
        curv_abs_output=frm_curv_abs,
        input1='@../' + name + '_coss_state/cosserat_state',
        input2='@../' + name + '_rigid_base/cosserat_base_mo',
        output='@FramesMO',
        debug=False,
        radius=re,
    )

    # ---- Visual model --------------------------------------------------------
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

    # frame_node returned so createScene can wire FramesMO into SSIM/BCM.
    return base_mo, coss_mo, tube_node, frame_node


# =============================================================================
#  CONTROLLER
# =============================================================================

class CTRController(Sofa.Core.Controller):
    """
    Two-phase actuation for the CTR.

    Phase 1 (ROTATE_STEPS steps) : rotate  Tube_3 around X
    Phase 2 (TRANS_STEPS  steps) : insert  Tube_3 along +X

    t2_x_offset stores Tube_3's retracted start position so that
    _set_translation_x() writes  absolute = rest_offset + delta  rather than
    just delta.  Without this, the first controller translation would snap
    Tube_3's base back to the world origin.
    """

    ROTATE_STEPS = 200
    TRANS_STEPS  = 200
    ROT_RATE_DEG = 1.0      # deg / step
    TRANS_RATE_M = 5e-4     # m   / step
    TIME_OUT     = 200

    def __init__(self, t1_base_mo, t2_base_mo, t2_x_offset=0.0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.t1_base_mo    = t1_base_mo
        self.t2_base_mo    = t2_base_mo
        self._t2_x0        = t2_x_offset   # Tube_3's resting base position [m]
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

    def _set_translation_x(self, mo, tx, x0=0.0):
        """Write absolute x = rest_offset + insertion_delta to the base MO."""
        with mo.position.writeable() as pos:
            p = list(pos[0])
            p[0] = x0 + tx
            pos[0] = p

    def onAnimateBeginEvent(self, event):
        self._step += 1
        n = self._step

        # if 1 <= n <= self.ROTATE_STEPS:
        #     self._total_angle += self._rot_rate_rad
        #     self._set_rotation_x(self.t2_base_mo, self._total_angle)
        #     if n % 50 == 0:
        #         print(f"[CTR | step {n:4d}] Phase 1 - "
        #               f"Tube_3 rotation: {math.degrees(self._total_angle):6.1f} deg")

        # self._set_translation_x(self.t1_base_mo, 0, x0=0.0)
        # self._set_translation_x(self.t2_base_mo, 0, x0=self._t2_x0)
        #
        # if self.TIME_OUT < n <= self.TRANS_STEPS + self.TIME_OUT:
        #     self._total_tx += self.TRANS_RATE_M
        #     # x0=self._t2_x0: insertion measured from Tube_3's retracted start.
        #     self._set_translation_x(self.t2_base_mo, self._total_tx,
        #                             x0=self._t2_x0)
        #     phase_step = n - self.ROTATE_STEPS
        #     if phase_step % 50 == 0:
        #         print(f"[CTR | step {n:4d}] Phase 2 - "
        #               f"Tube_3 insertion: {self._total_tx * 1e3:6.1f} mm")
        #
        # elif n == self.ROTATE_STEPS + self.TRANS_STEPS + 1:
        #     print("[CTR] Actuation complete -- free simulation continues.")


# =============================================================================
#  CREATE SCENE
# =============================================================================

def createScene(root_node):

    # ---- Required plugins ---------------------------------------------------
    # Classic collision pipeline plugins (CollisionDetection.Algorithm,
    # CollisionDetection.Intersection, Collision.Response.Contact,
    # Collision.Geometry, Topology.Container.Dynamic, MultiThreading)
    # are removed -- they are fully replaced by the SSIM custom pipeline.
    root_node.addObject('RequiredPlugin', pluginName=[
        'Cosserat',                                           # DiscreteCosseratMapping, BeamHookeLawForceField,
                                                              # SphereSweptIntersectionMethod,
                                                              # BeamContactMapping, ContactFeeder
        'Sofa.Component.AnimationLoop',                       # FreeMotionAnimationLoop
        'Sofa.Component.Constraint.Lagrangian.Correction',   # UncoupledConstraintCorrection
        'Sofa.Component.Constraint.Lagrangian.Model',        # UnilateralLagrangianConstraint
        'Sofa.Component.Constraint.Lagrangian.Solver',       # BlockGaussSeidelConstraintSolver
        'Sofa.Component.Constraint.Projective',              # PartialFixedProjectiveConstraint
        'Sofa.Component.LinearSolver.Direct',                # SparseLDLSolver
        'Sofa.Component.Mapping.Linear',                     # IdentityMapping
        'Sofa.Component.Mapping.NonLinear',                  # RigidMapping
        'Sofa.Component.Mass',                               # UniformMass
        'Sofa.Component.ODESolver.Backward',                 # EulerImplicitSolver
        'Sofa.Component.Setting',                            # BackgroundSetting
        'Sofa.Component.StateContainer',                     # MechanicalObject
        'Sofa.Component.Topology.Container.Constant',        # MeshTopology
        'Sofa.Component.Visual',                             # VisualStyle
        'Sofa.GL.Component.Rendering3D',                     # OglModel
        'Sofa.GUI.Component',                                # Camera
    ])

    root_node.gravity = [0., 0., 0.]
    root_node.dt      = 0.001

    root_node.addObject('DefaultVisualManagerLoop')
    root_node.addObject('FreeMotionAnimationLoop')
    root_node.addObject('BackgroundSetting', color=[1.0, 1.0, 1.0, 0])

    # BlockGaussSeidelConstraintSolver must sit at root so its VecIds are
    # visible to all MechanicalObjects during FreeMotionAnimationLoop's
    # global visitors.
    root_node.addObject('BlockGaussSeidelConstraintSolver',
                        name='ConstraintSolver',
                        tolerance=1e-5,
                        maxIterations=500)

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

    # ---- Compute Tube_3's retraction offset ---------------------------------
    # X = L_inner - L_outer = 0.33 - 0.17 = 0.16 m = 160 mm
    # Tube_3's rigid base sits at x = -0.16 m so at t=0 (both tubes straight)
    # the tips coincide at x = +0.17 m.  Tube_1 is entirely unmodified.
    x_t3 = compute_concentric_offset(T1_PARAMS, T2_PARAMS)   # -> -0.16 m

    # ---- Build tubes --------------------------------------------------------
    # add_cosserat_tube now returns frame_node as the 4th value so we can
    # wire FramesMO into SSIM and BCM in the contact_node below.
    t1_base_mo, _, _, t1_frame_node = add_cosserat_tube(root_node, T1_PARAMS)
    t2_base_mo, _, _, t2_frame_node = add_cosserat_tube(root_node, T2_PARAMS,
                                                        x_offset=x_t3)

    root_node.addObject(CTRController(
        name='CTRController',
        t1_base_mo=t1_base_mo,
        t2_base_mo=t2_base_mo,
        t2_x_offset=x_t3,
    ))

    # =========================================================================
    #  Custom contact pipeline
    # =========================================================================
    contact_node = root_node.addChild('contact_node')

    # ---- FramesMO handles ---------------------------------------------------
    # t{1,2}_frame_node is the '..._frames' node returned by add_cosserat_tube.
    # .FramesMO accesses the MechanicalObject<Rigid3d> inside it by name.
    t1_MO = t1_frame_node.FramesMO
    t3_MO = t2_frame_node.FramesMO

    # ---- SphereSweptIntersectionMethod (SSIM) --------------------------------
    #
    # CTR internal contact geometry
    # --------------------------------
    # Tube_3 (inner) runs coaxially INSIDE Tube_1 (outer).  Contact occurs
    # between the INNER WALL of Tube_1 (rin_1 = 13.5e-4 m) and the OUTER WALL
    # of Tube_3 (rex_3 = 4e-4 m).
    #
    # For internal contact the gap formula is:
    #   gap_internal = rin_1 - d_centreline - rex_3
    #      > 0  ->  clearance (no contact)
    #      = 0  ->  touching
    #      < 0  ->  penetrating
    #
    # This is the opposite sign of the external formula (d - r1 - r2).
    # SSIM must be set to internal-contact mode so it computes r1 - d - r2.
    # The radius arguments below are physically correct regardless of mode:
    #   radius1 = rin_1  (inner wall of outer tube -- the lumen surface)
    #   radius2 = rex_3  (outer wall of inner tube -- the tube surface)
    #
    # If '.position' suffix causes a link error, try beam1Frames=t1_MO.getLinkPath()
    contact_node.addObject(
        'SphereSweptIntersectionMethod',
        name='ssim',
        beam1Frames=t1_MO.getLinkPath() + '.position',
        beam2Frames=t3_MO.getLinkPath() + '.position',
        radius1=T1_PARAMS['rin'],      # inner wall of Tube_1 [m]
        radius2=T2_PARAMS['rex'],      # outer wall of Tube_3 [m]
        algorithmType=ALGORITHM,
    )

    # ---- Output MOs for BeamContactMapping (gap mode) -----------------------
    #
    # contactMO_ref -- fixed zero-reference MO (NOT an output of BCM).
    #   Pre-allocated with MAX_K = 420 zero DOFs.
    #   ContactFeeder.addContact(params, n, 0, k, k) reads ref[k] = [0,0,0]
    #   and gap[k] = delta[k].  As long as k < MAX_K the access is always valid.
    #   The reaction force -lambda*n applied by ULC to contactMO_ref is
    #   discarded (no solver, no mass); Newton's 3rd law is already embedded in
    #   BCM's antisymmetric Jacobian: d_delta/d_q = d_Pc_B/d_q2 - d_Pc_A/d_q1.
    #
    # contactMO_gap -- BCM's sole output (gap mode).
    #   BCM.apply() resizes this to K DOFs each step and writes:
    #     contactMO_gap[k] = delta[k] = Pc_B[k] - Pc_A[k]
    #   BCM.applyJT distributes contact forces antisymmetrically to both tubes:
    #     Tube_1 frame forces: -weight * lambda * n  (push inward)
    #     Tube_3 frame forces: +weight * lambda * n  (push outward)
    contactMO_ref = contact_node.addObject(
        'MechanicalObject', template='Vec3d',
        name='contactMO_ref',
        position=[[0., 0., 0.]] * MAX_K)   # pre-allocated zeros; never resized

    contactMO_gap = contact_node.addObject(
        'MechanicalObject', template='Vec3d',
        name='contactMO_gap')              # BCM resizes this each step to K DOFs

    # ---- BeamContactMapping (BCM) -------------------------------------------
    #
    # input1 = Tube_1 FramesMO,  input2 = Tube_3 FramesMO  (Rigid3d).
    # output = contactMO_gap                                (Vec3d, gap mode).
    #
    # contactSectionIds / curvilinearParams: Data-to-Data links to SSIM outputs.
    # Reading these triggers SSIM.doUpdate() lazily when frames are dirty.
    #
    # radius1 / radius2: same physical surfaces as passed to SSIM above.
    # mappingMode = 'gap': BCM writes delta[k] = Pc_B[k] - Pc_A[k] and
    # back-propagates forces to BOTH FramesMOs without a second output MO,
    # avoiding the SceneCheckMapping / MappingGraph diamond-graph conflict.
    contact_node.addObject(
        'BeamContactMapping',
        name='bcm',
        input1=t1_MO.getLinkPath(),
        input2=t3_MO.getLinkPath(),
        output=contactMO_gap.getLinkPath(),
        radius1=T1_PARAMS['rin'],          # inner wall of Tube_1
        radius2=T2_PARAMS['rex'],          # outer wall of Tube_3
        isAlgo2=ALGORITHM == 'ALGO_2',
        mappingMode='gap',
        contactSectionIds='@ssim.contactSectionIds',
        curvilinearParams='@ssim.curvilinearParams',
    )

    # ---- UnilateralLagrangianConstraint (ULC) --------------------------------
    #
    # object1 = contactMO_ref  (fixed zero, MAX_K zero DOFs)
    # object2 = contactMO_gap  (BCM sole output, K gap-vector DOFs)
    #
    # Gap violation at pair k:
    #   dfree = dot(gap_free[k] - ref[k], n) = dot(delta_free[k], n) = delta_n_free[k]
    #   dfree < 0  =>  penetration  =>  lambda >= 0 enforced
    #   dfree >= 0 =>  clearance    =>  lambda = 0 (no contact force)
    #
    # 'object1'/'object2' are the standard PairInteractionConstraint link aliases.
    # Use 'mstate1'/'mstate2' if your SOFA version uses those names instead.
    contact_node.addObject(
        'UnilateralLagrangianConstraint',
        template='Vec3d',
        name='ulc',
        object1=contactMO_ref.getLinkPath(),
        object2=contactMO_gap.getLinkPath(),
    )

    # ---- ContactFeeder -------------------------------------------------------
    #
    # surfacePoints1 / surfacePoints2: SSIM surface-point outputs (Pc_A, Pc_B).
    # distances:   SSIM distances output; distances[k][0] = signed gap delta_n.
    #   delta_n < 0               -> penetrating
    #   0 <= delta_n < ALARM_DIST -> near-contact -> constraint activated
    #   delta_n >= ALARM_DIST     -> far -> no constraint
    #
    # ALARM_DISTANCE = wall-to-wall clearance (9.5e-4 m) + 0.5e-4 m margin
    # so the feeder activates as soon as any eccentricity consumes the gap.
    #
    # mu=0.0: frictionless -- UnilateralConstraintResolution (1 DOF per contact).
    # Fires at AnimateBeginEvent, before FreeMotionAnimationLoop runs.
    # SSIM.doUpdate() is triggered here using the t_n frame positions.
    contact_node.addObject(
        'ContactFeeder',
        name='feeder',
        surfacePoints1='@ssim.surfacePoints1',
        surfacePoints2='@ssim.surfacePoints2',
        distances='@ssim.distances',
        constraint='@ulc',
        alarmDistance=ALARM_DISTANCE,
        mu=0.0,
    )

    # ---- Diagnostic print ---------------------------------------------------
    for p in (T1_PARAMS, T2_PARAMS):
        L_crv = p['crv_radius'] * math.radians(p['crv_angle_deg'])
        L_str = p['str_length'] - L_crv
        print(f"[CTR] {p['name']:8s}  "
              f"L={p['str_length']*1e3:.0f} mm  "
              f"L_str={L_str*1e3:.1f} mm  "
              f"L_crv={L_crv*1e3:.1f} mm  "
              f"kappa={1/p['crv_radius']:.2f} rad/m  "
              f"mass={tube_mass(p)*1e3:.2f} g")
    print(f"[CTR] ALGORITHM={ALGORITHM}  MAX_K={MAX_K}  "
          f"ALARM_DISTANCE={ALARM_DISTANCE*1e3:.2f} mm  "
          f"wall_gap={_GAP_WALL*1e3:.2f} mm")

    return root_node