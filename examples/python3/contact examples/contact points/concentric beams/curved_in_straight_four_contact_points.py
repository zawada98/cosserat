import math
import Sofa
import Sofa.Core

from init_monitoring import InitializationMonitor
from gui2 import CTRGuiBridgeStraightOuter
from live_monitor   import LiveContactMonitor
# =============================================================================
#  TUBE PHYSICAL PARAMETERS
# =============================================================================

T1_PARAMS = {
    'name':          'Outer_Tube',
    'tube_number':   1,
    'rex':           6e-3, # 6 mm
    'rin':           5e-3,  #5 mm
    'E':             6e10,
    'v':             0.33,
    'density':       6450,
    # Multiples of 3 put the S-shape contact extrema at actual frame/section
    # locations: s = 0, L/3, 2L/3, L.
    'nb_sections':   30,
    'nb_frames':     60,
    'color':         [0.75, 0.20, 0.75, 0.35],
}

T2_PARAMS = {
    'name':          'Arc',
    'tube_number':   2,
    'crv_radius':    30e-3,        # 3 cm
    'rex':           1.5e-3,        # outer radius [m]
    'rin':           8e-4,      # inner (lumen) radius [m]
    'E':             6e10,
    'v':             0.33,
    'density':       6450,
    # Multiples of 3 put the S-shape contact extrema at actual frame/section
    # locations: s = 0, L/3, 2L/3, L.
    'nb_sections':   30,
    'nb_frames':     60,
    'color':         [0.15, 0.50, 1.00, 1.0],
}

N_CIRCLE = 10   # points per cross-sectional ring for visual model
DEFAULT_NORMAL      = '0 -1 0'

# The relaxed S-shaped inner beam is tangent to the outer tube's inner wall at
# four extrema: endpoint, interior extremum, interior extremum, endpoint.
S_CONTACT_CLEARANCE = 0.0
S_SAMPLE_COUNT = 2000

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


_GAP_WALL      = T1_PARAMS['rin'] - T2_PARAMS['rex']
ALARM_DISTANCE = T2_PARAMS['rex']


# =============================================================================
#  GEOMETRY HELPERS
# =============================================================================
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

    TWO_PI = 2.0 * math.pi
    ring = [[0.0,
             r * math.cos(TWO_PI * k / n_sides),
             r * math.sin(TWO_PI * k / n_sides)] for k in range(n_sides)]
    return ring * n_frames



def build_hollow_tube_surface(rex, rin, n_frames_total, n_circle):
    if not (rex > rin > 0.0):
        raise ValueError(
            f"Hollow tube surface requires rex > rin > 0 "
            f"(got rex={rex}, rin={rin})."
        )

    N      = n_circle
    TWO_PI = 2.0 * math.pi

    # ---- Positions: per frame, outer ring then inner ring --------------------
    positions = []
    for _ in range(n_frames_total):
        for k in range(N):
            ang = TWO_PI * k / N
            positions.append([0.0,
                              rex * math.cos(ang),
                              rex * math.sin(ang)])
        for k in range(N):
            ang = TWO_PI * k / N
            positions.append([0.0,
                              rin * math.cos(ang),
                              rin * math.sin(ang)])

    # ---- Quads ---------------------------------------------------------------
    quads = []

    # Outer shell -- front face has radial-outward normal.
    # Frame i outer ring starts at index (i * 2N).
    for i in range(n_frames_total - 1):
        for j in range(N):
            a = i * 2 * N + j
            b = i * 2 * N + (j + 1) % N
            c = (i + 1) * 2 * N + (j + 1) % N
            d = (i + 1) * 2 * N + j
            quads.append([a, b, c, d])


    for i in range(n_frames_total - 1):
        for j in range(N):
            a = i * 2 * N + N + j
            b = i * 2 * N + N + (j + 1) % N
            c = (i + 1) * 2 * N + N + (j + 1) % N
            d = (i + 1) * 2 * N + N + j
            quads.append([a, d, c, b])

    rigid_idx = [f for f in range(n_frames_total) for _ in range(2 * N)]

    return positions, quads, rigid_idx


def add_tube_visual(frame_node, p, color=None, n_circle=N_CIRCLE):

    name = p['name']
    re   = p['rex']
    ri   = p['rin']
    nf   = p['nb_frames']
    n_frames_total = nf + 1

    # ---- Resolve color -------------------------------------------------------
    # A `None` color falls back to the tube's default color from the param
    # dict; a list/tuple is converted to OglModel's expected
    # whitespace-separated string format; a string is passed through verbatim.
    if color is None:
        color = p['color']

    if isinstance(color, str):
        color_str = color
    else:
        color_str = " ".join(str(v) for v in color)

    # ---- Hollow-tube surface geometry ----------------------------------------
    #                rigid_idx = [f for f in range(n_frames_total)
    #                             for _ in range(n_circle)]
    #                quads     = build_tube_quads(n_frames_total, n_circle)
    # The new helper builds outer + inner shells (open ends, no caps) so
    # the rendered geometry matches the physical tube wall: rex outside,
    # rin inside, with the lumen visible through the open ends.
    positions, quads, rigid_idx = build_hollow_tube_surface(
        re, ri, n_frames_total, n_circle
    )

    # ---- Visual node ---------------------------------------------------------
    visu_node = frame_node.addChild(name + '_visu')
    visu_node.addObject('MeshTopology', name='topo',
                        position=positions, quads=quads)
    visu_node.addObject('MechanicalObject', name='visMO',
                        template='Vec3d', position=positions)
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

    return visu_node


def calculate_x_offset(outer_tube, arc):
    L_outer = outer_tube['length']
    L_arc = arc['crv_radius'] * arc['angle']
    return (L_outer - L_arc)/2


# =============================================================================
#  GEOMETRY HELPERS
# =============================================================================

def compute_s_shape_in_straight_beam(outer_beam, beam,
                                     contact_clearance=0.0,
                                     sample_count=S_SAMPLE_COUNT):
    """
    Set the straight outer beam length and S-shaped inner-beam rest geometry.

    The S centerline is:
        y(x) = A * sin(3*pi*x/X)
    for x in [-X/2, X/2].

    This gives four wall contacts at:
        x = -X/2   y = +A   endpoint
        x = -X/6   y = -A   interior extremum
        x = +X/6   y = +A   interior extremum
        x = +X/2   y = -A   endpoint

    X is chosen so the peak curvature radius at those extrema equals
    beam['crv_radius'].
    """
    rin_outer = outer_beam['rin']
    rex_beam  = beam['rex']
    r_min     = beam['crv_radius']

    if contact_clearance < 0.0:
        raise ValueError("contact_clearance must be non-negative.")

    amplitude = rin_outer - rex_beam - contact_clearance
    if amplitude <= 0.0:
        raise ValueError(
            "S beam cannot fit inside the outer tube: "
            f"rin_outer={rin_outer}, rex_beam={rex_beam}, "
            f"contact_clearance={contact_clearance}."
        )

    if not (r_min > amplitude):
        raise ValueError(
            "S beam curvature radius must be larger than the contact amplitude "
            f"(got r_min={r_min}, amplitude={amplitude})."
        )

    wave_number = math.sqrt(1.0 / (amplitude * r_min))
    x_span = 3.0 * math.pi / wave_number

    samples = _sample_s_shape(amplitude, x_span, sample_count)
    total_length = samples[-1]['s']

    outer_beam['length'] = x_span
    outer_beam['y0'] = 0.0

    beam['length'] = total_length
    beam['s_shape'] = {
        'amplitude': amplitude,
        'x_span': x_span,
        'wave_number': wave_number,
        'samples': samples,
    }
    midpoint = _interp_s_shape(samples, 0.5 * total_length)
    qz_mid = math.sin(0.5 * midpoint['theta'])
    qw_mid = math.cos(0.5 * midpoint['theta'])
    beam['crest_pose'] = [
        midpoint['x'], midpoint['y'], 0.0,
        0.0, 0.0, qz_mid, qw_mid,
    ]

    print(
        "[Geometry] symmetric four-contact S placement: "
        f"outer L={outer_beam['length'] * 1e3:.3f} mm, "
        f"S length={beam['length'] * 1e3:.3f} mm, "
        f"amplitude={amplitude * 1e3:.3f} mm, "
        f"contact clearance={contact_clearance * 1e6:.1f} um"
    )


def _sample_s_shape(amplitude, x_span, sample_count):
    samples = []
    previous_x = previous_y = None
    s = 0.0

    for i in range(sample_count + 1):
        u = i / float(sample_count)
        x = -0.5 * x_span + u * x_span
        k = 3.0 * math.pi / x_span
        y = amplitude * math.sin(k * x)
        dy = amplitude * k * math.cos(k * x)
        ddy = -amplitude * k * k * math.sin(k * x)
        theta = math.atan(dy)
        curvature = ddy / ((1.0 + dy * dy) ** 1.5)

        if previous_x is not None:
            dx = x - previous_x
            dy_step = y - previous_y
            s += math.sqrt(dx * dx + dy_step * dy_step)

        samples.append({
            's': s,
            'x': x,
            'y': y,
            'theta': theta,
            'curvature': curvature,
        })
        previous_x, previous_y = x, y

    return samples


def _interp_s_shape(samples, s_target):
    if s_target <= 0.0:
        return samples[0]
    if s_target >= samples[-1]['s']:
        return samples[-1]

    lo = 0
    hi = len(samples) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if samples[mid]['s'] < s_target:
            lo = mid
        else:
            hi = mid

    a = samples[lo]
    b = samples[hi]
    span = b['s'] - a['s']
    t = 0.0 if span <= 1e-15 else (s_target - a['s']) / span

    return {
        's': s_target,
        'x': (1.0 - t) * a['x'] + t * b['x'],
        'y': (1.0 - t) * a['y'] + t * b['y'],
        'theta': (1.0 - t) * a['theta'] + t * b['theta'],
        'curvature': (1.0 - t) * a['curvature'] + t * b['curvature'],
    }


def tube_mass(p):                                                    
    """Hollow-cylinder mass = density * pi * (re^2 - ri^2) * length."""
    ri, re = p['rin'], p['rex']
    return p['density'] * math.pi * (re * re - ri * ri) * p['length']


# =============================================================================
#  COSSERAT TUBE GEOMETRY
# =============================================================================

def compute_tube_geometry(p, beam_type,
                          init_strategy='natural',
                          outer_params=None):
    """
    Build all geometry quantities consumed by add_cosserat_tube.

    Returns
    -------
    section_lengths : list[float]            len = nb_sections
    rest_states     : list[Vec3]             len = nb_sections
    init_states     : list[Vec3]             len = nb_sections
    sec_curv_abs    : list[float]            len = nb_sections + 1   (intrinsic, [0, L])
    frame_positions : list[Rigid3d]          len = nb_frames + 1     (WORLD coords)
    frm_curv_abs    : list[float]            len = nb_frames + 1     (intrinsic, [0, L])
    edge_indices    : list[[int,int]]        len = nb_frames

    Conventions
    -----------
    - Strain Vec3 is [k_torsion, k_y, k_z] in the BODY frame (BeamHookeLaw,
      Vec3d template). Pure planar bending in the local x-y plane uses k_z.
    - frame_positions[0] equals the rigid-base pose by construction; the
      caller writes that into the rigid_base MO.
    - 'straight': base at (-L/2, 0, 0), identity quat. Frames evenly spaced
      along world +X.
    - 'arc'    : base at the LEFT endpoint of an arc that lies on the
      circle x^2 + (y - y0)^2 = r_crv^2 (in the z=0 plane), symmetric
      about the world Y axis. Body +x at base points along the arc toward
      s = L. Requires outer_params['y0'].

    init_strategy
    -------------
      'natural'  : init_states == rest_states (system starts at rest)
      'straight' : init_states == 0 regardless of rest (system relaxes
                   from a straight initial shape toward rest_states)
    """
    if beam_type == 'straight':
        return _geometry_straight(p, init_strategy)
    elif beam_type == 'arc':
        return _geometry_arc(p, init_strategy, outer_params)
    elif beam_type == 's_shape':
        return _geometry_s_shape(p, init_strategy)
    else:
        raise ValueError(
            f"compute_tube_geometry: unknown beam_type={beam_type!r}; "
            f"expected 'straight', 'arc', or 's_shape'."
        )


def _geometry_straight(p, init_strategy):
    L      = p['length']
    nb_sec = p['nb_sections']
    nb_frm = p['nb_frames']      # number of frame intervals; nb_frm+1 frames

    # ---- Sections --------------------------------------------------------
    ls              = L / nb_sec
    section_lengths = [ls] * nb_sec
    sec_curv_abs    = [i * ls for i in range(nb_sec + 1)]
    sec_curv_abs[-1] = L                                  # snap the endpoint

    # ---- Strain (rest = init = 0) ---------------------------------------
    rest_states = [[0.0, 0.0, 0.0] for _ in range(nb_sec)]
    if init_strategy in ('natural', 'straight'):
        init_states = [[0.0, 0.0, 0.0] for _ in range(nb_sec)]
    else:
        raise ValueError(
            f"_geometry_straight: unsupported init_strategy={init_strategy!r}"
        )

    # ---- Frames ---------------------------------------------------------
    lf            = L / nb_frm
    frm_curv_abs  = [i * lf for i in range(nb_frm + 1)]
    frm_curv_abs[-1] = L

    x_base = -0.5 * L
    frame_positions = [
        [x_base + s, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        for s in frm_curv_abs
    ]

    edge_indices = [[i, i + 1] for i in range(nb_frm)]

    return (section_lengths, rest_states, init_states,
            sec_curv_abs, frame_positions, frm_curv_abs,
            edge_indices)


def _geometry_arc(p, init_strategy, outer_params):                        
    L      = p['length']
    r_crv  = p['crv_radius']
    nb_sec = p['nb_sections']
    nb_frm = p['nb_frames']

    if init_strategy not in ('natural', 'straight'):
        raise ValueError(
            f"_geometry_arc: unsupported init_strategy={init_strategy!r}"
        )

    # ---- Sections (intrinsic) -------------------------------------------
    ls               = L / nb_sec
    section_lengths  = [ls] * nb_sec
    sec_curv_abs     = [i * ls for i in range(nb_sec + 1)]
    sec_curv_abs[-1] = L

    # ---- Rest strain (intrinsic; same for both init strategies) --------
    kappa_z     = -1.0 / r_crv
    rest_states = [[0.0, 0.0, kappa_z] for _ in range(nb_sec)]

    # ---- Frame curvilinear abscissas (intrinsic) ------------------------
    lf               = L / nb_frm
    frm_curv_abs     = [i * lf for i in range(nb_frm + 1)]
    frm_curv_abs[-1] = L

    if init_strategy == 'natural':
        # Base at LEFT endpoint of the symmetric arc on the (0,y0)-centred
        # circle.  Rod starts in its natural curved shape (init = rest).
        if outer_params is None or 'y0' not in outer_params:              
            raise ValueError(
                "_geometry_arc: 'natural' requires outer_params['y0'] "
                "(set by compute_symmetric_arc_in_straight_beam)."
            )
        y0       = outer_params['y0']
        phi_half = 0.5 * (L / r_crv)

        x_base    = -r_crv * math.sin(phi_half)
        y_base    =  y0 + r_crv * math.cos(phi_half)
        beta      =  phi_half
        qz_b      = math.sin(0.5 * beta)
        qw_b      = math.cos(0.5 * beta)
        base_pose = [x_base, y_base, 0.0, 0.0, 0.0, qz_b, qw_b]

        init_states = [list(v) for v in rest_states]

        cb, sb          = math.cos(beta), math.sin(beta)
        frame_positions = []
        for s in frm_curv_abs:
            phi     = s / r_crv
            x_body  =  r_crv * math.sin(phi)
            y_body  = -r_crv * (1.0 - math.cos(phi))
            wx      = x_base + cb * x_body - sb * y_body
            wy      = y_base + sb * x_body + cb * y_body
            theta_w = beta - phi
            fqz     = math.sin(0.5 * theta_w)
            fqw     = math.cos(0.5 * theta_w)
            frame_positions.append([wx, wy, 0.0, 0.0, 0.0, fqz, fqw])

    else:
        # Start straight, centered on the straight outer beam centerline.
        # The centerlines are concentric at initialization: y=z=0 and the
        # inner straightened arc is centered about x=0.
        x_base = -0.5 * L
        y_base = 0.0
        base_pose = [x_base, y_base, 0.0, 0.0, 0.0, 0.0, 1.0]
        init_states = [[0.0, 0.0, 0.0] for _ in range(nb_sec)]

        frame_positions = [
            [x_base + s, y_base, 0.0, 0.0, 0.0, 0.0, 1.0]
            for s in frm_curv_abs
        ]

    assert frame_positions[0] == base_pose, (
        "_geometry_arc: frame_positions[0] != base_pose -- code bug."
    )

    edge_indices = [[i, i + 1] for i in range(nb_frm)]

    return (section_lengths, rest_states, init_states,
            sec_curv_abs, frame_positions, frm_curv_abs,
            edge_indices)


def _geometry_s_shape(p, init_strategy):
    L      = p['length']
    nb_sec = p['nb_sections']
    nb_frm = p['nb_frames']

    if init_strategy not in ('natural', 'straight'):
        raise ValueError(
            f"_geometry_s_shape: unsupported init_strategy={init_strategy!r}"
        )
    if 's_shape' not in p:
        raise ValueError(
            "_geometry_s_shape requires p['s_shape']; call "
            "compute_s_shape_in_straight_beam first."
        )

    samples = p['s_shape']['samples']

    ls = L / nb_sec
    section_lengths = [ls] * nb_sec
    sec_curv_abs = [i * ls for i in range(nb_sec + 1)]
    sec_curv_abs[-1] = L

    rest_states = []
    for i in range(nb_sec):
        s_mid = (i + 0.5) * ls
        sample = _interp_s_shape(samples, s_mid)
        rest_states.append([0.0, 0.0, sample['curvature']])

    if init_strategy == 'natural':
        init_states = [list(v) for v in rest_states]
    else:
        init_states = [[0.0, 0.0, 0.0] for _ in range(nb_sec)]

    lf = L / nb_frm
    frm_curv_abs = [i * lf for i in range(nb_frm + 1)]
    frm_curv_abs[-1] = L

    if init_strategy == 'natural':
        frame_positions = []
        for s in frm_curv_abs:
            sample = _interp_s_shape(samples, s)
            qz = math.sin(0.5 * sample['theta'])
            qw = math.cos(0.5 * sample['theta'])
            frame_positions.append([
                sample['x'], sample['y'], 0.0,
                0.0, 0.0, qz, qw,
            ])
    else:
        x_base = -0.5 * L
        frame_positions = [
            [x_base + s, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
            for s in frm_curv_abs
        ]

    edge_indices = [[i, i + 1] for i in range(nb_frm)]

    return (section_lengths, rest_states, init_states,
            sec_curv_abs, frame_positions, frm_curv_abs,
            edge_indices)

# =============================================================================
#  COSSERAT TUBE: SCENE-GRAPH BUILDER
# =============================================================================

def add_cosserat_tube(root_node, p, beam_type,
                      init_strategy='natural',
                      outer_params=None,
                      fixed_directions=None,
                      fix_strain=False):
    """
    Build a Cosserat tube under root_node.

    Scene graph:
        root_node/<name>                              (tube_node)
          SolverNode/                                 (solver_node)
            EulerImplicitSolver, SparseLDLSolver,
            GenericConstraintCorrection
            <name>_rigid_base/                        (rigid_base)
              MechanicalObject<Rigid3d>  (1 dof)      <- base_pose from geometry
              PartialFixedProjectiveConstraint        (free X-trans, X-rot)
            <name>_coss_state/                        (coss_state)
              MechanicalObject<Vec3d>    (nb_sec dofs)
              BeamHookeLawForceField                  (linear elastic, circular hollow)
            <name>_frames/                            (frame_node)
              MechanicalObject<Rigid3d>  (nb_frm+1 dofs)
              UniformMass
              DiscreteCosseratMapping  (input1=coss, input2=base, output=frames)

    Parameters
    ----------
    beam_type     : 'straight' or 'arc'
    init_strategy : 'natural' (init = rest) or 'straight' (init = zero)
    outer_params  : required for beam_type='arc' (provides 'y0').
    fixed_directions : list[int] of length 6, or None
                      If None (default), no PartialFixedProjectiveConstraint is added.
                      Otherwise, a list of 6 ints (1=fixed, 0=free) for [tx,ty,tz,rx,ry,rz].

    Returns
    -------
    base_mo, coss_mo, tube_node, frame_node, odesolver
    """
    name = p['name']

    (section_lengths, rest_states, init_states,
     sec_curv_abs, frame_positions, frm_curv_abs,
     _edge_indices) = compute_tube_geometry(
        p, beam_type,
        init_strategy=init_strategy,
        outer_params=outer_params,
    )

    # frame_positions[0] is the base pose in WORLD coords (set by the
    # geometry function for the chosen beam_type and init_strategy).
    base_pose = frame_positions[0]

    re, ri = p['rex'], p['rin']
    mass   = tube_mass(p)

    # ---- Solver scope ----------------------------------------------------
    tube_node   = root_node.addChild(name)
    solver_node = tube_node.addChild('SolverNode')

    odesolver = solver_node.addObject(
        'EulerImplicitSolver',
        name='odesolver',
        rayleighStiffness=0.2,
        rayleighMass=0.1,
        firstOrder=False,
    )
    solver_node.addObject(
        'SparseLDLSolver',
        name='Solver',
        template='CompressedRowSparseMatrixd',
    )
    solver_node.addObject('GenericConstraintCorrection')

    # ---- Rigid base ------------------------------------------------------
    rigid_base = solver_node.addChild(name + '_rigid_base')
    base_mo = rigid_base.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='cosserat_base_mo',
        position=[base_pose],
        showObject=True,
        showObjectScale=0.001,
    )
    # Free X-translation (insertion) and X-rotation (axial twist) only.
    if fixed_directions is not None:
        rigid_base.addObject(
            'PartialFixedProjectiveConstraint',
            name='proximal_bc',
            fixedDirections=list(fixed_directions),
            indices=[0],
        )

    # ---- Cosserat strain state ------------------------------------------
    coss_state = solver_node.addChild(name + '_coss_state')
    coss_mo = coss_state.addObject(
        'MechanicalObject',
        template='Vec3d',
        name='cosserat_state',
        position=init_states,
        rest_position=rest_states,
    )
    if fix_strain:
        coss_state.addObject(
            'FixedProjectiveConstraint',
            name='fixed_strain',
            indices=list(range(len(init_states))),
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

    # ---- Output frames (single parent: SolverNode) ----------------------
    # DiscreteCosseratMapping reaches both inputs by relative sibling paths;
    # do NOT addChild this node from rigid_base or coss_state.
    frame_node = solver_node.addChild(name + '_frames')
    frame_node.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='FramesMO',
        position=frame_positions,
        showObject=True,
        showObjectScale=0.001,
    )
    frame_node.addObject('UniformMass', name='mass', totalMass=mass)
    frame_node.addObject(
        'DiscreteCosseratMapping',
        name='cosseratMapping',
        curv_abs_input=sec_curv_abs,
        curv_abs_output=frm_curv_abs,
        input1='@../' + name + '_coss_state/cosserat_state',
        input2='@../' + name + '_rigid_base/cosserat_base_mo',
        output='@FramesMO',
        debug=False,
        radius=re,
    )

    return base_mo, coss_mo, tube_node, frame_node, odesolver


# =============================================================================
#  CONTROLLER
# =============================================================================
class CTRController(Sofa.Core.Controller):
    """
    GUI-driven midpoint control-point actuation for the inner S-shaped beam.

    The inner tube starts straight and relaxes toward its curved rest shape
    during initialization. During that phase the controller does not write the
    inner base, so the Cosserat solve can find the natural symmetric placement.
    When the initialization monitor switches to control mode, the controller
    captures that settled midpoint-control pose and then applies the five GUI DOFs: world
    tx/ty/tz, world-Z yaw, and local-X twist.
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
        self.root_node  = root_node
        self.t1_base_mo = t1_base_mo
        self.t2_base_mo = t2_base_mo
        self.t2_control_mo = t2_control_mo
        self.gui        = gui_bridge

        # Outer: pin to construction-time pose forever.
        self._t1_rest_pose = list(t1_base_mo.position.value[0])

        self._t2_initial_pose = list(t2_base_mo.position.value[0])
        self._t2_control_initial_pose = list(t2_control_mo.position.value[0])
        self._t2_control_rest_pose = None

        self._t2_tx = 0.0
        self._t2_ty = 0.0
        self._t2_tz = 0.0
        self._t2_rx = 0.0
        self._t2_rz = 0.0
        try:
            dt0 = float(self.root_node.dt.value)
        except Exception:
            dt0 = 1e-4
        self._dt_current = dt0
        self._dt_target  = dt0

        self._step       = 0
        self._prev_phase = 'waiting'

    # ------------------------------------------------------------------
    def onAnimateBeginEvent(self, event):
        self._step += 1
        snap  = self.gui.snapshot()
        phase = snap['phase']

        # 1) Outer tube: pin every step regardless of phase.
        self._write_pose(self.t1_base_mo, self._t1_rest_pose)

        # 2) Before initialization starts, hold the straight concentric pose.
        # During initialization, do not write the inner base: the straight arc
        # must relax under its own rest curvature and contact constraints.
        if phase == 'waiting':
            self._write_pose(self.t2_base_mo, self._t2_initial_pose)
            self._write_pose(self.t2_control_mo, self._t2_control_initial_pose)
            self._prev_phase = phase
            return

        if phase == 'initializing':
            self._write_pose(self.t2_control_mo, self._t2_control_initial_pose)
            self._prev_phase = phase
            return

        if phase != 'control':
            self._prev_phase = phase
            return

        if self._prev_phase != 'control':
            self._t2_control_rest_pose = list(self.t2_control_mo.position.value[0])
            x, y, z, qx, qy, qz, qw = self._t2_control_rest_pose
            print(f"[CTRController] init -> control at step {self._step}; "
                  f"captured inner-tube midpoint control pose: "
                  f"pos=({x:+.4f}, {y:+.4f}, {z:+.4f}) "
                  f"quat=({qx:+.4f}, {qy:+.4f}, {qz:+.4f}, {qw:+.4f})")
        self._prev_phase = phase

        # 3) Control phase: dt ramp + DOF advancement + write crest control point.
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
        self._t2_ty = self._step_toward(self._t2_ty,
                                        float(snap['t2_ty_target_m']),  max_t)
        self._t2_tz = self._step_toward(self._t2_tz,
                                        float(snap['t2_tz_target_m']),  max_t)
        self._t2_rx = self._step_toward(self._t2_rx,
                                        float(snap['t2_rx_target_rad']), max_r)
        self._t2_rz = self._step_toward(self._t2_rz,
                                        float(snap['t2_rz_target_rad']), max_r)

        self._write_pose(self.t2_control_mo, self._compose_inner_pose())

    def _compose_inner_pose(self):
        pose = list(self._t2_control_rest_pose)
        pose[0] += self._t2_tx
        pose[1] += self._t2_ty
        pose[2] += self._t2_tz

        q_rest = tuple(self._t2_control_rest_pose[3:7])
        q_new = self._quat_mul(
            self._quat_z(self._t2_rz),
            self._quat_mul(q_rest, self._quat_x(self._t2_rx)),
        )
        pose[3:7] = list(self._quat_normalize(q_new))
        return pose

    def _disable_init_anchor(self):
        if self.init_anchor is None:
            return
        for data_name in ('stiffness', 'angularStiffness'):
            try:
                setattr(self.init_anchor, data_name, 0.0)
            except Exception:
                try:
                    self.init_anchor.findData(data_name).value = 0.0
                except Exception:
                    pass
        self.init_anchor = None

        # ------------------------------------------------------------------
    def _unused_write_kinematic_dofs(self, mo, new_x): 
        raise RuntimeError("Deprecated path: use _compose_inner_pose().")
        """
        Write only kinematic DOFs:
          - tx  : insertion axis (written explicitly)
          - Rx  : axial twist (composed onto CURRENT physics quaternion)
        ty, tz, Rz are NOT written — physics and projective constraint own them.
        """
        with mo.position.writeable() as p:
            current = list(p[0])
            current[0] = new_x  # tx: insertion

            # ---- Rx-only quaternion update --------------------------------
            # Read the CURRENT physics quaternion (includes physics-determined Rz).
            q_current = tuple(current[3:7])

            # Remove the Rx we applied last step to isolate the physics Rz:
            q_rx_prev_inv = self._quat_conj(  
                self._quat_x(self._t2_rx_written))  
            q_rz_physics = self._quat_normalize(  
                self._quat_mul(q_current, q_rx_prev_inv))  

            # Reapply the newly commanded cumulative Rx on top of physics Rz:
            q_rx_new = self._quat_x(self._t2_rx)  
            q_new = self._quat_normalize(  
                self._quat_mul(q_rz_physics, q_rx_new))  

            current[3:7] = list(q_new)  
            p[0] = current

        self._t2_rx_written = self._t2_rx

    @staticmethod
    def _quat_conj(q):
        return (-q[0], -q[1], -q[2], q[3])

    @staticmethod
    def _write_pose(mo, pose):
        with mo.position.writeable() as p:
            p[0] = list(pose)

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
    def _quat_slerp(a, b, alpha):
        ax, ay, az, aw = CTRController._quat_normalize(a)
        bx, by, bz, bw = CTRController._quat_normalize(b)
        dot = ax*bx + ay*by + az*bz + aw*bw
        if dot < 0.0:
            bx, by, bz, bw = -bx, -by, -bz, -bw
            dot = -dot
        if dot > 0.9995:
            return CTRController._quat_normalize((
                ax + alpha * (bx - ax),
                ay + alpha * (by - ay),
                az + alpha * (bz - az),
                aw + alpha * (bw - aw),
            ))
        theta_0 = math.acos(max(-1.0, min(1.0, dot)))
        theta = theta_0 * alpha
        sin_theta = math.sin(theta)
        sin_theta_0 = math.sin(theta_0)
        s0 = math.cos(theta) - dot * sin_theta / sin_theta_0
        s1 = sin_theta / sin_theta_0
        return (
            s0 * ax + s1 * bx,
            s0 * ay + s1 * by,
            s0 * az + s1 * bz,
            s0 * aw + s1 * bw,
        )

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
        'Sofa.Component.SolidMechanics.Spring',              # RestShapeSpringsForceField
        'Sofa.Component.StateContainer',                     # MechanicalObject
        'Sofa.Component.Topology.Container.Constant',        # MeshTopology
        'Sofa.Component.Visual',                             # VisualStyle
        'Sofa.GL.Component.Rendering3D',                     # OglModel
        'Sofa.GUI.Component',                                # Camera
    ])

    root_node.gravity = [0., 0., 0.]
    root_node.dt      = 1e-3

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

    compute_s_shape_in_straight_beam(
        T1_PARAMS,
        T2_PARAMS,
        contact_clearance=S_CONTACT_CLEARANCE,
    )

    # ---- Build tubes --------------------------------------------------------
    t1_base_mo, _, _, t1_frame_node, t1_solver = add_cosserat_tube(
        root_node, T1_PARAMS, beam_type='straight',
        init_strategy='natural',
        fixed_directions=[1, 1, 1, 1, 1, 1],
        fix_strain=True,
    )
    t2_base_mo, _, _, t2_frame_node, t2_solver = add_cosserat_tube(
        root_node, T2_PARAMS, beam_type='s_shape',
        init_strategy='straight',
        outer_params=T1_PARAMS,
        fixed_directions=None,
    )

    # GUI control point attached midway along the inner S. During
    # initialization it also removes the free global translation mode while the
    # tube bends from straight to its natural S shape.
    t2_mid_frame = T2_PARAMS['nb_frames'] // 2
    control_node = root_node.addChild('S_Mid_ControlPoint')
    t2_control_mo = control_node.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='controlPointMO',
        position=[T2_PARAMS['crest_pose']],
        showObject=True,
        showObjectScale=0.005,
    )
    t2_frame_node.addObject(
        'RestShapeSpringsForceField',
        name='midpoint_control_spring',
        points=[t2_mid_frame],
        external_rest_shape=t2_control_mo.getLinkPath(),
        external_points=[0],
        mstate='@FramesMO',
        stiffness=1e7,
        angularStiffness=1e7,
        template='Rigid3d',
        activeDirections=[1, 1, 1, 1, 1, 1, 1],
    )
    # ---- Visual models ------------------------------------------------------

    add_tube_visual(t1_frame_node, T1_PARAMS, color=T1_PARAMS['color'])  # outer
    add_tube_visual(t2_frame_node, T2_PARAMS, color=T2_PARAMS['color'])  # inner

    # ---- GUI bridge ---------------------------------------------------------
    gui_bridge = CTRGuiBridgeStraightOuter(
        root_node=root_node,
        max_tx_m=0.04,  # tune for your scene
        max_ty_m=0.01,
        max_tz_m=0.01,
        max_rx_deg=180.0,
        max_rz_deg=15.0,
        init_dt=float(root_node.dt.value),
        default_control_dt=1e-3,
        dt_min=1e-6, dt_max=1e-1,
    )

    root_node.addObject(CTRController(
        name='CTRController',
        root_node=root_node,
        t1_base_mo=t1_base_mo,
        t2_base_mo=t2_base_mo,
        t2_control_mo=t2_control_mo,
        gui_bridge=gui_bridge,
    ))

    intersection_node = root_node.addChild('IntersectionNode')

    # ---- FramesMO handles ---------------------------------------------------
    t1_MO = t1_frame_node.FramesMO
    t2_MO = t2_frame_node.FramesMO


    ssim = intersection_node.addObject(
        'SphereSweptIntersectionMethod',
        name='ssim',
        beam1Frames=t1_MO.getLinkPath() + '.position',
        beam2Frames=t2_MO.getLinkPath() + '.position',
        beam1Velocities=t1_MO.getLinkPath() + '.velocity',
        beam2Velocities=t2_MO.getLinkPath() + '.velocity',
        radius1=T1_PARAMS['rex'],  # inner wall of Tube_1 [m]
        radius2=T2_PARAMS['rex'],
        innerRadius1=T1_PARAMS['rin'],
        innerRadius2=T2_PARAMS['rin'],
        contactConfiguration = "nested",
        defaultNormal = DEFAULT_NORMAL,
        broadPhaseMarginFactor = 2
    )

    contact_output = t1_frame_node.addChild('contactOutput')
    t2_frame_node.addChild(contact_output)

    contactMO = contact_output.addObject(
        'MechanicalObject', template='Vec3d',
        name='contactMO_gap',
        position=[[0., 0., 0.]] * 2 * MAX_K)

    bcm = contact_output.addObject(
        'BeamContactMapping',
        name='bcm',
        input1=t1_MO.getLinkPath(),
        input2=t2_MO.getLinkPath(),
        output=contactMO.getLinkPath(),
        ssim=ssim.getLinkPath(),
        mappingMode='contactPoints'
    )

    contact_output.addObject(
        'ContactPointsUnilateralConstraint',
        name='cpuc',
        mu=0,
        contactTriads=bcm.getLinkPath() + '.contactTriads',
        gapSign=bcm.getLinkPath() + '.gapSign'
    )

    t2_curv_abs_frames = list(
        t2_frame_node.cosseratMapping.curv_abs_output.value
    )

    intersection_node.addObject(InitializationMonitor(
        name='InitMonitor',
        t1_MO=t1_MO,
        t2_MO=t2_MO,
        bcm=bcm,
        contact_mo=contactMO,
        t2_frame_curv_abs=t2_curv_abs_frames,
        vel_threshold=10e-3,
        quiet_steps=100,
        warmup_steps=500,
        require_contact=True,
        min_contact_pairs=4,
        contact_gap_threshold=1e-5,
        log_every=200,
        auto_open=True,
        png_path=None,
        on_init_complete=gui_bridge.signal_init_complete,
    ))

    intersection_node.addObject(LiveContactMonitor(
        name='LiveMonitor',
        t2_MO=t2_MO,
        bcm=bcm,
        contact_mo=contactMO,
        t2_frame_curv_abs=t2_curv_abs_frames,
        bridge=gui_bridge,
        every_n_steps=20,
        contact_constraint=cpuc,
        force_unit_scale=1.0,  # SI scene units: kg*m/s^2 = N
    ))
    # -------------------------------------------------------------------------

    root_node.animate = False

    return root_node
