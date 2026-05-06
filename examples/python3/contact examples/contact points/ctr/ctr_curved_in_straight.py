# -*- coding: utf-8 -*-
"""
ctr_curved_in_straight.py
=========================
SOFA / Cosserat Plugin -- nested CTR with a STRAIGHT outer beam and a
PRE-CURVED inner arc.  Heavily based on simulate_ctr_two_tubes.py; only
the geometry, the initial-state strategy, the rigid-base BCs, the origin
shift, and the 5-DOF actuation block differ.

Geometry summary (all SI units)
-------------------------------
  Outer (Tube_outer) : STRAIGHT, L = 0.17 m, OD = 3.0 mm,  ID = 2.7 mm.
                       Stationary -- all 6 base DOFs pinned.

  Inner (Tube_inner) : PRE-CURVED ARC about z-axis (curvature in xy plane).
                       L = 0.22 m, OD = 0.8 mm, ID = 0.54 mm.
                       Curve direction is NEGATIVE about +z, so the natural
                       shape opens upward (C is the lowest point of the U).
                       Initial configuration: STRAIGHTENED, concentric.

Origin convention  (DIFFERS from ctr_two_tubes.py)
-------------------------------------------------
  The world origin (0,0,0) is at the GEOMETRIC CENTER of the outer beam,
  midway along its longitudinal axis.

    Outer rigid base @ x = -L_outer/2 = -0.085 m
    Outer tip        @ x = +L_outer/2 = +0.085 m

    Inner rigid base @ x = +L_outer/2 - L_inner = -0.135 m
    Inner tip (when straightened) @ x = +L_outer/2 = +0.085 m

  At t=0 the inner tube is STRAIGHTENED concentric with the outer beam.
  When the user presses Initialize the inner tube relaxes from straight to
  its natural curved shape, accumulating contact with the outer beam's
  inner wall.  C (the lowest y point of the natural U-curve) is at
  approximately x = midpoint-of-curved-portion, y = -sagitta.

5-DOF actuation at the inner tube's proximal extremity
------------------------------------------------------
  Once equilibrium is detected (InitializationMonitor fires), the GUI
  exposes 5 absolute-target sliders driving the inner tube's rigid base:

    tx [m]               translation along world +X
    ty [m]               translation along world +Y
    tz [m]               translation along world +Z
    rx [deg]             rotation about local +X (TWIST -- beam axis)
    rz [deg]             rotation about world +Z (YAW   -- bend in xy)

  ROTATION ABOUT y IS NOT EXPOSED.  The PartialFixedProjective constraint
  fixedDirections=[0,0,0,0,1,0] enforces zero rotation about local y.

  Compose the rigid base orientation as q = q_z(rz) * q_x(rx)  (intrinsic
  yaw-then-twist in active form).  Components in SOFA's (qx,qy,qz,qw)
  layout (alpha = rx, beta = rz):

      qx =  cos(beta/2) * sin(alpha/2)
      qy =  sin(alpha/2) * sin(beta/2)
      qz =  cos(alpha/2) * sin(beta/2)
      qw =  cos(alpha/2) * cos(beta/2)

  All five DOFs use ABSOLUTE TARGETS (not rates).  The controller rate-
  limits each axis independently per step:
    - translation_step_m (shared) caps |tx,ty,tz| step size (m/step)
    - rotation_step_rad  (shared) caps |rx,rz|    step size (rad/step)
  Per-step (NOT per-sec-sim) is the right unit because real-time
  perception of motion is roughly dt-independent -- see the long comment
  on shared['translation_step_m'] in gui_5dof.py for the rationale.

GEOMETRIC INFEASIBILITY OF THE "2-CONTACT-AT-RIMS, NONE AT C" REGIME
-------------------------------------------------------------------
The constraints originally proposed for this scene
    C1:  r_c + r_ex_c < y_O + r_in       (C does not touch lower wall)
    C2:  (x_O - L/2)^2 + y_O^2 < r_c^2   (arc circle covers the rim)
were intended to engineer a configuration with exactly two contact
points at the outer beam's extremities and no contact at C.  With
x_O = 0 and a symmetric U-arc inside a straight cylinder of constant
clearance c = r_in - r_ex_c, this regime is GEOMETRICALLY UNREACHABLE:
the natural arc dips MOST at C (its lowest y point), so if the dip at
C is less than c there is no contact anywhere, and if it exceeds c
then C is the FIRST point to touch -- not the rims.

Substituting C1 (corrected, wall-aware) into C2 yields (L/2)^2 <= 0,
strictly infeasible for L > 0.  See verify_geometry() below for the
formal check that runs at scene-build time and prints a regime
classification.

What this scene actually delivers
---------------------------------
With the parameters chosen below (r_c = 2.0 m, theta_arc = 5 deg ->
sagitta ~ 1.9 mm, wall clearance = 0.95 mm), the inner tube ends up in
the "mild distributed contact" regime: contact along the lower wall in
a region around C, distributed continuously rather than as two clean
rim points.  The geometry can be retuned for milder or no contact by
increasing crv_radius or decreasing crv_angle_deg.

NODE TOPOLOGY -- identical to ctr_two_tubes.py (single-parent frame_node).
CUSTOM CONTACT PIPELINE -- identical (SSIM + BCM + CPULC).
"""

import math
import Sofa
import Sofa.Core

from init_monitoring import InitializationMonitor
from gui2       import CTRGuiBridgeStraightOuter   # modified  was: from gui import CTRGuiBridge
from live_monitor   import LiveContactMonitor
# =============================================================================
#  TUBE PHYSICAL PARAMETERS
# =============================================================================


T1_PARAMS = {
    'name':          'Tube_outer',
    'tube_number':   1,
    'str_length':    0.17,                   # total arc length [m]
    'crv_radius':    1.0,                    # modified  was: 0.12  -- placeholder, unused (theta=0)
    'crv_angle_deg': 0.0,                    # modified  was: 60.0  -- STRAIGHT
    'rex':           110e-4,                  # outer (external) radius [m]
    'rin':           100e-4,                # inner (lumen) radius   [m]
    'E':             6e10,
    'v':             0.33,
    'density':       6450,
    'nb_sections':   50,
    'nb_frames':     100,
    'color':         [0.75, 0.20, 0.75, 0.35],
}


T2_PARAMS = {
    'name':          'Tube_inner',
    'tube_number':   2,
    'str_length':    0.040 + 0.2 * math.radians(60.0), #0.2 * math.radians(60.0)
    'crv_radius':    0.2,
    'crv_angle_deg': 60.0,
    'rex':           4e-4,
    'rin':           2.7e-4,
    'E':             6e10,
    'v':             0.33,
    'density':       6450,
    'nb_sections':   50,
    'nb_frames':     100,
    'color':         [0.15, 0.50, 1.00, 1.0],
    'crv_sign':      -1,
    'curve_placement': 'centered',
}

N_CIRCLE       = 10
DEFAULT_NORMAL = '0 1 0'

# =============================================================================
#  CONTACT PIPELINE PARAMETERS  (same as ctr_two_tubes.py)
# =============================================================================
ALGORITHM = "ALGO_1"

MAX_K = max(
    T1_PARAMS['nb_sections'] * T2_PARAMS['nb_sections'],
    (T1_PARAMS['nb_frames'] + 1) * T2_PARAMS['nb_sections'],
)

_GAP_WALL      = T1_PARAMS['rin'] - T2_PARAMS['rex']
ALARM_DISTANCE = T2_PARAMS['rex']


def verify_geometry(p_outer, p_inner, x_O=0.0):
    """
    Evaluate the two geometric constraints originally proposed for the
    "2 contact points at extremities, none at C" configuration, and
    classify the regime the current parameters fall into.

    The natural arc circle has center (x_O, y_O, 0) and radius r_c.  The
    arc's lowest point C sits at (x_O, y_O - r_c, 0).  For the natural
    arc to dip below the outer beam centerline at C we need y_O < r_c
    (so y_C < 0).  The user's original constraints, in the notation of
    this scene (r_ex -> r_ex_inner = T2_PARAMS['rex']):

        C1 (user):  r_c + r_ex_inner < y_O + r_in_outer
                    -> y_C + r_ex_inner > -r_in_outer
                    -> the inner tube body at C clears the lower wall

        C2 (user):  (x_O - L/2)^2 + y_O^2 < r_c^2
                    -> the arc circle covers the rim point (L/2, 0)
                    -> arc passes below y=0 at x=+L/2

    The CORRECTED version of C2 -- which actually expresses "the natural
    arc penetrates the lower wall at the rim", i.e. y_arc(L/2) < -c with
    c = r_in_outer - r_ex_inner the wall clearance -- is

        C2*:  (L/2)^2 + (y_O + r_in_outer - r_ex_inner)^2 <= r_c^2

    Substituting C1 into C2* yields (L/2)^2 <= 0, which is infeasible
    for any L > 0.  This is the geometric reason the
    "2-contact-at-rims, none at C" regime is unreachable for a
    symmetric U-arc inside a straight cylinder of constant clearance:
    a symmetric arc dips MOST at C, so contact starts there before it
    can reach the rims.

    Regime classification (all derived from the natural arc, before
    any constraint forces are applied):
        sagitta = r_c - sqrt(r_c^2 - (L/2)^2)        [arc dip at x=+/-L/2]
        s_C     = r_c - y_O                          [arc dip at C]
        c       = r_in_outer - r_ex_inner            [clearance]

      no-contact      : s_C < c                  (everything above lower wall)
      contact-at-C    : s_C >= c, no contact at extremities
      mild-distributed: contact in a centered region around C
      infeasible-rims : the user's intended regime -- not reachable here
    """
    L     = p_outer['str_length']
    r_in  = p_outer['rin']
    r_ex  = p_inner['rex']
    r_c   = p_inner['crv_radius']
    theta = math.radians(p_inner['crv_angle_deg'])
    c     = r_in - r_ex                                  # wall clearance [m]

    # The simulation places the natural arc with its midpoint at depth
    # = sagitta below the centerline (the curved portion's start tangent
    # is along +X, so the arc dips symmetrically about its arc midpoint).
    # For the *mathematical* analysis below we adopt the user's framing:
    # arc center at (x_O, y_O), radius r_c, with y_O treated as a free
    # parameter chosen to satisfy C1.
    sagitta_full = r_c * (1.0 - math.cos(theta / 2.0))   # arc dip over full chord


    y_O = r_c - 0.999 * c

    # ---- C1 (user's original) -----------------------------------------------
    # r_c + r_ex < y_O + r_in   <=>   y_O > r_c + r_ex - r_in = r_c - c
    c1_lhs = r_c + r_ex
    c1_rhs = y_O + r_in
    c1_ok  = c1_lhs < c1_rhs

    # ---- C2 (user's original) -----------------------------------------------
    c2u_lhs = (x_O - 0.5 * L) ** 2 + y_O ** 2
    c2u_rhs = r_c ** 2
    c2u_ok  = c2u_lhs < c2u_rhs

    # ---- C2* (corrected, wall-aware) ----------------------------------------
    # arc passes below the lower wall at x = +/- L/2:
    #   y_arc(L/2) = y_O - sqrt(r_c^2 - (L/2)^2) < -c
    #   <=>  sqrt(r_c^2 - (L/2)^2) > y_O + c
    #   <=>  r_c^2 > (L/2)^2 + (y_O + c)^2  (assuming y_O + c > 0)
    c2c_lhs = (0.5 * L) ** 2 + (y_O + c) ** 2
    c2c_rhs = r_c ** 2
    c2c_ok  = c2c_lhs <= c2c_rhs

    # ---- Sagitta vs clearance regime ----------------------------------------
    # The arc dip at C (relative to the centerline) is r_c - y_O = 0.999*c.
    # The arc dip at x = +/- L/2 is r_c - sqrt(r_c^2 - (L/2)^2).  This is
    # only meaningful if the chord covers x = +/- L/2, i.e., the curved
    # portion is long enough -- here we check it geometrically.
    L_crv_chord_half = r_c * math.sin(theta / 2.0)
    chord_covers_rim = L_crv_chord_half >= 0.5 * L
    if chord_covers_rim:
        sag_at_rim = r_c - math.sqrt(max(0.0, r_c ** 2 - (0.5 * L) ** 2))
    else:
        sag_at_rim = float('inf')        # arc doesn't reach the rim

    # Regime
    if sagitta_full < c:
        regime = "no-contact (natural sagitta < clearance)"
    elif chord_covers_rim and sag_at_rim < c:
        regime = "contact-near-C-only (arc rises back above wall at rim)"
    elif chord_covers_rim and sag_at_rim >= c:
        regime = "mild-distributed (contact extends to the rim)"
    else:
        regime = "contact-at-C-only (curved chord does not reach rim)"

    print("=" * 76)
    print("[CTR-geometry] Constraint check for 2-contact regime")
    print("=" * 76)
    print(f"  L       = {L * 1e3:.3f} mm   (outer beam length)")
    print(f"  r_in    = {r_in * 1e3:.3f} mm   (outer beam inner radius)")
    print(f"  r_ex    = {r_ex * 1e3:.3f} mm   (inner beam outer radius)")
    print(f"  c       = {c * 1e3:.3f} mm   (wall clearance r_in - r_ex)")
    print(f"  r_c     = {r_c * 1e3:.3f} mm   (curvature radius of inner arc)")
    print(f"  theta   = {math.degrees(theta):.3f} deg "
          f"({theta:.5f} rad, arc length = {r_c * theta * 1e3:.3f} mm)")
    print(f"  x_O,y_O = ({x_O * 1e3:.3f}, {y_O * 1e3:.3f}) mm "
          f"(arc center, picked at C1 boundary)")
    print(f"  C       = ({x_O * 1e3:.3f}, {(y_O - r_c) * 1e3:.3f}) mm "
          f"(natural arc lowest point)")
    print()
    print(f"  C1 user:  r_c + r_ex < y_O + r_in")
    print(f"            {c1_lhs * 1e3:.5f} < {c1_rhs * 1e3:.5f}  ->  "
          f"{'PASS' if c1_ok else 'FAIL'}")
    print(f"  C2 user:  (x_O - L/2)^2 + y_O^2 < r_c^2")
    print(f"            {c2u_lhs:.6e} < {c2u_rhs:.6e}    ->  "
          f"{'PASS' if c2u_ok else 'FAIL'}")
    print(f"  C2* wall: (L/2)^2 + (y_O + c)^2 <= r_c^2  "
          f"(arc penetrates wall at rim)")
    print(f"            {c2c_lhs:.6e} <= {c2c_rhs:.6e}   ->  "
          f"{'PASS' if c2c_ok else 'FAIL'}")
    print()
    print(f"  Sagitta over full chord  = {sagitta_full * 1e3:.4f} mm")
    if chord_covers_rim:
        print(f"  Sagitta at rim (x=L/2)   = {sag_at_rim * 1e3:.4f} mm   "
              f"(chord 2*r_c*sin(theta/2) = {2 * L_crv_chord_half * 1e3:.3f} "
              f"mm covers L = {L * 1e3:.3f} mm)")
    else:
        print(f"  Curved chord half-length = {L_crv_chord_half * 1e3:.4f} mm "
              f"(< L/2 = {0.5 * L * 1e3:.4f} mm; arc does not reach rim)")
    print(f"  Sagitta / clearance      = {sagitta_full / c:.3f}")
    print(f"  Regime                   = {regime}")
    print()
    print(f"  [Note] The 'rims-only, no contact at C' regime that motivated")
    print(f"  the original C1+C2 constraints is geometrically unreachable")
    print(f"  for a symmetric arc inside a straight tube of constant")
    print(f"  clearance: substituting C1 into C2* gives (L/2)^2 <= 0.  See")
    print(f"  the module docstring for the full derivation.")
    print()

    # ---- Spatial layout report (ADDED) --------------------------------------
    # Where the proximal base, curve midpoint, curve end, and distal tip
    # actually land in world coordinates depends on curve_placement and
    # the inner_base_x picked by compute_inner_offset.  Print both so the
    # user can see how moving the placement parameter changes the
    # geometry.
    L_inner   = p_inner['str_length']
    placement = p_inner.get('curve_placement', 'distal')
    L_crv     = r_c * theta
    L_str     = max(0.0, L_inner - L_crv)
    sign      = p_inner.get('crv_sign', +1)
    # With base tangent +X (no rotation), the curve's tangent rotates by
    # sign*theta over the curve's arc length.  Curve position vs arc
    # length s in (1) on a circle of radius r_c centered at
    # (x_curve_start, y_circ) where y_circ = -sign * r_c:
    #   x(s) = x_curve_start + (r_c / sign) * sin(sign * s / r_c)
    #        = x_curve_start + r_c * sin(s / r_c)        (sin is odd; sign*sin(sign*x)=sin(x))
    # Wait that's not right either.  Cleanest form:
    #   x(s) = x_curve_start + (1/kappa) * sin(kappa * s),  kappa = sign / r_c
    #   y(s) =                 (1/kappa) * (1 - cos(kappa * s))
    # i.e., y is on the SAME side as kappa's sign, and the straight is
    # tangent at the ends.  For sign=-1, y is negative; for sign=+1, y is
    # positive.
    end_phi   = sign * theta                          # tangent angle at curve end

    if placement == 'centered':
        inner_base_x = -0.5 * L_inner                          # see compute_inner_offset
        x_curve_start = inner_base_x + 0.5 * L_str
    else:  # 'distal' (legacy) or any other placement collapses to this
        inner_base_x = 0.5 * L - L_inner
        x_curve_start = inner_base_x + L_str

    # Curve midpoint (arc length L_crv/2)
    x_curve_mid   = x_curve_start + (r_c / sign) * math.sin(0.5 * sign * theta)
    y_curve_mid   = (r_c / sign) * (1.0 - math.cos(0.5 * sign * theta))
    # Curve end (arc length L_crv)
    x_curve_end   = x_curve_start + (r_c / sign) * math.sin(sign * theta)
    y_curve_end   = (r_c / sign) * (1.0 - math.cos(sign * theta))

    if placement == 'centered':
        # Distal straight at angle end_phi, length L_str/2:
        x_tip = x_curve_end + 0.5 * L_str * math.cos(end_phi)
        y_tip = y_curve_end + 0.5 * L_str * math.sin(end_phi)
    else:
        x_tip = x_curve_end                                    # no distal straight
        y_tip = y_curve_end

    print(f"  Inner-tube layout:  curve_placement = {placement!r}")
    print(f"    base       @ x = {inner_base_x   * 1e3:+8.2f} mm,  "
          f"y = {0.0:+7.3f} mm")
    print(f"    curve start@ x = {x_curve_start  * 1e3:+8.2f} mm,  "
          f"y = {0.0:+7.3f} mm")
    print(f"    curve mid  @ x = {x_curve_mid    * 1e3:+8.2f} mm,  "
          f"y = {y_curve_mid * 1e3:+7.3f} mm     "
          f"<-- arc-length midpoint of the curve")
    print(f"    curve end  @ x = {x_curve_end    * 1e3:+8.2f} mm,  "
          f"y = {y_curve_end * 1e3:+7.3f} mm")
    print(f"    tip        @ x = {x_tip          * 1e3:+8.2f} mm,  "
          f"y = {y_tip * 1e3:+7.3f} mm     "
          f"<-- LOWEST y of natural shape (sign = {sign:+d})")
    print(f"    Outer beam @ x in [{-0.5 * L * 1e3:+.2f}, "
          f"{+0.5 * L * 1e3:+.2f}] mm")
    if placement == 'centered':
        eps_x = abs(x_curve_mid)
        print(f"    Curve midpoint x-offset from origin = {eps_x * 1e3:.4f} mm  "
              f"(O((theta)^3) error from 'midway' ideal)")
    print("=" * 76)


# =============================================================================
#  GEOMETRY HELPERS
# =============================================================================

def compute_inner_offset(p_outer, p_inner):
    """
    Return the absolute world-x positions of the two rigid bases under the
    new origin convention (world (0,0,0) at the GEOMETRIC CENTER of the
    outer beam).

    The placement of the inner tube depends on p_inner['curve_placement']:

      'distal' (legacy)
        Inner tube placed so its tip coincides with the outer tip when
        STRAIGHTENED:
            outer base @ x = -L_outer/2
            outer tip  @ x = +L_outer/2
            inner base @ x = +L_outer/2 - L_inner    (always negative)
            inner tip  @ x = +L_outer/2              (when straightened)
        The inner tube protrudes only on the proximal (left) side.

      'centered' (Reading B)
        Inner tube placed SYMMETRICALLY about the outer beam's center, so
        the curved portion (which sits at the inner tube's arc-length
        midpoint by construction; see compute_tube_geometry) lands AT
        x = 0 -- both in the STRAIGHTENED initial state and (to within a
        microscopic O(theta_arc^3) offset) in the natural shape:
            outer base @ x = -L_outer/2
            outer tip  @ x = +L_outer/2
            inner base @ x = -L_inner/2
            inner tip  @ x = +L_inner/2              (when straightened)
        The inner tube protrudes equally on both sides of the outer beam
        by (L_inner - L_outer) / 2.

    Parameters
    ----------
    p_outer, p_inner : dict   Outer/inner tube param dicts.

    Returns
    -------
    (outer_base_x, inner_base_x) : (float, float)
        Absolute x-positions [m] for the two rigid bases.

    Raises
    ------
    ValueError  If the inner tube is not strictly longer than the outer tube.
    """
    L_outer = p_outer['str_length']
    L_inner = p_inner['str_length']
    placement = p_inner.get('curve_placement', 'distal')      # [ADDED]

    if L_inner <= L_outer:
        raise ValueError(
            f"Inner tube '{p_inner['name']}' (L = {L_inner * 1e3:.1f} mm) must "
            f"be strictly longer than outer tube '{p_outer['name']}' "
            f"(L = {L_outer * 1e3:.1f} mm) for valid concentric placement."
        )

    outer_base_x = -0.5 * L_outer

    if placement == 'centered':                               # [ADDED]
        # Inner tube symmetric about origin -- protrudes equally on both
        # sides.  Curve midpoint (in arc length) is at x = 0 because the
        # curve sits at the inner tube's arc-length midpoint (compute_tube
        # _geometry's 'centered' layout) and the tube is straight at t=0.
        inner_base_x = -0.5 * L_inner
        protrusion = 0.5 * (L_inner - L_outer)
        print(
            f"[CTR] Origin at center of outer beam; placement = 'centered'.  "
            f"Outer base @ x = {outer_base_x * 1e3:.1f} mm,  "
            f"inner base @ x = {inner_base_x * 1e3:.1f} mm.  "
            f"Inner tube protrudes equally on both sides by "
            f"{protrusion * 1e3:.1f} mm."
        )
    elif placement == 'distal':
        # Legacy: tip-coincidence on the right.
        inner_base_x = 0.5 * L_outer - L_inner
        print(
            f"[CTR] Origin at center of outer beam; placement = 'distal'.  "
            f"Outer base @ x = {outer_base_x * 1e3:.1f} mm,  "
            f"inner base @ x = {inner_base_x * 1e3:.1f} mm.  "
            f"Tips coincide at x = +{0.5 * L_outer * 1e3:.1f} mm when "
            f"inner is straightened."
        )
    else:
        raise ValueError(
            f"Unknown curve_placement: {placement!r}.  "
            f"Expected 'distal' or 'centered'."
        )

    return outer_base_x, inner_base_x


def integrate_frame_positions(section_lengths, strains, frm_curv_abs,
                              x_offset=0.0):
    """
    Integrate piecewise-constant Cosserat strains along the arc-length to
    produce Rigid3d frame positions in WORLD coordinates, starting from
    rigid base (x_offset, 0, 0) with identity orientation.

    Strain convention follows the rest of the file: only the third
    component k_z is used (bending around local Z axis -> planar curve in
    the world XY plane, since base orientation is identity and twist is
    zero).  For a constant strain k_z over a section of length l:

      orientation rotates by (k_z * l) around Z
      local-frame translation:
        k_z != 0 :   ( sin(k_z*l)/k_z , (1 - cos(k_z*l))/k_z , 0 )
        k_z == 0 :   ( l ,             0 ,                    0 )

    [SAME AS ctr_two_tubes.py -- the helper handles negative kappa
    correctly because both sin and (1 - cos)/kappa are odd / even
    in kappa as needed.]
    """
    sec_ends = [0.0]
    for sl in section_lengths:
        sec_ends.append(sec_ends[-1] + sl)
    L_total = sec_ends[-1]

    frame_positions = []
    for s_target in frm_curv_abs:
        s_clip = min(s_target, L_total)

        x_w, y_w, z_w = x_offset, 0.0, 0.0
        theta = 0.0
        done = False

        for i, sl in enumerate(section_lengths):
            s_start = sec_ends[i]
            s_end   = sec_ends[i + 1]
            kappa   = strains[i][2]

            if s_clip <= s_end + 1e-12:
                ds = max(0.0, s_clip - s_start)
                if abs(kappa) < 1e-12:
                    dx_l, dy_l = ds, 0.0
                else:
                    dx_l = math.sin(kappa * ds) / kappa
                    dy_l = (1.0 - math.cos(kappa * ds)) / kappa

                ct, st = math.cos(theta), math.sin(theta)
                x_w += ct * dx_l - st * dy_l
                y_w += st * dx_l + ct * dy_l
                theta_at = theta + kappa * ds

                qz = math.sin(theta_at * 0.5)
                qw = math.cos(theta_at * 0.5)
                frame_positions.append([x_w, y_w, z_w, 0.0, 0.0, qz, qw])
                done = True
                break

            if abs(kappa) < 1e-12:
                dx_l, dy_l = sl, 0.0
            else:
                dx_l = math.sin(kappa * sl) / kappa
                dy_l = (1.0 - math.cos(kappa * sl)) / kappa

            ct, st = math.cos(theta), math.sin(theta)
            x_w += ct * dx_l - st * dy_l
            y_w += st * dx_l + ct * dy_l
            theta += kappa * sl

        if not done:
            qz = math.sin(theta * 0.5)
            qw = math.cos(theta * 0.5)
            frame_positions.append([x_w, y_w, z_w, 0.0, 0.0, qz, qw])

    return frame_positions


def compute_tube_geometry(p, x_offset=0.0,
                         init_strategy='natural',
                         outer_params=None):
    """
    Build Cosserat discretization for a tube.

    Layout: [BASE] --- straight (L_str) ---+--- curved arc (L_crv) --- [TIP]

    Modified from ctr_two_tubes.py:
    -------------------------------
    * Reads p.get('crv_sign', +1) and applies it to kappa so the curved
      portion can curve about local +z OR local -z.  +1 reproduces the
      old behaviour (arc opens "downward" in the local frame, i.e. the
      tangent rotates CCW about +z).  -1 makes the arc open upward
      (tangent rotates CW about +z), which is what the new scene wants.
    * crv_angle_deg = 0 produces n_str = ns and n_crv = 0 (a pure
      straight beam) -- this is exactly what we want for T1_PARAMS.

    init_strategy
    -------------
      'straight'         : init_states = 0 everywhere (tube starts perfectly
                           straight, snaps to rest_states on the first
                           integration step).  For the inner tube in this
                           scene, this is the canonical choice -- it places
                           the inner tube concentric with the (straight)
                           outer beam at t=0.
      'natural'          : init_states = rest_states (no initial release of
                           strain energy).  For T1_PARAMS this collapses to
                           init_states = 0 since rest_states is all zero
                           when crv_angle_deg = 0.
      'conform_to_outer' : Inner-tube CTR case (same as ctr_two_tubes.py).
                           If the OUTER tube is straight (theta_o = 0), this
                           reduces to init_states = 0 too -- equivalent to
                           'straight' here.

    See ctr_two_tubes.py compute_tube_geometry docstring for the full
    background on the conform_to_outer arc-length-weighted averaging.
    """
    L     = p['str_length']
    R     = p['crv_radius']
    theta = math.radians(p['crv_angle_deg'])
    ns    = p['nb_sections']
    nf    = p['nb_frames']
    sign  = p.get('crv_sign', +1)        # [ADDED] curvature sign about local +z
    placement = p.get('curve_placement', 'distal')   # [ADDED] 'distal' | 'centered'

    L_crv = R * theta
    L_str = max(0.0, L - L_crv)
    kappa = sign * 1.0 / R               # modified  was: kappa = 1.0 / R

    if L_str < 1e-12:
        n_str, n_crv = 0, ns
    elif L_crv < 1e-12:
        # [ADDED] Pure-straight tube branch (theta = 0).  Allocate every
        # section to the straight portion; rest_states will be all zero.
        n_str, n_crv = ns, 0
    else:
        frac_str = L_str / L
        n_str    = max(1, round(frac_str * ns))
        n_crv    = max(1, ns - n_str)
        n_str    = ns - n_crv

    ls = L_str / n_str if n_str > 0 else 0.0
    lc = L_crv / n_crv if n_crv > 0 else 0.0

    # [ADDED] Layout selection.  The two layouts share section_lengths
    # totals but differ in WHERE the curved sections are placed:
    #
    #   'distal'   -- [straight L_str ─────][curved L_crv ─────]
    #                 base                                       tip
    #                 Legacy layout.  Lowest y of the natural shape is at
    #                 the distal tip (since the curve's tangent rotates
    #                 monotonically from 0 at the start to -theta_arc at
    #                 the end, and the tube continues to dip throughout).
    #
    #   'centered' -- [straight L_str/2 ─][curved L_crv ─][straight L_str/2 ─]
    #                 The curve sits at the ARC-LENGTH midpoint of the
    #                 inner tube.  When inner_base_x = -L_inner/2, the
    #                 curve midpoint (in arc length) lands at x = 0 in the
    #                 STRAIGHT initial state AND at x ~ 0 in the natural
    #                 shape (the small offset = r_c*sin(theta/2) - L_crv/2
    #                 is on the order of (theta_arc)^3 * r_c / 24, which
    #                 is microscopic for theta_arc < 10 deg).
    #
    # IMPORTANT geometric note for 'centered':
    # The lowest y of the natural shape is NOT at the curve midpoint --
    # it is at the distal tip, because the curve's tangent rotates
    # monotonically from 0 (at the start) through -theta_arc/2 (at the
    # midpoint) to -theta_arc (at the end), and the distal straight then
    # continues at -theta_arc below horizontal.  To get a TRUE symmetric
    # U with C at the curve midpoint as the lowest y, the rigid base must
    # be ROTATED by +theta_arc/2 about +z so the tangent enters the
    # curve at +theta_arc/2 and exits at -theta_arc/2.  That rotation
    # makes the STRAIGHT initial state non-concentric with the outer
    # beam, which conflicts with init_strategy='straight'.  The
    # 'centered' layout here is the most-symmetric thing achievable
    # while keeping a horizontal base tangent.  See verify_geometry()
    # for the spatial coordinates of the curve midpoint and the tip
    # under both layouts.
    if placement == 'centered' and n_str > 0 and n_crv > 0:
        n_str_a = n_str // 2                         # proximal half
        n_str_b = n_str - n_str_a                    # distal half (catches odd ns)
        section_lengths = ([ls] * n_str_a
                           + [lc] * n_crv
                           + [ls] * n_str_b)
        rest_states     = ([[0., 0., 0.]]   * n_str_a
                           + [[0., 0., kappa]] * n_crv
                           + [[0., 0., 0.]]   * n_str_b)
    elif placement == 'distal' or n_str == 0 or n_crv == 0:
        # Legacy layout (also the only sensible thing when there are no
        # straight sections to split or no curved sections at all).
        section_lengths = [ls] * n_str + [lc] * n_crv
        rest_states     = [[0., 0., 0.]] * n_str + [[0., 0., kappa]] * n_crv
    else:
        raise ValueError(
            f"Unknown curve_placement: {placement!r}.  "
            f"Expected 'distal' or 'centered'."
        )

    if init_strategy == 'straight':
        init_states = [[0., 0., 0.]] * ns

    elif init_strategy == 'natural':
        init_states = [list(s) for s in rest_states]

    elif init_strategy == 'conform_to_outer':
        if outer_params is None:
            raise ValueError(
                "init_strategy='conform_to_outer' requires outer_params "
                "(the param dict of the enclosing outer tube)."
            )

        L_o     = outer_params['str_length']
        R_o     = outer_params['crv_radius']
        theta_o = math.radians(outer_params['crv_angle_deg'])
        sign_o  = outer_params.get('crv_sign', +1)        # [ADDED] for symmetry
        L_crv_o = R_o * theta_o
        L_str_o = max(0.0, L_o - L_crv_o)
        kappa_o = sign_o * 1.0 / R_o                      # modified  was: kappa_o = 1.0 / R_o

        # x_offset interpretation: in this scene the outer base is at
        # x = -L_outer/2 (NOT 0).  The inner base is at x = -L_outer/2 - X
        # where X is the inner tube's arc-length coordinate at which it
        # crosses the outer base.  X = -(x_offset_inner - x_offset_outer)
        # is what the conform_to_outer logic actually needs -- a *relative*
        # offset between the two bases.  Because we only support straight-
        # outer here (theta_o = 0 -> kappa_o = 0 -> _outer_bend_angle == 0
        # for any (a,b)), this branch reduces to init_states = 0 regardless
        # of how X is computed.  We compute X anyway for parity with the
        # original implementation.
        outer_base_x = -0.5 * L_o
        X = -(x_offset - outer_base_x)

        def _outer_bend_angle(a, b):
            a = max(0.0, min(L_o, a))
            b = max(0.0, min(L_o, b))
            if b <= L_str_o:
                return 0.0
            if a >= L_str_o:
                return kappa_o * (b - a)
            return kappa_o * (b - L_str_o)

        init_states = []
        s_run = 0.0
        for sl in section_lengths:
            s_in1 = s_run      - X
            s_in2 = s_run + sl - X
            angle = _outer_bend_angle(s_in1, s_in2)
            avg_kappa = angle / sl if sl > 0 else 0.0
            init_states.append([0.0, 0.0, avg_kappa])
            s_run += sl
    else:
        raise ValueError(f"Unknown init_strategy: {init_strategy!r}")

    sec_curv_abs = [0.0]
    s = 0.0
    for sl in section_lengths:
        s += sl
        sec_curv_abs.append(round(s, 10))

    lf           = L / nf
    frm_curv_abs = [round(i * lf, 10) for i in range(nf + 1)]

    frame_positions = integrate_frame_positions(
        section_lengths, init_states, frm_curv_abs, x_offset
    )

    edge_indices = [[i, i + 1] for i in range(nf)]

    return (section_lengths, rest_states, init_states,
            sec_curv_abs, frame_positions, frm_curv_abs,
            edge_indices)


def tube_mass(p):
    ri, re = p['rin'], p['rex']
    return p['density'] * math.pi * (re**2 - ri**2) * p['str_length']


# =============================================================================
#  TUBE VISUAL MODEL  (HOLLOW cylindrical surface -- IDENTICAL to
#  ctr_two_tubes.py.  Re-imported below to avoid duplicating ~150 lines
#  of unchanged code.)
# =============================================================================

def build_hollow_tube_surface(rex, rin, n_frames_total, n_circle):
    """[Identical to ctr_two_tubes.build_hollow_tube_surface -- see that
    module's docstring for full details.]"""
    if not (rex > rin > 0.0):
        raise ValueError(
            f"Hollow tube surface requires rex > rin > 0 "
            f"(got rex={rex}, rin={rin})."
        )

    N      = n_circle
    TWO_PI = 2.0 * math.pi

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

    quads = []
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
    """[Identical to ctr_two_tubes.add_tube_visual.]"""
    name = p['name']
    re   = p['rex']
    ri   = p['rin']
    nf   = p['nb_frames']
    n_frames_total = nf + 1

    if color is None:
        color = p['color']
    if isinstance(color, str):
        color_str = color
    else:
        color_str = " ".join(str(v) for v in color)

    positions, quads, rigid_idx = build_hollow_tube_surface(
        re, ri, n_frames_total, n_circle
    )

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


# =============================================================================
#  ADD ONE COSSERAT TUBE
# =============================================================================

def add_cosserat_tube(root_node, p, x_offset=0.0,
                     init_strategy='natural',
                     outer_params=None,
                     fixed_directions=(0, 1, 1, 0, 1, 1)):           # [ADDED]
    """
    Build the full Cosserat beam hierarchy for one tube.

    Modified from ctr_two_tubes.py:
    -------------------------------
    * Added `fixed_directions` kwarg to control PartialFixedProjective-
      Constraint per tube.  In this scene:
         outer (Tube_outer) : (1,1,1,1,1,1) -- all 6 DOFs pinned
         inner (Tube_inner) : (0,0,0,0,1,0) -- ry pinned, all others free

    Returns: (base_mo, coss_mo, tube_node, frame_node, odesolver) -- same
    as ctr_two_tubes.py.
    """
    name = p['name']

    (section_lengths, rest_states, init_states,
     sec_curv_abs, frame_positions, frm_curv_abs,
     _edge_indices) = compute_tube_geometry(p, x_offset,
                                            init_strategy=init_strategy,
                                            outer_params=outer_params)

    re, ri = p['rex'], p['rin']
    nf     = p['nb_frames']
    mass   = tube_mass(p)

    I_sec      = math.pi / 4.0 * (re ** 4 - ri ** 4)
    L_avg      = sum(section_lengths) / len(section_lengths)
    compliance = L_avg / (p['E'] * I_sec)                    # [unused -- see note]

    tube_node   = root_node.addChild(name)
    solver_node = tube_node.addChild('SolverNode')
    odesolver   = solver_node.addObject('EulerImplicitSolver',
                          name='odesolver',
                          rayleighStiffness=0.2,
                          rayleighMass=0.1,
                          firstOrder=False)
    solver_node.addObject('SparseLDLSolver',
                          name='Solver',
                          template='CompressedRowSparseMatrixd')
    # Same as ctr_two_tubes.py -- using GenericConstraintCorrection.  See
    # the module-level note in that file on the BeamHookeLawForceField
    # addKToMatrix issue; the existing scene runs with this and we
    # preserve the pattern for parity with the user's known-good config.
    solver_node.addObject('GenericConstraintCorrection')

    # ---- Rigid base ---------------------------------------------------------
    rigid_base = solver_node.addChild(name + '_rigid_base')
    base_mo = rigid_base.addObject(
        'MechanicalObject',
        template='Rigid3d',
        name='cosserat_base_mo',
        position=[[x_offset, 0., 0., 0., 0., 0., 1.]],
        showObject=True,
        showObjectScale=0.001,
    )
    rigid_base.addObject('PartialFixedProjectiveConstraint',
                         name='proximal_bc',
                         fixedDirections=list(fixed_directions),
                         indices=[0])

    # ---- Cosserat strain state ----------------------------------------------
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
#  CONTROLLER  (modified for 5-DOF inner actuation, fixed outer)
# =============================================================================

class CTRController(Sofa.Core.Controller):
    """
    GUI-driven actuation for the curved-in-straight CTR.

    Differences from ctr_two_tubes.CTRController
    --------------------------------------------
    * The outer tube is fully pinned (no actuation).  We do NOT write its
      pose every step -- the PartialFixedProjective constraint with all-
      ones fixedDirections holds it in place.
    * The inner tube exposes 5 absolute-target DOFs:
        tx, ty, tz   [m]      -- translation of the proximal base
        rx           [rad]    -- twist about local +X
        rz           [rad]    -- yaw   about world  +Z
      Rotation about y is NOT exposed (fixedDirections=[0,0,0,0,1,0]
      pins it).
    * Each axis is rate-limited per step:
        translation_step_m   (shared)  -> caps |tx|,|ty|,|tz| step
        rotation_step_rad    (shared)  -> caps |rx|,|rz|     step
      Both read live from gui_5dof.snapshot() every step.
    * Quaternion composition: q = q_z(rz) * q_x(rx)  (intrinsic
      yaw-then-twist -- yaw rotates the whole frame, twist about the
      yawed local x).  See the module docstring for the qx,qy,qz,qw
      formulae.

    Phase semantics: same as ctr_two_tubes.CTRController.  dt ramping is
    unchanged.  Inner-tube pose is rewritten every step so the rigid base
    stays pinned at the cumulative target while the inner tube relaxes.
    """

    DT_RAMP_PER_STEP = 1.02

    def __init__(self,
                 root_node,
                 t2_base_mo,
                 t2_x_offset,
                 gui_bridge,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.root_node = root_node
        self.t2_base_mo = t2_base_mo
        self.t2_x0     = float(t2_x_offset)
        self.gui       = gui_bridge

        # Cumulative pose, tracked across steps and never reset.  All five
        # DOFs are absolute-target driven; no rate integration needed
        # because the per-step cap is applied directly when chasing the
        # target each step (see _step_toward).
        self._t2_tx_m   = 0.0       # cumulative tx in [m]
        self._t2_ty_m   = 0.0
        self._t2_tz_m   = 0.0
        self._t2_rx_rad = 0.0       # cumulative twist (about local +x)
        self._t2_rz_rad = 0.0       # cumulative yaw   (about world  +z)

        # ---- dt ramp state (identical to ctr_two_tubes.py) ----
        try:
            dt0 = float(self.root_node.dt.value)
        except Exception:
            dt0 = 1e-4
        self._dt_current = dt0
        self._dt_target  = dt0

        self._step = 0
        self._bc_relaxed = False

    # ------------------------------------------------------------------
    def onAnimateBeginEvent(self, event):
        self._step += 1
        snap  = self.gui.snapshot()
        phase = snap['phase']

        if phase == 'control' and not self._bc_relaxed:
            bc = self.t2_base_mo.getContext().proximal_bc
            with bc.fixedDirections.writeable() as fd:
                fd[:] = [0, 0, 0, 0, 1, 0]  # only ry pinned, rest controller-driven
            self._bc_relaxed = True
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

            # 2) Translation: rate-limit toward each slider target.
            max_tstep = float(snap['translation_step_m'])
            self._t2_tx_m = self._step_toward(
                self._t2_tx_m, snap['t2_tx_target_m'], max_tstep)
            self._t2_ty_m = self._step_toward(
                self._t2_ty_m, snap['t2_ty_target_m'], max_tstep)
            self._t2_tz_m = self._step_toward(
                self._t2_tz_m, snap['t2_tz_target_m'], max_tstep)

            # 3) Rotation: rate-limit toward each slider target.  Both
            #    rx (twist) and rz (yaw) are absolute targets in
            #    radians.  Per-step cap in [rad/step] is read live from
            #    the GUI; same per-step rationale as translation_step_m.
            max_rstep = float(snap['rotation_step_rad'])
            self._t2_rx_rad = self._step_toward(
                self._t2_rx_rad, snap['t2_rx_target_rad'], max_rstep)
            self._t2_rz_rad = self._step_toward(
                self._t2_rz_rad, snap['t2_rz_target_rad'], max_rstep)

        # 'initializing' (and 'waiting' if we ever get here) -> hold
        # cumulative state where it is.  Pose is written every step
        # regardless to keep the rigid base pinned at the initial
        # offset while strain DOFs of the inner tube relax.
        if phase == 'control':
            self._set_inner_pose(self.t2_base_mo,
                                 self._t2_tx_m, self._t2_ty_m, self._t2_tz_m,
                                 self._t2_rx_rad, self._t2_rz_rad,
                                 x0=self.t2_x0)

    # ------------------------------------------------------------------
    @staticmethod
    def _step_toward(current, target, max_step):
        """Move `current` toward `target` by at most `max_step`.  Same as
        ctr_two_tubes.py."""
        delta = target - current
        if delta >  max_step: return current + max_step
        if delta < -max_step: return current - max_step
        return target

    @staticmethod
    def _set_inner_pose(mo, tx, ty, tz, rx, rz, x0=0.0):
        """
        Write the inner tube's rigid base pose:
          x_world      = x0 + tx           (x0 = initial offset, tx = cumulative)
          y_world      = ty
          z_world      = tz
          orientation  = q_z(rz) * q_x(rx) (intrinsic yaw-then-twist)

        Quaternion components (alpha = rx, beta = rz):
          qx =  cos(beta/2) * sin(alpha/2)
          qy =  sin(alpha/2) * sin(beta/2)
          qz =  cos(alpha/2) * sin(beta/2)
          qw =  cos(alpha/2) * cos(beta/2)

        SOFA Rigid3d layout: [x, y, z, qx, qy, qz, qw].  When alpha = 0
        (no twist) and beta = 0 (no yaw) this reduces to (0,0,0,0,0,0,1)
        -- the identity orientation.  qy is non-zero only when BOTH rx
        and rz are non-zero (it is the geometric coupling between
        intrinsic-x and global-z rotations); it does NOT mean there is
        an independent rotation about y -- the Euler decomposition has
        zero y-angle, which is the constraint the user requested.
        """
        ha = 0.5 * rx
        hb = 0.5 * rz
        sa, ca = math.sin(ha), math.cos(ha)
        sb, cb = math.sin(hb), math.cos(hb)

        qx =  cb * sa
        qy =  sa * sb
        qz =  ca * sb
        qw =  ca * cb

        with mo.position.writeable() as pos:
            p = list(pos[0])
            p[0] = x0 + tx
            p[1] = ty
            p[2] = tz
            p[3] = qx
            p[4] = qy
            p[5] = qz
            p[6] = qw
            pos[0] = p


# =============================================================================
#  CREATE SCENE
# =============================================================================

def createScene(root_node):

    # ---- Required plugins  (same as ctr_two_tubes.py) ---------------------
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
    root_node.dt      = 1e-3

    root_node.addObject('DefaultVisualManagerLoop')
    root_node.addObject('FreeMotionAnimationLoop')
    root_node.addObject('BackgroundSetting', color=[1.0, 1.0, 1.0, 0])

    root_node.addObject('BlockGaussSeidelConstraintSolver',
                        name='ConstraintSolver',
                        tolerance=1e-2,
                        maxIterations=500)

    # Camera looking at the origin (center of outer beam).  Adjusted from
    # ctr_two_tubes.py since the world origin moved.
    root_node.addObject('Camera',
                        position=[0.4, -0.3, 0.3],     # modified  was: [0.5, -0.3, 0.3]
                        lookAt=[0.0, 0.0, 0.0])         # modified  was: [0.08, 0.0, 0.0]

    root_node.addObject('VisualStyle',
                        displayFlags='showVisualModels hideBehaviorModels '
                                     'hideCollisionModels '
                                     'hideBoundingCollisionModels '
                                     'hideForceFields '
                                     'hideInteractionForceFields '
                                     'hideWireframe '
                                     'hideMechanicalMappings')

    # ---- Geometry verification (ADDED) -------------------------------------
    # Print the regime the chosen parameters fall into BEFORE building
    # any SOFA objects, so the user sees this even if construction
    # subsequently fails.
    verify_geometry(T1_PARAMS, T2_PARAMS, x_O=0.0)

    # ---- Compute origin-shifted base positions -----------------------------
    outer_base_x, inner_base_x = compute_inner_offset(T1_PARAMS, T2_PARAMS)

    # ---- Build tubes -------------------------------------------------------
    # Outer tube: STRAIGHT, base at x=-L_outer/2, ALL DOFs fixed.
    # init_strategy='natural' here is equivalent to 'straight' since
    # crv_angle_deg=0 forces rest_states = 0 (see compute_tube_geometry).
    t1_base_mo, _, _, t1_frame_node, t1_solver = add_cosserat_tube(
        root_node, T1_PARAMS,
        x_offset=outer_base_x,
        init_strategy='natural',
        fixed_directions=(1, 1, 1, 1, 1, 1),
    )

    # Inner tube: PRE-CURVED, base at x=L_outer/2 - L_inner, starts STRAIGHT
    # concentric with outer.  PartialFixedProjective fixes only ry; tx, ty,
    # tz, rx, rz are driven by the CTRController.
    t2_base_mo, _, _, t2_frame_node, t2_solver = add_cosserat_tube(
        root_node, T2_PARAMS,
        x_offset=inner_base_x,
        init_strategy='straight',
        outer_params=T1_PARAMS,
        fixed_directions=(1, 1, 1, 1, 1, 1),
    )

    # ---- Visual models -----------------------------------------------------
    add_tube_visual(t1_frame_node, T1_PARAMS, color=T1_PARAMS['color'])
    add_tube_visual(t2_frame_node, T2_PARAMS, color=T2_PARAMS['color'])

    # ---- GUI bridge --------------------------------------------------------
    # 5-DOF inner-tube actuation, no outer-tube controls.  Slider ranges
    # are conservative defaults; tweak in createScene() per experiment.
    gui_bridge = CTRGuiBridgeStraightOuter(
        root_node=root_node,
        # Inner-tube absolute-target ranges (sliders run -max..+max):
        max_tx_m   = 0.04,                # +/- 4 cm  insertion/retraction
        max_ty_m   = 0.005,               # +/- 5 mm  lateral
        max_tz_m   = 0.005,               # +/- 5 mm  lateral
        max_rx_deg = 180.0,               # +/- 180 deg twist
        max_rz_deg = 15.0,                # +/- 15 deg yaw  (keeps tube inside outer beam)
        # Per-step rate limits (defaults match ctr_two_tubes.py for tx):
        default_trans_step_um = 50.0,
        default_rot_step_deg  = 0.05,
        init_dt               = 1e-4,
        default_control_dt    = 1e-3,
        dt_min                = 1e-6,
        dt_max                = 1e-1,
    )

    # ---- Controller --------------------------------------------------------
    root_node.addObject(CTRController(
        name='CTRController',
        root_node=root_node,
        t2_base_mo=t2_base_mo,
        t2_x_offset=inner_base_x,
        gui_bridge=gui_bridge,
    ))

    # ---- Custom contact pipeline (SSIM + BCM + CPULC, identical wiring)----
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
        contactConfiguration="nested",
        defaultNormal=DEFAULT_NORMAL,
        broadPhaseMarginFactor=1.5,
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
        mappingMode='contactPoints',
    )

    # contact_output.addObject(
    #     'ContactPointsUnilateralConstraint',
    #     name='cpuc',
    #     mu=0,
    #     contactTriads=bcm.getLinkPath() + '.contactTriads',
    #     gapSign=bcm.getLinkPath() + '.gapSign',
    # )

    # ---- Initialization & live monitors  (REUSED AS-IS) -------------------
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
    ))

    # ---- Pause until Initialize is clicked ---------------------------------
    root_node.animate = False