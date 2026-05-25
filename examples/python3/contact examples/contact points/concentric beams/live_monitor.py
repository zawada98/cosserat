# -*- coding: utf-8 -*-
"""
live_monitor.py
================================================================================
LiveContactMonitor — continuous diagnostic Sofa.Core.Controller for two-tube
CTR scenes using the SSIM + BeamContactMapping + ContactPointsUnilateralConstraint
pipeline.

PURPOSE
-------
Companion to InitializationMonitor.  Where InitializationMonitor fires ONCE at
the end of relaxation and renders a static PNG, this monitor pushes a fresh
(abscissae, gaps) frame to a CTRGuiBridge every `every_n_steps` simulation
steps.  The bridge's Tk thread renders the data into a matplotlib
FigureCanvasTkAgg embedded in a separate Toplevel window — a live plot.

THREADING
---------
This class lives on the SOFA thread.  It performs NO Tk / matplotlib calls
whatsoever.  Its only cross-thread interaction is one lock-guarded write into
the bridge via `bridge.push_contact_profile(step, abscissae_m, gaps_m)`.  Tk
widget mutation happens entirely on the GUI thread.  This is the same pattern
that InitializationMonitor.on_init_complete uses for signal_init_complete.

REQUIREMENTS
------------
- BeamContactMapping must be in 'contactPoints' mode (interleaved Pc_A/Pc_B in
  contact_mo.position; Pc_B at index 2k+1).  In 'gap' mode the contactMO holds
  δ = Pc_B − Pc_A and Pc_B is not directly available, so abscissa projection
  cannot be reconstructed without re-reading SSIM outputs separately.  A
  warning is emitted at construction if the mode is not 'contactPoints'.

INVALID-PAIR FILTER
-------------------
BCM pre-allocates `distances` to MAX_K slots and writes the sentinel
Vec3(kInvalidGap, 0, 0) with kInvalidGap = 1e9 in unused slots.  Any pair k
with δn ≥ INVALID_GAP_THRESHOLD = 1e8 is dropped.  When no valid pairs remain,
an empty frame is still pushed so the live plot can clear stale markers
rather than freeze on the previous frame's data.

USAGE
-----
    from live_monitor import LiveContactMonitor

    intersection_node.addObject(LiveContactMonitor(
        name='LiveMonitor',
        t2_MO              = t2_frame_node.FramesMO,
        bcm                = bcm,
        contact_mo         = contactMO,
        t2_frame_curv_abs  = list(t2_frame_node.cosseratMapping.curv_abs_output.value),
        bridge             = gui_bridge,
        every_n_steps      = 20,
    ))
"""

import numpy as np
import Sofa
import Sofa.Core


# Sentinel from BeamContactMapping.cpp (kInvalidGap = 1e9) marks unused slots.
INVALID_GAP_THRESHOLD = 1e8


class LiveContactMonitor(Sofa.Core.Controller):
    """
    Continuous-time pusher of (curvilinear-abscissa, normal-gap) frames into
    a CTRGuiBridge for live plotting.  See module docstring for full details.
    """

    # ------------------------------------------------------------------
    #  Construction
    # ------------------------------------------------------------------
    def __init__(self,
                 t2_MO,                       # Rigid3d FramesMO of Tube_3 (inner)
                 bcm,                         # BeamContactMapping object
                 contact_mo,                  # BCM output MO (Vec3d, contactPoints mode)
                 t2_frame_curv_abs,           # frame curv abs of Tube_3 (intrinsic)
                 bridge,                      # CTRGuiBridge instance
                 every_n_steps=20,            # throttle: push every N steps
                 constraint_solver=None,       # optional BGS solver for contact forces
                 contact_constraint=None,      # optional CPULC/CPULC2 object
                 force_unit_scale=1.0,         # model force unit -> newtons
                 force_conversion='auto',      # auto | impulse | lambda
                 *args, **kwargs):
        super().__init__(*args, **kwargs)

        # ---- External handles ----
        self.t2_MO       = t2_MO
        self.bcm         = bcm
        self.contact_mo  = contact_mo
        self.t2_curv_abs = np.asarray(t2_frame_curv_abs, dtype=float)
        self.bridge      = bridge
        self.constraint_solver = constraint_solver
        self.contact_constraint = contact_constraint
        self.force_unit_scale = float(force_unit_scale)
        self.force_conversion = str(force_conversion).lower()

        # ---- Tunables ----
        # Throttle.  At dt = 1e-4 s and 100 sim-steps / real-second, every=20
        # gives a 5 Hz push rate -- well below the GUI's 50 ms (=20 Hz)
        # repaint cadence, so the plot never skips a frame, and we don't
        # waste CPU pushing data the GUI couldn't display anyway.
        self.every = max(1, int(every_n_steps))

        # ---- Internal state ----
        self._step = 0
        self._force_warning_emitted = False
        self._force_conversion_warning_emitted = False
        self._force_conversion_mode = None
        self._empty_reason_counts = {}

        # ---- Sanity: contactPoints mode required ----
        # In 'gap' mode the contactMO holds δ = Pc_B − Pc_A, and the world
        # position of Pc_B alone is not directly recoverable from one MO.
        # The projection step would silently produce wrong abscissae.
        try:
            mode = str(self.bcm.mappingMode.value)
            if mode != 'contactPoints':
                Sofa.msg_warning(
                    'LiveContactMonitor',
                    f"BeamContactMapping is in '{mode}' mode; this monitor "
                    "expects 'contactPoints' mode (Pc_A/Pc_B interleaved). "
                    "Pc_B abscissa computation will be wrong.")
        except Exception:
            pass  # non-fatal

    # ------------------------------------------------------------------
    #  Per-step hook
    # ------------------------------------------------------------------
    def onAnimateEndEvent(self, event):
        self._step += 1

        # Throttle.  Cheap modulo guard before any data read.
        if self._step % self.every:
            return

        # ---- Read BCM distances ----
        dists = np.asarray(self.bcm.distances.value)        # (Kbuf, 3)
        if dists.size == 0:
            valid_k = self._active_contact_indices_from_constraint()
            gaps = np.zeros(valid_k.shape, dtype=float)
            if valid_k.size == 0:
                self._push_empty("empty_distances")
                return

        # Component 0 of d_distances is δn (normal gap).
        delta_n = dists[:, 0]

        # Filter out unused slots (kInvalidGap = 1e9 sentinel).
        valid_k = np.flatnonzero(delta_n < INVALID_GAP_THRESHOLD)
        if valid_k.size == 0:
            valid_k = self._active_contact_indices_from_constraint()
            gaps = np.zeros(valid_k.shape, dtype=float)
            if valid_k.size == 0:
                self._push_empty("no_valid_distances")
                return
        else:
            gaps = delta_n[valid_k]

        # ---- Read contact MO (contactPoints mode: Pc_A,Pc_B interleaved) ----
        pos = np.asarray(self.contact_mo.position.value)    # (2*MAX_K, 3)
        idx_B = 2 * valid_k + 1                             # Pc_B for each k

        # Defensive bound check: a wiring-time size mismatch (MAX_K wrong)
        # would crash here.  Skip rather than crash; the user sees no live
        # data but the simulation keeps running.
        if idx_B.size and idx_B.max() >= pos.shape[0]:
            self._note_empty_reason("contact_mo_too_small")
            return

        pcB  = pos[idx_B]                                    # (K, 3)
        forces_n = self._read_contact_forces(valid_k)

        # ---- Project Pc_B onto Tube_3 centerline polyline ----
        frames = np.asarray(self.t2_MO.position.value)            # (Nf, 7)
        centers = frames[:, :3]
        quats = frames[:, 3:7]
        if centers.shape[0] != self.t2_curv_abs.shape[0]:
            # Frame count mismatch (curv_abs was captured before the mapping
            # finished its first apply, etc.).  Skip this frame.
            self._note_empty_reason("frame_count_mismatch")
            return

        projections = [
            self._project_to_polyline_with_frame(
                p, centers, quats, self.t2_curv_abs)
            for p in pcB
        ]
        absc = np.array([proj[0] for proj in projections])
        angles = np.array([proj[1] for proj in projections])

        # Sort by abscissa for a clean line plot and consistent x-axis order.
        order = np.argsort(absc)
        self._push_contact_profile(absc[order], gaps[order], forces_n[order])
        self._push_surface_profile(absc[order], angles[order], gaps[order])

    # ------------------------------------------------------------------
    #  Helpers
    # ------------------------------------------------------------------
    def _push_empty(self, reason="empty"):
        """Send an empty frame so the live plot can clear stale markers."""
        self._note_empty_reason(reason)
        self._push_contact_profile(
            np.empty(0, dtype=float),
            np.empty(0, dtype=float),
            np.empty(0, dtype=float))
        self._push_surface_profile(
            np.empty(0, dtype=float),
            np.empty(0, dtype=float),
            np.empty(0, dtype=float))

    def _note_empty_reason(self, reason):
        count = self._empty_reason_counts.get(reason, 0) + 1
        self._empty_reason_counts[reason] = count
        if count in (1, 10, 50) or count % 200 == 0:
            try:
                dists_size = len(self.bcm.distances.value)
            except Exception:
                dists_size = -1
            try:
                triads_size = len(self.bcm.contactTriads.value)
            except Exception:
                triads_size = -1
            try:
                pos_size = len(self.contact_mo.position.value)
            except Exception:
                pos_size = -1
            print(
                f"[LiveContactMonitor] step={self._step} pushed empty "
                f"({reason}); counts={self._empty_reason_counts} "
                f"dists={dists_size} triads={triads_size} contactMO={pos_size}"
            )

    def _push_contact_profile(self, abscissae, gaps, forces):
        """
        Push the longitudinal live profile.  New GUI bridges accept forces as
        a fourth argument; older/custom bridges still accept the original
        three-argument signature.
        """
        try:
            self.bridge.push_contact_profile(
                self._step, abscissae, gaps, forces)
        except TypeError:
            self.bridge.push_contact_profile(self._step, abscissae, gaps)

    def _push_surface_profile(self, abscissae, angles, gaps):
        """
        Optional extension point for GUI bridges that can display where Pc_B
        lies around the inner tube's local cross-section.
        """
        push = getattr(self.bridge, 'push_contact_surface_profile', None)
        if push is not None:
            push(self._step, abscissae, angles, gaps)

    def _active_contact_indices_from_constraint(self):
        """
        Fallback when BCM.distances is empty/stale but CPULC2 still exposes the
        active BCM pair ids used by the constraint rows.
        """
        constraint = self._get_contact_constraint()
        if constraint is None:
            return np.empty(0, dtype=int)

        try:
            ids = np.asarray(
                constraint.activeContactPairIndices.value, dtype=int)
        except Exception:
            return np.empty(0, dtype=int)

        if ids.size == 0:
            return ids

        try:
            mo_size = len(self.contact_mo.position.value)
        except Exception:
            return np.empty(0, dtype=int)

        return ids[(ids >= 0) & (2 * ids + 1 < mo_size)]

    def _read_contact_forces(self, valid_k):
        """
        Return normal contact reactions aligned with `valid_k`.

        Newer ContactPointsUnilateralConstraint variants expose the solved
        normal lambda together with the original BCM pair index k.  That is
        the only exact way to align forces after activation filtering skips
        some geometric pairs.  For older plugin binaries, fall back to the
        previous solver-vector path but only for active rows inferred from the
        constraint activation tolerance.
        """
        indexed = self._read_indexed_contact_impulses(valid_k)
        if indexed is not None:
            return indexed

        solver = self._get_constraint_solver()
        if solver is None:
            return np.full(valid_k.shape, np.nan, dtype=float)

        try:
            solver.computeConstraintForces.value = True
        except Exception:
            try:
                solver.findData('computeConstraintForces').value = True
            except Exception:
                if not self._force_warning_emitted:
                    Sofa.msg_warning(
                        'LiveContactMonitor',
                        "Could not enable computeConstraintForces on the "
                        "constraint solver; force subplot will stay empty.")
                    self._force_warning_emitted = True
                return np.full(valid_k.shape, np.nan, dtype=float)

        try:
            raw = np.asarray(solver.constraintForces.value, dtype=float).reshape(-1)
        except Exception:
            return np.full(valid_k.shape, np.nan, dtype=float)

        if raw.size == 0:
            return np.full(valid_k.shape, np.nan, dtype=float)

        active_mask = self._active_mask_from_constraint(valid_k)
        forces = np.zeros(valid_k.shape, dtype=float)
        active_pos = np.flatnonzero(active_mask)
        count = active_pos.size
        if count == 0:
            return forces

        if raw.size >= count:
            lambdas = raw[:count]
        else:
            lambdas = np.full(count, np.nan, dtype=float)
            lambdas[:raw.size] = raw

        forces[active_pos] = self._lambdas_to_contact_forces(lambdas)
        return forces

    def _read_indexed_contact_impulses(self, valid_k):
        constraint = self._get_contact_constraint()
        if constraint is None:
            return None

        try:
            pair_ids = np.asarray(
                constraint.activeContactPairIndices.value, dtype=int).reshape(-1)
            lambdas = np.asarray(
                constraint.normalContactImpulses.value, dtype=float).reshape(-1)
        except Exception:
            return None

        if pair_ids.size == 0 or lambdas.size == 0:
            return np.zeros(valid_k.shape, dtype=float)

        n = min(pair_ids.size, lambdas.size)
        by_pair = {int(pair_ids[i]): float(lambdas[i]) for i in range(n)}
        out = np.zeros(valid_k.shape, dtype=float)
        for i, k in enumerate(valid_k):
            out[i] = by_pair.get(int(k), 0.0)
        return self._lambdas_to_contact_forces(out)

    def _active_mask_from_constraint(self, valid_k):
        """
        Best-effort fallback for older binaries that do not expose per-pair
        impulses.  It prevents inactive geometric pairs from consuming entries
        from solver.constraintForces, but cannot recover global cId offsets.
        """
        constraint = self._get_contact_constraint()
        if constraint is None:
            return np.ones(valid_k.shape, dtype=bool)

        try:
            activation = float(constraint.activationTolerance.value)
        except Exception:
            return np.ones(valid_k.shape, dtype=bool)

        try:
            dists = np.asarray(self.bcm.distances.value)
            gaps = dists[valid_k, 0]
        except Exception:
            return np.ones(valid_k.shape, dtype=bool)

        return gaps <= activation

    def _lambdas_to_contact_forces(self, lambdas):
        mode = self._resolved_force_conversion_mode()
        if mode == 'lambda':
            return lambdas * self.force_unit_scale

        dt = 1.0
        try:
            dt = float(self.getContext().getRoot().dt.value)
        except Exception:
            try:
                dt = float(self.getContext().dt.value)
            except Exception:
                pass
        if abs(dt) < 1e-20:
            dt = 1.0
        return (lambdas / dt) * self.force_unit_scale

    def _resolved_force_conversion_mode(self):
        """
        Return how the solved contact lambda should be displayed as force.

        In second-order dynamics the lambda exposed by CPULC2 is impulse-like,
        so the plotted average force is lambda / dt.  In first-order
        quasi-static Euler, the same Data field name is misleading: the
        position correction solve stores a force/reaction-like lambda already,
        and dividing by dt inflates the display by 1 / dt.
        """
        if self.force_conversion in ('impulse', 'impulses', 'lambda_over_dt'):
            return 'impulse'
        if self.force_conversion in ('lambda', 'force', 'raw'):
            return 'lambda'
        if self.force_conversion != 'auto':
            if not self._force_conversion_warning_emitted:
                Sofa.msg_warning(
                    'LiveContactMonitor',
                    f"Unknown force_conversion={self.force_conversion!r}; "
                    "falling back to auto.")
                self._force_conversion_warning_emitted = True

        if self._force_conversion_mode is None:
            self._force_conversion_mode = (
                'lambda' if self._uses_first_order_ode_solver()
                else 'impulse'
            )
        return self._force_conversion_mode

    def _uses_first_order_ode_solver(self):
        """
        Best-effort scene inspection.  The contact monitor is often attached
        outside the tube solver node, so search from root and accept the first
        ODE solver exposing a firstOrder Data field.  Existing non-first-order
        scenes keep the historical impulse/dt conversion.
        """
        contexts = []
        for obj in (self.t2_MO, self.contact_mo, self.bcm, self):
            try:
                ctx = obj.getContext()
                if ctx is not None and ctx not in contexts:
                    contexts.append(ctx)
            except Exception:
                pass
        try:
            root = self.getContext().getRoot()
            if root is not None and root not in contexts:
                contexts.append(root)
        except Exception:
            pass

        for ctx in contexts:
            solver = self._find_object_with_data_in_tree(
                ctx, ('firstOrder',), ('odesolver', 'EulerImplicitSolver'))
            if solver is None:
                continue
            try:
                return bool(solver.firstOrder.value)
            except Exception:
                try:
                    return bool(solver.findData('firstOrder').value)
                except Exception:
                    continue
        return False

    @staticmethod
    def _find_object_with_data_in_tree(node, data_names, preferred_names=()):
        for name in preferred_names:
            try:
                obj = node.getObject(name)
                if obj is not None:
                    for data_name in data_names:
                        try:
                            if obj.findData(data_name) is not None:
                                return obj
                        except Exception:
                            pass
            except Exception:
                pass
            try:
                obj = getattr(node, name)
                if obj is not None:
                    for data_name in data_names:
                        try:
                            if obj.findData(data_name) is not None:
                                return obj
                        except Exception:
                            pass
            except Exception:
                pass

        try:
            objects = list(node.objects)
        except Exception:
            objects = []
        for obj in objects:
            for data_name in data_names:
                try:
                    if obj.findData(data_name) is not None:
                        return obj
                except Exception:
                    pass

        children = []
        for attr in ('children', 'getChildren'):
            try:
                value = getattr(node, attr)
                children = list(value() if callable(value) else value)
                break
            except Exception:
                pass

        for child in children:
            found = LiveContactMonitor._find_object_with_data_in_tree(
                child, data_names, preferred_names)
            if found is not None:
                return found
        return None

    def _get_contact_constraint(self):
        if self.contact_constraint is not None:
            return self.contact_constraint

        # Most scenes put CPULC in the same node as the BCM/contact MO, but
        # a few attach the monitor elsewhere.  Search both local context and
        # root, and accept either the old or new class name/object name.
        candidates = (
            'cpuc',
            'cpulc',
            'ContactPointsUnilateralConstraint',
            'ContactPointsUnilateralConstraint2',
        )
        contexts = []
        for obj in (self.contact_mo, self.bcm, self):
            try:
                ctx = obj.getContext()
                if ctx is not None and ctx not in contexts:
                    contexts.append(ctx)
            except Exception:
                pass
        try:
            root = self.getContext().getRoot()
            if root is not None and root not in contexts:
                contexts.append(root)
        except Exception:
            pass

        for ctx in contexts:
            found = self._find_object_in_tree(ctx, candidates)
            if found is not None:
                self.contact_constraint = found
                return found
        return None

    @staticmethod
    def _find_object_in_tree(node, names):
        for name in names:
            try:
                obj = node.getObject(name)
                if obj is not None:
                    return obj
            except Exception:
                pass
            try:
                obj = getattr(node, name)
                if obj is not None:
                    return obj
            except Exception:
                pass

        children = []
        for attr in ('children', 'getChildren'):
            try:
                value = getattr(node, attr)
                children = list(value() if callable(value) else value)
                break
            except Exception:
                pass

        for child in children:
            found = LiveContactMonitor._find_object_in_tree(child, names)
            if found is not None:
                return found
        return None

    def _get_constraint_solver(self):
        if self.constraint_solver is not None:
            return self.constraint_solver

        try:
            root = self.getContext().getRoot()
        except Exception:
            return None

        for name in (
            'ConstraintSolver',
            'BlockGaussSeidelConstraintSolver',
            'constraintSolver',
        ):
            try:
                solver = getattr(root, name)
                if solver is not None:
                    self.constraint_solver = solver
                    return solver
            except Exception:
                pass
            try:
                solver = root.getObject(name)
                if solver is not None:
                    self.constraint_solver = solver
                    return solver
            except Exception:
                pass

        return None

    @staticmethod
    def _project_to_polyline(point, centers, abscissae):
        """
        Project `point` ∈ R³ onto the polyline whose vertices are `centers`,
        return the curvilinear abscissa of the projection.

        Per segment [C_i, C_{i+1}] with intrinsic abscissae [s_i, s_{i+1}]:
            t   = clamp( (P - C_i) · (C_{i+1} - C_i) / ‖seg‖² , 0, 1 )
            d²  = ‖ P - (C_i + t·seg) ‖²
        Pick the segment with smallest d²; return  s_i + t · (s_{i+1} - s_i).

        Identical to InitializationMonitor._project_to_polyline -- duplicated
        here to keep this module self-contained (no cross-imports between
        diagnostic monitors).
        """
        best_d2  = np.inf
        best_abs = float(abscissae[0])
        N = len(centers)
        for i in range(N - 1):
            C0 = centers[i]
            C1 = centers[i + 1]
            seg = C1 - C0
            seg_l2 = float(seg @ seg)
            if seg_l2 < 1e-20:
                continue                                # degenerate segment
            t = float((point - C0) @ seg) / seg_l2
            t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
            proj = C0 + t * seg
            d2   = float((point - proj) @ (point - proj))
            if d2 < best_d2:
                best_d2  = d2
                best_abs = float(abscissae[i] + t * (abscissae[i + 1] - abscissae[i]))
        return best_abs

    @staticmethod
    def _project_to_polyline_with_frame(point, centers, quats, abscissae):
        """
        Project `point` onto the centerline polyline and return:

            (s, theta)

        where `s` is the curvilinear abscissa of the projection and `theta`
        is the angle of the contact point around the local cross-section:

            theta = atan2(local_z, local_y)

        The local frame is interpolated on the closest segment by nlerp of
        neighboring Rigid3d quaternions, then inverted to express
        (point - projected_center) in local tube coordinates.
        """
        best_d2 = np.inf
        best = None
        N = len(centers)
        for i in range(N - 1):
            C0 = centers[i]
            C1 = centers[i + 1]
            seg = C1 - C0
            seg_l2 = float(seg @ seg)
            if seg_l2 < 1e-20:
                continue
            t = float((point - C0) @ seg) / seg_l2
            t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
            proj = C0 + t * seg
            d2 = float((point - proj) @ (point - proj))
            if d2 < best_d2:
                q = LiveContactMonitor._quat_nlerp(quats[i], quats[i + 1], t)
                s = float(abscissae[i] + t * (abscissae[i + 1] - abscissae[i]))
                best_d2 = d2
                best = (s, proj, q)

        if best is None:
            return float(abscissae[0]), 0.0

        s, proj, q = best
        local = LiveContactMonitor._quat_rotate(
            LiveContactMonitor._quat_conj(q), point - proj)
        theta = float(np.arctan2(local[2], local[1]))
        return s, theta

    @staticmethod
    def _quat_nlerp(q0, q1, t):
        q0 = np.asarray(q0, dtype=float)
        q1 = np.asarray(q1, dtype=float)
        if float(q0 @ q1) < 0.0:
            q1 = -q1
        q = (1.0 - t) * q0 + t * q1
        n = float(np.linalg.norm(q))
        if n < 1e-24:
            return np.array([0.0, 0.0, 0.0, 1.0], dtype=float)
        return q / n

    @staticmethod
    def _quat_conj(q):
        return np.array([-q[0], -q[1], -q[2], q[3]], dtype=float)

    @staticmethod
    def _quat_rotate(q, v):
        qx, qy, qz, qw = q
        vx, vy, vz = v
        # Quaternion-vector rotation using q * [v,0] * q^-1.
        tx = 2.0 * (qy * vz - qz * vy)
        ty = 2.0 * (qz * vx - qx * vz)
        tz = 2.0 * (qx * vy - qy * vx)
        return np.array([
            vx + qw * tx + (qy * tz - qz * ty),
            vy + qw * ty + (qz * tx - qx * tz),
            vz + qw * tz + (qx * ty - qy * tx),
        ], dtype=float)
