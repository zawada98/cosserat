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
                 *args, **kwargs):
        super().__init__(*args, **kwargs)

        # ---- External handles ----
        self.t2_MO       = t2_MO
        self.bcm         = bcm
        self.contact_mo  = contact_mo
        self.t2_curv_abs = np.asarray(t2_frame_curv_abs, dtype=float)
        self.bridge      = bridge

        # ---- Tunables ----
        # Throttle.  At dt = 1e-4 s and 100 sim-steps / real-second, every=20
        # gives a 5 Hz push rate -- well below the GUI's 50 ms (=20 Hz)
        # repaint cadence, so the plot never skips a frame, and we don't
        # waste CPU pushing data the GUI couldn't display anyway.
        self.every = max(1, int(every_n_steps))

        # ---- Internal state ----
        self._step = 0

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
            # Brand-new step before SSIM/BCM has run.  Push empty so the
            # live plot can clear stale data instead of locking on the
            # previous frame.
            self._push_empty()
            return

        # Component 0 of d_distances is δn (normal gap).
        delta_n = dists[:, 0]

        # Filter out unused slots (kInvalidGap = 1e9 sentinel).
        valid_k = np.flatnonzero(delta_n < INVALID_GAP_THRESHOLD)
        if valid_k.size == 0:
            self._push_empty()
            return

        # ---- Read contact MO (contactPoints mode: Pc_A,Pc_B interleaved) ----
        pos = np.asarray(self.contact_mo.position.value)    # (2*MAX_K, 3)
        idx_B = 2 * valid_k + 1                             # Pc_B for each k

        # Defensive bound check: a wiring-time size mismatch (MAX_K wrong)
        # would crash here.  Skip rather than crash; the user sees no live
        # data but the simulation keeps running.
        if idx_B.size and idx_B.max() >= pos.shape[0]:
            return

        pcB  = pos[idx_B]                                    # (K, 3)
        gaps = delta_n[valid_k]                              # (K,)

        # ---- Project Pc_B onto Tube_3 centerline polyline ----
        centers = np.asarray(self.t2_MO.position.value)[:, :3]   # (Nf, 3)
        if centers.shape[0] != self.t2_curv_abs.shape[0]:
            # Frame count mismatch (curv_abs was captured before the mapping
            # finished its first apply, etc.).  Skip this frame.
            return

        absc = np.array([
            self._project_to_polyline(p, centers, self.t2_curv_abs)
            for p in pcB
        ])

        # Sort by abscissa for a clean line plot and consistent x-axis order.
        order = np.argsort(absc)
        self.bridge.push_contact_profile(
            self._step, absc[order], gaps[order])

    # ------------------------------------------------------------------
    #  Helpers
    # ------------------------------------------------------------------
    def _push_empty(self):
        """Send an empty frame so the live plot can clear stale markers."""
        self.bridge.push_contact_profile(
            self._step, np.empty(0, dtype=float), np.empty(0, dtype=float))

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