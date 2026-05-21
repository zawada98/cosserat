# -*- coding: utf-8 -*-
"""
init_monitor.py
================================================================================
InitializationMonitor — diagnostic Sofa.Core.Controller for two-tube CTR scenes
that use the SSIM + BeamContactMapping + ContactPointsUnilateralConstraint
pipeline.

PURPOSE
-------
Detects the END of the initialization phase (the inner tube relaxing from its
conform-to-outer initial state into contact with the outer tube) and pops up
ONCE a 2D plot of the normal gap δn at every contact pair k, plotted against
the curvilinear abscissa of Pc_B[k] along the inner tube (Tube_3).

This is purely diagnostic.  It does not modify any DOFs, does not push
constraints, and does not interact with the solver — it only READS state.

EQUILIBRIUM CRITERION
---------------------
Trips when the maximum LINEAR velocity magnitude across the FramesMO of BOTH
tubes stays below VEL_THRESHOLD for QUIET_STEPS consecutive steps, after a
WARMUP_STEPS warm-up.  Velocity is preferred over position deltas because:
  - SOFA already publishes FramesMO.velocity directly (no finite differencing).
  - It's the same signal as a position delta divided by dt, just one step less
    arithmetic.

Linear-only norm (the ωx,ωy,ωz angular components are ignored): translational
quietness is what governs whether interpenetration is still evolving.

CURVILINEAR ABSCISSA OF Pc_B
----------------------------
Neither SSIM nor BCM expose contact section IDs (i, j) as a SOFA Data field.
SSIM exposes only `curvilinearParams` (normalised {s1*, s2*} ∈ [0,1] WITHIN an
unknown segment), so we cannot reconstruct the absolute abscissa from s2* alone.

Workaround: project Pc_B[k] (read from contactMO.position[2k+1] in
'contactPoints' mode) onto Tube_3's centerline polyline {C_0, C_1, …, C_Nf}
made of FramesMO frame centers.  For each segment [C_i, C_{i+1}] the closest-
point parameter t ∈ [0,1] is computed in closed form, and the absolute abscissa
is interpolated from the intrinsic frame curvilinear-abscissa table
(curv_abs_output, the same vector you passed to DiscreteCosseratMapping).

INVALID-PAIR FILTER
-------------------
BeamContactMapping pre-allocates BCM.distances to MAX_K slots and writes a
sentinel Vec3(kInvalidGap, 0, 0) with kInvalidGap = 1e9 in unused slots
(see BeamContactMapping.cpp:212, 366, 400, 424, 444).  We drop any pair k
whose δn = distances[k][0] ≥ INVALID_GAP_THRESHOLD = 1e8.

USAGE
-----
    monitor = InitializationMonitor(
        t1_MO              = t1_frame_node.FramesMO,
        t2_MO              = t2_frame_node.FramesMO,
        bcm                = bcm,
        contact_mo         = contactMO,
        t2_frame_curv_abs  = list(t2_frame_node.cosseratMapping.curv_abs_output.value),
        vel_threshold      = 1e-3,    # m/s — tune for your scene
        quiet_steps        = 100,     # consecutive sub-threshold steps
        warmup_steps       = 500,     # ignore equilibrium check before this step
        log_every          = 200,     # print max-vel every N steps (0 = silent)
        auto_open          = True,    # open PNG in OS image viewer when fired
        png_path           = None,    # default: ./init_phase_gap_profile.png
    )
    intersection_node.addObject(monitor)


WHY NOT A LIVE matplotlib POPUP?
--------------------------------
Inside runSofa, calling plt.show()/plt.ion() from onAnimateEndEvent crashes
with:
    Fatal Python error: PyEval_RestoreThread: the function must be called
    with the GIL held … (PyInit__tkinter / TclNRRunCallbacks / GLFW main loop)
SOFA's GLFW main loop owns the main thread + GIL state in a way that
matplotlib's Tk and Qt backends cannot share.  The robust workaround is the
headless 'Agg' backend → save PNG → os.startfile/xdg-open/open.  No event
loop, no GIL conflict.  See the comment block above _plot() for full details.
"""

import os
import numpy as np
import Sofa
import Sofa.Core


# Sentinel from BeamContactMapping.cpp (kInvalidGap = 1e9) marks unused slots.
INVALID_GAP_THRESHOLD = 1e8


class InitializationMonitor(Sofa.Core.Controller):
    """
    One-shot detector + visualizer of the post-init equilibrium gap profile.
    See module docstring for design rationale.
    """

    # ------------------------------------------------------------------
    #  Construction
    # ------------------------------------------------------------------
    def __init__(self,
                 t1_MO,                       # [ADDED] Rigid3d FramesMO of Tube_1
                 t2_MO,                       # [ADDED] Rigid3d FramesMO of Tube_3
                 bcm,                         # [ADDED] BeamContactMapping object
                 contact_mo,                  # [ADDED] BCM output MO (Vec3d, contactPoints mode)
                 t2_frame_curv_abs,           # [ADDED] frame curv abs of Tube_3 (intrinsic)
                 contact_constraint=None,      # optional CPULC/CPULC2 object for contact fallback
                 vel_threshold=1e-3,          # [ADDED] m/s
                 quiet_steps=100,             # [ADDED] consecutive sub-threshold steps
                 warmup_steps=500,            # [ADDED] step before which check is suppressed
                 require_contact=False,        # wait for absolute contact/gap criterion before quiet test
                 min_contact_pairs=1,          # minimum valid BCM pairs when require_contact=True
                 contact_gap_threshold=1e-4,   # [m] considered close/contact-ready
                 log_every=200,               # [ADDED] periodic max-vel report (0 = silent)
                 auto_open=True,              # [ADDED] open PNG in OS image viewer
                 png_path=None,               # [ADDED] explicit PNG path
                 on_init_complete=None,       # [MODIFIED] callable() invoked once when init settles
                 *args, **kwargs):
        super().__init__(*args, **kwargs)

        # ---- External handles ----
        self.t1_MO       = t1_MO
        self.t2_MO       = t2_MO
        self.bcm         = bcm
        self.contact_mo  = contact_mo
        self.t2_curv_abs = np.asarray(t2_frame_curv_abs, dtype=float)
        self.contact_constraint = contact_constraint

        # ---- Tunables ----
        self.vel_threshold = float(vel_threshold)
        self.quiet_steps   = int(quiet_steps)
        self.warmup_steps  = int(warmup_steps)
        self.require_contact = bool(require_contact)
        self.min_contact_pairs = int(min_contact_pairs)
        self.contact_gap_threshold = float(contact_gap_threshold)
        self.log_every     = int(log_every)
        self.auto_open     = bool(auto_open)
        self.png_path      = png_path or os.path.abspath('init_phase_gap_profile.png')

        # ---- External listener hook ----
        # [MODIFIED] One-shot callback invoked the moment equilibrium is detected.
        # Used by ctr_gui.CTRGuiBridge.signal_init_complete to enable the control
        # panel.  Called BEFORE the PNG plot is rendered so the user gets control
        # responsiveness immediately rather than waiting for matplotlib to save.
        # Any callable taking no arguments is acceptable; exceptions are caught
        # and logged so a faulty listener cannot bring down the simulation.
        self.on_init_complete = on_init_complete

        # ---- Internal state ----
        self._step  = 0
        self._quiet = 0
        self._fired = False

        # Sanity: monitor depends on contactPoints mapping mode (interleaved Pc_A/Pc_B).
        # gap mode would put a single Vec3 (δn,δt1,δt2) per pair — unsupported here.
        try:
            mode = str(self.bcm.mappingMode.value)
            if mode != 'contactPoints':
                Sofa.msg_warning(
                    'InitializationMonitor',
                    f"BeamContactMapping is in '{mode}' mode; this monitor expects "
                    "'contactPoints' mode (Pc_A/Pc_B interleaved in contactMO). "
                    "Pc_B abscissa computation will be wrong.")
        except Exception:
            pass  # non-fatal — Data field naming may vary across builds

    # ------------------------------------------------------------------
    #  Per-step hook (end-of-step: BCM.apply has already written distances)
    # ------------------------------------------------------------------
    def onAnimateEndEvent(self, event):
        if self._fired:
            return

        self._step += 1

        # Periodic progress so the user can pick threshold/quiet from data.
        v_max = self._max_linear_speed()
        contact_count, min_gap = self._contact_progress()
        contact_ok = (
            (not self.require_contact) or
            (contact_count >= self.min_contact_pairs and
             min_gap <= self.contact_gap_threshold)
        )
        if self.log_every and (self._step % self.log_every == 0):
            gap_text = "none" if contact_count == 0 else f"{min_gap:.3e} m"
            print(f"[InitMonitor | step {self._step:6d}] "
                  f"max |v_lin| = {v_max:.3e} m/s   "
                  f"valid_pairs = {contact_count:4d}   "
                  f"min_gap = {gap_text}   "
                  f"contact_ok = {contact_ok}   "
                  f"quiet_count = {self._quiet}/{self.quiet_steps}")

        if self._step < self.warmup_steps:
            return

        if contact_ok and v_max < self.vel_threshold:
            self._quiet += 1
        else:
            self._quiet = 0

        if self._quiet >= self.quiet_steps:
            self._fired = True
            self._on_initialization_complete(v_max)

    # ------------------------------------------------------------------
    #  Equilibrium signal
    # ------------------------------------------------------------------
    def _max_linear_speed(self):
        """
        Max |v_linear| across FramesMO of both tubes.

        Rigid3d Deriv layout: [vx, vy, vz, ωx, ωy, ωz].  We threshold on the
        first three components only (translational quietness governs
        interpenetration).  Angular spin doesn't change the gap.
        """
        v1 = np.asarray(self.t1_MO.velocity.value)   # (N1, 6)
        v2 = np.asarray(self.t2_MO.velocity.value)   # (N2, 6)

        m1 = float(np.linalg.norm(v1[:, :3], axis=1).max()) if v1.size else 0.0
        m2 = float(np.linalg.norm(v2[:, :3], axis=1).max()) if v2.size else 0.0
        return max(m1, m2)

    def _contact_progress(self):
        """
        Return (valid_pair_count, min_gap_m) from BCM.distances.

        This is an absolute state check, not a rate-of-change check.  With a
        tiny dt, slow gap evolution cannot satisfy this criterion early; the
        monitor waits until the gap itself is below contact_gap_threshold.
        """
        try:
            dists = np.asarray(self.bcm.distances.value)
        except Exception:
            return self._constraint_contact_progress()

        if dists.size == 0:
            return self._constraint_contact_progress()

        delta_n = dists[:, 0]
        valid = delta_n < INVALID_GAP_THRESHOLD
        if not np.any(valid):
            return self._constraint_contact_progress()

        active = delta_n[valid]
        return int(active.size), float(active.min())

    def _constraint_contact_progress(self):
        """
        Fallback for builds/steps where BCM public distance buffers are empty
        or stale while CPULC/CPULC2 has already built active contact rows.
        """
        if self.contact_constraint is None:
            return 0, np.inf

        try:
            impulses = list(self.contact_constraint.normalContactImpulses.value)
        except Exception:
            impulses = []

        if impulses:
            # Treat active impulses as "contact exists".  The min-gap value is
            # only used against contact_gap_threshold, so 0 is the conservative
            # touching value.
            return len(impulses), 0.0

        try:
            ids = list(self.contact_constraint.activeContactPairIndices.value)
        except Exception:
            ids = []

        if ids:
            return len(ids), 0.0

        return 0, np.inf

    # ------------------------------------------------------------------
    #  One-shot fire path
    # ------------------------------------------------------------------
    def _on_initialization_complete(self, v_at_trip):
        """Read SSIM/BCM outputs, compute abscissae, plot."""

        dists = np.asarray(self.bcm.distances.value)            # (Kbuf, 3)
        pos   = np.asarray(self.contact_mo.position.value)      # (2*MAX_K, 3)

        if dists.size == 0:
            print(f"[InitMonitor] step {self._step}: equilibrium reached, "
                  f"but bcm.distances is empty.  Nothing to plot.")
            return

        # Component 0 of d_distances is δn (normal gap).
        delta_n = dists[:, 0]

        # Filter out unused slots (kInvalidGap = 1e9 sentinel).
        valid_k = np.flatnonzero(delta_n < INVALID_GAP_THRESHOLD)
        if valid_k.size == 0:
            print(f"[InitMonitor] step {self._step}: equilibrium reached, "
                  f"but no valid contact pairs (all slots are kInvalidGap).")
            return

        # In contactPoints mode: contactMO[2k] = Pc_A,  contactMO[2k+1] = Pc_B.
        idx_B = 2 * valid_k + 1
        # Defensive bound check (size mismatch would be a wiring bug).
        if idx_B.max() >= pos.shape[0]:
            print(f"[InitMonitor] step {self._step}: contactMO size "
                  f"({pos.shape[0]}) too small for K={valid_k.size}. "
                  f"Increase MAX_K in the scene.")
            return

        pcB  = pos[idx_B]            # (K, 3)  world coordinates
        gaps = delta_n[valid_k]      # (K,)    normal gap [m]

        # Tube_3 frame centers (world coords) for projection.
        t2_centers = np.asarray(self.t2_MO.position.value)[:, :3]   # (Nf, 3)
        if t2_centers.shape[0] != self.t2_curv_abs.shape[0]:
            print(f"[InitMonitor] WARNING: Tube_3 frame count "
                  f"({t2_centers.shape[0]}) != curv_abs length "
                  f"({self.t2_curv_abs.shape[0]}). Plot may be inconsistent.")

        # Project each Pc_B onto Tube_3's centerline polyline → absolute abscissa.
        abscissae = np.array([
            self._project_to_polyline(p, t2_centers, self.t2_curv_abs)
            for p in pcB
        ])

        # Sort by abscissa for a clean line plot.
        order     = np.argsort(abscissae)
        abscissae = abscissae[order]
        gaps      = gaps[order]

        # ---- Console summary ----
        n_pen = int((gaps < 0).sum())
        print("=" * 72)
        print(f"[InitMonitor] INITIALIZATION COMPLETE at step {self._step}")
        print(f"              max |v_lin| at trip   = {v_at_trip:.3e} m/s")
        print(f"              valid contact pairs   = {valid_k.size}")
        print(f"              gap range             = "
              f"[{gaps.min():.3e}, {gaps.max():.3e}] m")
        print(f"              interpenetrating pairs= {n_pen}"
              + (f"   (worst δn = {gaps.min():.3e} m)" if n_pen else ""))
        print("=" * 72)

        # ---- Signal external listener (e.g. GUI bridge) ----
        # [MODIFIED] Fire the callback BEFORE the plot.  Setting a flag in the
        # GUI bridge is essentially free; rendering the PNG can take several
        # seconds.  We want the user's control panel to come alive the instant
        # equilibrium is detected, not after matplotlib finishes saving.
        if self.on_init_complete is not None:
            try:
                self.on_init_complete()
            except Exception as e:
                print(f"[InitMonitor] on_init_complete callback failed: {e!r}")

        # ---- Plot ----
        try:
            self._plot(abscissae, gaps, n_pen)
        except Exception as e:
            print(f"[InitMonitor] Plot failed: {e!r}")

    # ------------------------------------------------------------------
    #  Geometry helper
    # ------------------------------------------------------------------
    @staticmethod
    def _project_to_polyline(point, centers, abscissae):
        """
        Project `point` ∈ R³ onto the polyline whose vertices are `centers`,
        return the curvilinear abscissa of the projection.

        Per segment [C_i, C_{i+1}] with intrinsic abscissae [s_i, s_{i+1}]:
            t   = clamp( (P - C_i) · (C_{i+1} - C_i) / ‖seg‖² , 0, 1 )
            d²  = ‖ P - (C_i + t·seg) ‖²
        Pick the segment with smallest d²; return  s_i + t · (s_{i+1} - s_i).
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

    # ------------------------------------------------------------------
    #  Plotting  ── headless render + OS image viewer
    # ------------------------------------------------------------------
    #
    # WHY NOT A LIVE matplotlib POPUP?
    #
    # SOFA's runSofa runs the GLFW main loop on the main thread and owns the
    # GIL state in a way that conflicts with matplotlib's GUI backends (TkAgg,
    # QtAgg, …).  Calling plt.show()/plt.ion() from inside onAnimateEndEvent
    # crashes runSofa with:
    #
    #   Fatal Python error: PyEval_RestoreThread: the function must be called
    #   with the GIL held, after Python initialization and before Python
    #   finalization, but the GIL is released
    #     Current thread … (most recent call first):
    #       <no Python frame>
    #     PyInit__tkinter
    #     TclNRRunCallbacks
    #     …
    #     sofaglfw::SofaGLFWBaseGUI::runLoop
    #
    # Tk (and Qt) initialise their event loop assuming Python's main thread
    # state is available; under runSofa it is not, because the GLFW callback
    # context releases the GIL before the matplotlib backend tries to grab it.
    #
    # Robust workaround: render with the headless 'Agg' backend (no GUI, no
    # event loop), save a PNG, and let the OS image viewer (Windows Photo,
    # Preview, xdg-open) display it.  The simulation thread never touches a
    # GUI toolkit, so there is no GIL/thread conflict.
    # ------------------------------------------------------------------
    def _plot(self, abscissae, gaps, n_pen):
        """
        Render the gap profile to PNG with Agg, then auto-open in OS viewer.
        """
        import matplotlib
        import matplotlib.pyplot as plt
        plt.switch_backend('Agg')   # idempotent; overrides any GUI backend

        s_mm = abscissae * 1e3        # m -> mm  (tube length ≈ 330 mm)
        g_um = gaps      * 1e6        # m -> µm  (gap scale ≈ 100 µm)

        fig, ax = plt.subplots(figsize=(9.0, 4.5))

        # Reference line at δn = 0 (the wall-to-wall touching threshold).
        ax.axhline(0.0, color='k', linewidth=0.8, alpha=0.7)

        # Shade interpenetration region (δn < 0).
        if (g_um < 0).any():
            ax.fill_between(s_mm, g_um, 0.0,
                            where=(g_um < 0), interpolate=True,
                            color='red', alpha=0.18,
                            label=f'interpenetration ({n_pen} pair(s))')

        # Connecting line + per-pair markers, colored by gap sign.
        ax.plot(s_mm, g_um, '-', color='gray', linewidth=0.8, alpha=0.6)
        gmax = max(1e-12, float(np.abs(g_um).max()))
        sc = ax.scatter(s_mm, g_um, s=22, c=g_um, cmap='RdYlGn',
                        vmin=-gmax, vmax=gmax,
                        edgecolors='k', linewidths=0.4, zorder=3)

        cb = fig.colorbar(sc, ax=ax, pad=0.015)
        cb.set_label(r'$\delta_n$ [µm]')

        # ── Annotate interpenetrating points with their δn value ────────────
        # Why: when positive gaps reach hundreds of µm and penetration depths
        # are only a few µm, the penetrating markers sit invisibly on the y=0
        # axis.  The shaded fill_between underneath them is also too thin to
        # see.  Explicit text labels solve the visibility issue without
        # distorting the y-axis scale.
        #
        # Layout: a single label row near 30% of the positive y-range (or at
        # least 50 µm above the baseline), with a small horizontal stagger so
        # adjacent labels don't collide.  Each label has an arrow pointing
        # down to its data point.
        neg_mask = g_um < 0
        if neg_mask.any():
            # Re-mark penetrating points with a thicker red edge so the
            # *position* is locatable even if the value sits at y≈0.
            ax.scatter(s_mm[neg_mask], g_um[neg_mask],
                       s=70, facecolors='none',
                       edgecolors='red', linewidths=1.6, zorder=4)

            ymax_pos = float(g_um[g_um > 0].max()) if (g_um > 0).any() else 1.0
            label_y  = max(ymax_pos * 0.30, 50.0)   # at least 50 µm up

            neg_s = s_mm[neg_mask]
            neg_g = g_um[neg_mask]
            # Stagger label heights when several abscissae land close together,
            # so the boxes don't pile on top of each other.
            for i, (s_p, g_p) in enumerate(zip(neg_s, neg_g)):
                ax.annotate(
                    rf'$\delta_n$ = {g_p:+.2f} µm',
                    xy=(s_p, g_p),
                    xytext=(s_p, label_y * (1.0 + 0.18 * (i % 3))),
                    fontsize=8.5,
                    ha='center', va='bottom',
                    color='darkred',
                    bbox=dict(boxstyle='round,pad=0.3',
                              facecolor='white',
                              edgecolor='red',
                              alpha=0.9),
                    arrowprops=dict(arrowstyle='->',
                                    color='red',
                                    lw=1.0,
                                    shrinkA=2,
                                    shrinkB=3),
                    zorder=5,
                )

        ax.set_xlabel(r'Curvilinear abscissa of $P_{c,B}$ along Tube_3 [mm]')
        ax.set_ylabel(r'Normal gap $\delta_n$ [µm]')
        ax.set_title(f'Initialization-phase contact profile  '
                     f'(step {self._step},  K = {len(gaps)},  '
                     f'penetrating = {n_pen})')
        ax.grid(True, alpha=0.3)
        if (g_um < 0).any():
            ax.legend(loc='best', framealpha=0.85)

        fig.tight_layout()

        # ---- Save ----
        try:
            fig.savefig(self.png_path, dpi=150)
            print(f"[InitMonitor] Saved PNG -> {self.png_path}")
        except Exception as e:
            print(f"[InitMonitor] PNG save failed: {e!r}")
            plt.close(fig)
            return
        plt.close(fig)   # release figure resources

        # ---- Open in default image viewer ----
        if self.auto_open:
            self._open_in_system_viewer(self.png_path)

    # ------------------------------------------------------------------
    @staticmethod
    def _open_in_system_viewer(path):
        """
        Open `path` in the OS default image viewer.
        Windows: os.startfile (uses file association).
        macOS:   `open <path>`.
        Linux:   `xdg-open <path>`.

        All errors are swallowed and reported on stdout — failing to open the
        viewer must never crash the simulation.
        """
        import sys, subprocess
        try:
            if sys.platform.startswith('win'):
                os.startfile(path)                                # type: ignore[attr-defined]
            elif sys.platform == 'darwin':
                subprocess.Popen(['open', path])
            else:
                subprocess.Popen(['xdg-open', path])
        except Exception as e:
            print(f"[InitMonitor] Could not open viewer for {path}: {e!r}")
