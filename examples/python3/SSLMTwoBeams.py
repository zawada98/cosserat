# -*- coding: utf-8 -*-
"""
SphereSweptIntersectionMethod – demonstration scene
====================================================

Beam 1 – COMPLETELY FIXED  (acts as a rigid obstacle along +X)
Beam 2 – PARALLEL to Beam 1, CLAMPED at its base, free tip falls under gravity.

Beam 2 starts GAP_Z mm directly above Beam 1.  Gravity bends it downward
until the SSIM component detects contact along their shared length.

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
    runSofa scene_SSIM_TwoBeams.py
"""

import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")          # non-interactive backend – safe inside runSofa
import matplotlib.pyplot as plt
from matplotlib import cm

import Sofa
import Sofa.Core

# ──────────────────────────────────────────────────────────────────────────────
#  Scene parameters
# ──────────────────────────────────────────────────────────────────────────────
BEAM_LENGTH      = 120.0    # [mm]  total length of each beam
NB_SECTIONS      = 6        # number of Cosserat sections per beam
NB_FRAMES        = 12       # number of output Rigid3d frames per beam
RADIUS           = 1.0      # [mm]  cross-section radius
YOUNG_MODULUS    = 3.0e6    # [Pa]
POISSON_RATIO    = 0.49
STIFFNESS        = 1.0e8    # base-clamp stiffness  (Beam 2 clamped end)
ALGORITHM        = "ALGO_1" # "ALGO_1" (segment-seg) or "ALGO_2" (node-seg NR)
DT               = 0.01     # [s]   time step
MAX_STEPS        = 500      # stop automatically after this many steps
                             # set to 0 to disable auto-stop

# Output folder = same directory as this script
SCENE_DIR = os.path.dirname(os.path.abspath(__file__))

# ──────────────────────────────────────────────────────────────────────────────
#  Geometry
# ──────────────────────────────────────────────────────────────────────────────
#
#  Beam 1 – along global +X, centred at Y=0, Z=0  (fully fixed obstacle)
#    base at [0, 0, 0], quaternion [0,0,0,1]
#
#  Beam 2 – along global +Y, base at [L/2, 0, GAP_Z]
#    (so it crosses Beam 1 near the mid-span of both)
#    quaternion for 90° rotation about Z:  [0, 0, sin45°, cos45°]
#    With gravity in -Z the free end bends down and contacts Beam 1.
#
# Initial vertical gap between the two parallel beam axes
GAP_Z = 5.0    # [mm]   must be > 2*RADIUS to start without interpenetration


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
    return frames, curv_abs


def add_cosserat_beam(parent_node, name, base_pos, base_quat,
                      nb_sections, nb_frames, length, radius,
                      young_modulus, poisson_ratio, stiffness,
                      fully_fixed=False):
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
    frames,           curv_out = _make_frame_params(nb_frames, length)

    beam_node       = parent_node.addChild(name)
    rigid_base_node = beam_node.addChild('rigidBase')

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
        # ── Clamp the base with a stiff spring (allows small elastic reaction) ─
        rigid_base_node.addObject(
            'RestShapeSpringsForceField', name='clamp',
            stiffness=stiffness, angularStiffness=stiffness,
            external_points=0, points=0, template='Rigid3d')

    coord_node = beam_node.addChild('cosseratCoordinate')
    coord_mo   = coord_node.addObject(
        'MechanicalObject', template='Vec3d',
        name='cosserat_state', position=sec_pos)
    coord_node.addObject(
        'BeamHookeLawForceField',
        crossSectionShape='circular',
        length=sec_len, radius=radius,
        youngModulus=young_modulus,
        poissonRatio=poisson_ratio)

    if fully_fixed:
        # ── Fix all cosserat strains (zero = straight reference config) ────
        all_section_ids = list(range(nb_sections))
        coord_node.addObject(
            'FixedProjectiveConstraint', name='fixStrains',
            indices=all_section_ids)

    frame_node = rigid_base_node.addChild('mappedFrames')
    coord_node.addChild(frame_node)

    frames_mo = frame_node.addObject(
        'MechanicalObject', template='Rigid3d', name='FramesMO',
        position=frames, showObject=True, showObjectScale=2.0)

    if not fully_fixed:
        # Only the free beam needs mass; the fixed beam has none.
        frame_node.addObject('UniformMass', totalMass=0.1)

    frame_node.addObject(
        'DiscreteCosseratMapping', name='cosseratMapping',
        curv_abs_input=curv_in, curv_abs_output=curv_out,
        input1=coord_node.cosserat_state.getLinkPath(),
        input2=base_mo.getLinkPath(),
        output=frames_mo.getLinkPath(), debug=False)

    return frame_node


# ──────────────────────────────────────────────────────────────────────────────
#  SSIMPlotter  –  data collection + 3-D plots
# ──────────────────────────────────────────────────────────────────────────────

class SSIMPlotter(Sofa.Core.Controller):
    """
    Collects (t, s1*, s2*, δ) from the SSIM component every step,
    then generates and saves two 3-D plots on demand.

    Key bindings (inside the runSofa viewer)
    -----------------------------------------
      [Q]  →  save plots NOW and pause the simulation
      [P]  →  save plots NOW (simulation keeps running)

    Automatic stop
    --------------
      When step count reaches MAX_STEPS the simulation is paused and
      plots are saved (only if MAX_STEPS > 0).
    """

    def __init__(self, root_node, ssim_component,
                 beam1_frames_mo, beam2_frames_mo,
                 beam_length, dt, max_steps,
                 algorithm,
                 output_dir=SCENE_DIR,
                 print_every=50,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.root        = root_node
        self.ssim        = ssim_component
        self.beam1_mo    = beam1_frames_mo   # FramesMO – read live positions
        self.beam2_mo    = beam2_frames_mo
        self.algorithm   = algorithm         # label used in filenames + plots
        self.L           = beam_length
        self.dt          = dt
        self.max_steps   = max_steps if max_steps else None
        self.output_dir  = output_dir
        self.print_every = print_every

        self.step   = 0
        self._saved = False

        self._seg_len = beam_length / NB_FRAMES

        # ── History buffers (3-D plots) ───────────────────────────────────────
        self._times  = []
        self._s1_abs = []
        self._s2_abs = []
        self._deltas = []

        # ── Snapshot buffers (last-step 2-D plot) ─────────────────────────────
        # Overwritten every step; only the final state is retained.
        self._last_beam1_pos     = None   # np.ndarray (N_frames+1, 7) Rigid3d
        self._last_beam2_pos     = None
        self._last_contacts_snap = []     # list of (s1_abs, s2_abs, delta)
        self._last_t             = 0.0

        # ── Timing buffers ─────────────────────────────────────────────────────
        # _t_step_start is set in onAnimateBeginEvent just before the physics
        # solve; the elapsed wall-clock duration is stored in _step_times.
        self._t_step_start = None          # time.perf_counter() snapshot
        self._step_times   = []            # per-step wall-clock duration [s]

        print(f"\n[SSIMPlotter] ── Initialised ──────────────────────────────")
        print(f"[SSIMPlotter]   Algorithm     : {self.algorithm}")
        print(f"[SSIMPlotter]   Output folder : {self.output_dir}")
        print(f"[SSIMPlotter]   Auto-stop at  : step {self.max_steps}")
        print(f"[SSIMPlotter]   [Q] = save + stop    [P] = save only")
        print(f"[SSIMPlotter] ─────────────────────────────────────────────\n")

    def onAnimateBeginEvent(self, _event):
        """Record wall-clock time just before SOFA begins the physics step."""
        self._t_step_start = time.perf_counter()

    def onAnimateEndEvent(self, _event):
        self.step += 1
        t      = self.step * self.dt
        params = list(self.ssim.curvilinearParams.value)
        dists  = list(self.ssim.distances.value)
        ids    = list(self.ssim.contactSectionIds.value)

        # ── History (3-D plots) ───────────────────────────────────────────────
        for sv, d, ij in zip(params, dists, ids):
            s1_norm, s2_norm = float(sv[0]), float(sv[1])
            s1_abs = (int(ij[0]) + s1_norm) * self._seg_len
            s2_abs = (int(ij[1]) + s2_norm) * self._seg_len
            self._times .append(t)
            self._s1_abs.append(s1_abs)
            self._s2_abs.append(s2_abs)
            self._deltas.append(float(d))

        # ── Snapshot (last-step plot) – overwrite every step ─────────────────
        try:
            self._last_beam1_pos = np.array(self.beam1_mo.position.value)
            self._last_beam2_pos = np.array(self.beam2_mo.position.value)
        except Exception:
            pass
        self._last_contacts_snap = [
            ((int(ij[0]) + float(sv[0])) * self._seg_len,
             (int(ij[1]) + float(sv[1])) * self._seg_len,
             float(d))
            for sv, d, ij in zip(params, dists, ids)
        ]
        self._last_t = t

        # ── Timing ───────────────────────────────────────────────────────────
        if self._t_step_start is not None:
            self._step_times.append(time.perf_counter() - self._t_step_start)

        if self.step % self.print_every == 0:
            n     = len(dists)
            min_d = min(dists) if n else float('nan')
            # Running timing stats over the last print_every steps
            recent = self._step_times[-self.print_every:] if self._step_times else []
            t_mean = 1e6 * (sum(recent) / len(recent)) if recent else float('nan')
            t_max  = 1e6 * max(recent)                  if recent else float('nan')
            print(f"[SSIMPlotter] t={t:.3f}s  step={self.step:>5d}  "
                  f"contacts={n}  min_δ={min_d:+.4f} mm  "
                  f"step_time mean={t_mean:.1f} µs  max={t_max:.1f} µs  "
                  f"[{self.algorithm}]")

        if self.max_steps and self.step >= self.max_steps:
            print(f"\n[SSIMPlotter] MAX_STEPS={self.max_steps} reached. "
                  "Saving plots and pausing.")
            self._save_plots()
            self._stop_simulation()

    def onKeypressedEvent(self, event):
        key = event.get('key', '')
        if key in ('Q', 'q'):
            print("\n[SSIMPlotter] [Q] – saving plots and stopping.")
            self._save_plots()
            self._stop_simulation()
        elif key in ('P', 'p'):
            print("\n[SSIMPlotter] [P] – saving intermediate plots.")
            self._save_plots(tag=f"_step{self.step:05d}")

    def __del__(self):
        if not self._saved and len(self._times) > 1:
            print("[SSIMPlotter] Saving plots on exit (destructor fallback).")
            try:
                self._save_plots(tag="_on_exit")
            except Exception as exc:
                print(f"[SSIMPlotter] WARNING – could not save: {exc}")

    def _stop_simulation(self):
        try:
            self.root.animate.value = False
        except Exception:
            pass

    def _save_plots(self, tag=""):
        if self._saved and not tag:
            print("[SSIMPlotter] Plots already saved – skipping.")
            return
        if len(self._times) < 2:
            print("[SSIMPlotter] Not enough data yet.")
            return

        t  = np.asarray(self._times,  dtype=float)
        s1 = np.asarray(self._s1_abs, dtype=float)
        s2 = np.asarray(self._s2_abs, dtype=float)
        d  = np.asarray(self._deltas, dtype=float)

        os.makedirs(self.output_dir, exist_ok=True)

        algo = self.algorithm   # short prefix for all filenames

        self._plot_3d(
            s_vals  = s1,
            t_vals  = t,
            d_vals  = d,
            s_label = r"$s_1^*$  [mm]  (curvilinear abscissa – Beam 1, fixed)",
            title   = rf"Interpenetration $\delta$ vs $s_1$ and time  [{algo}]",
            fname   = os.path.join(self.output_dir,
                                   f"{algo}_delta_vs_s1_time{tag}.png"),
        )

        self._plot_3d(
            s_vals  = s2,
            t_vals  = t,
            d_vals  = d,
            s_label = r"$s_2^*$  [mm]  (curvilinear abscissa – Beam 2, free tip)",
            title   = rf"Interpenetration $\delta$ vs $s_2$ and time  [{algo}]",
            fname   = os.path.join(self.output_dir,
                                   f"{algo}_delta_vs_s2_time{tag}.png"),
        )

        self._save_snapshot_plot(tag=tag)
        self._save_timing_plot(tag=tag)

        if not tag:
            self._saved = True


    def _save_snapshot_plot(self, tag=""):
        """
        2-D snapshot of both beam centrelines at the **last simulated step**.

        Saved as  snapshot_last_step{tag}.png

        Layout (two stacked panels)
        ---------------------------
        Top panel – X-Z side view
            • Beam 1 centreline: straight horizontal line at Z = 0 (blue).
            • Beam 2 centreline: deformed shape read from FramesMO positions
              (orange).
            • ±radius shading around each centreline shows the physical
              cross-section footprint.
            • Each detected contact pair is marked with a filled dot on both
              centrelines, colour-coded by gap (red = penetration, green = gap).
            • A vertical double-headed arrow annotates the gap magnitude δ
              between the two contact points.

        Bottom panel – gap profile δ(s) at the last step
            • X-axis: curvilinear abscissa s along the beam [mm].
            • Y-axis: signed gap δ [mm] (negative = interpenetration).
            • Horizontal dashed line at δ = 0 for reference.
            • Markers colour-coded red / green by sign of δ.
        """

        if self._last_beam1_pos is None or self._last_beam2_pos is None:
            print("[SSIMPlotter] Snapshot skipped – no frame data yet.")
            return

        p1 = self._last_beam1_pos  # (N+1, 7) Rigid3d
        p2 = self._last_beam2_pos
        contacts = self._last_contacts_snap
        t_label = f"t = {self._last_t:.3f} s"

        # ── Build the curvilinear abscissa axis for the frames ────────────────
        n_frames = len(p1)  # NB_FRAMES + 1
        frame_len = self.L / (n_frames - 1)  # e.g. 120/12 = 10 mm
        s_frames = np.arange(n_frames) * frame_len  # [0, 10, …, 120]

        z1 = p1[:, 2]  # Z of each Beam-1 frame
        z2 = p2[:, 2]  # Z of each Beam-2 frame

        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, figsize=(12, 8),
            gridspec_kw={"height_ratios": [3, 1.6]})
        fig.suptitle(
            f"Beam snapshot at last step  ({t_label})",
            fontsize=12, fontweight="bold")

        # ── Top panel: curvilinear-abscissa / Z side view ────────────────────
        ax_top.plot(s_frames, z1, color="#2c6fad", lw=2.5,
                    label="Beam 1 (fixed)", zorder=3)
        ax_top.fill_between(s_frames, z1 - RADIUS, z1 + RADIUS,
                            color="#2c6fad", alpha=0.18, zorder=2)

        ax_top.plot(s_frames, z2, color="#e07b39", lw=2.5,
                    label="Beam 2 (cantilever)", zorder=3)
        ax_top.fill_between(s_frames, z2 - RADIUS, z2 + RADIUS,
                            color="#e07b39", alpha=0.18, zorder=2)

        # Contact annotations
        for s1_abs, s2_abs, delta in contacts:
            # Interpolate Z using curvilinear abscissa — both axes now consistent
            z1_c = float(np.interp(s1_abs, s_frames, z1))
            z2_c = float(np.interp(s2_abs, s_frames, z2))
            colour = "#d62728" if delta < 0 else "#2ca02c"

            ax_top.plot(s1_abs, z1_c, 'o', color=colour, ms=7, zorder=5)
            ax_top.plot(s2_abs, z2_c, 'o', color=colour, ms=7, zorder=5)

            x_mid = (s1_abs + s2_abs) / 2.0
            ax_top.annotate(
                "", xy=(s2_abs, z2_c), xytext=(s1_abs, z1_c),
                arrowprops=dict(arrowstyle="<->", color=colour, lw=1.5,
                                shrinkA=0, shrinkB=0),
                zorder=6)
            ax_top.text(
                (s1_abs + s2_abs) / 2.0 + 0.5,
                (z1_c + z2_c) / 2.0,
                f"δ={delta:+.3f}", fontsize=7, color=colour, va="center")

        ax_top.set_xlim(0, self.L)  # enforce [0, L] regardless of contact positions
        ax_top.set_xlabel("Curvilinear abscissa  s  [mm]", fontsize=9)


        ax_top.set_xlabel("Curvilinear abscissa  s  [mm]", fontsize=9)
        ax_top.set_ylabel("Z position  [mm]", fontsize=9)
        ax_top.set_title("Side view (X-Z plane)", fontsize=10)
        ax_top.legend(fontsize=9, loc="upper right")
        ax_top.axhline(0, color="grey", lw=0.6, ls="--", zorder=1)
        ax_top.grid(True, lw=0.4, alpha=0.5)

        # ── Bottom panel: gap profile δ(s) ───────────────────────────────────
        if contacts:
            s_vals  = np.array([(c[0] + c[1]) / 2.0 for c in contacts])
            d_vals  = np.array([c[2] for c in contacts])
            colours = ["#d62728" if dv < 0 else "#2ca02c" for dv in d_vals]
            sort_idx = np.argsort(s_vals)

            ax_bot.plot(s_vals[sort_idx], d_vals[sort_idx],
                        color="steelblue", lw=1.5, zorder=2)
            ax_bot.scatter(s_vals, d_vals, c=colours, s=40, zorder=3)
            ax_bot.axhline(0, color="grey", lw=1.2, ls="--", zorder=1,
                           label=r"$\delta = 0$")
            ax_bot.fill_between(s_vals[sort_idx], d_vals[sort_idx], 0,
                                where=d_vals[sort_idx] < 0,
                                color="#d62728", alpha=0.18,
                                label="penetration")
            ax_bot.fill_between(s_vals[sort_idx], d_vals[sort_idx], 0,
                                where=d_vals[sort_idx] >= 0,
                                color="#2ca02c", alpha=0.12,
                                label="gap")
            ax_bot.legend(fontsize=8, loc="upper right")
        else:
            ax_bot.text(0.5, 0.5, "No contact detected at this step",
                        ha="center", va="center",
                        transform=ax_bot.transAxes, fontsize=10, color="grey")

        ax_bot.set_xlabel("Curvilinear abscissa  s  [mm]", fontsize=9)
        ax_bot.set_ylabel(r"$\delta$  [mm]", fontsize=9)
        ax_bot.set_title(r"Gap profile $\delta(s)$ at last step", fontsize=10)
        ax_bot.grid(True, lw=0.4, alpha=0.5)

        plt.tight_layout()
        fname = os.path.join(self.output_dir, f"{self.algorithm}_snapshot_last_step{tag}.png")
        plt.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"[SSIMPlotter] ✓ Saved: {fname}")


    def _save_timing_plot(self, tag=""):
        """
        Plot the per-step wall-clock time of the full animation step.

        Because the SSIM intersection query is the dominant computation in
        this scene, the step time is a faithful proxy for the cost of the
        chosen algorithm.  Run the scene once with ALGO_1 and once with
        ALGO_2 to obtain two files:

            ALGO_1_ssim_timing.png
            ALGO_2_ssim_timing.png

        and compare them directly.

        Layout
        ------
        Main axes
            • Grey markers: raw per-step wall-clock time [µs].
            • Coloured line: rolling mean (window = min(20, n_steps)).
            • Horizontal dashed line: overall mean.
        Text box (upper right)
            Mean, median, std-dev, min, max, total cumulative time, and
            the algorithm label – everything needed for a quick comparison.
        """
        if len(self._step_times) < 2:
            print("[SSIMPlotter] Timing plot skipped – not enough data.")
            return

        dt_us   = np.array(self._step_times) * 1e6   # convert s → µs
        steps   = np.arange(1, len(dt_us) + 1)

        mean_v  = dt_us.mean()
        med_v   = np.median(dt_us)
        std_v   = dt_us.std()
        min_v   = dt_us.min()
        max_v   = dt_us.max()
        total_s = dt_us.sum() / 1e6                   # back to seconds

        # Rolling mean
        window   = min(20, len(dt_us))
        roll_mean = np.convolve(dt_us,
                                np.ones(window) / window,
                                mode='valid')
        roll_steps = steps[window - 1:]

        fig, ax = plt.subplots(figsize=(11, 5))

        ax.plot(steps, dt_us, color="lightgrey", lw=0.8,
                marker=".", ms=3, zorder=2, label="per-step time")
        ax.plot(roll_steps, roll_mean,
                color="#1f77b4", lw=2.0, zorder=3,
                label=f"rolling mean (w={window})")
        ax.axhline(mean_v, color="#d62728", lw=1.4, ls="--", zorder=4,
                   label=f"overall mean  {mean_v:.1f} µs")

        # Stats box
        stats_text = (
            f"Algorithm : {self.algorithm}\n"
            f"Steps     : {len(dt_us)}\n"
            f"Mean      : {mean_v:.2f} µs\n"
            f"Median    : {med_v:.2f} µs\n"
            f"Std dev   : {std_v:.2f} µs\n"
            f"Min       : {min_v:.2f} µs\n"
            f"Max       : {max_v:.2f} µs\n"
            f"Total     : {total_s:.4f} s"
        )
        ax.text(0.98, 0.97, stats_text,
                transform=ax.transAxes,
                fontsize=8, family="monospace",
                va="top", ha="right",
                bbox=dict(boxstyle="round,pad=0.5",
                          facecolor="white", edgecolor="grey",
                          alpha=0.85))

        ax.set_xlabel("Simulation step", fontsize=10)
        ax.set_ylabel("Wall-clock time per step  [µs]", fontsize=10)
        ax.set_title(
            f"SSIM step timing – {self.algorithm}",
            fontsize=12, fontweight="bold")
        ax.legend(fontsize=9, loc="upper left")
        ax.grid(True, lw=0.4, alpha=0.5)
        ax.set_xlim(left=1)
        ax.set_ylim(bottom=0)

        plt.tight_layout()
        fname = os.path.join(self.output_dir,
                             f"{self.algorithm}_ssim_timing{tag}.png")
        plt.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"[SSIMPlotter] ✓ Saved: {fname}")

        # Also print final summary to terminal
        print(f"\n[SSIMPlotter] ── Timing summary  [{self.algorithm}] ────────")
        print(f"[SSIMPlotter]   steps  : {len(dt_us)}")
        print(f"[SSIMPlotter]   mean   : {mean_v:.2f} µs")
        print(f"[SSIMPlotter]   median : {med_v:.2f} µs")
        print(f"[SSIMPlotter]   std    : {std_v:.2f} µs")
        print(f"[SSIMPlotter]   min    : {min_v:.2f} µs")
        print(f"[SSIMPlotter]   max    : {max_v:.2f} µs")
        print(f"[SSIMPlotter]   total  : {total_s:.4f} s")
        print(f"[SSIMPlotter] ─────────────────────────────────────────────\n")

    @staticmethod
    def _plot_3d(s_vals, t_vals, d_vals, s_label, title, fname):
        fig = plt.figure(figsize=(11, 7))
        ax  = fig.add_subplot(111, projection='3d')

        sc = ax.scatter(
            s_vals, t_vals, d_vals,
            c=d_vals, cmap='RdYlGn',
            vmin=d_vals.min(), vmax=max(d_vals.max(), 0.01),
            s=20, alpha=0.9, depthshade=True, zorder=5,
        )

        try:
            from scipy.interpolate import griddata
            s_unique = np.unique(np.round(s_vals, 2))
            t_unique = np.unique(np.round(t_vals, 3))
            if len(s_unique) >= 3 and len(t_unique) >= 3:
                n_s = min(50, len(s_unique))
                n_t = min(50, len(t_unique))
                si  = np.linspace(s_vals.min(), s_vals.max(), n_s)
                ti  = np.linspace(t_vals.min(), t_vals.max(), n_t)
                Si, Ti = np.meshgrid(si, ti)
                Di = griddata((s_vals, t_vals), d_vals,
                              (Si, Ti), method='linear')
                mask = ~np.isnan(Di)
                if mask.any():
                    ax.plot_surface(Si, Ti, Di,
                                    cmap='RdYlGn', alpha=0.22,
                                    linewidth=0, antialiased=True, zorder=2)
            else:
                sort_idx = np.argsort(t_vals)
                ax.plot(s_vals[sort_idx], t_vals[sort_idx], d_vals[sort_idx],
                        color='steelblue', linewidth=1.5, alpha=0.7,
                        label='contact trace', zorder=3)
        except (ImportError, Exception):
            pass

        s_edge = np.array([s_vals.min(), s_vals.max()])
        t_edge = np.array([t_vals.min(), t_vals.max()])
        SE, TE = np.meshgrid(s_edge, t_edge)
        ax.plot_surface(SE, TE, np.zeros_like(SE),
                        color='grey', alpha=0.18, zorder=1)
        ax.text(s_vals.mean(), t_vals.max(), 0.,
                r"$\delta = 0$", color='grey', fontsize=8)

        cbar = fig.colorbar(sc, ax=ax, shrink=0.45, pad=0.12)
        cbar.set_label(r"$\delta$ [mm]   (negative = penetration)", fontsize=9)

        ax.set_xlabel(s_label,               fontsize=9, labelpad=10)
        ax.set_ylabel("Simulation time [s]", fontsize=9, labelpad=10)
        ax.set_zlabel(r"$\delta$ [mm]",      fontsize=9, labelpad=10)
        ax.set_title(title, fontsize=11, pad=16)
        ax.view_init(elev=28, azim=-55)

        plt.tight_layout()
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"[SSIMPlotter] ✓ Saved: {fname}")


# ──────────────────────────────────────────────────────────────────────────────
#  SSIMMonitor  –  lightweight terminal logger
# ──────────────────────────────────────────────────────────────────────────────

class SSIMMonitor(Sofa.Core.Controller):
    """Prints a compact summary every `print_every` steps."""

    def __init__(self, ssim_component, print_every=50, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ssim  = ssim_component
        self.every = print_every
        self.step  = 0

    def onAnimateEndEvent(self, _event):
        self.step += 1
        if self.step % self.every != 0:
            return
        dists  = list(self.ssim.distances.value)
        params = list(self.ssim.curvilinearParams.value)
        ids    = list(self.ssim.contactSectionIds.value)
        if not dists:
            return
        print(f"\n── SSIM monitor  step {self.step} ───────────────────────")
        for k, (d, sv, ij) in enumerate(zip(dists, params, ids)):
            tag = "PENETRATING" if d < 0 else "gap"
            print(f"  #{k}: sections({ij[0]},{ij[1]})  "
                  f"s1*={sv[0]:.3f}  s2*={sv[1]:.3f}  "
                  f"δ={d:+.4f} mm  [{tag}]")


# ──────────────────────────────────────────────────────────────────────────────
#  createScene – SOFA entry point
# ──────────────────────────────────────────────────────────────────────────────

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

    Beam 2 (cantilever falling under gravity)
    -----------------------------------------
    • Lies along the global +X axis – PARALLEL to Beam 1.
    • Base clamped at [0, 0, GAP_Z] (directly above Beam 1's base), same
      orientation quaternion [0,0,0,1].
    • Clamped at its base by a RestShapeSpringsForceField.
    • Carries a UniformMass; gravity (−Z) bends it downward until it contacts
      Beam 1 along their shared length.
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
    ])
    root_node.addObject('DefaultVisualManagerLoop')
    root_node.addObject('DefaultAnimationLoop')
    root_node.addObject('BackgroundSetting', color=[0.05, 0.05, 0.12, 1.0])

    # Camera positioned to see both beams from a 3/4-angle view
    root_node.addObject('Camera',
                        position=[-40, -80, 120],
                        lookAt=[50, 0, 0])
    root_node.addObject('VisualStyle',
                        displayFlags='showVisualModels showBehaviorModels '
                                     'hideCollisionModels '
                                     'hideBoundingCollisionModels '
                                     'hideForceFields '
                                     'hideInteractionForceFields '
                                     'hideWireframe '
                                     'showMechanicalMappings')

    # ── Shared solver ─────────────────────────────────────────────────────────
    solver_node = root_node.addChild('solverNode')
    solver_node.addObject('EulerImplicitSolver',
                          rayleighStiffness=0.2, rayleighMass=0.1)
    solver_node.addObject('SparseLDLSolver', name='solver',
                          template='CompressedRowSparseMatrixd')
    solver_node.addObject('GenericConstraintCorrection')

    # ── Beam 1 – horizontal along +X, COMPLETELY FIXED ───────────────────────
    #
    #   Axis:  X ∈ [0, L]   Y = 0   Z = 0
    #   This beam never deforms or translates; it is a rigid obstacle.
    #
    beam1_frames = add_cosserat_beam(
        solver_node, 'Beam1_Fixed',
        base_pos   = [0., 0., 0.],
        base_quat  = [0., 0., 0., 1.],
        nb_sections = NB_SECTIONS,
        nb_frames   = NB_FRAMES,
        length      = BEAM_LENGTH,
        radius      = RADIUS,
        young_modulus  = YOUNG_MODULUS,
        poisson_ratio  = POISSON_RATIO,
        stiffness      = STIFFNESS,
        fully_fixed    = True,        # ← freeze every DOF
    )

    # ── Beam 2 – along +X (parallel to Beam 1), CLAMPED at base, free end falls ─
    #
    #   Same orientation as Beam 1 (quaternion [0,0,0,1]).
    #   Base at [0, 0, GAP_Z] so both beams share the same X axis range and
    #   Beam 2 starts GAP_Z mm directly above Beam 1.
    #
    #   Gravity acts in -Z: the free end bends downward until it contacts Beam 1.
    #   The base end (X = 0) is held by the RestShapeSpringsForceField clamp.
    #
    beam2_frames = add_cosserat_beam(
        solver_node, 'Beam2_Cantilever',
        base_pos   = [0., 0., GAP_Z],
        base_quat  = [0., 0., 0., 1.],            # same orientation as Beam 1
        nb_sections = NB_SECTIONS,
        nb_frames   = NB_FRAMES,
        length      = BEAM_LENGTH,
        radius      = RADIUS,
        young_modulus  = YOUNG_MODULUS,
        poisson_ratio  = POISSON_RATIO,
        stiffness      = STIFFNESS,
        fully_fixed    = False,       # ← free to deform; clamped at base only
    )

    # ── SSIM component ────────────────────────────────────────────────────────
    ssim = solver_node.addObject(
        'SphereSweptIntersectionMethod',
        name            = 'ssim',
        beam1Frames     = beam1_frames.FramesMO.position.getLinkPath(),
        beam2Frames     = beam2_frames.FramesMO.position.getLinkPath(),
        radius1         = RADIUS,
        radius2         = RADIUS,
        algorithmType   = ALGORITHM,
        maxNRIterations = 20,
        nrTolerance     = 1e-12,
    )

    # ── Attach both controllers ───────────────────────────────────────────────
    root_node.addObject(SSIMMonitor(
        ssim_component = ssim,
        print_every    = 50,
        name           = 'monitor',
    ))

    root_node.addObject(SSIMPlotter(
        root_node       = root_node,
        ssim_component  = ssim,
        beam1_frames_mo = beam1_frames.FramesMO,
        beam2_frames_mo = beam2_frames.FramesMO,
        algorithm       = ALGORITHM,
        beam_length     = BEAM_LENGTH,
        dt              = DT,
        max_steps       = MAX_STEPS,
        output_dir      = SCENE_DIR,
        print_every     = 50,
        name            = 'plotter',
    ))

    return root_node