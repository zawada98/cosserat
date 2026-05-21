#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_comparison.py
================================================================================
Post-processor: load ./comparison_cpulc.npz and ./comparison_feeder.npz, emit
overlay PNGs per snapshot, an animation scrubbing through snapshots, and
scalar comparison plots.

Outputs (under --out-dir, default ./plots/):
    overlay_snap_<NN>.png      one per snapshot, both modes overlaid
    overlay_animation.mp4      (or .gif if ffmpeg unavailable)
    scalar_max_pen.png         max |delta_n| vs snapshot index
    scalar_n_active.png        active pair count vs snapshot index
    scalar_iters.png           solver iterations vs simulation step
    scalar_residual.png        solver residual vs simulation step
    scalar_walltime.png        wall-clock per step (rolling median)
    summary.txt                aggregate totals

Usage:
    python plot_comparison.py
    python plot_comparison.py --cpulc path/to/comparison_cpulc.npz \
                              --feeder path/to/comparison_feeder.npz \
                              --out-dir ./plots
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter


COLOR_CPULC = '#1f77b4'   # tab:blue
COLOR_FEEDER = '#d62728'  # tab:red


# =============================================================================
def load_run(path: Path) -> dict:
    z = np.load(path, allow_pickle=True)
    return {k: z[k] for k in z.files}


def s_to_mm(arr):
    return np.asarray(arr, dtype=float) * 1e3


def gap_to_mm(arr):
    return np.asarray(arr, dtype=float) * 1e3


# =============================================================================
def plot_snapshot_overlay(out_dir: Path, A: dict, B: dict, snap_idx: int):
    """
    Two-panel snapshot PNG:
      top    -- delta_n(s) for both modes, with distinct line styles so
                even bit-identical curves are visually distinguishable.
      bottom -- absolute difference |gap_cpulc - gap_feeder|(s), log scale.
                If the curves agree to floating-point noise, this panel
                will show very small values, making it obvious that the
                modes are physically equivalent.
    """
    has_a = snap_idx < len(A['snap_s_B']) and len(A['snap_s_B'][snap_idx])
    has_b = snap_idx < len(B['snap_s_B']) and len(B['snap_s_B'][snap_idx])

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8, 6.5), dpi=110, sharex=True,
        gridspec_kw={'height_ratios': [2.2, 1.0]},
    )

    # ---- Top: overlay --------------------------------------------------
    if has_a:
        s = s_to_mm(A['snap_s_B'][snap_idx])
        g = gap_to_mm(A['snap_gap_n'][snap_idx])
        # cpulc: solid line, filled circles
        ax_top.plot(s, g, '-', color=COLOR_CPULC, linewidth=1.4,
                    label='cpulc', zorder=2)
        ax_top.plot(s, g, 'o', color=COLOR_CPULC, markersize=4,
                    zorder=3)
    if has_b:
        s = s_to_mm(B['snap_s_B'][snap_idx])
        g = gap_to_mm(B['snap_gap_n'][snap_idx])
        # feeder: dashed line, OPEN markers, slightly transparent so
        # an underlying cpulc curve still shows through.
        ax_top.plot(s, g, '--', color=COLOR_FEEDER, linewidth=1.4,
                    alpha=0.85, label='feeder', zorder=4)
        ax_top.plot(s, g, 'o', mfc='none', mec=COLOR_FEEDER,
                    markersize=6, markeredgewidth=1.2, alpha=0.9,
                    zorder=5)

    ax_top.axhline(0.0, color='black', linewidth=0.5, alpha=0.5)
    pos_mm  = float(A['snap_t2_pos_m'][snap_idx]) * 1e3
    ang_deg = float(np.degrees(A['snap_t2_angle_rad'][snap_idx]))
    phase   = str(A['snap_phase'][snap_idx])
    ax_top.set_title(f"Snapshot {snap_idx}  [{phase}]  "
                     f"t2_translation={pos_mm:.2f} mm, "
                     f"t2_rotation={ang_deg:.2f} deg")
    ax_top.set_ylabel(r'normal gap  $\delta_n$  [mm]')
    ax_top.legend(loc='best')
    ax_top.grid(True, alpha=0.3)

    # ---- Bottom: |difference| on a log scale --------------------------
    # Only meaningful if both runs have the same number of contact pairs
    # at this snapshot AND the abscissae match (which they should -- same
    # geometry, same schedule).  If pair counts differ we just blank the
    # panel and annotate.
    if has_a and has_b:
        sa = np.asarray(A['snap_s_B'][snap_idx], dtype=float)
        ga = np.asarray(A['snap_gap_n'][snap_idx], dtype=float)
        sb = np.asarray(B['snap_s_B'][snap_idx], dtype=float)
        gb = np.asarray(B['snap_gap_n'][snap_idx], dtype=float)
        if sa.size == sb.size:
            diff_mm = np.abs(ga - gb) * 1e3                # m -> mm
            # Replace exact zeros with a tiny positive value so log scale
            # doesn't drop them.  Floor at the smallest positive nonzero
            # value among diff_mm, or 1e-18 mm if all are zero.
            nonzero = diff_mm[diff_mm > 0]
            floor = float(nonzero.min()) if nonzero.size else 1e-18
            diff_mm_plot = np.where(diff_mm > 0, diff_mm, floor)
            ax_bot.semilogy(s_to_mm(sa), diff_mm_plot, '-', color='#444444',
                            linewidth=1.2, marker='o', markersize=3)
            max_diff_mm = float(diff_mm.max())
            ax_bot.set_title(
                f"|cpulc - feeder|   max = {max_diff_mm:.3e} mm   "
                f"({'bit-identical' if max_diff_mm == 0 else 'fp-noise level' if max_diff_mm < 1e-9 else 'measurable'})",
                fontsize=10,
            )
        else:
            ax_bot.text(0.5, 0.5,
                        f"pair counts differ: cpulc={sa.size}, feeder={sb.size}",
                        ha='center', va='center', transform=ax_bot.transAxes)
    else:
        ax_bot.text(0.5, 0.5, 'one or both runs empty for this snapshot',
                    ha='center', va='center', transform=ax_bot.transAxes)

    ax_bot.set_xlabel(r'curvilinear abscissa  $s_B$  on Tube_3  [mm]')
    ax_bot.set_ylabel(r'$|\Delta \delta_n|$  [mm]')
    ax_bot.grid(True, which='both', alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / f"overlay_snap_{snap_idx:02d}.png")
    plt.close(fig)


# =============================================================================
def make_animation(out_dir: Path, A: dict, B: dict):
    """Animation scrubbing through all snapshots.  mp4 if ffmpeg present,
    else gif via Pillow."""
    n_snaps = min(len(A['snap_s_B']), len(B['snap_s_B']))
    if n_snaps == 0:
        print("[plot_comparison] no snapshots -- skipping animation.")
        return

    # X/Y bounds from union over all snapshots
    all_s, all_g = [], []
    for run in (A, B):
        for s, g in zip(run['snap_s_B'], run['snap_gap_n']):
            if len(s):
                all_s.append(np.asarray(s, dtype=float))
                all_g.append(np.asarray(g, dtype=float))
    if not all_s:
        print("[plot_comparison] all snapshots empty -- skipping animation.")
        return
    all_s = np.concatenate(all_s) * 1e3
    all_g = np.concatenate(all_g) * 1e3
    s_lo, s_hi = float(all_s.min()), float(all_s.max())
    g_lo, g_hi = float(all_g.min()), float(all_g.max())
    pad_g = 0.05 * (g_hi - g_lo + 1e-12)

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=110)
    # cpulc: solid + filled markers
    line_a, = ax.plot([], [], '-', color=COLOR_CPULC, linewidth=1.4,
                      label='cpulc', zorder=2)
    pts_a,  = ax.plot([], [], 'o', color=COLOR_CPULC, markersize=4,
                      zorder=3)
    # feeder: dashed + open markers, slightly transparent so cpulc shows through
    line_b, = ax.plot([], [], '--', color=COLOR_FEEDER, linewidth=1.4,
                      alpha=0.85, label='feeder', zorder=4)
    pts_b,  = ax.plot([], [], 'o', mfc='none', mec=COLOR_FEEDER,
                      markersize=6, markeredgewidth=1.2, alpha=0.9,
                      zorder=5)
    ax.axhline(0.0, color='black', linewidth=0.5, alpha=0.5)
    ax.set_xlim(s_lo, s_hi)
    ax.set_ylim(g_lo - pad_g, g_hi + pad_g)
    ax.set_xlabel(r'curvilinear abscissa  $s_B$  [mm]')
    ax.set_ylabel(r'normal gap  $\delta_n$  [mm]')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    title = ax.set_title('')

    def update(idx):
        for run, line, pts in [(A, line_a, pts_a), (B, line_b, pts_b)]:
            if idx < len(run['snap_s_B']):
                s = s_to_mm(run['snap_s_B'][idx])
                g = gap_to_mm(run['snap_gap_n'][idx])
                line.set_data(s, g)
                pts.set_data(s, g)
            else:
                line.set_data([], [])
                pts.set_data([], [])
        pos_mm = float(A['snap_t2_pos_m'][idx]) * 1e3
        ang_deg = float(np.degrees(A['snap_t2_angle_rad'][idx]))
        phase = str(A['snap_phase'][idx])
        title.set_text(f"Snap {idx}  [{phase}]  "
                       f"tx={pos_mm:.2f} mm, theta={ang_deg:.2f} deg")
        return line_a, pts_a, line_b, pts_b, title

    anim = FuncAnimation(fig, update, frames=n_snaps, interval=400, blit=False)

    mp4 = out_dir / 'overlay_animation.mp4'
    gif = out_dir / 'overlay_animation.gif'
    try:
        anim.save(mp4, writer=FFMpegWriter(fps=2.5, bitrate=2400))
        print(f"[plot_comparison] wrote {mp4}")
    except Exception as e:
        print(f"[plot_comparison] FFMpeg unavailable ({e}); falling back to GIF")
        try:
            anim.save(gif, writer=PillowWriter(fps=2))
            print(f"[plot_comparison] wrote {gif}")
        except Exception as e2:
            print(f"[plot_comparison] GIF write also failed: {e2}")
    plt.close(fig)


# =============================================================================
def _phase_info(run: dict) -> dict:
    """
    Reconstruct phase boundaries (in step-index space) from the schedule
    constants stored in the .npz, so we can shade the per-step plots and
    mark snapshot positions consistently across runs.

    Returns a dict with:
      init_end          step index where init phase ends
      trans_end         step index where translation phase ends
      rot_end           step index where rotation phase ends (== total)
      trans_cycle_len   length of one (actuate + hold) translation cycle
      rot_cycle_len     length of one (actuate + hold) rotation cycle
      snap_step_idx     step indices of all snapshots (length 1+N_t+N_r)
    """
    init   = int(run['sched_init_steps'])
    hold   = int(run['sched_hold_steps'])
    n_t    = int(run['sched_trans_n_increments'])
    sp_t   = int(run['sched_trans_steps_per_incr'])
    n_r    = int(run['sched_rot_n_increments'])
    sp_r   = int(run['sched_rot_steps_per_incr'])

    trans_cycle = sp_t + hold
    rot_cycle   = sp_r + hold
    init_end    = init
    trans_end   = init_end  + n_t * trans_cycle
    rot_end     = trans_end + n_r * rot_cycle

    return {
        'init_end':        init_end,
        'trans_end':       trans_end,
        'rot_end':         rot_end,
        'trans_cycle_len': trans_cycle,
        'rot_cycle_len':   rot_cycle,
        'snap_step_idx':   np.asarray(run['snap_step_idx'], dtype=np.int64),
    }


# Colors for phase shading -- light enough not to drown the data lines.
_PHASE_COLORS = {
    'init':  '#e8e8e8',   # light gray
    'trans': '#dceaf7',   # pale blue
    'rot':   '#fce5d4',   # pale orange
}


def _shade_phases(ax, run: dict, x_max: int | None = None,
                  draw_snapshots: bool = True):
    """
    Shade the init / trans / rot regions of the per-step x-axis on `ax`,
    and optionally draw thin vertical lines at every snapshot step.

    Call AFTER all data lines are plotted on `ax` so the shading sits
    underneath; we set zorder=0 explicitly.
    """
    info = _phase_info(run)
    if x_max is None:
        x_max = info['rot_end']

    # Three colored bands
    ax.axvspan(0,                    info['init_end'],  facecolor=_PHASE_COLORS['init'],
               alpha=0.55, zorder=0, lw=0)
    ax.axvspan(info['init_end'],     info['trans_end'], facecolor=_PHASE_COLORS['trans'],
               alpha=0.55, zorder=0, lw=0)
    ax.axvspan(info['trans_end'],    min(info['rot_end'], x_max),
               facecolor=_PHASE_COLORS['rot'], alpha=0.55, zorder=0, lw=0)

    # Snapshot tick lines
    if draw_snapshots:
        for s in info['snap_step_idx']:
            ax.axvline(int(s), color='black', alpha=0.18, linewidth=0.5,
                       zorder=0.5)


# =============================================================================
def plot_actuation_timeline(out_dir: Path, A: dict, B: dict):
    """
    Two-panel timeline of t2_pos (top) and t2_angle (bottom) vs simulation
    step, with phase shading.  Sanity check that both runs followed the
    same actuation schedule.  If the curves visibly disagree, something
    is wrong with the determinism between modes.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5.5), dpi=110, sharex=True)

    for run, color, label in [(A, COLOR_CPULC,  'cpulc'),
                              (B, COLOR_FEEDER, 'feeder')]:
        x = np.asarray(run['step_idx'], dtype=np.int64)
        ax1.plot(x, np.asarray(run['t2_pos_m']) * 1e3, '-',
                 color=color, label=label, linewidth=1.0)
        ax2.plot(x, np.degrees(np.asarray(run['t2_angle_rad'])), '-',
                 color=color, label=label, linewidth=1.0)

    _shade_phases(ax1, A)
    _shade_phases(ax2, A)

    ax1.set_ylabel('Tube_3 translation  [mm]')
    ax2.set_ylabel('Tube_3 rotation  [deg]')
    ax2.set_xlabel('simulation step')

    for ax in (ax1, ax2):
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)

    ax1.set_title("Actuation timeline (sanity check: both modes should match exactly)")
    fig.tight_layout()
    fig.savefig(out_dir / 'actuation_timeline.png')
    plt.close(fig)


# =============================================================================
def plot_scalar_per_snap(out_dir: Path, A: dict, B: dict, key: str,
                         ylabel: str, fname: str, scale: float = 1.0,
                         abs_val: bool = False):
    """Generic per-snapshot scalar overlay."""
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=110)
    for run, color, label in [(A, COLOR_CPULC,  'cpulc'),
                              (B, COLOR_FEEDER, 'feeder')]:
        if key == 'max_pen_per_snap':
            # Recompute per-snapshot max pen from gap arrays
            vals = []
            for g in run['snap_gap_n']:
                if len(g):
                    g = np.asarray(g, dtype=float)
                    vals.append(max(0.0, float(-g.min())) if g.min() < 0 else 0.0)
                else:
                    vals.append(0.0)
            y = np.asarray(vals) * scale
        elif key == 'n_active_per_snap':
            y = np.asarray([len(g) for g in run['snap_gap_n']], dtype=float)
        else:
            y = np.asarray(run[key], dtype=float) * scale
        if abs_val:
            y = np.abs(y)
        x = np.arange(len(y))
        ax.plot(x, y, '-o', color=color, label=label, markersize=4,
                linewidth=1.3)
    ax.set_xlabel('snapshot index')
    ax.set_ylabel(ylabel)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / fname)
    plt.close(fig)


def plot_scalar_per_step(out_dir: Path, A: dict, B: dict, key: str,
                         ylabel: str, fname: str,
                         rolling_window: int = 0,
                         shade_phases: bool = True,
                         log_y: bool = False):
    """Per-step scalar overlay (long arrays).  Optionally shades phase
    regions and marks snapshot step indices."""
    fig, ax = plt.subplots(figsize=(10, 4.5), dpi=110)
    for run, color, label in [(A, COLOR_CPULC,  'cpulc'),
                              (B, COLOR_FEEDER, 'feeder')]:
        y = np.asarray(run[key], dtype=float)
        x = np.arange(y.size)
        if rolling_window > 1 and y.size > rolling_window:
            # Rolling median to suppress per-step jitter
            kernel = rolling_window
            y_pad = np.pad(y, kernel // 2, mode='edge')
            y_smooth = np.array([
                np.nanmedian(y_pad[i:i + kernel]) for i in range(y.size)
            ])
            ax.plot(x, y_smooth, '-', color=color, label=f"{label} (median, w={kernel})",
                    linewidth=1.2, zorder=3)
        else:
            ax.plot(x, y, '-', color=color, label=label, linewidth=0.8,
                    alpha=0.85, zorder=3)
    if shade_phases:
        _shade_phases(ax, A)
    if log_y:
        ax.set_yscale('symlog', linthresh=1e-9)
    ax.set_xlabel('simulation step')
    ax.set_ylabel(ylabel)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / fname)
    plt.close(fig)


# =============================================================================
def write_summary(out_dir: Path, A: dict, B: dict):
    lines = ['# Comparison summary', '']
    for run, label in [(A, 'cpulc'), (B, 'feeder')]:
        n_iters = np.asarray(run['n_iters'], dtype=float)
        max_pen = np.asarray(run['max_pen'], dtype=float)
        wall    = np.asarray(run['wall_dt_s'], dtype=float)
        nact    = np.asarray(run['n_active'], dtype=float)

        def _safe(fn, x):
            x = x[np.isfinite(x)]
            return float(fn(x)) if x.size else float('nan')

        lines += [
            f"## mode = {label}",
            f"  total steps recorded   : {len(run['step_idx'])}",
            f"  total wall-clock       : {_safe(np.nansum, wall):.2f} s",
            f"  mean iters/step        : {_safe(np.nanmean, n_iters):.2f}",
            f"  median iters/step      : {_safe(np.nanmedian, n_iters):.2f}",
            f"  max iters/step         : {_safe(np.nanmax, n_iters):.2f}",
            f"  mean active pairs/step : {_safe(np.nanmean, nact):.2f}",
            f"  max |delta_n| over run : {_safe(np.nanmax, max_pen):.3e} m",
            f"  mean wall_dt/step      : {_safe(np.nanmean, wall)*1e3:.3f} ms",
            f"  median wall_dt/step    : {_safe(np.nanmedian, wall)*1e3:.3f} ms",
            '',
        ]
    (out_dir / 'summary.txt').write_text('\n'.join(lines))
    print(f"[plot_comparison] wrote {out_dir/'summary.txt'}")


# =============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cpulc',  default='comparison_cpulc.npz')
    ap.add_argument('--feeder', default='comparison_feeder.npz')
    ap.add_argument('--out-dir', default='plots')
    ap.add_argument('--rolling', type=int, default=200,
                    help="Rolling median window for per-step plots.")
    args = ap.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    A = load_run(Path(args.cpulc))
    B = load_run(Path(args.feeder))

    # ---- Per-snapshot overlay PNGs ----
    n_snaps = min(len(A['snap_s_B']), len(B['snap_s_B']))
    print(f"[plot_comparison] writing {n_snaps} snapshot overlays...")
    for k in range(n_snaps):
        plot_snapshot_overlay(out_dir, A, B, k)

    # ---- Animation scrubbing through snapshots ----
    make_animation(out_dir, A, B)

    # ---- Actuation timeline (sanity check: both runs should match) ----
    plot_actuation_timeline(out_dir, A, B)

    # ---- Scalar comparisons -- snapshot-based ----
    plot_scalar_per_snap(
        out_dir, A, B, 'max_pen_per_snap',
        r'max penetration $|\delta_n^-|$  [mm]', 'scalar_max_pen_per_snap.png',
        scale=1e3,
    )
    plot_scalar_per_snap(
        out_dir, A, B, 'n_active_per_snap',
        'active pair count', 'scalar_n_active_per_snap.png',
    )

    # ---- Scalar comparisons -- per-step (every step recorded) ----
    plot_scalar_per_step(
        out_dir, A, B, 'max_pen', r'max penetration $|\delta_n^-|$  [m]',
        'scalar_max_pen_per_step.png', rolling_window=0,
        log_y=True,
    )
    plot_scalar_per_step(
        out_dir, A, B, 'n_active', 'active pair count',
        'scalar_n_active_per_step.png', rolling_window=0,
    )
    plot_scalar_per_step(
        out_dir, A, B, 'n_iters', 'GS solver iterations',
        'scalar_iters.png', rolling_window=args.rolling,
    )
    plot_scalar_per_step(
        out_dir, A, B, 'residual', 'GS solver residual',
        'scalar_residual.png', rolling_window=args.rolling,
    )
    plot_scalar_per_step(
        out_dir, A, B, 'wall_dt_s', 'wall-clock per step  [s]',
        'scalar_walltime.png', rolling_window=args.rolling,
    )

    # ---- Summary ----
    write_summary(out_dir, A, B)
    print(f"[plot_comparison] all outputs in {out_dir}")


if __name__ == '__main__':
    main()