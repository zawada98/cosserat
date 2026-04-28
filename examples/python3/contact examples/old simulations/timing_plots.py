# -*- coding: utf-8 -*-
"""
plot_detection_times.py
=======================
Standalone script — run independently from any SOFA scene.

Reads two timing log files produced by TimingLogger and plots both
detection-time series on the same figure for direct comparison.

Usage
-----
    python plot_detection_times.py

Output
------
    detection_time_comparison.png   (saved next to this script)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# ── File paths ────────────────────────────────────────────────────────────────
DATA_DIR   = r"/"

CLASSIC_FILE = os.path.join(DATA_DIR, "classic_pipeline_detection_times.txt")
SSIM_FILE    = os.path.join(DATA_DIR, "ssim_pipeline_detection_times.txt")

OUTPUT_FILE  = os.path.join(DATA_DIR, "detection_time_comparison.png")


# ── Parser ────────────────────────────────────────────────────────────────────
def load_timing_file(path: str):
    """
    Parse a TimingLogger .txt file.

    Returns
    -------
    steps : np.ndarray  (int)
    times : np.ndarray  (float, milliseconds)
    """
    steps, times = [], []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("step"):          # header row
                continue
            parts = line.split(",")
            if len(parts) == 2:
                steps.append(int(parts[0]))
                times.append(float(parts[1]))
    return np.array(steps), np.array(times)


# ── Load data ─────────────────────────────────────────────────────────────────
steps_c, times_c = load_timing_file(CLASSIC_FILE)
steps_s, times_s = load_timing_file(SSIM_FILE)

# ── Rolling mean helper ───────────────────────────────────────────────────────
def rolling_mean(arr: np.ndarray, window: int = 20) -> np.ndarray:
    """Simple centred rolling mean; edges use available data only."""
    out = np.empty_like(arr)
    half = window // 2
    for i in range(len(arr)):
        lo = max(0, i - half)
        hi = min(len(arr), i + half + 1)
        out[i] = arr[lo:hi].mean()
    return out


rm_c = rolling_mean(times_c, window=20)
rm_s = rolling_mean(times_s, window=20)

# ── Summary stats ─────────────────────────────────────────────────────────────
def stats(arr):
    return arr.mean(), arr.min(), arr.max(), arr.std()

mean_c, min_c, max_c, std_c = stats(times_c)
mean_s, min_s, max_s, std_s = stats(times_s)

# ── Colours & style ───────────────────────────────────────────────────────────
CLR_CLASSIC = "#E05C2A"          # burnt orange
CLR_SSIM    = "#2A7AE0"          # steel blue
ALPHA_RAW   = 0.25
ALPHA_FILL  = 0.10

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.grid":         True,
    "grid.color":        "#DDDDDD",
    "grid.linewidth":    0.6,
    "grid.linestyle":    "--",
    "figure.facecolor":  "white",
    "axes.facecolor":    "#FAFAFA",
})

# ── Figure layout: 2 rows ─────────────────────────────────────────────────────
#   Row 0 : main time-series comparison
#   Row 1 : two side-by-side histograms
fig, axes = plt.subplots(
    2, 2,
    figsize=(14, 8),
    gridspec_kw={"height_ratios": [2.5, 1], "hspace": 0.45, "wspace": 0.3},
)
ax_main = plt.subplot2grid((2, 2), (0, 0), colspan=2, fig=fig)
ax_hist_c = axes[1, 0]
ax_hist_s = axes[1, 1]

# Remove the auto-created top-row axes (replaced by ax_main spanning both cols)
axes[0, 0].remove()
axes[0, 1].remove()

# ── Row 0: main time-series ───────────────────────────────────────────────────
# Raw scatter (very light)
ax_main.scatter(steps_c, times_c, color=CLR_CLASSIC, s=4,
                alpha=ALPHA_RAW, zorder=2, label="_nolegend_")
ax_main.scatter(steps_s, times_s, color=CLR_SSIM,    s=4,
                alpha=ALPHA_RAW, zorder=2, label="_nolegend_")

# ±1 std fill around rolling mean
ax_main.fill_between(steps_c,
                     np.clip(rm_c - std_c, 0, None), rm_c + std_c,
                     color=CLR_CLASSIC, alpha=ALPHA_FILL, zorder=1)
ax_main.fill_between(steps_s,
                     np.clip(rm_s - std_s, 0, None), rm_s + std_s,
                     color=CLR_SSIM, alpha=ALPHA_FILL, zorder=1)

# Rolling mean lines
ax_main.plot(steps_c, rm_c, color=CLR_CLASSIC, lw=2.0, zorder=3,
             label=f"Classic pipeline  (mean {mean_c:.3f} ms)")
ax_main.plot(steps_s, rm_s, color=CLR_SSIM,    lw=2.0, zorder=3,
             label=f"SSIM pipeline     (mean {mean_s:.3f} ms)")

# Mean horizontal dashed lines
ax_main.axhline(mean_c, color=CLR_CLASSIC, lw=1.0, ls=":", alpha=0.8)
ax_main.axhline(mean_s, color=CLR_SSIM,    lw=1.0, ls=":", alpha=0.8)

ax_main.set_xlabel("Simulation step", fontsize=12)
ax_main.set_ylabel("Detection time (ms)", fontsize=12)
ax_main.set_title("Contact detection time: Classic pipeline vs. SSIM pipeline",
                  fontsize=13, fontweight="bold", pad=12)
ax_main.set_xlim(left=1)
ax_main.set_ylim(bottom=0)
ax_main.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
ax_main.tick_params(axis="both", which="both", direction="in")

# Custom legend (raw dots + mean line explained)
legend_handles = [
    Line2D([0], [0], color=CLR_CLASSIC, lw=2.5,
           label=f"Classic  — mean {mean_c:.3f} ms, max {max_c:.3f} ms"),
    Line2D([0], [0], color=CLR_SSIM,    lw=2.5,
           label=f"SSIM     — mean {mean_s:.3f} ms, max {max_s:.3f} ms"),
    Line2D([0], [0], color="grey", lw=1.0, ls="--", alpha=0.5,
           label="20-step rolling mean  (shaded band = ±1 std)"),
]
ax_main.legend(handles=legend_handles, framealpha=0.92,
               loc="upper right", fontsize=10)

# Speedup annotation
speedup = mean_c / mean_s if mean_s > 0 else float("inf")
ax_main.annotate(
    f"SSIM is {speedup:.1f}× faster on average",
    xy=(0.02, 0.93), xycoords="axes fraction",
    fontsize=11, color="#333333",
    bbox=dict(boxstyle="round,pad=0.35", fc="white",
              ec="#CCCCCC", lw=0.8, alpha=0.9),
)

# ── Row 1: histograms ─────────────────────────────────────────────────────────
bins = np.linspace(0, max(times_c.max(), times_s.max()) * 1.05, 45)

for ax, times, clr, label, mean_v, std_v in [
    (ax_hist_c, times_c, CLR_CLASSIC, "Classic pipeline", mean_c, std_c),
    (ax_hist_s, times_s, CLR_SSIM,    "SSIM pipeline",    mean_s, std_s),
]:
    ax.hist(times, bins=bins, color=clr, alpha=0.75, edgecolor="white",
            linewidth=0.4)
    ax.axvline(mean_v, color=clr, lw=1.8, ls="--",
               label=f"mean = {mean_v:.3f} ms")
    ax.axvline(mean_v + std_v, color=clr, lw=1.0, ls=":",
               label=f"±1σ  = {std_v:.3f} ms")
    ax.axvline(max(0, mean_v - std_v), color=clr, lw=1.0, ls=":")
    ax.set_xlabel("Detection time (ms)", fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title(label, fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.9)
    ax.tick_params(axis="both", direction="in")
    ax.set_xlim(left=0)

# ── Save ──────────────────────────────────────────────────────────────────────
fig.tight_layout(rect=[0, 0, 1, 1])
plt.savefig(OUTPUT_FILE, dpi=150, bbox_inches="tight")
print(f"Plot saved → {OUTPUT_FILE}")
plt.show()