#!/usr/bin/env python3
"""FSI1 SU2 + CalculiX histories (histories/level{1,2,4}.csv) with the Featflow 7+0 values.
Usage: plot_histories.py OUT.png   (run from paper-data/fsi1_su2_calculix)"""
import csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FEATFLOW = {"ux_m": 2.270493e-5, "uy_m": 8.208773e-4, "drag_N_per_m": 14.29426, "lift_N_per_m": 0.7637460}
RUNS = [("level1", "1x", "#2a78d6", "-"), ("level2", "2x", "#eb6834", "--"), ("level4", "4x", "#1baf7a", "-.")]
PANELS = [("ux_m", "$u_x(A)$ (m)"), ("uy_m", "$u_y(A)$ (m)"), ("drag_N_per_m", "drag (N/m)"), ("lift_N_per_m", "lift (N/m)")]
data = {r: list(csv.DictReader(open(f"histories/{r}.csv"))) for r, *_ in RUNS}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#8a8a85", "axes.labelcolor": "#2b2b29",
                     "xtick.color": "#5c5c58", "ytick.color": "#5c5c58"})
fig, axs = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
for ax, (key, lab) in zip(axs.flat, PANELS):
    for r, name, col, ls in RUNS:
        ax.plot([float(x["time_s"]) for x in data[r]], [float(x[key]) for x in data[r]], color=col, ls=ls, lw=1.5, label=f"SU2 + CalculiX {name}")
    ax.axhline(FEATFLOW[key], color="#2b2b29", lw=1.0, ls=":", label="Featflow 7+0")
    ax.set_ylabel(lab); ax.grid(True, color="#e4e3df", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
for ax in axs[1]: ax.set_xlabel("t (s)")
# the coupled start (t = 2 s) jumps from the rigid-plate state; zoom on the approach to the steady value
for ax, (key, _) in zip(axs.flat, PANELS):
    vals = [float(x[key]) for r, *_ in RUNS for x in data[r] if float(x["time_s"]) > 8]
    lo, hi = min(vals + [FEATFLOW[key]]), max(vals + [FEATFLOW[key]]); pad = 0.15 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
axs[0, 0].legend(loc="lower right", frameon=False, fontsize=8)
fig.suptitle("Hron–Turek FSI1, SU2 (incompressible) + CalculiX via preCICE: coupled from t = 2 s (y-range zoomed on t > 8 s)", fontsize=9)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=150)
