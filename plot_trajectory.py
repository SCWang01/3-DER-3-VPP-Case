# plot_trajectory.py — Python equivalent of Plot_Profit_Taking_Ratio.m
# Reads a CSV produced by Sub_VPP.py (CLI mode) and plots profit-taking ratios.
#
# Usage:
#   python plot_trajectory.py gs_trajectory_eps0.1.csv

import sys
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams["font.family"] = "Calibri"
plt.rcParams["axes.unicode_minus"] = False

# ── Load data ──────────────────────────────────────────────────────────────
csv_path = sys.argv[1] if len(sys.argv) > 1 else "gs_trajectory_eps0.1.csv"

iters, x1s, x2s, x3s = [], [], [], []
with open(csv_path, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        iters.append(int(row["iter"]))
        x1s.append(float(row["x_1"]))
        x2s.append(float(row["x_2"]))
        x3s.append(float(row["x_3"]))

# ── Plot ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 4))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

ax.plot(iters, x1s, "-o", linewidth=2.2, markersize=5, label="VPP 1")
ax.plot(iters, x2s, "-s", linewidth=2.2, markersize=5, label="VPP 2")
ax.plot(iters, x3s, "-^", linewidth=2.2, markersize=5, label="VPP 3")

ax.set_xlabel("Iteration count", fontsize=13)
ax.set_ylabel("Profit-taking ratio", fontsize=13)
ax.tick_params(labelsize=11)
ax.set_ylim(-0.05, 1.05)
ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
ax.legend(fontsize=12, loc="best")
ax.grid(True, linestyle="--", alpha=0.4)
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
out_path = csv_path.replace(".csv", ".png")
fig.savefig(out_path, dpi=300)
print(f"Saved -> {out_path}")
