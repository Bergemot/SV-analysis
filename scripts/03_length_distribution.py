#!/usr/bin/env python3
"""
SV Length Distribution by Type
Input:  /data/liujt/data/cover2.vcf
Filter: SUPP >= 2
Output:
  01_basic_stats/length_distribution/reproduced/结构变异-类型长度分布.png   (KDE plot)
  01_basic_stats/length_distribution/reproduced/结构变异-类型数量.png        (bar chart)
"""

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde

VCF      = "/data/liujt/data/cover2.vcf"
OUT_KDE  = "/data/liujt/SV/01_basic_stats/length_distribution/reproduced/结构变异-类型长度分布.png"
OUT_BAR  = "/data/liujt/SV/01_basic_stats/length_distribution/reproduced/结构变异-类型数量.png"
MIN_SUPP = 2

FILL_C = {"DEL": "#aec6cf", "INV": "#f4a8a8", "DUP": "#b3ddf0", "INS": "#f4dca8"}
LINE_C = {"DEL": "#1f3b6e", "INV": "#8b1a1a", "DUP": "#1a6b8b", "INS": "#8b6a1a"}
ORDER  = ["DEL", "INV", "DUP", "INS"]

# ── 1. Parse ─────────────────────────────────────────────────────────────────
lengths = {t: [] for t in ORDER}

with open(VCF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        info = line.split("\t")[7]
        supp = int(re.search(r"SUPP=(\d+)", info).group(1))
        if supp < MIN_SUPP:
            continue
        m_t = re.search(r"SVTYPE=(\w+)", info)
        m_l = re.search(r"SVLEN=(-?\d+)", info)
        if not m_t or not m_l:
            continue
        t = m_t.group(1)
        l = abs(int(m_l.group(1)))
        if t in lengths and l > 0:
            lengths[t].append(l)

counts = {t: len(v) for t, v in lengths.items()}
total  = sum(counts.values())
print("Counts:", counts, "Total:", total)

# ── 2. KDE density plot ───────────────────────────────────────────────────────
x_grid  = np.logspace(1, 8, 3000)
x_log10 = np.log10(x_grid)
bw_fixed = 0.20

fig, ax = plt.subplots(figsize=(10, 6))
for svtype in ORDER:
    arr = np.log10(np.array(lengths[svtype]))
    kde = gaussian_kde(arr, bw_method=bw_fixed)
    y   = kde(x_log10) * (counts[svtype] / total)
    ax.fill_between(x_grid, y, alpha=0.35, color=FILL_C[svtype])
    ax.plot(x_grid, y, color=LINE_C[svtype], linewidth=1.0)

ax.set_xscale("log")
ax.set_xlim(1, 1e8)
ax.set_xlabel("SV Length (bp, log10 scale)", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title("Structural Variant Length Distribution by Type",
             fontsize=14, fontweight="bold")
legend_patches = [
    mpatches.Patch(fc=FILL_C[t], ec=LINE_C[t], lw=1.5, label=t)
    for t in ORDER
]
ax.legend(handles=legend_patches, fontsize=11, loc="upper right")
plt.tight_layout()
plt.savefig(OUT_KDE, dpi=150)
print(f"Saved → {OUT_KDE}")
plt.close()

# ── 3. Bar chart of counts (log10 y-axis) ────────────────────────────────────
bar_order  = ["DEL", "DUP", "INS", "INV"]
bar_labels = [f"{t}s" for t in bar_order]
bar_vals   = [counts[t] for t in bar_order]

fig2, ax2 = plt.subplots(figsize=(5, 5))
bars = ax2.bar(
    bar_labels, bar_vals,
    color=[FILL_C[t] for t in bar_order],
    edgecolor=[LINE_C[t] for t in bar_order],
    linewidth=1.5, width=0.55,
)
ax2.set_yscale("log")
ax2.set_ylabel("SV Count (log10)", fontsize=11)
ax2.set_xlabel("SV Type", fontsize=11)
ax2.set_title("Counts of Structural Variant Types",
              fontsize=12, fontweight="bold")

for bar, val in zip(bars, bar_vals):
    ax2.text(bar.get_x() + bar.get_width() / 2,
             bar.get_height() * 1.4,
             f"{val:,}", ha="center", va="bottom", fontsize=9.5)

ax2.set_ylim(top=max(bar_vals) * 8)
ax2.yaxis.grid(True, which="both", linestyle="-", linewidth=0.4, alpha=0.5)
ax2.set_axisbelow(True)
plt.tight_layout()
plt.savefig(OUT_BAR, dpi=150)
print(f"Saved → {OUT_BAR}")
plt.close()
