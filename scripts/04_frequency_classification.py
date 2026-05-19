#!/usr/bin/env python3
"""
SV Frequency Classification: Shared / Major / Polymorphic / Singleton
Uses ALL records from cover2.vcf (including SUPP=1 singletons) to show
the full landscape; SUPP>=2 records follow the delly.vcf classification.

Output:
  01_basic_stats/frequency_classification/reproduced/结构变异-类型统计.png   (stacked bar)
  01_basic_stats/frequency_classification/reproduced/结构变异-频率统计.png   (pie chart)
"""

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict

VCF       = "/data/liujt/data/cover2.vcf"
OUT_BAR   = "/data/liujt/SV/01_basic_stats/frequency_classification/reproduced/结构变异-类型统计.png"
OUT_PIE   = "/data/liujt/SV/01_basic_stats/frequency_classification/reproduced/结构变异-频率统计.png"

# Number of samples in VCF
N_SAMPLES = 177

CATEGORIES = ["Singleton SV", "Polymorphic SV", "Major SV", "Shared SV"]
COLORS = {
    "Singleton SV":   "#f4a8a8",
    "Polymorphic SV": "#f4dca8",
    "Major SV":       "#aec6cf",
    "Shared SV":      "#1f3b6e",
}
SVTYPES = ["DEL", "DUP", "INS", "INV"]

# ── 1. Parse all records ──────────────────────────────────────────────────────
counts = defaultdict(lambda: defaultdict(int))   # counts[svtype][category]

with open(VCF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        info = line.split("\t")[7]
        supp_m = re.search(r"SUPP=(\d+)", info)
        type_m = re.search(r"SVTYPE=(\w+)", info)
        if not supp_m or not type_m:
            continue
        supp   = int(supp_m.group(1))
        svtype = type_m.group(1)
        if svtype not in SVTYPES:
            continue

        n = N_SAMPLES
        if supp == 1:
            cat = "Singleton SV"
        elif supp == n:
            cat = "Shared SV"
        elif supp >= n / 2:
            cat = "Major SV"
        else:
            cat = "Polymorphic SV"

        counts[svtype][cat] += 1

# ── 2. Stacked bar chart ─────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

bottoms = [0] * len(SVTYPES)
for cat in CATEGORIES:
    vals = [counts[t][cat] for t in SVTYPES]
    ax.bar(SVTYPES, vals, bottom=bottoms,
           color=COLORS[cat], edgecolor="white", linewidth=0.5,
           label=cat, width=0.5)
    bottoms = [b + v for b, v in zip(bottoms, vals)]

ax.set_xlabel("SVTYPE", fontsize=12)
ax.set_ylabel("SV Count", fontsize=12)
ax.set_title("SV Category Composition by SVTYPE", fontsize=14, fontweight="bold")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
legend = ax.legend(title="Category", fontsize=9, title_fontsize=9,
                   loc="upper right", framealpha=0.9)
legend.get_frame().set_edgecolor("grey")
plt.tight_layout()
plt.savefig(OUT_BAR, dpi=150)
print(f"Saved → {OUT_BAR}")
plt.close()

# ── 3. Pie chart (overall) ────────────────────────────────────────────────────
totals = {cat: sum(counts[t][cat] for t in SVTYPES) for cat in CATEGORIES}
total_all = sum(totals.values())

labels = [f"{cat}\n{totals[cat]:,} ({totals[cat]/total_all*100:.1f}%)" for cat in CATEGORIES]
colors  = [COLORS[cat] for cat in CATEGORIES]

fig2, ax2 = plt.subplots(figsize=(7, 7))
wedges, _ = ax2.pie(
    [totals[c] for c in CATEGORIES],
    colors=colors,
    startangle=90,
    wedgeprops={"edgecolor": "white", "linewidth": 1.5},
)
# Add labels with leader lines
for i, (wedge, label) in enumerate(zip(wedges, labels)):
    angle = (wedge.theta1 + wedge.theta2) / 2
    x = 1.25 * np.cos(np.deg2rad(angle))
    y = 1.25 * np.sin(np.deg2rad(angle))
    ax2.annotate(label,
                 xy=(0.7 * np.cos(np.deg2rad(angle)),
                     0.7 * np.sin(np.deg2rad(angle))),
                 xytext=(x, y),
                 ha="center", va="center", fontsize=10,
                 arrowprops=dict(arrowstyle="-", color="grey", lw=0.8))

ax2.set_title("SV Category Distribution", fontsize=16, fontweight="bold", pad=20)
plt.tight_layout()
plt.savefig(OUT_PIE, dpi=150, bbox_inches="tight")
print(f"Saved → {OUT_PIE}")
plt.close()

# ── 4. Summary ────────────────────────────────────────────────────────────────
print("\n=== Frequency Classification Summary ===")
print(f"{'Category':<20} {'Count':>8} {'Pct':>6}")
for cat in CATEGORIES:
    n = totals[cat]
    print(f"{cat:<20} {n:>8,} {n/total_all*100:>5.1f}%")
print(f"{'Total':<20} {total_all:>8,}")
