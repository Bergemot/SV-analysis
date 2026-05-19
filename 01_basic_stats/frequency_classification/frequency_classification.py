#!/usr/bin/env python3
"""
SV Frequency Classification: Shared / Major / Polymorphic / Singleton
Input:  /data/liujt/data/cover2.vcf  (SURVIVOR merged, 177 samples)
Output: sv_category_stacked_bar.png, sv_category_pie.png  (same directory as script)

NOTE: This script uses the COMPLETE unfiltered dataset (all SVs including SUPP=1
singletons) to provide a full overview of SV frequency distribution. The total SV
count (~176,000) is therefore higher than the SUPP>=2 filtered set (~139,000) used
in other analyses (variant_type_stats, length_distribution, growth_curve).
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict

HERE      = os.path.dirname(os.path.abspath(__file__))
VCF       = "/data/liujt/data/cover2.vcf"
OUT_BAR   = os.path.join(HERE, "sv_category_stacked_bar.png")
OUT_PIE   = os.path.join(HERE, "sv_category_pie.png")
N_SAMPLES = 177

CATEGORIES = ["Singleton SV", "Polymorphic SV", "Major SV", "Shared SV"]
COLORS = {
    "Shared SV":      "#34679a",
    "Major SV":       "#87b9d6",
    "Polymorphic SV": "#f5b785",
    "Singleton SV":   "#e27861",
}
SVTYPES = ["DEL", "DUP", "INS", "INV"]

# ── 1. Parse ─────────────────────────────────────────────────────────────────
counts = defaultdict(lambda: defaultdict(int))

with open(VCF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        info   = line.split("\t")[7]
        supp_m = re.search(r"SUPP=(\d+)", info)
        type_m = re.search(r"SVTYPE=(\w+)", info)
        if not supp_m or not type_m:
            continue
        supp   = int(supp_m.group(1))
        svtype = type_m.group(1)
        if svtype not in SVTYPES:
            continue
        if supp == 1:
            cat = "Singleton SV"
        elif supp == N_SAMPLES:
            cat = "Shared SV"
        elif supp >= N_SAMPLES / 2:
            cat = "Major SV"
        else:
            cat = "Polymorphic SV"
        counts[svtype][cat] += 1

# ── 2. Stacked bar ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))
bottoms = [0] * len(SVTYPES)
for cat in CATEGORIES:
    vals = [counts[t][cat] for t in SVTYPES]
    ax.bar(SVTYPES, vals, bottom=bottoms,
           color=COLORS[cat], edgecolor="white", linewidth=0.5,
           label=cat, width=0.6)
    bottoms = [b + v for b, v in zip(bottoms, vals)]

ax.set_xlabel("SVTYPE", fontsize=12)
ax.set_ylabel("SV Count", fontsize=12)
ax.set_title("SV Category Composition by SVTYPE", fontsize=14, fontweight="bold")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
legend = ax.legend(title="Category", fontsize=10, title_fontsize=10,
                   loc="upper right", framealpha=0.9)
legend.get_frame().set_edgecolor("grey")
plt.tight_layout()
plt.savefig(OUT_BAR, dpi=300)
print(f"Saved → {OUT_BAR}")
plt.close()

# ── 3. Pie chart ──────────────────────────────────────────────────────────────
totals    = {cat: sum(counts[t][cat] for t in SVTYPES) for cat in CATEGORIES}
total_all = sum(totals.values())

labels = [
    f"{cat}\n{totals[cat]:,} ({totals[cat]/total_all*100:.1f}%)"
    for cat in CATEGORIES
]
colors = [COLORS[cat] for cat in CATEGORIES]

fig2, ax2 = plt.subplots(figsize=(7, 7))
wedges, _ = ax2.pie(
    [totals[c] for c in CATEGORIES],
    colors=colors,
    startangle=90,
    counterclock=False,
    wedgeprops={"width": 0.6, "edgecolor": "white", "linewidth": 1.5},
)
for wedge, label in zip(wedges, labels):
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
plt.savefig(OUT_PIE, dpi=300, bbox_inches="tight")
print(f"Saved → {OUT_PIE}")
plt.close()

# ── 4. Summary ────────────────────────────────────────────────────────────────
print("\n=== Frequency Classification Summary ===")
for cat in CATEGORIES:
    n = totals[cat]
    print(f"{cat:<20} {n:>8,}  ({n/total_all*100:.1f}%)")
print(f"{'Total':<20} {total_all:>8,}")
