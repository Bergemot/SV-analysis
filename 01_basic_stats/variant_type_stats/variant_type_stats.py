#!/usr/bin/env python3
"""
SV Type Statistics: table + bar chart
Input:  /data/liujt/data/cover2.vcf  (SURVIVOR merged, 177 samples)
Filter: SUPP >= 2
Output: sv_type_table.png, sv_type_barplot.png  (same directory as script)
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

HERE     = os.path.dirname(os.path.abspath(__file__))
VCF      = "/data/liujt/data/cover2.vcf"
OUT_TBL  = os.path.join(HERE, "sv_type_table.png")
OUT_BAR  = os.path.join(HERE, "sv_type_barplot.png")
MIN_SUPP = 2

COLOR = {
    "DEL": "#34679a",
    "DUP": "#87b9d6",
    "INS": "#f5b785",
    "INV": "#e27861",
}
ORDER = ["DEL", "DUP", "INS", "INV"]

# ── 1. Parse ─────────────────────────────────────────────────────────────────
lengths = {t: [] for t in ORDER}

with open(VCF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        info  = line.split("\t")[7]
        supp  = int(re.search(r"SUPP=(\d+)", info).group(1))
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

# ── 2. Compute stats ──────────────────────────────────────────────────────────
rows = []
for t in ORDER:
    arr = np.array(lengths[t])
    rows.append([
        f"{t}s",
        len(arr),
        round(float(arr.mean()), 1),
        int(arr.sum()),
        round(float(np.percentile(arr, 2.5)), 1),
        round(float(np.percentile(arr, 97.5)), 1),
    ])

for r in rows:
    print("\t".join(str(x) for x in r))

# ── 3. Table PNG ──────────────────────────────────────────────────────────────
col_labels = ["Type", "Count", "Mean Length (bp)", "Total Length (bp)", "2.5% Length", "97.5% Length"]
cell_text  = [
    [r[0], f"{r[1]:,}", f"{r[2]:,}", f"{r[3]:,}", f"{r[4]:,}", f"{r[5]:,}"]
    for r in rows
]

fig, ax = plt.subplots(figsize=(10, 2.4))
ax.axis("off")
tbl = ax.table(cellText=cell_text, colLabels=col_labels, loc="center", cellLoc="center")
tbl.auto_set_font_size(False)
tbl.set_fontsize(11)
tbl.scale(1, 1.7)
for j in range(len(col_labels)):
    tbl[0, j].set_facecolor("#e8eaf0")
    tbl[0, j].set_text_props(weight="bold")
for i, t in enumerate(ORDER):
    tbl[i + 1, 0].set_facecolor(mcolors.to_rgba(COLOR[t], 0.25))
plt.tight_layout()
plt.savefig(OUT_TBL, dpi=300, bbox_inches="tight")
print(f"Saved → {OUT_TBL}")
plt.close()

# ── 4. Bar chart (log10 y-axis) ───────────────────────────────────────────────
bar_vals   = [len(lengths[t]) for t in ORDER]
bar_labels = [f"{t}s" for t in ORDER]
face_colors = [mcolors.to_rgba(COLOR[t], 0.2) for t in ORDER]
edge_colors = [COLOR[t] for t in ORDER]

fig2, ax2 = plt.subplots(figsize=(5, 5))
bars = ax2.bar(bar_labels, bar_vals, color=face_colors,
               edgecolor=edge_colors, linewidth=2, width=0.55)
ax2.set_yscale("log")
ax2.set_ylabel("SV Count (log10)", fontsize=12)
ax2.set_xlabel("SV Type", fontsize=12)
ax2.set_title("Counts of Structural Variant Types", fontsize=12, fontweight="bold", pad=10)

for bar, val in zip(bars, bar_vals):
    ax2.text(bar.get_x() + bar.get_width() / 2,
             val * 1.6, f"{val:,}",
             ha="center", va="bottom", fontsize=11)

ax2.set_ylim(top=max(bar_vals) * 10)
ax2.yaxis.grid(True, which="both", linestyle="--", linewidth=0.4, alpha=0.5)
ax2.set_axisbelow(True)
ax2.spines[["top", "right"]].set_visible(False)
plt.tight_layout(pad=2.0)
plt.savefig(OUT_BAR, dpi=300)
print(f"Saved → {OUT_BAR}")
plt.close()
