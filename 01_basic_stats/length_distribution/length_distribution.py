#!/usr/bin/env python3
"""
SV Length Distribution by Type
Input:  /data/liujt/data/cover2.vcf  (SURVIVOR merged, 177 samples)
Filter: SUPP >= 2, SVLEN < 1,000,000
Output: sv_length_kde.png, sv_type_count_bar.png  (same directory as script)
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import seaborn as sns
import pandas as pd

HERE     = os.path.dirname(os.path.abspath(__file__))
VCF      = "/data/liujt/data/cover2.vcf"
OUT_KDE  = os.path.join(HERE, "sv_length_kde.png")
OUT_BAR  = os.path.join(HERE, "sv_type_count_bar.png")
MIN_SUPP = 2

COLOR = {
    "DEL": "#34679a",
    "DUP": "#87b9d6",
    "INS": "#f5b785",
    "INV": "#e27861",
}
ORDER = ["DEL", "DUP", "INS", "INV"]

# ── 1. Parse ─────────────────────────────────────────────────────────────────
records = {"Variant Type": [], "Length": []}

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
        if t in COLOR and 0 < l < 1_000_000:
            records["Variant Type"].append(t)
            records["Length"].append(l)

df = pd.DataFrame(records)
df["log_Length"] = np.log10(df["Length"])
counts = df["Variant Type"].value_counts()
print("Counts:", counts.to_dict())

# ── 2. KDE plot ───────────────────────────────────────────────────────────────
g = sns.displot(
    data=df,
    x="log_Length",
    hue="Variant Type",
    hue_order=ORDER,
    kind="kde",
    fill=True,
    alpha=0.2,
    height=6,
    aspect=1.4,
    bw_adjust=1.0,
    palette=COLOR,
)

xtick_vals   = [1, 2, 3, 4, 5, 6]
xtick_labels = ["10", "100", "1K", "10K", "100K", "1M"]
plt.xticks(xtick_vals, xtick_labels)
plt.xlabel("SV Length (bp, log10 scale)", fontsize=12)
plt.title("Structural Variant Length Distribution by Type", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT_KDE, dpi=300)
print(f"Saved → {OUT_KDE}")
plt.close()

# ── 3. Bar chart (log10 y-axis) ───────────────────────────────────────────────
bar_vals   = [counts.get(t, 0) for t in ORDER]
bar_labels = [f"{t}s" for t in ORDER]
face_colors = [mcolors.to_rgba(COLOR[t], 0.2) for t in ORDER]
edge_colors = [COLOR[t] for t in ORDER]

fig, ax = plt.subplots(figsize=(5, 5))
bars = ax.bar(bar_labels, bar_vals, color=face_colors,
              edgecolor=edge_colors, linewidth=2, width=0.55)
ax.set_yscale("log")
ax.set_ylabel("SV Count (log10)", fontsize=12)
ax.set_xlabel("SV Type", fontsize=12)
ax.set_title("Counts of Structural Variant Types", fontsize=12, fontweight="bold", pad=10)

for bar, val in zip(bars, bar_vals):
    ax.text(bar.get_x() + bar.get_width() / 2,
            val * 1.6, f"{val:,}",
            ha="center", va="bottom", fontsize=11)

ax.set_ylim(top=max(bar_vals) * 10)
ax.yaxis.grid(False)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout(pad=2.0)
plt.savefig(OUT_BAR, dpi=300)
print(f"Saved → {OUT_BAR}")
plt.close()
