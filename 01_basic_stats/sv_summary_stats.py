#!/usr/bin/env python3
"""
SV summary statistics table: count, mean/total/2.5%/97.5% length per SV type.
Input:  /data/liujt/data/cover2.vcf  (all 176,001 SVs, no filter)
Output: sv_summary_table.tsv, sv_summary_table.png
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE    = os.path.dirname(os.path.abspath(__file__))
VCF     = "/data/liujt/data/cover2.vcf"
OUT_TSV = os.path.join(HERE, "sv_summary_table.tsv")
OUT_PNG = os.path.join(HERE, "sv_summary_table.png")
ORDER   = ["DEL", "DUP", "INS", "INV"]

# ── 1. Parse VCF ──────────────────────────────────────────────────────────────
lengths = {t: [] for t in ORDER}

with open(VCF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        info  = line.split("\t")[7]
        m_t   = re.search(r"SVTYPE=(\w+)", info)
        m_l   = re.search(r"SVLEN=(-?\d+)", info)
        if not m_t or not m_l:
            continue
        t = m_t.group(1)
        l = abs(int(m_l.group(1)))
        if t in lengths:
            lengths[t].append(l)

# ── 2. Compute stats ──────────────────────────────────────────────────────────
rows = []
for t in ORDER:
    arr = np.array(lengths[t])
    rows.append({
        "Type":             f"{t}s",
        "Count":            len(arr),
        "Mean_length(bp)":  round(float(arr.mean()), 1),
        "Total_length(bp)": int(arr.sum()),
        "P2.5_length(bp)":  round(float(np.percentile(arr, 2.5)), 1),
        "P97.5_length(bp)": round(float(np.percentile(arr, 97.5)), 1),
    })

# Add total row
all_lens = np.concatenate([lengths[t] for t in ORDER])
rows.append({
    "Type":             "All",
    "Count":            len(all_lens),
    "Mean_length(bp)":  round(float(all_lens.mean()), 1),
    "Total_length(bp)": int(all_lens.sum()),
    "P2.5_length(bp)":  round(float(np.percentile(all_lens, 2.5)), 1),
    "P97.5_length(bp)": round(float(np.percentile(all_lens, 97.5)), 1),
})

df = pd.DataFrame(rows)

# ── 3. Save TSV ───────────────────────────────────────────────────────────────
df.to_csv(OUT_TSV, sep="\t", index=False)
print(df.to_string(index=False))
print(f"\nSaved → {OUT_TSV}")

# ── 4. Save PNG table ─────────────────────────────────────────────────────────
col_labels = ["Type", "Count", "Mean Length (bp)", "Total Length (bp)", "2.5% Length (bp)", "97.5% Length (bp)"]
cell_text  = [
    [r["Type"],
     f"{r['Count']:,}",
     f"{r['Mean_length(bp)']:,}",
     f"{r['Total_length(bp)']:,}",
     f"{r['P2.5_length(bp)']:,}",
     f"{r['P97.5_length(bp)']:,}"]
    for _, r in df.iterrows()
]

fig, ax = plt.subplots(figsize=(11, 2.8))
ax.axis("off")
tbl = ax.table(cellText=cell_text, colLabels=col_labels, loc="center", cellLoc="center")
tbl.auto_set_font_size(False)
tbl.set_fontsize(11)
tbl.scale(1, 1.8)

for j in range(len(col_labels)):
    tbl[0, j].set_facecolor("#e8eaf0")
    tbl[0, j].set_text_props(weight="bold")

# Highlight total row
for j in range(len(col_labels)):
    tbl[len(rows), j].set_facecolor("#f5f5f5")
    tbl[len(rows), j].set_text_props(style="italic")

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
print(f"Saved → {OUT_PNG}")
