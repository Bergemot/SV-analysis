#!/usr/bin/env python3
"""
SV Type Statistics Table
Input:  /data/liujt/data/cover2.vcf  (SURVIVOR merged, 177 samples)
Filter: SUPP >= 2
Output: 01_basic_stats/variant_type_stats/reproduced/count.png
"""

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VCF     = "/data/liujt/data/cover2.vcf"
OUT     = "/data/liujt/SV/01_basic_stats/variant_type_stats/reproduced/count.png"
MIN_SUPP = 2

# ── 1. Parse ─────────────────────────────────────────────────────────────────
lengths = {"DEL": [], "DUP": [], "INS": [], "INV": []}

with open(VCF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        info = line.split("\t")[7]
        supp = int(re.search(r"SUPP=(\d+)", info).group(1))
        if supp < MIN_SUPP:
            continue
        svtype_m = re.search(r"SVTYPE=(\w+)", info)
        svlen_m  = re.search(r"SVLEN=(-?\d+)", info)
        if not svtype_m or not svlen_m:
            continue
        svtype = svtype_m.group(1)
        svlen  = abs(int(svlen_m.group(1)))
        if svtype in lengths:
            lengths[svtype].append(svlen)

# ── 2. Compute stats ──────────────────────────────────────────────────────────
ORDER = ["INS", "DEL", "DUP", "INV"]
rows  = []
for t in ORDER:
    arr = np.array(lengths[t])
    rows.append([
        f"{t}s",
        len(arr),
        round(arr.mean(), 5),
        int(arr.sum()),
        round(float(np.percentile(arr, 2.5)), 3),
        round(float(np.percentile(arr, 97.5)), 3),
    ])

for r in rows:
    print("\t".join(str(x) for x in r))

# ── 3. Render as table PNG ────────────────────────────────────────────────────
col_labels = ["变异类型", "数量", "平均长度", "总长度", "2.5% 长度", "97.5% 长度"]
cell_text  = [
    [str(r[0]), f"{r[1]:,}", f"{r[2]}", f"{r[3]:,}", f"{r[4]}", f"{r[5]}"]
    for r in rows
]

fig, ax = plt.subplots(figsize=(9, 2.2))
ax.axis("off")
tbl = ax.table(
    cellText=cell_text,
    colLabels=col_labels,
    loc="center",
    cellLoc="center",
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(11)
tbl.scale(1, 1.6)

# Style header
for j in range(len(col_labels)):
    tbl[0, j].set_facecolor("#f0f0f0")
    tbl[0, j].set_text_props(weight="bold")

plt.tight_layout()
plt.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"Saved → {OUT}")
