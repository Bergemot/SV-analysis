#!/usr/bin/env python3
"""
SV Discovery Growth Curve
Input:  /data/liujt/data/cover2.vcf  (SURVIVOR merged, all 177 samples)
Filter: SUPP >= 2
Output: 01_basic_stats/growth_curve/reproduced/结构变异-增长曲线.png
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict

VCF = "/data/liujt/data/cover2.vcf"
OUT = "/data/liujt/SV/01_basic_stats/growth_curve/reproduced/结构变异-增长曲线.png"
N_PERM = 100
MIN_SUPP = 2

# ── 1. Parse SUPP_VEC from VCF ──────────────────────────────────────────────
print("Parsing VCF …")
supp_vecs = []   # list of list[int], one per passing SV
n_samples = None

with open(VCF) as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            n_samples = len(line.strip().split("\t")) - 9
            continue
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        info = fields[7]
        # SUPP
        supp = int([x for x in info.split(";") if x.startswith("SUPP=")][0].split("=")[1])
        if supp < MIN_SUPP:
            continue
        # SUPP_VEC
        vec_str = [x for x in info.split(";") if x.startswith("SUPP_VEC=")][0].split("=")[1]
        supp_vecs.append([int(c) for c in vec_str])

n_sv   = len(supp_vecs)
print(f"  {n_sv} SVs (SUPP>={MIN_SUPP}), {n_samples} samples")

# Convert to numpy array: shape (n_sv, n_samples), dtype uint8
mat = np.array(supp_vecs, dtype=np.uint8)   # 1 = sample has this SV

# ── 2. Rarefaction ───────────────────────────────────────────────────────────
print(f"Running {N_PERM} permutations …")
sample_sizes = list(range(1, n_samples + 1))
counts = np.zeros((N_PERM, n_samples), dtype=np.int32)

rng = np.random.default_rng(42)
for p in range(N_PERM):
    order = rng.permutation(n_samples)
    observed = np.zeros(n_sv, dtype=bool)
    for k, idx in enumerate(order):
        observed |= mat[:, idx].astype(bool)
        counts[p, k] = observed.sum()
    if (p + 1) % 10 == 0:
        print(f"  perm {p+1}/{N_PERM}", flush=True)

mean_sv = counts.mean(axis=0)
std_sv  = counts.std(axis=0)

# ── 3. Plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))

ax.fill_between(sample_sizes, mean_sv - std_sv, mean_sv + std_sv,
                color="#aec6cf", alpha=0.5, label="±1 SD")
ax.plot(sample_sizes, mean_sv, color="#1f3b6e", linewidth=1.8, label="Mean SV count")

ax.set_xlabel("Number of Samples", fontsize=12)
ax.set_ylabel("Number of Non-redundant SVs", fontsize=12)
ax.set_title("Growth Curve of Structural Variations (SVs)", fontsize=14, fontweight="bold")
ax.legend(fontsize=10)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax.grid(axis="y", linestyle="--", linewidth=0.5, color="grey", alpha=0.5)

plt.tight_layout()
plt.savefig(OUT, dpi=150)
print(f"Saved → {OUT}")
