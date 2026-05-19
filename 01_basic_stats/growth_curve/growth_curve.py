#!/usr/bin/env python3
"""
SV Discovery Growth Curve
Input:  /data/liujt/data/cover2.vcf  (SURVIVOR merged, 177 samples)
Filter: none (all 176,001 SVs)
Output: sv_growth_curve.png, sv_growth_mean_std.csv  (same directory as script)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE     = os.path.dirname(os.path.abspath(__file__))
VCF      = "/data/liujt/data/cover2.vcf"
OUT_PNG  = os.path.join(HERE, "sv_growth_curve.png")
OUT_CSV  = os.path.join(HERE, "sv_growth_mean_std.csv")
N_PERM   = 100

# ── 1. Parse SUPP_VEC ────────────────────────────────────────────────────────
print("Parsing VCF ...")
supp_vecs = []
n_samples = None

with open(VCF) as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            n_samples = len(line.strip().split("\t")) - 9
            continue
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        info   = fields[7]
        vec_str = next(x for x in info.split(";") if x.startswith("SUPP_VEC=")).split("=")[1]
        supp_vecs.append([int(c) for c in vec_str])

mat = np.array(supp_vecs, dtype=np.uint8)
n_sv = len(supp_vecs)
print(f"  {n_sv} SVs (no filter), {n_samples} samples")

# ── 2. Rarefaction ───────────────────────────────────────────────────────────
print(f"Running {N_PERM} permutations ...")
x      = np.arange(1, n_samples + 1)
counts = np.zeros((N_PERM, n_samples), dtype=np.int32)
rng    = np.random.default_rng(42)

for p in range(N_PERM):
    order    = rng.permutation(n_samples)
    observed = np.zeros(n_sv, dtype=bool)
    for k, idx in enumerate(order):
        observed |= mat[:, idx].astype(bool)
        counts[p, k] = observed.sum()
    if (p + 1) % 10 == 0:
        print(f"  perm {p+1}/{N_PERM}", flush=True)

mean_sv = counts.mean(axis=0)
std_sv  = counts.std(axis=0)

# ── 3. Save CSV ──────────────────────────────────────────────────────────────
pd.DataFrame({
    "sample_number":  x,
    "mean_sv_count":  mean_sv,
    "std_sv_count":   std_sv,
}).to_csv(OUT_CSV, index=False)
print(f"Saved → {OUT_CSV}")

# ── 4. Plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
ax.fill_between(x, mean_sv - std_sv, mean_sv + std_sv,
                color="#87b9d6", alpha=0.4, label="±1 SD")
ax.plot(x, mean_sv, color="#34679a", linewidth=2, label="Mean SV count")

ax.set_xlabel("Number of Samples", fontsize=12)
ax.set_ylabel("Number of Non-redundant SVs", fontsize=12)
ax.set_title("Growth Curve of Structural Variations (SVs)", fontsize=14, fontweight="bold")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
ax.grid(axis="y", linestyle="--", linewidth=0.5, color="grey", alpha=0.5)
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
print(f"Saved → {OUT_PNG}")
