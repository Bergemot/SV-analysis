#!/usr/bin/env python3
"""
SV summary statistics: count, avg/total/2.5%/97.5% length per SV type.
Uses the ANNOVAR-annotated set (SUPP>=2 filtered, ~139k SVs).

Input:
  SV_annovar.output2.variant_function  — annotation labels (col 0) + coordinates
  sv_type.txt                          — chr, start, ID, SVTYPE  (key: chr+start)
Output:
  sv_summary_table.tsv
"""

import numpy as np
import pandas as pd

VARFUNC  = "/data/liujt/SV/02_ANNOVAR_annotation/SV_annovar.output2.variant_function"
SVTYPE_F = "/data/liujt/SV/02_ANNOVAR_annotation/sv_type.txt"
OUT      = "/data/liujt/SV/01_basic_stats/sv_summary_table.tsv"

# ── genomic region classification ──────────────────────────────────────────
EXONIC     = {"exonic", "ncRNA_exonic", "splicing", "ncRNA_splicing"}
REGULATORY = {"UTR3", "UTR5", "upstream", "downstream",
              "upstream;downstream", "UTR5;UTR3"}

def classify(annot):
    if annot in EXONIC:
        return "Exonic"
    if annot in REGULATORY:
        return "Regulatory"
    return "Noncoding"

# ── load SVTYPE lookup (chr+start → svtype) ────────────────────────────────
svtype_map = {}
with open(SVTYPE_F) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        key = (parts[0], parts[1])   # (chr, start)  0-based in VCF/sv_type.txt
        svtype_map[key] = parts[3]   # SVTYPE

print(f"SVTYPE entries loaded: {len(svtype_map):,}")

# ── load annotation ────────────────────────────────────────────────────────
rows = []
missing = 0
with open(VARFUNC) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        annot  = parts[0]
        chrom  = parts[2]
        start  = parts[3]
        end    = int(parts[4])
        length = end - int(start)
        region = classify(annot)

        # Some entries match directly; others need start-1 (ANNOVAR 1-based shift)
        svtype = svtype_map.get((chrom, start)) or \
                 svtype_map.get((chrom, str(int(start) - 1)))
        if svtype is None:
            missing += 1
            continue

        rows.append({"svtype": svtype, "length": length, "region": region})

print(f"Records matched: {len(rows):,}  |  unmatched: {missing:,}")
df = pd.DataFrame(rows)

# ── compute stats per group ────────────────────────────────────────────────
def stats(sub):
    L   = sub["length"]
    n   = len(sub)
    reg = sub["region"].value_counts()
    exo = int(reg.get("Exonic", 0))
    rg  = int(reg.get("Regulatory", 0))
    nc  = int(reg.get("Noncoding", 0))
    return pd.Series({
        "Count":             n,
        "Mean_length(bp)":   round(L.mean(), 1),
        "Total_length(bp)":  int(L.sum()),
        "P2.5_length(bp)":   int(np.percentile(L, 2.5)),
        "P97.5_length(bp)":  int(np.percentile(L, 97.5)),
        "Exonic_n":          exo,
        "Exonic_%":          round(exo / n * 100, 1),
        "Regulatory_n":      rg,
        "Regulatory_%":      round(rg / n * 100, 1),
        "Noncoding_n":       nc,
        "Noncoding_%":       round(nc / n * 100, 1),
    })

types  = [t for t in ["DEL", "DUP", "INS", "INV"] if t in df["svtype"].values]
result = df.groupby("svtype").apply(stats).loc[types]
overall = stats(df).rename("All")
result  = pd.concat([result, overall.to_frame().T])
result.index.name = "SV_type"

result.to_csv(OUT, sep="\t")
print(result.to_string())
print(f"\nSaved → {OUT}")
