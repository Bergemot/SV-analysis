#!/usr/bin/env python3
"""
SV Genomic Distribution Statistics
Classifies SVs into three categories based on ANNOVAR annotation:
  - Exonic:     exonic, ncRNA_exonic, splicing, ncRNA_splicing
  - Regulatory: UTR3, UTR5, upstream, downstream, upstream;downstream, UTR5;UTR3
  - Noncoding:  intronic, intergenic, ncRNA_intronic, and all others

Input:
  SV_annovar.output2.variant_function  (ANNOVAR output)
  sv_type.txt                          (chr, start, ID, SVTYPE)
Output:
  sv_genomic_distribution.tsv
"""

import os
import numpy as np
import pandas as pd

HERE     = os.path.dirname(os.path.abspath(__file__))
VARFUNC  = "/data/liujt/SV/02_ANNOVAR_annotation/SV_annovar.output2.variant_function"
SVTYPE_F = "/data/liujt/SV/02_ANNOVAR_annotation/sv_type.txt"
OUT_TSV  = os.path.join(HERE, "sv_genomic_distribution.tsv")

EXONIC     = {"exonic", "ncRNA_exonic", "splicing", "ncRNA_splicing"}
REGULATORY = {"UTR3", "UTR5", "upstream", "downstream",
              "upstream;downstream", "UTR5;UTR3"}

def classify(annot):
    if annot in EXONIC:
        return "Exonic"
    if annot in REGULATORY:
        return "Regulatory"
    return "Noncoding"

# ── 1. Load SVTYPE lookup ─────────────────────────────────────────────────────
svtype_map = {}
with open(SVTYPE_F) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        svtype_map[(p[0], p[1])] = p[3]

# ── 2. Parse annotation ───────────────────────────────────────────────────────
rows = []
with open(VARFUNC) as fh:
    for line in fh:
        p      = line.rstrip("\n").split("\t")
        annot  = p[0]
        chrom  = p[2]
        start  = p[3]
        region = classify(annot)
        svtype = (svtype_map.get((chrom, start)) or
                  svtype_map.get((chrom, str(int(start) - 1))))
        if svtype is None:
            continue
        rows.append({"svtype": svtype, "region": region})

df = pd.DataFrame(rows)
print(f"Total records: {len(df):,}")

# ── 3. Compute stats ──────────────────────────────────────────────────────────
ORDER      = ["DEL", "DUP", "INS", "INV"]
CATEGORIES = ["Exonic", "Regulatory", "Noncoding"]
results    = []

for svtype in ORDER + ["All"]:
    sub = df if svtype == "All" else df[df["svtype"] == svtype]
    n   = len(sub)
    row = {"SV_type": svtype, "Total": n}
    for cat in CATEGORIES:
        cnt = (sub["region"] == cat).sum()
        row[f"{cat}_n"]   = int(cnt)
        row[f"{cat}_%"]   = round(cnt / n * 100, 1) if n > 0 else 0.0
    results.append(row)

result_df = pd.DataFrame(results)
col_order = ["SV_type", "Total",
             "Exonic_n", "Exonic_%",
             "Regulatory_n", "Regulatory_%",
             "Noncoding_n", "Noncoding_%"]
result_df = result_df[col_order]

result_df.to_csv(OUT_TSV, sep="\t", index=False)
print(result_df.to_string(index=False))
print(f"\nSaved → {OUT_TSV}")
