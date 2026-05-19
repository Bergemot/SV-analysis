# SV Analysis — *Atrijuglans hetaohei*

Structural variation (SV) analysis pipeline for 177 whole-genome resequencing samples of the walnut shoot moth (*Atrijuglans hetaohei*). SVs were called with Delly, Lumpy, and Manta, then merged with SURVIVOR.

---

## Dataset

| Item | Detail |
|---|---|
| Species | *Atrijuglans hetaohei* |
| Samples | 177 individuals |
| Reference genome | Feng et al., 2025 |
| SV callers | Delly · Lumpy · Manta |
| Merge tool | SURVIVOR 1.0.7 |
| Merged VCF | `cover2.vcf` |
| Filter | SUPP ≥ 2 (singletons excluded unless noted) |

**SV counts after filtering (SUPP ≥ 2)**

| Type | Count |
|---|---|
| DEL | 127,187 |
| DUP | 7,402 |
| INS | 2,557 |
| INV | 2,146 |
| **Total** | **139,292** |

---

## Repository structure

```
SV-analysis/
├── scripts/                          # Early-version scripts (reference)
│   ├── 01_growth_curve.py
│   ├── 02_variant_type_stats.py
│   ├── 03_length_distribution.py
│   └── 04_frequency_classification.py
│
└── 01_basic_stats/
    ├── growth_curve/
    │   └── growth_curve.py           # SV accumulation curve (100 permutations)
    ├── variant_type_stats/
    │   └── variant_type_stats.py     # Per-type count table + bar chart
    ├── length_distribution/
    │   └── length_distribution.py    # KDE length distribution + count bar chart
    └── frequency_classification/
        └── frequency_classification.py  # Shared/Major/Polymorphic/Singleton
```

---

## Analyses

### 01 · Basic statistics ✅
| Script | Output |
|---|---|
| `growth_curve.py` | `sv_growth_curve.png`, `sv_growth_mean_std.csv` |
| `variant_type_stats.py` | `sv_type_table.png`, `sv_type_barplot.png` |
| `length_distribution.py` | `sv_length_kde.png`, `sv_type_count_bar.png` |
| `frequency_classification.py` | `sv_category_stacked_bar.png`, `sv_category_pie.png` |

Frequency categories are defined as:

| Category | Criterion |
|---|---|
| Shared SV | Present in all 177 samples |
| Major SV | Present in ≥ 50% of samples |
| Polymorphic SV | Present in > 1 sample |
| Singleton SV | Present in exactly 1 sample |

### 02 · ANNOVAR annotation
Functional annotation of SVs against the reference genome.

### 03 · Genomic distribution
SVs classified into Exonic, Regulatory (UTR / upstream / downstream), and Noncoding (intronic / intergenic) regions.

### 04 · Repeat sequences & formation mechanisms
Overlap analysis with repeat elements; SVs classified into VNTR, NAHR, TE-mediated, and NHR formation mechanisms.

### 05 · Population structure
ADMIXTURE and PCA based on SV genotypes across four populations.

### 06 · Population differentiation (F~ST~)
Pairwise F~ST~ calculated among the four populations.

### 07 · SV density & nucleotide diversity (π)
Per-chromosome distribution of SV density and π.

### 08 · Selection sweep
Identification of selection signals and functional enrichment of candidate regions.

### 09 · GEA
Genotype–environment association analysis linking SV variation to climate variables.

---

## Color palette

```python
color_dict = {
    "DEL": "#34679a",
    "DUP": "#87b9d6",
    "INS": "#f5b785",
    "INV": "#e27861",
}
```

---

## Dependencies

```
python >= 3.8
numpy · pandas · matplotlib · seaborn · scipy
```
