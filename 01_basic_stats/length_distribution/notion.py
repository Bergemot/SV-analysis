import pysam
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# 你的 VCF 文件路径
vcf_file = "merged.vcf"

variants = {"Variant Type": [], "Length": []}

with pysam.VariantFile(vcf_file) as vcf:
    for record in vcf:
        svtype = record.info.get("SVTYPE", None)
        svlen = record.info.get("SVLEN", None)
        if isinstance(svlen, list):
            svlen = svlen[0]
        if svtype and svlen and svlen < 1000000:
            variants["Variant Type"].append(svtype)
            variants["Length"].append(abs(svlen))

df = pd.DataFrame(variants)

# 统一Variant Type格式，确保和color_dict对应
df['Variant Type'] = df['Variant Type'].map({
    "INS": "INS",
    "DEL": "DEL",
    "DUP": "DUP",
    "INV": "INV"
})

df['log_Length'] = np.log10(df['Length'].replace(0, np.nan))

color_dict = {
    "DEL": "#374f99",
    "DUP": "#6ea4c4",
    "INS": "#f0b169",
    "INV": "#df5b3e"
}

# 画分布图，线和填充同色，fill=True + alpha
g = sns.displot(
    data=df,
    x="log_Length",
    hue="Variant Type",
    kind="kde",
    fill=True,
    alpha=0.2,          # 透明度，填充色半透明
    height=6,
    aspect=1.3,
    bw_adjust=1.0,
    palette=color_dict
)
# g._legend.remove()

xtick_values = [0, 1, 2, 3, 4, 5, 6, 7, 8]
xtick_labels = ["1", "10", "100", "1K", "10K", "100K", "1M", "10M", "100M"]
plt.xticks(xtick_values, xtick_labels)

plt.xlabel("SV Length (bp, log10 scale)")
plt.title("Structural Variant Length Distribution by Type")
plt.tight_layout()
plt.savefig("sv_kde_displot.png", dpi=300)
plt.show()
