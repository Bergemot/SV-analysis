<aside>
💡

色卡：

```python
color_dict = {
    "DEL": "#34679a",
    "DUP": "#87b9d6",
    "INS": "#f5b785",
    "INV": "#e27861"
}
```

</aside>

### 1. 核密度长度分布图

- python length.py
    
    ```python
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
    
    ```
    

### 2. 数量柱状图

- python count.py
    
    ```python
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # 创建数据
    data = {
        "Type": ["DELs", "DUPs", "INSs", "INVs"],
        "Count": [126713, 7252, 3283, 2592]
    }
    df = pd.DataFrame(data)
    
    # 统一颜色（Hex 颜色码）
    color_dict = {
        "DELs": "#374f99",
        "DUPs": "#6ea4c4",
        "INSs": "#f0b169",
        "INVs": "#df5b3e"
    }
    palette = [color_dict[svtype] for svtype in df["Type"]]
    
    # 设置主题
    sns.set_theme(style="whitegrid")
    
    # 画图
    plt.figure(figsize=(5,5))
    ax = sns.barplot(data=df, x="Type", y="Count", fill=True)
    
    # 设置边框颜色和填充透明度
    for i, patch in enumerate(ax.patches):
        svtype = df["Type"].iloc[i]
        # 设置边框颜色为原色（不透明）
        patch.set_edgecolor(color_dict[svtype])
        patch.set_linewidth(2)
        # 设置填充颜色为带透明度的颜色
        import matplotlib.colors as mcolors
        face_color = mcolors.to_rgba(color_dict[svtype], alpha=0.2)
        patch.set_facecolor(face_color)
        # 不要设置整体alpha，这样边框保持不透明
        # patch.set_alpha(0.2)  # 注释掉这行
    
    for index, row in df.iterrows():
        plt.text(index, row["Count"] * 1.1,  # 上移一点
                 f'{row["Count"]:,}',
                 ha='center', va='bottom',
                 fontsize=12)
    
    # 标签和样式
    plt.title("Counts of Structural Variant Types", fontsize=12, pad=10) 
    plt.xlabel("SV Type")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.yscale("log")
    plt.ylabel("SV Count (log10)", fontsize=12)
    plt.tight_layout(pad=2.0)
    plt.margins(y=0.1)
    sns.despine()
    
    # 保存图像
    plt.savefig("sv_type_barplot_colors.png", dpi=300)
    plt.show()
    ```
    

### 3. 分类

- Method
    
    
    | 类型 | 判定条件 | 颜色 |
    | --- | --- | --- |
    | **共享 SV (shared SV)** | 在全部的样本中支持 | #34679a |
    | **主要 SV (major SV)** | 在 ≥50% 的样本中支持 | #87b9d6 |
    | **多态 SV (polymorphic SV)** | 在 >1 个样本中支持 | #f5b785 |
    | **个体 SV (Singleton SV)** | 仅 1 个样本支持 | #e27861 |
    - 1. 从 VCF 中提取 GT 信息，生成PAV矩阵  `python vcf2pav.py`
        
        ```python
        # This script is used to convert vcf file to SV's PAV file
        # 6.26.20205 Fosquer.
        
        import csv
        
        vcf_file = "merged.vcf"
        sample_file = "sample_name.txt"
        output_file = "SV_PAV_matrix.tsv"
        
        # 读取样本名列表
        with open(sample_file) as f:
            sample_names = [line.strip() for line in f if line.strip()]
        
        n_samples = len(sample_names)
        print(f"Detected {n_samples} samples from sample_name.txt")
        
        # 处理VCF文件
        matrix = []
        with open(vcf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                chrom = fields[0]
                sv_id = fields[2]
                info = fields[7]
        
                # 提取SUPP_VEC
                supp_vec = None
                for item in info.split(';'):
                    if item.startswith("SUPP_VEC="):
                        supp_vec = item.replace("SUPP_VEC=", "")
                        break
        
                if supp_vec is None:
                    continue  # 没有SUPP_VEC字段，跳过
        
                if len(supp_vec) != n_samples:
                    print(f"[WARNING] SV {sv_id}: SUPP_VEC 长度 {len(supp_vec)} 与样本数 {n_samples} 不一致，跳过")
                    continue
        
                matrix.append([chrom, sv_id] + list(supp_vec))
        
        # 写出PAV矩阵（加上染色体列）
        with open(output_file, 'w', newline='') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(["CHROM", "SV_ID"] + sample_names)
            writer.writerows(matrix)
        
        print(f"✅ 输出完成：{output_file}")
        
        ```
        
        输出 TSV 示例如下：
        
        ```
        Chr01	DEL00000000	1	0	0	0	1	0	1	1	1	0
        Chr01	DEL00000004	0	0	0	0	0	0	0	0	0	0
        ```
        
    - 2. summary VCF信息，生成SV摘要 `python summarize_PAV.py`
        
        <aside>
        👉🏻
        
        包括：CHROM   SV_ID   SVTYPE  SVLEN   SUPP_NUM    *CATEGORY*
        
        </aside>
        
        ```python
        # This script is used to summary important SV info into tsv file
        # 6.26.2025 Fosquer.
        
        import csv
        
        vcf_file = "merged.vcf"
        sample_file = "sample_name.txt"
        output_file = "SV_summary.tsv"
        
        # 读取样本名
        with open(sample_file) as f:
            sample_names = [line.strip() for line in f if line.strip()]
        n_samples = len(sample_names)
        
        print(f"Detected {n_samples} samples from sample_name.txt")
        
        def classify_supp(n):
            if n == n_samples:
                return "Shared SV"
            elif n >= n_samples / 2:
                return "Major SV"
            elif n > 1:
                return "Polymorphic SV"
            else:
                return "Singleton SV"
        
        summary = []
        
        with open(vcf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                chrom = fields[0]
                sv_id = fields[2]
                info = fields[7]
        
                # 提取字段
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        info_dict[k] = v
        
                supp_vec = info_dict.get("SUPP_VEC", "")
                if len(supp_vec) != n_samples:
                    print(f"[WARNING] {sv_id}: SUPP_VEC 长度不符，跳过")
                    continue
        
                svtype = info_dict.get("SVTYPE", "NA")
                svlen = info_dict.get("SVLEN", "NA")
                try:
                    svlen = abs(int(svlen))
                except:
                    svlen = "NA"
        
                supp_num = supp_vec.count("1")
                category = classify_supp(supp_num)
        
                summary.append([chrom, sv_id, svtype, svlen, supp_num, category])
        
        # 写出结果
        with open(output_file, "w", newline='') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(["CHROM", "SV_ID", "SVTYPE", "SVLEN", "SUPP_NUM", "CATEGORY"])
            writer.writerows(summary)
        
        print(f"✅ SV 分类统计完成，共 {len(summary)} 条记录，输出至 {output_file}")
        
        ```
        
    - 3. 分类堆叠柱状图 `python cate_bardraw.py`
        
        ```python
        import pandas as pd
        import matplotlib.pyplot as plt
        
        # 读取数据
        df = pd.read_csv("SV_summary.tsv", sep='\t')
        
        # 统计每种 SVTYPE 下的各类 CATEGORY 数量
        count_df = df.groupby(["SVTYPE", "CATEGORY"]).size().unstack(fill_value=0)
        
        # 自定义堆叠顺序和颜色
        desired_order = ["Singleton SV", "Polymorphic SV", "Major SV", "Shared SV"]
        color_map = {
            "Singleton SV": "#e27861",
            "Polymorphic SV": "#f5b785",
            "Major SV": "#87b9d6",
            "Shared SV": "#34679a"
        }
        
        # 只保留实际存在的列并排序
        columns_in_data = [c for c in desired_order if c in count_df.columns]
        count_df = count_df[columns_in_data]
        colors = [color_map[c] for c in columns_in_data]
        
        # 绘图（柱子粗一点）
        fig, ax = plt.subplots(figsize=(8, 6))
        bar_width = 0.6  # 默认是 0.5，可以调大一点
        count_df.plot(
            kind='bar',
            stacked=True,
            ax=ax,
            width=bar_width,
            color=colors,
            edgecolor='black'
        )
        
        # 图形美化
        ax.set_ylabel("SV Count", fontsize=12)
        ax.set_xlabel("SVTYPE", fontsize=12)
        ax.set_title("SV Category Composition by SVTYPE", fontsize=14)
        ax.set_xticklabels(count_df.index, rotation=0)
        
        # 图例放在图内部右上角
        ax.legend(
            title="Category",
            loc='upper right',
            bbox_to_anchor=(0.98, 0.98),
            borderaxespad=0.5,
            frameon=True
        )
        
        plt.tight_layout()
        plt.savefig("SVTYPE_category_stacked_bar.png", dpi=300)
        plt.close()
        
        print("✅ 图已保存为 SVTYPE_category_stacked_bar.png")
        ```
        
    - 4. 分类饼图 `python cate_piedraw.py`
        
        ```python
        import pandas as pd
        import matplotlib.pyplot as plt
        
        # 读取数据
        df = pd.read_csv("SV_summary.tsv", sep='\t')
        
        # 汇总各类数量
        cat_counts = df["CATEGORY"].value_counts().reindex([
            "Shared SV", "Major SV", "Polymorphic SV", "Singleton SV"
        ], fill_value=0)
        
        # 颜色（与你的设定保持一致）
        color_map = {
            "Shared SV": "#34679a",
            "Major SV": "#87b9d6",
            "Polymorphic SV": "#f5b785",
            "Singleton SV": "#e27861"
        }
        colors = [color_map[c] for c in cat_counts.index]
        
        # 设置标签：类别（数量，占比%）
        labels = [
            f"{cat}\n{count} ({count / cat_counts.sum() * 100:.1f}%)"
            for cat, count in zip(cat_counts.index, cat_counts.values)
        ]
        
        # 绘制饼图
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.pie(
            cat_counts,
            labels=labels,
            colors=colors,
            startangle=90,
            counterclock=False,
            wedgeprops=dict(width=0.6, edgecolor='white')
        )
        
        plt.title("SV Category Distribution", fontsize=14)
        plt.tight_layout()
        plt.savefig("SV_Category_Pie.png", dpi=300)
        plt.close()
        
        print("✅ 饼图完成，已保存为 SV_Category_Pie.png")
        
        ```