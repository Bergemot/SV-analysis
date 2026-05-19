### 1. 提取每个样本的基因型数据（GT）

- 生成 `sv_genotypes.tsv` 文件

```python
# 提取 GT 信息，每个变异一行，样本是列
bcftools query -f '%CHROM\t%POS[\t%GT]\n' merged.vcf > sv_genotypes.tsv
```

### 2. 将基因型的tsv文件转换成SV的PAV矩阵文件

- 生成 `sv_pav_matrix.csv` 文件
- python tsv2pav.py
    
    ```python
    import pandas as pd
    
    # 读取文件（注意：没有表头）
    df = pd.read_csv("sv_genotypes.tsv", sep="\t", header=None)
    
    # 添加 SV ID，例如 chr1_12345
    df['sv_id'] = df[0].astype(str) + "_" + df[1].astype(str)
    
    # 提取 GT 信息列（从第3列开始到倒数第2列）
    gt_df = df.iloc[:, 2:-1]
    gt_df.index = df['sv_id']
    
    # 将 GT 转换为 PAV（0, 1, NA）
    def gt_to_pav(gt):
        if gt in ["1/1", "0/1", "1/0"]:  # 只要携带 SV，就算存在
            return 1
        elif gt == "0/0":
            return 0
        else:
            return pd.NA  # 缺失或其他情况
    
    pav_df = gt_df.applymap(gt_to_pav)
    
    # 可选：为样本加列名（如果你有一个样本名列表）
    # pav_df.columns = ["Sample1", "Sample2", ..., "SampleN"]
    
    # 保存输出
    pav_df.to_csv("sv_pav_matrix.csv")
    
    ```
    

### 3. 生成增长曲线

- 生成 `sv_growth_curve.png` 和 `sv_growth_mean_std.csv` 文件
- python curve.py
    
    ```python
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    # === 参数 ===
    input_file = "sv_pav_matrix.csv"  # 输入文件
    n_repeats = 100  # 随机重复次数
    
    # === 读取 PAV 数据 ===
    df = pd.read_csv(input_file, index_col=0)
    n_samples = df.shape[1]
    
    # === 存储每次的增长曲线 ===
    growth_curves = np.zeros((n_repeats, n_samples))
    
    for i in range(n_repeats):
        shuffled_samples = np.random.permutation(df.columns)
        cumulative_counts = []
    
        for j in range(n_samples):
            current_samples = shuffled_samples[:j+1]
            subset = df[current_samples]
            present_sv = subset.sum(axis=1) > 0
            cumulative_counts.append(present_sv.sum())
        
        growth_curves[i, :] = cumulative_counts
    
    # === 计算平均值和标准差 ===
    mean_curve = growth_curves.mean(axis=0)
    std_curve = growth_curves.std(axis=0)
    x = np.arange(1, n_samples + 1)
    
    # === 保存为 CSV 文件 ===
    
    # 1. 平均值 + 标准差
    summary_df = pd.DataFrame({
        'sample_number': x,
        'mean_sv_count': mean_curve,
        'std_sv_count': std_curve
    })
    summary_df.to_csv("sv_growth_mean_std.csv", index=False)
    
    # 2. 所有重复的增长曲线
    growth_df = pd.DataFrame(growth_curves, columns=[f"sample_{i}" for i in x])
    growth_df.to_csv("sv_growth_repeats.csv", index=False)
    
    # === 绘图 ===
    plt.figure(figsize=(8,6))
    plt.plot(x, mean_curve, color='blue', label='Mean SV count')
    plt.fill_between(x, mean_curve - std_curve, mean_curve + std_curve, alpha=0.3, color='blue', label='±1 SD')
    plt.xlabel("Number of Samples")
    plt.ylabel("Cumulative Non-redundant SVs")
    plt.title("Structural Variant Accumulation Curve")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig("sv_growth_curve.png", dpi=300)
    plt.show()
    ```
    

### 4. 绘图

- 生成 `sv_growth_curve_debug.png` 文件
- python draw.py
    
    ```python
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    # === 读取保存好的平均值 + 标准差文件 ===
    df = pd.read_csv("sv_growth_mean_std.csv")
    
    x = df['sample_number']
    mean_curve = df['mean_sv_count']
    std_curve = df['std_sv_count']
    
    # === 绘图 ===
    plt.figure(figsize=(10,6))
    plt.plot(x, mean_curve, color='#374f99', label='Mean SV count')
    plt.fill_between(x, mean_curve - std_curve, mean_curve + std_curve, alpha=0.3, color='#6ea3c4', label='±1 SD')
    plt.xlabel("Number of Samples", fontsize=12)
    plt.ylabel("Number of Non-redundant SVs", fontsize=12)
    plt.title("Growth Curve of Structural Variations (SVs)", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5, axis='y')
    plt.legend()
    plt.tight_layout()
    plt.savefig("sv_growth_curve_debug.png", dpi=300)
    plt.show()
    
    ```