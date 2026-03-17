# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
from mofapy2.run.entry_point import entry_point
import seaborn as sns
import matplotlib.pyplot as plt

# ==================== 1. 设置路径 ====================
base_path = r"C:\Users\ASUS\Desktop\Gut"

gut_file = os.path.join(base_path, "Gut microbiota.xlsx")
met_neg_file = os.path.join(base_path, "Metabolites quantification profiles in negative ion modes.xlsx")
met_pos_file = os.path.join(base_path, "Metabolites quantification profiles in positive ion modes.xlsx")

# ==================== 2. 读取菌群数据 ====================
print("正在读取菌群数据...")
gut = pd.read_excel(gut_file, sheet_name="Sheet1", header=None)
# 第一列是物种名，从第2行开始是数据
taxa_names = gut.iloc[1:, 0].values
samples = gut.iloc[0, 1:].values  # 样本名
microbiome_raw = gut.iloc[1:, 1:].astype(float).values

# 替换0为最小非零值的1/10（避免log(0)）
microbiome_raw[microbiome_raw == 0] = np.nan
min_val = np.nanmin(microbiome_raw) * 0.1
microbiome_raw = np.nan_to_num(microbiome_raw, nan=min_val)

# CLR 变换（必须！）
def clr_transform(x):
    logx = np.log(x)
    return logx - np.mean(logx)

microbiome_clr = np.apply_along_axis(clr_transform, 0, microbiome_raw)
micro_clr_df = pd.DataFrame(microbiome_clr, index=taxa_names, columns=samples).T

# 挑选丰度前150的属（避免特征太多）
top_genera = micro_clr_df.abs().mean().sort_values(ascending=False).head(150).index
micro_clr_df = micro_clr_df[top_genera]

print(f"菌群数据处理完成：{micro_clr_df.shape}")  # 应该是 (50, 150)

# ==================== 3. 读取并合并代谢物数据 ====================
print("正在读取代谢物数据...")
met_neg = pd.read_excel(met_neg_file, sheet_name="Sheet1")
met_pos = pd.read_excel(met_pos_file, sheet_name="Sheet1")

# 合并正负离子模式（以 Name 去重）
met_neg = met_neg[["Name"] + list(samples)]
met_pos = met_pos[["Name"] + list(samples)]
met_all = pd.concat([met_neg, met_pos], axis=0).drop_duplicates(subset="Name")

# 转置 + 去掉全为0或缺失>50%的代谢物
met_values = met_all[list(samples)].T
met_values = met_values.loc[:, met_values.isna().mean() < 0.5]
met_values = met_values.loc[:, (met_values != 0).mean() > 0.1]
met_values = met_values.fillna(met_values.min().min() / 2)

# Log2 + Pareto scaling
met_log = np.log2(met_values)
met_scaled = (met_log - met_log.mean()) / np.sqrt(met_log.std())

met_df = pd.DataFrame(met_scaled.values, index=samples, columns=met_values.columns)

print(f"代谢物数据处理完成：{met_df.shape}")

# ==================== 4. 确保样本顺序一致 ====================
common_samples = sorted(set(micro_clr_df.index) & set(met_df.index))
micro_clr_df = micro_clr_df.loc[common_samples]
met_df = met_df.loc[common_samples]

print(f"最终样本数：{len(common_samples)}")

# ==================== 5. 运行 MOFA+（小样本最优参数） ====================
print("开始运行 MOFA+...")
ent = entry_point()

ent.set_data_matrix(
    data=[met_df.values, micro_clr_df.values],
    views_names=["Metabolites", "Microbiome"],
    samples_names=common_samples,
    features_names=[met_df.columns.tolist(), micro_clr_df.columns.tolist()]
)

ent.set_model_options(factors=10, spikeslab_weights=True, ard_factors=True, ard_weights=True)
ent.set_train_options(iter=1500, convergence_mode="medium", seed=42, verbose=False)

ent.build()
ent.run()

# 保存模型
os.makedirs("MOFA_result", exist_ok=True)
ent.save("MOFA_result/mofa_model.hdf5")

# ==================== 6. 出 Figure 1 四联图 ====================
print("正在生成 Figure 1...")

# A. 方差解释图
ent.plot_r2_variance_explained(file="MOFA_result/Fig1A_variance_explained.pdf")

# B. Factor 1 vs Factor 2（用你的分组上色）
factors = ent.get_factors()
factor_df = pd.DataFrame({
    "Factor1": factors["Factor1"][:, 0],
    "Factor2": factors["Factor2"][:, 1],
    "Group": ["RCD" if "NCD" in s else "Normal" for s in common_samples]
}, index=common_samples)

plt.figure(figsize=(6,5))
sns.scatterplot(data=factor_df, x="Factor1", y="Factor2", hue="Group", s=100, palette=["#E74C3C", "#3498DB"])
plt.legend(title="")
plt.title("Factor 1 vs Factor 2")
plt.savefig("MOFA_result/Fig1B_scatter.pdf", dpi=300, bbox_inches='tight')

# C. Top weights 热图（关键！）
ent.plot_top_weights(
    factors=[0,1,2,3],
    view="all",
    n_features=12,
    file="MOFA_result/Fig1C_top_weights.pdf"
)

# D. 权重热图（手动画更清楚）
weights_met = ent.get_weights(view="Metabolites", as_data_frame=True)
weights_mic = ent.get_weights(view="Microbiome", as_data_frame=True)

print("\n=== MOFA+ 运行完成！===")
print("结果保存在：MOFA_result 文件夹")
print(f"代谢物特征数：{met_df.shape[1]}")
print(f"菌群属数：{micro_clr_df.shape[1]}")
print("请查看 Fig1C_top_weights.pdf —— PEA 和 Ruminococcus 应该排在前面！")