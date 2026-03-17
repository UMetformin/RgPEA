# -*- coding: utf-8 -*-
# 文件名: scheme_B_final.py
# 路径保存到: C:\Users\ASUS\Desktop\Gut\scheme_B_final.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
import shap
import os

plt.rcParams['font.size'] = 12
sns.set(style="whitegrid", font_scale=1.2)

# ==================== 1. 路径和文件 ====================
base_path = r"C:\Users\ASUS\Desktop\Gut"
gut_file = os.path.join(base_path, "Gut microbiota.xlsx")
met_neg_file = os.path.join(base_path, "Metabolites quantification profiles in negative ion modes.xlsx")
met_pos_file = os.path.join(base_path, "Metabolites quantification profiles in positive ion modes.xlsx")

# ==================== 2. 读取并合并代谢组（正负离子去重） ====================
print("正在读取并合并代谢组...")
met_neg = pd.read_excel(met_neg_file)
met_pos = pd.read_excel(met_pos_file)

# 用 Name 列，缺失的用 Compound_ID 补
met_neg['feat'] = met_neg['Name'].fillna(met_neg['Compound_ID'])
met_pos['feat'] = met_pos['Name'].fillna(met_pos['Compound_ID'])
met_neg = met_neg.set_index('feat')
met_pos = met_pos.set_index('feat')

# 去重 + 转置
met_all = pd.concat([met_neg, met_pos]).loc[:, ~pd.concat([met_neg, met_pos]).columns.duplicated()]
met = met_all.T

# 清洗
met = met.loc[:, (met > 0).mean() > 0.1]                    # 至少10%样本检出
met = met.loc[:, met.isna().mean() < 0.5]                  # 缺失率<50%
met = met.fillna(met.min().min() / 2)

# log2 + Pareto scaling
met_log = np.log2(met.replace(0, np.nan)).fillna(0)
met_final = (met_log - met_log.mean()) / np.sqrt(met_log.std())

print(f"代谢物最终特征数: {met_final.shape[1]}")

# ==================== 3. 读取菌群并处理到属水平 ====================
print("正在处理菌群数据...")
gut = pd.read_excel(gut_file, header=None)
samples = gut.iloc[0, 1:].tolist()
taxa = gut.iloc[1:, 0].astype(str)

# 提取属名
def get_genus(t):
    if 'g__' in t:
        return t.split('g__')[1].split(';')[0]
    else:
        return "Unknown"

genera = [get_genus(t) for t in taxa]
micro_raw = gut.iloc[1:, 1:].astype(float).values
micro_df = pd.DataFrame(micro_raw, index=genera, columns=samples).T
micro_df = micro_df.groupby(micro_df.index.name, axis=1).sum()  # 相同属合并

# CLR 变换
micro_clr = micro_df + 1e-6
micro_clr = np.log(micro_clr) - np.log(micro_clr).mean()
micro_final = pd.DataFrame(micro_clr, index=micro_df.index, columns=micro_df.columns)

print(f"菌群属水平特征数: {micro_final.shape[1]}")

# ==================== 4. 对齐样本 + 标签 ====================
common_samples = sorted(set(met_final.index) & set(micro_final.index))
X_met = met_final.loc[common_samples]
X_mic = micro_final.loc[common_samples]
X = pd.concat([X_met.add_prefix('Met_'), X_mic.add_prefix('Mic_')], axis=1)

y = [1 if 'NCD' in s else 0 for s in common_samples]
group_names = ['NCD' if 'NCD' in s else 'Normal' for s in common_samples]

print(f"最终样本数: {len(common_samples)} (NCD={sum(y)}, Normal={len(y)-sum(y)})")

# ==================== 5. Figure 1A：火山图 ====================
print("正在生成火山图...")
pvals = []
log2fc = []
names = []

for col in X.columns:
    g1 = X.loc[y==1, col]
    g2 = X.loc[y==0, col]
    stat, p = mannwhitneyu(g1, g2, alternative='two-sided')
    fc = g1.mean() - g2.mean()
    pvals.append(p)
    log2fc.append(fc)
    names.append(col)

volcano = pd.DataFrame({'Feature': names, 'log2FC': log2fc, '-log10p': -np.log10(pvals)})
volcano['Significant'] = (volcano['-log10p'] > 1.3) & (abs(volcano['log2FC']) > 0.5)

# 高亮 PEA 和 Ruminococcus
volcano['Label'] = ''
pea_idx = volcano['Feature'].str.contains('phenylethylamine|phenethylamine', case=False, na=False)
rum_idx = volcano['Feature'].str.contains('Ruminococcus', case=False, na=False)
volcano.loc[pea_idx, 'Label'] = 'Phenylethylamine'
volcano.loc[rum_idx, 'Label'] = 'Ruminococcus'

plt.figure(figsize=(8,6))
sns.scatterplot(data=volcano, x='log2FC', y='-log10p', hue='Significant', palette={True:'red', False:'grey'}, s=60)
for i, row in volcano[volcano['Label']!=''].iterrows():
    plt.text(row['log2FC']+0.05, row['-log10p']+0.1, row['Label'], fontsize=14, fontweight='bold')
plt.axhline(1.3, color='black', linestyle='--')
plt.axvline(0.5, color='black', linestyle='--')
plt.axvline(-0.5, color='black', linestyle='--')
plt.xlabel('Log2 Fold Change (NCD vs Normal)')
plt.ylabel('-Log10 p-value')
plt.title('Volcano Plot of Multi-omics Features')
plt.legend(title='Significant', loc='upper left')
plt.tight_layout()
plt.savefig(r"C:\Users\ASUS\Desktop\Gut\Fig1A_volcano.pdf", dpi=300)
plt.savefig(r"C:\Users\ASUS\Desktop\Gut\Fig1A_volcano.png", dpi=300)
plt.show()

# ==================== 6. Figure 1B：LASSO 特征选择 ====================
print("正在运行 LASSO...")
lasso = LogisticRegressionCV(
    Cs=10, penalty='l1', solver='liblinear', cv=StratifiedKFold(5),
    max_iter=10000, random_state=42, scoring='roc_auc'
)
lasso.fit(X, y)
selected = np.where(lasso.coef_[0] != 0)[0]
selected_features = X.columns[selected]
print(f"LASSO 选出 {len(selected_features)} 个特征")

# 画系数路径
plt.figure(figsize=(8,5))
plt.barh(range(len(selected_features)), lasso.coef_[0][selected])
plt.yticks(range(len(selected_features)), selected_features)
plt.xlabel('LASSO Coefficient')
plt.title('LASSO Selected Features')
plt.tight_layout()
plt.savefig(r"C:\Users\ASUS\Desktop\Gut\Fig1B_LASSO.pdf", dpi=300)
plt.show()

# ==================== 7. Figure 1C：XGBoost + SHAP ====================
print("正在运行 XGBoost + SHAP...")
model = GradientBoostingClassifier(
    n_estimators=500, learning_rate=0.01, max_depth=3, subsample=0.8, random_state=42
)
model.fit(X, y)

explainer = shap.Explainer(model, X)
shap_values = explainer(X)

plt.figure(figsize=(8,6))
shap.summary_plot(shap_values, X, plot_type="bar", max_display=15, show=False)
plt.title("Top 15 Features by Mean |SHAP value|")
plt.tight_layout()
plt.savefig(r"C:\Users\ASUS\Desktop\Gut\Fig1C_SHAP_bar.pdf", dpi=300, bbox_inches='tight')
plt.show()

# 详细SHAP图
plt.figure(figsize=(8,6))
shap.summary_plot(shap_values, X, max_display=20, show=False)
plt.savefig(r"C:\Users\ASUS\Desktop\Gut\Fig1C_SHAP_beeswarm.pdf", dpi=300, bbox_inches='tight')
plt.show()

print("\n所有图已保存到桌面 Gut 文件夹！")
print("Fig1A: 火山图（PEA和Ruminococcus高亮）")
print("Fig1B: LASSO选中的特征")
print("Fig1C: SHAP值排名（大概率PEA和Ruminococcus在前5）")