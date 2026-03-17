import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
sns.set_style('white')

print("Loading results...")
with open(r"C:\Users\ASUS\Desktop\Gut\mofa_results.pkl", 'rb') as f:
    mofa_results = pickle.load(f)

with open(r"C:\Users\ASUS\Desktop\Gut\xgboost_results.pkl", 'rb') as f:
    xgb_results = pickle.load(f)

Z = mofa_results['Z']
W_microbiome = mofa_results['W_microbiome']
W_metabolome = mofa_results['W_metabolome']
var_df = mofa_results['variance_explained']
selected_factors = mofa_results['selected_factors']
feature_names_microbiome = mofa_results['feature_names_microbiome']
feature_names_metabolome = mofa_results['feature_names_metabolome']
ncd_status = mofa_results['ncd_status']
ncd_binary = mofa_results['ncd_binary']

importance_df = xgb_results['importance_df']
mean_auc = xgb_results['mean_auc']
ci_lower = xgb_results['ci_lower']
ci_upper = xgb_results['ci_upper']

# Create figure
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35,
                      left=0.08, right=0.95, top=0.95, bottom=0.05)

# ==================== Panel A: Variance Explained ====================
ax_a = fig.add_subplot(gs[0, :2])

# Prepare data for selected factors
selected_var_df = var_df.iloc[selected_factors].copy()
selected_var_df['Factor'] = [f'Factor {i+1}' for i in range(len(selected_factors))]

x = np.arange(len(selected_var_df))
width = 0.35

bars1 = ax_a.bar(x - width/2, selected_var_df['Metabolome'], width,
                 label='Metabolome', color='#E64B35', alpha=0.8, edgecolor='black', linewidth=1)
bars2 = ax_a.bar(x + width/2, selected_var_df['Microbiome'], width,
                 label='Microbiome', color='#4DBBD5', alpha=0.8, edgecolor='black', linewidth=1)

ax_a.set_xlabel('Latent Factor', fontsize=12, fontweight='bold')
ax_a.set_ylabel('Variance Explained (%)', fontsize=12, fontweight='bold')
ax_a.set_title('A. Variance explained by latent factors', fontsize=13, fontweight='bold', loc='left')
ax_a.set_xticks(x)
ax_a.set_xticklabels(selected_var_df['Factor'], rotation=0)
ax_a.legend(frameon=False, fontsize=11, loc='upper right')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.grid(False)

# Add value labels on bars
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        if height > 0.5:
            ax_a.text(bar.get_x() + bar.get_width()/2., height,
                     f'{height:.1f}',
                     ha='center', va='bottom', fontsize=8)

# ==================== Panel B: Factor Scatter Plot ====================
ax_b = fig.add_subplot(gs[0, 2])

# Use Factor 1 and Factor 2 (index 0 and 1)
factor1_idx = 0
factor2_idx = 1

colors = ['#00A087' if status == 'Normal' else '#E64B35' for status in ncd_status]

for i, (status, color) in enumerate(zip(ncd_status, colors)):
    if status == 'Normal':
        label = 'Normal' if i == 0 or ncd_status[i-1] != 'Normal' else None
        marker = 'o'
    else:
        label = 'NCD' if i == 0 or ncd_status[i-1] != 'NCD' else None
        marker = 's'

    ax_b.scatter(Z[i, factor1_idx], Z[i, factor2_idx],
                c=color, s=80, alpha=0.7, edgecolors='black',
                linewidth=0.5, label=label, marker=marker)

# Statistical test
factor1_normal = Z[np.array(ncd_binary) == 0, factor1_idx]
factor1_ncd = Z[np.array(ncd_binary) == 1, factor1_idx]
stat, pval = mannwhitneyu(factor1_normal, factor1_ncd)

ax_b.set_xlabel(f'Factor {factor1_idx+1}', fontsize=11, fontweight='bold')
ax_b.set_ylabel(f'Factor {factor2_idx+1}', fontsize=11, fontweight='bold')
ax_b.set_title('B. Factor separation by NCD status', fontsize=13, fontweight='bold', loc='left')

# Format p-value
if pval < 0.001:
    pval_text = 'p < 0.001'
elif pval < 0.01:
    pval_text = f'p = {pval:.3f}'
else:
    pval_text = f'p = {pval:.2f}'

ax_b.text(0.05, 0.95, f'Wilcoxon {pval_text}',
         transform=ax_b.transAxes, fontsize=9,
         verticalalignment='top', bbox=dict(boxstyle='round',
         facecolor='white', alpha=0.8, edgecolor='gray'))

ax_b.legend(frameon=False, fontsize=10, loc='lower right')
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.grid(False)

# ==================== Panel C: Top Features Heatmap ====================
ax_c = fig.add_subplot(gs[1:, :])

# Select top 4 factors for display
n_factors_display = min(4, len(selected_factors))
top_n_features = 10

# Combine weights from both views for top factors
all_weights = []
all_features = []
all_views = []

for i in range(n_factors_display):
    # Metabolome weights
    metab_weights = W_metabolome[:, i]
    metab_abs = np.abs(metab_weights)
    top_metab_idx = np.argsort(metab_abs)[-top_n_features:][::-1]

    for idx in top_metab_idx:
        all_weights.append(metab_weights[idx])
        all_features.append(feature_names_metabolome[idx])
        all_views.append(f'Factor {selected_factors[i]+1}')

    # Microbiome weights
    micro_weights = W_microbiome[:, i]
    micro_abs = np.abs(micro_weights)
    top_micro_idx = np.argsort(micro_abs)[-top_n_features:][::-1]

    for idx in top_micro_idx:
        all_weights.append(micro_weights[idx])
        all_features.append(feature_names_microbiome[idx])
        all_views.append(f'Factor {selected_factors[i]+1}')

# Create DataFrame
heatmap_df = pd.DataFrame({
    'Feature': all_features,
    'Factor': all_views,
    'Weight': all_weights
})

# Pivot for heatmap
heatmap_pivot = heatmap_df.pivot_table(index='Feature', columns='Factor',
                                       values='Weight', fill_value=0)

# Sort by absolute weight in Factor 1
heatmap_pivot['abs_sum'] = heatmap_pivot.abs().sum(axis=1)
heatmap_pivot = heatmap_pivot.sort_values('abs_sum', ascending=False)
heatmap_pivot = heatmap_pivot.drop('abs_sum', axis=1)

# Take top 15 features
heatmap_pivot = heatmap_pivot.head(15)

# Shorten feature names for display
short_names = []
for feat in heatmap_pivot.index:
    if len(feat) > 40:
        short_names.append(feat[:37] + '...')
    else:
        short_names.append(feat)

heatmap_pivot.index = short_names

# Plot heatmap
sns.heatmap(heatmap_pivot, cmap='RdBu_r', center=0,
           cbar_kws={'label': 'Feature Weight', 'shrink': 0.6},
           linewidths=0.5, linecolor='gray',
           ax=ax_c, vmin=-0.3, vmax=0.3,
           annot=False, fmt='.2f')

ax_c.set_title('C. Top weighted features in latent factors',
              fontsize=13, fontweight='bold', loc='left', pad=15)
ax_c.set_xlabel('Latent Factor', fontsize=12, fontweight='bold')
ax_c.set_ylabel('Feature', fontsize=12, fontweight='bold')
ax_c.set_xticklabels(ax_c.get_xticklabels(), rotation=0)
ax_c.set_yticklabels(ax_c.get_yticklabels(), rotation=0, fontsize=9)

# Highlight phenylethylamine and Ruminococcus if present
for i, label in enumerate(heatmap_pivot.index):
    if 'phenylethylamine' in label.lower() or 'pea' in label.lower():
        ax_c.get_yticklabels()[i].set_color('#E64B35')
        ax_c.get_yticklabels()[i].set_weight('bold')
    if 'Ruminococcus' in label:
        ax_c.get_yticklabels()[i].set_color('#4DBBD5')
        ax_c.get_yticklabels()[i].set_weight('bold')

# ==================== Panel D: Permutation Importance ====================
# Create inset axis for panel D
ax_d = fig.add_axes([0.65, 0.08, 0.28, 0.22])

importance_sorted = importance_df.sort_values('Importance', ascending=True)

bars = ax_d.barh(range(len(importance_sorted)), importance_sorted['Importance'],
                 color='#3C5488', alpha=0.8, edgecolor='black', linewidth=1)

ax_d.set_yticks(range(len(importance_sorted)))
ax_d.set_yticklabels(importance_sorted['Factor'], fontsize=9)
ax_d.set_xlabel('Permutation Importance', fontsize=11, fontweight='bold')
ax_d.set_title(f'D. XGBoost Feature Importance\nAUC = {mean_auc:.2f} (95% CI: {ci_lower:.2f}-{ci_upper:.2f})',
              fontsize=11, fontweight='bold', loc='left')
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)
ax_d.grid(False)

# Highlight top 2 factors
for i, bar in enumerate(bars):
    if i >= len(bars) - 2:
        bar.set_color('#E64B35')

# Save figure
plt.savefig(r"C:\Users\ASUS\Desktop\Gut\Figure1_MOFA_Analysis.png",
           dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(r"C:\Users\ASUS\Desktop\Gut\Figure1_MOFA_Analysis.pdf",
           bbox_inches='tight', facecolor='white')

print("\nFigure 1 saved!")
print("Files created:")
print("  - Figure1_MOFA_Analysis.png")
print("  - Figure1_MOFA_Analysis.pdf")

plt.show()
