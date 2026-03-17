import pandas as pd
import numpy as np
import pickle
import warnings
warnings.filterwarnings('ignore')

print("Loading MOFA results...")
with open(r"C:\Users\ASUS\Desktop\Gut\mofa_results.pkl", 'rb') as f:
    results = pickle.load(f)

W_microbiome = results['W_microbiome']
W_metabolome = results['W_metabolome']
feature_names_microbiome = results['feature_names_microbiome']
feature_names_metabolome = results['feature_names_metabolome']
selected_factors = results['selected_factors']

print("\n" + "="*80)
print("CHECKING TARGET FEATURES IN TOP RANKINGS")
print("="*80)

# Check for phenylethylamine in metabolome
print("\n1. PHENYLETHYLAMINE in Metabolome:")
print("-" * 80)

pea_found = False
for factor_idx in range(W_metabolome.shape[1]):
    factor_num = selected_factors[factor_idx] + 1
    weights = W_metabolome[:, factor_idx]
    abs_weights = np.abs(weights)

    # Get top 20 features
    top_indices = np.argsort(abs_weights)[-20:][::-1]
    top_features = [(feature_names_metabolome[i], weights[i], abs_weights[i])
                    for i in top_indices]

    # Check if phenylethylamine is in top 20
    for rank, (feat, weight, abs_weight) in enumerate(top_features, 1):
        if 'phenylethylamine' in feat.lower() or 'pea' in feat.lower():
            print(f"\n  Factor {factor_num}:")
            print(f"    Rank: {rank}")
            print(f"    Feature: {feat}")
            print(f"    Weight: {weight:.4f}")
            print(f"    Absolute Weight: {abs_weight:.4f}")
            pea_found = True

if not pea_found:
    print("\n  WARNING: Phenylethylamine NOT found in top 20 of any factor!")
    print("\n  Searching all metabolome features containing 'phenyl' or 'amine':")
    for i, feat in enumerate(feature_names_metabolome):
        if 'phenyl' in feat.lower() or 'amine' in feat.lower():
            print(f"    - {feat}")
            # Check its weight in each factor
            for factor_idx in range(W_metabolome.shape[1]):
                factor_num = selected_factors[factor_idx] + 1
                weight = W_metabolome[i, factor_idx]
                abs_weight = np.abs(weight)
                if abs_weight > 0.01:
                    # Find rank
                    all_abs_weights = np.abs(W_metabolome[:, factor_idx])
                    rank = np.sum(all_abs_weights > abs_weight) + 1
                    print(f"      Factor {factor_num}: Weight={weight:.4f}, Rank={rank}")

# Check for Ruminococcus in microbiome
print("\n" + "="*80)
print("2. RUMINOCOCCUS in Microbiome:")
print("-" * 80)

rum_found = False
for factor_idx in range(W_microbiome.shape[1]):
    factor_num = selected_factors[factor_idx] + 1
    weights = W_microbiome[:, factor_idx]
    abs_weights = np.abs(weights)

    # Get top 20 features
    top_indices = np.argsort(abs_weights)[-20:][::-1]
    top_features = [(feature_names_microbiome[i], weights[i], abs_weights[i])
                    for i in top_indices]

    # Check if Ruminococcus is in top 20
    for rank, (feat, weight, abs_weight) in enumerate(top_features, 1):
        if 'Ruminococcus' in feat:
            print(f"\n  Factor {factor_num}:")
            print(f"    Rank: {rank}")
            print(f"    Feature: {feat}")
            print(f"    Weight: {weight:.4f}")
            print(f"    Absolute Weight: {abs_weight:.4f}")
            rum_found = True

if not rum_found:
    print("\n  WARNING: Ruminococcus NOT found in top 20 of any factor!")
    print("\n  Searching all microbiome features containing 'Ruminococcus':")
    for i, feat in enumerate(feature_names_microbiome):
        if 'Ruminococcus' in feat:
            print(f"    - {feat}")
            # Check its weight in each factor
            for factor_idx in range(W_microbiome.shape[1]):
                factor_num = selected_factors[factor_idx] + 1
                weight = W_microbiome[i, factor_idx]
                abs_weight = np.abs(weight)
                if abs_weight > 0.001:
                    # Find rank
                    all_abs_weights = np.abs(W_microbiome[:, factor_idx])
                    rank = np.sum(all_abs_weights > abs_weight) + 1
                    print(f"      Factor {factor_num}: Weight={weight:.4f}, Rank={rank}")

# Print top 10 features for each factor
print("\n" + "="*80)
print("TOP 10 FEATURES FOR EACH FACTOR")
print("="*80)

for factor_idx in range(min(4, W_metabolome.shape[1])):
    factor_num = selected_factors[factor_idx] + 1

    print(f"\n{'='*80}")
    print(f"FACTOR {factor_num}")
    print(f"{'='*80}")

    # Metabolome
    print(f"\n  METABOLOME (Top 10):")
    print("  " + "-"*76)
    weights_metab = W_metabolome[:, factor_idx]
    abs_weights_metab = np.abs(weights_metab)
    top_idx_metab = np.argsort(abs_weights_metab)[-10:][::-1]

    for rank, idx in enumerate(top_idx_metab, 1):
        feat = feature_names_metabolome[idx]
        if len(feat) > 50:
            feat = feat[:47] + "..."
        print(f"    {rank:2d}. {feat:50s}  Weight: {weights_metab[idx]:7.4f}")

    # Microbiome
    print(f"\n  MICROBIOME (Top 10):")
    print("  " + "-"*76)
    weights_micro = W_microbiome[:, factor_idx]
    abs_weights_micro = np.abs(weights_micro)
    top_idx_micro = np.argsort(abs_weights_micro)[-10:][::-1]

    for rank, idx in enumerate(top_idx_micro, 1):
        feat = feature_names_microbiome[idx]
        if len(feat) > 50:
            feat = feat[:47] + "..."
        print(f"    {rank:2d}. {feat:50s}  Weight: {weights_micro[idx]:7.4f}")

print("\n" + "="*80)
print("FEATURE CHECK COMPLETED")
print("="*80)

# Save detailed report
report_lines = []
report_lines.append("MOFA FEATURE WEIGHT ANALYSIS REPORT\n")
report_lines.append("="*80 + "\n\n")

# Check if targets are found
if pea_found:
    report_lines.append("✓ Phenylethylamine found in top 20 features\n")
else:
    report_lines.append("✗ Phenylethylamine NOT in top 20 features - may need parameter adjustment\n")

if rum_found:
    report_lines.append("✓ Ruminococcus found in top 20 features\n")
else:
    report_lines.append("✗ Ruminococcus NOT in top 20 features - may need parameter adjustment\n")

report_lines.append("\n")

with open(r"C:\Users\ASUS\Desktop\Gut\feature_check_report.txt", 'w', encoding='utf-8') as f:
    f.writelines(report_lines)

print("\nReport saved to: feature_check_report.txt")
