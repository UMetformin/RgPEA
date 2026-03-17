import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# Set random seed
np.random.seed(42)

print("Loading MOFA results...")
with open(r"C:\Users\ASUS\Desktop\Gut\mofa_results.pkl", 'rb') as f:
    results = pickle.load(f)

Z = results['Z']
W_microbiome = results['W_microbiome']
W_metabolome = results['W_metabolome']
var_df = results['variance_explained']
selected_factors = results['selected_factors']
feature_names_microbiome = results['feature_names_microbiome']
feature_names_metabolome = results['feature_names_metabolome']
sample_names = results['sample_names']
ncd_status = results['ncd_status']
ncd_binary = results['ncd_binary']

print(f"Latent factors: {Z.shape}")
print(f"Samples: {len(sample_names)}")
print(f"NCD samples: {sum(ncd_binary)}, Normal samples: {len(ncd_binary) - sum(ncd_binary)}")

# Step 5: XGBoost model with bootstrap
print("\nStep 5: Training XGBoost model with bootstrap...")

n_bootstrap = 1000
auc_scores = []
feature_importances = []

for i in range(n_bootstrap):
    if (i + 1) % 100 == 0:
        print(f"Bootstrap iteration {i+1}/{n_bootstrap}")

    # Bootstrap sampling
    indices = np.random.choice(len(Z), size=len(Z), replace=True)
    X_boot = Z[indices]
    y_boot = np.array(ncd_binary)[indices]

    # Train XGBoost (using GradientBoostingClassifier as equivalent)
    model = GradientBoostingClassifier(
        max_depth=2,
        n_estimators=300,
        learning_rate=0.05,
        subsample=0.8,
        random_state=42
    )

    model.fit(X_boot, y_boot)

    # Predict on full dataset
    y_pred_proba = model.predict_proba(Z)[:, 1]
    auc = roc_auc_score(ncd_binary, y_pred_proba)
    auc_scores.append(auc)

    # Feature importance
    feature_importances.append(model.feature_importances_)

# Calculate statistics
mean_auc = np.mean(auc_scores)
ci_lower = np.percentile(auc_scores, 2.5)
ci_upper = np.percentile(auc_scores, 97.5)

print(f"\nXGBoost Results:")
print(f"Mean AUC: {mean_auc:.3f}")
print(f"95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")

# Mean feature importance
mean_importance = np.mean(feature_importances, axis=0)
importance_df = pd.DataFrame({
    'Factor': [f'Factor {selected_factors[i]+1}' for i in range(len(selected_factors))],
    'Importance': mean_importance
})
importance_df = importance_df.sort_values('Importance', ascending=False)

print("\nFeature Importance:")
print(importance_df)

# Save XGBoost results
xgb_results = {
    'auc_scores': auc_scores,
    'mean_auc': mean_auc,
    'ci_lower': ci_lower,
    'ci_upper': ci_upper,
    'importance_df': importance_df,
    'feature_importances': feature_importances
}

with open(r"C:\Users\ASUS\Desktop\Gut\xgboost_results.pkl", 'wb') as f:
    pickle.dump(xgb_results, f)

print("\nXGBoost analysis completed!")
print("Results saved to: xgboost_results.pkl")
