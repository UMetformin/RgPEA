import pandas as pd
import numpy as np

# Load and process data
microbiome = pd.read_csv(r"C:\Users\ASUS\Desktop\Gut\Gut microbiota.csv", index_col=0, encoding='utf-8-sig')

def extract_genus(idx):
    if pd.isna(idx):
        return 'Unknown'
    idx_str = str(idx)
    if '|g__' in idx_str:
        return idx_str.split('|')[-2].replace('g__', '')
    else:
        return idx_str

microbiome.index = [extract_genus(idx) for idx in microbiome.index]
microbiome = microbiome.groupby(microbiome.index).sum()
microbiome_t = microbiome.T
microbiome_t.index = microbiome_t.index.str.replace(' ', '_')

# Process metabolomics
metab_neg_raw = pd.read_excel(r"C:\Users\ASUS\Desktop\Gut\Metabolites quantification profiles in negative ion modes.xlsx")
metab_pos_raw = pd.read_excel(r"C:\Users\ASUS\Desktop\Gut\Metabolites quantification profiles in positive ion modes.xlsx")

metab_neg = metab_neg_raw.copy()
metab_neg['CompoundName'] = metab_neg['Name'].fillna(metab_neg['Compound_ID'])
metab_neg = metab_neg.set_index('CompoundName')
metab_neg_numeric = metab_neg.select_dtypes(include=[np.number])
metab_neg_t = metab_neg_numeric.T

metab_pos = metab_pos_raw.copy()
metab_pos['CompoundName'] = metab_pos['Name'].fillna(metab_pos['Compound_ID'])
metab_pos = metab_pos.set_index('CompoundName')
metab_pos_numeric = metab_pos.select_dtypes(include=[np.number])
metab_pos_t = metab_pos_numeric.T

metabolome_combined = pd.concat([metab_neg_t, metab_pos_t], axis=1)
metabolome_dedup = metabolome_combined.loc[:, ~metabolome_combined.columns.duplicated(keep='first')]
metabolome_dedup.index = metabolome_dedup.index.str.replace(' ', '_')

# Get common samples
common_samples = list(set(microbiome_t.index) & set(metabolome_dedup.index))
common_samples.sort()

microbiome_final = microbiome_t.loc[common_samples]
metabolome_final = metabolome_dedup.loc[common_samples]

# Check for duplicates
micro_features = set(microbiome_final.columns)
metab_features = set(metabolome_final.columns)

print(f"Microbiome features: {len(micro_features)}")
print(f"Metabolome features: {len(metab_features)}")

# Find overlapping feature names
overlapping = micro_features & metab_features

if overlapping:
    print(f"\nOVERLAPPING FEATURE NAMES FOUND: {len(overlapping)}")
    print("\nOverlapping names:")
    for name in sorted(overlapping):
        print(f"  - {name}")
else:
    print("\nNo overlapping feature names found!")

# Check for duplicates within each view
micro_dups = [item for item, count in pd.Series(list(micro_features)).value_counts().items() if count > 1]
metab_dups = [item for item, count in pd.Series(list(metab_features)).value_counts().items() if count > 1]

if micro_dups:
    print(f"\nDuplicates within microbiome: {micro_dups}")
else:
    print("\nNo duplicates within microbiome")

if metab_dups:
    print(f"\nDuplicates within metabolome: {metab_dups}")
else:
    print("\nNo duplicates within metabolome")

# When adding prefixes
micro_with_prefix = ['Microbiome_' + str(f) for f in microbiome_final.columns]
metab_with_prefix = ['Metabolome_' + str(f) for f in metabolome_final.columns]

all_features = micro_with_prefix + metab_with_prefix
print(f"\nTotal features with prefixes: {len(all_features)}")
print(f"Unique features with prefixes: {len(set(all_features))}")

if len(all_features) != len(set(all_features)):
    print("\nDUPLICATES FOUND EVEN WITH PREFIXES!")
    from collections import Counter
    counts = Counter(all_features)
    dups = [item for item, count in counts.items() if count > 1]
    print(f"Duplicated entries: {dups[:10]}")
