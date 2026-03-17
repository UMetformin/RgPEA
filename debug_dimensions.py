import pandas as pd
import numpy as np

# Quick check of data dimensions
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

print("Microbiome data:")
print(f"Original microbiome shape: {microbiome.shape} (features x samples)")
print(f"Transposed microbiome shape: {microbiome_t.shape} (samples x features)")
print(f"Number of samples: {len(microbiome_t)}")
print(f"Number of genera: {len(microbiome_t.columns)}")
print(f"\nFirst few sample names: {list(microbiome_t.index[:5])}")

metab_neg_raw = pd.read_excel(r"C:\Users\ASUS\Desktop\Gut\Metabolites quantification profiles in negative ion modes.xlsx")
metab_pos_raw = pd.read_excel(r"C:\Users\ASUS\Desktop\Gut\Metabolites quantification profiles in positive ion modes.xlsx")

print("\n" + "="*80)
print("Metabolome negative:")
print(f"Shape: {metab_neg_raw.shape}")
print(f"Columns: {list(metab_neg_raw.columns[:10])}")

print("\nMetabolome positive:")
print(f"Shape: {metab_pos_raw.shape}")
print(f"Columns: {list(metab_pos_raw.columns[:10])}")

# Check sample columns
sample_cols_neg = [col for col in metab_neg_raw.columns if 'Normal' in str(col) or 'NCD' in str(col)]
print(f"\nNumber of sample columns in negative mode: {len(sample_cols_neg)}")
print(f"Sample columns: {sample_cols_neg[:5]}")
