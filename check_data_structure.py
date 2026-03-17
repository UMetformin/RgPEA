import pandas as pd
import numpy as np

print("Checking data structure...\n")

# Load metabolomics data
print("="*80)
print("NEGATIVE ION MODE METABOLITES")
print("="*80)
metab_neg = pd.read_excel(r"C:\Users\ASUS\Desktop\Gut\Metabolites quantification profiles in negative ion modes.xlsx")
print(f"\nShape: {metab_neg.shape}")
print(f"\nFirst few rows:")
print(metab_neg.head())
print(f"\nColumn names (first 10):")
print(list(metab_neg.columns[:10]))
print(f"\nData types:")
print(metab_neg.dtypes.head(10))

print("\n" + "="*80)
print("POSITIVE ION MODE METABOLITES")
print("="*80)
metab_pos = pd.read_excel(r"C:\Users\ASUS\Desktop\Gut\Metabolites quantification profiles in positive ion modes.xlsx")
print(f"\nShape: {metab_pos.shape}")
print(f"\nFirst few rows:")
print(metab_pos.head())
print(f"\nColumn names (first 10):")
print(list(metab_pos.columns[:10]))

# Search for phenylethylamine
print("\n" + "="*80)
print("SEARCHING FOR PHENYLETHYLAMINE")
print("="*80)

all_metab_cols = list(metab_neg.columns) + list(metab_pos.columns)
pea_matches = []

for col in all_metab_cols:
    col_str = str(col).lower()
    if 'phenyl' in col_str or 'ethylamine' in col_str or 'pea' in col_str:
        pea_matches.append(col)

print(f"\nMatches found: {len(pea_matches)}")
if pea_matches:
    for match in pea_matches[:20]:
        print(f"  - {match}")
else:
    print("\nNo exact matches. Showing sample metabolite names:")
    print("\nFrom negative mode (first 20):")
    for col in list(metab_neg.columns[1:21]):
        print(f"  - {col}")
    print("\nFrom positive mode (first 20):")
    for col in list(metab_pos.columns[1:21]):
        print(f"  - {col}")

# Check if first column is identifier
print("\n" + "="*80)
print("CHECKING INDEX COLUMN")
print("="*80)
print(f"\nFirst column name (negative): {metab_neg.columns[0]}")
print(f"First column values (negative):")
print(metab_neg.iloc[:5, 0])

print(f"\nFirst column name (positive): {metab_pos.columns[0]}")
print(f"First column values (positive):")
print(metab_pos.iloc[:5, 0])
