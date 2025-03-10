import pandas as pd
import os
from tabulate import tabulate

# 📂 File Paths
raw_file = "raw.csv"
snp_file = "data/uniq_snips.csv"

# 📌 Load raw DNA data
if not os.path.exists(raw_file):
    raise FileNotFoundError(f"File not found: {raw_file}")

df = pd.read_csv(raw_file, sep=',', dtype=str, comment='#', low_memory=False)

# 🔍 Normalize Column Names
df.columns = df.columns.str.strip().str.lower()

# ✅ FIX: Rename 'result' to 'genotype'
df.rename(columns={'result': 'genotype'}, inplace=True)

# Ensure required columns exist
required_cols = {'rsid', 'chromosome', 'position', 'genotype'}
if not required_cols.issubset(df.columns):
    raise KeyError(f"Missing columns in raw.csv. Found: {df.columns}, Expected: {required_cols}")

# 📌 Load SNPedia Database
if not os.path.exists(snp_file):
    raise FileNotFoundError(f"File not found: {snp_file}")

snp_df = pd.read_csv(snp_file)

# 🔍 Normalize Column Names
snp_df.columns = snp_df.columns.str.strip().str.lower()
print(f"Columns in snp_df after renaming: {snp_df.columns}")

# ✅ Fix Genotype Extraction from SNPedia
snp_df['genotype'] = snp_df['rsid'].str.extract(r'\(([AGCT]);([AGCT])\)')[0] + \
                      snp_df['rsid'].str.extract(r'\(([AGCT]);([AGCT])\)')[1]

# ✅ Ensure rsid and genotype formats match exactly
df['rsid'] = df['rsid'].str.lower().str.strip()
snp_df['rsid'] = snp_df['rsid'].str.replace(r'\(.*\)', '', regex=True).str.lower().str.strip()


# ✅ Convert Genotype Format ("A;A" → "AA") to Match SNPedia
df['genotype'] = df['genotype'].str.replace(';', '', regex=True).str.upper().str.strip()

snp_df['genotype'] = snp_df['genotype'].str.upper().str.strip()

# 🛠 Debugging: Print sample data to check why merge is failing
print("\n🔍 First few rows of raw.csv (Your DNA Data - After Fix):")
print(df[['rsid', 'genotype']].head())

print("\n🔍 First few rows of SNPedia Database:")
print(snp_df[['rsid', 'genotype']].head())

print("\n🔍 Unique Genotypes in raw.csv (After Fix):")
print(df['genotype'].unique())

print("\n🔍 Unique Genotypes in SNPedia (uniq_snips.csv):")
print(snp_df['genotype'].unique())

# ✅ Merge raw DNA data with SNPedia using rsid and genotype
new_df = snp_df.merge(df, how='inner', on=['rsid', 'genotype'])
print("\n✅ Rows After Merge:", len(new_df))

# 🚦 Check for Good and Bad Genes
gene_type = "Good"  # Change to "Good" if needed

filtered_genes = new_df[new_df['repute'].str.strip().str.lower() == gene_type.lower()]

# 🔗 Generate SNPedia URLs
base_url = 'https://www.snpedia.com/index.php/'
gene_urls = [f"{base_url}{rsid}" for rsid in filtered_genes['rsid']]

# 📌 Print Output
if filtered_genes.empty:
    print(f"No '{gene_type}' genes found. Double-check rsid and genotype matching.")
else:
    for url in gene_urls:
        print(url)

    print("\n" + "="*40 + "\n")
    print(tabulate(filtered_genes[['rsid', 'genotype', 'repute', 'summary']], headers='keys', tablefmt='psql'))
