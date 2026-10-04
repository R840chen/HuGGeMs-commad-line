#!/usr/bin/env python3

import pandas as pd

# -----------------------------
# Parameters
# -----------------------------
input_file = "metaphlan_species_abundance_matrix-new.tsv"
output_file = "metaphlan_species_abundance_matrix.prevalence5pct-new.tsv"
prevalence_threshold = 0.05   # 5%

# -----------------------------
# Load data
# -----------------------------
print("Loading matrix...")
df = pd.read_csv(input_file, sep="\t", index_col=0)

n_samples = df.shape[1]
min_samples = int(n_samples * prevalence_threshold)

print(f"Total samples: {n_samples}")
print(f"Prevalence threshold: {prevalence_threshold*100:.1f}% "
      f"({min_samples} samples)")

# -----------------------------
# Prevalence filtering
# -----------------------------
print("Calculating prevalence...")
prevalence = (df > 0).sum(axis=1)

filtered_df = df.loc[prevalence >= min_samples]

# -----------------------------
# Output
# -----------------------------
filtered_df.to_csv(output_file, sep="\t")

print("Filtering completed.")
print(f"Original species: {df.shape[0]}")
print(f"Remaining species: {filtered_df.shape[0]}")
print(f"Saved: {output_file}")
