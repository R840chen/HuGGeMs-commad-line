#!/usr/bin/env python3

import pandas as pd
import numpy as np

# -----------------------------
# Parameters
# -----------------------------
input_file = "metaphlan_species_abundance_matrix.prevalence5pct.filtered-new.tsv"
output_file = "metaphlan_species_abundance_matrix.prevalence5pct.filtered.CLR-new.tsv"
pseudocount = 1e-6

# -----------------------------
# Load data
# -----------------------------
print("Loading abundance matrix...")
df = pd.read_csv(input_file, sep="\t", index_col=0)

print(f"Matrix shape: {df.shape}")

# -----------------------------
# CLR transformation
# -----------------------------
print("Applying CLR transformation...")

# Add pseudocount
df_pc = df + pseudocount

# Log transform
log_df = np.log(df_pc)

# Subtract column-wise geometric mean (i.e. mean of logs)
clr_df = log_df.sub(log_df.mean(axis=0), axis=1)

# -----------------------------
# Output
# -----------------------------
clr_df.to_csv(output_file, sep="\t")

print("CLR transformation completed.")
print(f"Saved: {output_file}")
