#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re

# ==============================
# Input / Output files
# ==============================
matrix_file = "metaphlan_species_abundance_matrix.prevalence5pct-new.tsv"
sample_info_file = "2174-sample-info.tsv"
output_file = "metaphlan_species_abundance_matrix.prevalence5pct.filtered-new.tsv"

# ==============================
# 1. Load MetaPhlAn abundance matrix
# Rows: species
# Columns: samples
# ==============================
df = pd.read_csv(matrix_file, sep="\t", index_col=0)

original_columns = df.columns.tolist()

# ==============================
# 2. Remove MetaPhlAn sample name prefix
# Example:
# mpa4-lab-downlaod-sample-result-DRR127476_1
# -> DRR127476_1
# ==============================
prefix = "mpa4-lab-downlaod-sample-result-"
clean_columns = [col.replace(prefix, "") for col in original_columns]
df.columns = clean_columns

# ==============================
# 3. Load sample metadata
# SampleName column contains IDs like: DRR127476
# ==============================
sample_info = pd.read_csv(sample_info_file, sep="\t")
valid_sample_ids = set(sample_info["SampleName"].astype(str))

# ==============================
# 4. Keep only samples listed in metadata
# Match rule:
# DRR127476_1 / DRR127476_2  -> DRR127476
# ==============================
def keep_sample(column_name):
    """
    column_name format: DRR127476_1 or DRR127476_2
    """
    base_id = re.sub(r"_[12]$", "", column_name)
    return base_id in valid_sample_ids

kept_columns = [col for col in df.columns if keep_sample(col)]

df_filtered = df[kept_columns]

# ==============================
# 5. Write filtered matrix
# ==============================
df_filtered.to_csv(output_file, sep="\t")

# ==============================
# 6. Summary
# ==============================
print(f"Original number of samples: {len(original_columns)}")
print(f"Samples after prefix removal: {len(clean_columns)}")
print(f"Samples kept after metadata filtering: {len(kept_columns)}")
print(f"Filtered matrix written to: {output_file}")
