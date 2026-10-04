#!/usr/bin/env python3

import os
import glob
import pandas as pd

# -----------------------------
# Configuration
# -----------------------------

# Input directory (Result from the previous step)
INPUT_DIR = "/ME5012-Vol01/chenc/real-metagenomes-validate/all-lab-real-metagenomes-analysis/GMB-all-metaphlan-result/all-GMB-result-20260124/2-matrices-file"

# Output directory (Where modified files will be saved)
OUTPUT_DIR = "/ME5012-Vol01/chenc/real-metagenomes-validate/all-lab-real-metagenomes-analysis/GMB-all-metaphlan-result/all-GMB-result-20260124/3-modified-matrices-file/0.1-prevalence"

# Threshold for filtering (0.05 = 5%, 0.1 = 10%)
# MaAsLin 3 typically recommends at least 0.1 (10%) to reduce sparsity/noise.
PREVALENCE_THRESHOLD = 0.1 

def process_matrix(file_path, output_dir, threshold):
    """
    Reads a Species x Samples matrix, filters by prevalence, 
    transposes it to Samples x Species (MaAsLin style), and saves it.
    """
    filename = os.path.basename(file_path)
    print(f"Processing: {filename}")

    try:
        # 1. Load Data
        # Assuming current format is Index=Species, Columns=Samples
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        
        # Check if the dataframe is empty
        if df.empty:
            print(f"  [WARN] File is empty. Skipping.")
            return

        n_species_initial = df.shape[0]
        n_samples = df.shape[1]
        
        # Calculate minimum sample count required
        min_samples = int(n_samples * threshold)

        print(f"  - Input shape: {df.shape} (Species x Samples)")
        print(f"  - Prevalence cutoff: {threshold*100}% ({min_samples} samples)")

        # 2. Prevalence Filtering
        # Logic: Calculate how many samples have abundance > 0 for each species (row)
        prevalence_counts = (df > 0).sum(axis=1)
        
        # Keep rows (species) that meet the threshold
        filtered_df = df.loc[prevalence_counts >= min_samples]
        
        n_species_final = filtered_df.shape[0]
        print(f"  - Filtering: {n_species_initial} -> {n_species_final} species retained")

        if n_species_final == 0:
            print(f"  [WARN] No species passed the filter! Skipping output.")
            return

        # 3. Transpose for MaAsLin 3
        # Current: Rows=Species, Cols=Samples
        # Target:  Rows=Samples, Cols=Species
        final_df = filtered_df.T
        
        # Ensure the index name is typically 'Sample' or 'ID' for MaAsLin
        final_df.index.name = "SampleID"

        # 4. Save Output
        # Construct new filename
        base_name = os.path.splitext(filename)[0] # remove .tsv
        new_filename = f"{base_name}.prev{int(threshold*100)}pct.maaslin_ready.tsv"
        output_path = os.path.join(output_dir, new_filename)

        final_df.to_csv(output_path, sep="\t")
        print(f"  - Saved to: {output_path}")
        print(f"  - Final shape: {final_df.shape} (Samples x Species)")

    except Exception as e:
        print(f"  [ERROR] Failed to process {filename}: {e}")

def main():
    # Ensure output directory exists
    if not os.path.exists(OUTPUT_DIR):
        print(f"Creating output directory: {OUTPUT_DIR}")
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Find all TSV files in the input directory
    # Adjust pattern if your files end with .txt or .csv
    input_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*_species_abundance_matrix.tsv")))

    if not input_files:
        print(f"No matching files found in {INPUT_DIR}")
        return

    print(f"Found {len(input_files)} matrix files to process.\n")
    print("-" * 50)

    for file_path in input_files:
        process_matrix(file_path, OUTPUT_DIR, PREVALENCE_THRESHOLD)
        print("-" * 50)

    print("\nBatch processing completed.")

if __name__ == "__main__":
    main()