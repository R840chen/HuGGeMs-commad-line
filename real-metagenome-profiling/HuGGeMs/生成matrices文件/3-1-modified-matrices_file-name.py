#!/usr/bin/env python3

import os
import glob
import pandas as pd

# -----------------------------
# Configuration
# -----------------------------

# Input directory (The output from the previous Step 3)
INPUT_DIR = "/ME5012-Vol01/chenc/real-metagenomes-validate/all-lab-real-metagenomes-analysis/GMB-all-metaphlan-result/all-GMB-result-20260124/3-modified-matrices-file"

# Output directory (Where cleaned files will be saved)
OUTPUT_DIR = "/ME5012-Vol01/chenc/real-metagenomes-validate/all-lab-real-metagenomes-analysis/GMB-all-metaphlan-result/all-GMB-result-20260124/4-modified-txt-file-name"

# The prefix string to remove from sample names
PREFIX_TO_REMOVE = "mpa4-lab-downlaod-sample-result-"

def clean_and_rename(file_path, output_dir):
    """
    Reads a TSV, removes specific prefix from the index (sample names),
    renames the index header to 'ID', and saves the file.
    """
    filename = os.path.basename(file_path)
    
    try:
        # 1. Read the TSV file
        # Assuming the first column is the Sample Index
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        
        # 2. Clean the Sample Names (Index)
        # Using string replace to remove the prefix
        if df.index.dtype == 'object':
            df.index = df.index.str.replace(PREFIX_TO_REMOVE, "", regex=False)
            # Also remove .txt if it exists in the sample name
            df.index = df.index.str.replace(".txt", "", regex=False)
        
        # 3. Rename the Index Header to "ID"
        df.index.name = "ID"
        
        # 4. Save to new directory
        output_path = os.path.join(output_dir, filename)
        df.to_csv(output_path, sep="\t")
        
        print(f"[DONE] Processed: {filename}")
        # Optional: Print first few index names to verify
        # print(f"       First sample: {df.index[0]}")

    except Exception as e:
        print(f"[ERROR] Failed to process {filename}: {e}")

def main():
    # Ensure output directory exists
    if not os.path.exists(OUTPUT_DIR):
        print(f"Creating output directory: {OUTPUT_DIR}")
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Find all TSV files in the input directory
    input_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.tsv")))

    if not input_files:
        print(f"No TSV files found in {INPUT_DIR}")
        return

    print(f"Found {len(input_files)} files. Starting cleanup...\n")

    for file_path in input_files:
        clean_and_rename(file_path, OUTPUT_DIR)

    print("\nAll files processed and saved to 'modified-txt-file-name'.")

if __name__ == "__main__":
    main()