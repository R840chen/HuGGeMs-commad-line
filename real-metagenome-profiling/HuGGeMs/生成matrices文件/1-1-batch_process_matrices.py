#!/usr/bin/env python3

import os
import glob
import pandas as pd
import shutil

# --- CONFIGURATION ---

# 1. Base directory containing all Project folders
BASE_DIR = "/ME5012-Vol01/chenc/real-metagenomes-validate/all-lab-real-metagenomes-analysis/GMB-all-metaphlan-result/all-GMB-result-20260124/separate-project-txt-file/left"

# 2. Central directory to collect all result matrices
CENTRAL_OUTPUT_DIR = "/ME5012-Vol01/chenc/real-metagenomes-validate/all-lab-real-metagenomes-analysis/GMB-all-metaphlan-result/all-GMB-result-20260124/2-matrices-file"

# ---------------------

def process_single_project(project_path, project_name):
    """
    Process a single project folder: 
    1. Read txt files and generate an abundance matrix.
    2. Save the matrix to the project folder.
    3. Save a copy of the matrix to the central output folder.
    """
    
    # Find txt files matching the specific naming pattern
    txt_pattern = os.path.join(project_path, "mpa4-lab-downlaod-sample-result-*.txt")
    txt_files = sorted(glob.glob(txt_pattern))

    if not txt_files:
        print(f"  [SKIP] No matching txt files found in: {project_name}")
        return

    print(f"  Processing {len(txt_files)} files in: {project_name}")

    species_abundance = {}

    for txt_file in txt_files:
        # Get Sample name from filename
        sample_name = os.path.basename(txt_file).replace(".txt", "")
        species_abundance[sample_name] = {}

        try:
            with open(txt_file, "r", errors="ignore") as f:
                for line in f:
                    line = line.strip()

                    # Skip empty lines or headers
                    if not line or line.startswith("#"):
                        continue

                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue

                    clade = parts[0]
                    try:
                        abundance = float(parts[2])
                    except ValueError:
                        continue

                    # ? ONLY use species level (s__)
                    if "|s__" not in clade:
                        continue

                    # Extract species name
                    species = clade.split("|s__")[-1]
                    
                    # Remove possible trailing info after species name
                    species = species.split("|")[0]

                    if species == "" or species.lower().startswith("unclassified"):
                        continue

                    species_abundance[sample_name][species] = abundance
        
        except Exception as e:
            print(f"    [ERROR] Reading file {txt_file}: {e}")

    # Build DataFrame
    if not species_abundance:
        print(f"  [WARN] No valid species data extracted for {project_name}")
        return

    df = pd.DataFrame.from_dict(species_abundance, orient="columns")
    df.index.name = "Species"
    df = df.fillna(0)

    # --- SAVE OUTPUTS ---
    
    output_filename = f"{project_name}_species_abundance_matrix.tsv"

    # Location 1: Inside the Project Folder (Local)
    local_output_path = os.path.join(project_path, output_filename)
    df.to_csv(local_output_path, sep="\t")
    
    # Location 2: Inside the Central Folder (Aggregation)
    central_output_path = os.path.join(CENTRAL_OUTPUT_DIR, output_filename)
    df.to_csv(central_output_path, sep="\t")

    print(f"  [DONE] Saved local:   {local_output_path}")
    print(f"  [DONE] Saved central: {central_output_path}")


def main():
    # Check if base directory exists
    if not os.path.exists(BASE_DIR):
        print(f"Error: Base directory not found: {BASE_DIR}")
        return

    # Create the central output directory if it doesn't exist
    if not os.path.exists(CENTRAL_OUTPUT_DIR):
        print(f"Creating central output directory: {CENTRAL_OUTPUT_DIR}")
        os.makedirs(CENTRAL_OUTPUT_DIR, exist_ok=True)

    print(f"Scanning directories in: {BASE_DIR}\n")

    # Iterate through each subdirectory in BASE_DIR
    subdirs = sorted(os.listdir(BASE_DIR))
    
    for item in subdirs:
        full_path = os.path.join(BASE_DIR, item)
        
        # Only process directories
        if os.path.isdir(full_path):
            print(f"-> Found Project: {item}")
            process_single_project(full_path, item)
            print("-" * 50)

    print("\nAll projects processed and aggregated.")

if __name__ == "__main__":
    main()