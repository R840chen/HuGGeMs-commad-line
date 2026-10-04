#!/usr/bin/env python3

import os
import glob
import pandas as pd

INPUT_DIR = "./all-txt-result"
OUTPUT_FILE = "metaphlan_species_abundance_matrix-new.tsv"

species_abundance = {}

txt_files = sorted(
    glob.glob(os.path.join(INPUT_DIR, "mpa4-lab-downlaod-sample-result-*.txt"))
)

print(f"Found {len(txt_files)} MetaPhlAn files")

for txt_file in txt_files:
    sample_name = os.path.basename(txt_file).replace(".txt", "")
    species_abundance[sample_name] = {}

    with open(txt_file, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()

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

            species = clade.split("|s__")[-1]

            # remove possible trailing info after species name
            species = species.split("|")[0]

            if species == "" or species.lower().startswith("unclassified"):
                continue

            species_abundance[sample_name][species] = abundance

print("Building species abundance matrix...")

df = pd.DataFrame.from_dict(species_abundance, orient="columns")
df.index.name = "Species"
df = df.fillna(0)

df.to_csv(OUTPUT_FILE, sep="\t")

print("Saved:", OUTPUT_FILE)
print("Matrix shape:", df.shape)
