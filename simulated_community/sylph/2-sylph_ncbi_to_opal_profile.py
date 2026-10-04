#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert Sylph MetaPhlAn-style output with NCBI taxids
into an OPAL-compatible taxonomy profile.

ONLY species-level entries are kept.
"""

import sys
from collections import defaultdict

if len(sys.argv) != 4:
    sys.exit(
        "Usage: python sylph_ncbi_to_opal_taxonomy_profile.py "
        "<input.sylph.ncbi.tsv> <sample_id> <output.opal.tsv>"
    )

input_file = sys.argv[1]
sample_id = sys.argv[2]
output_file = sys.argv[3]

rank_map = {
    "s": "species"
}

# (taxid, rank, taxpath, taxpathsn) -> abundance
profile = defaultdict(float)

with open(input_file, encoding="utf-8") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(header)}

    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        clade = fields[idx["clade_name"]]
        abundance = float(fields[idx["relative_abundance"]])
        taxid = fields[idx["taxid"]]

        if not taxid.isdigit():
            continue

        parts = clade.split("|")

        for part in parts:
            if not part.startswith("s__"):
                continue

            species_name = part.split("__", 1)[1]

            profile[(taxid, "species", taxid, species_name)] += abundance

# Write OPAL file
with open(output_file, "w", encoding="utf-8") as out:
    out.write(f"@SampleID:{sample_id}\n")
    out.write("@Version:0.9.3\n")
    out.write("@Ranks:species\n")
    out.write("@TaxonomyID:NCBI\n")
    out.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")

    for (taxid, rank, taxpath, taxpathsn), abund in sorted(
        profile.items(), key=lambda x: -x[1]
    ):
        out.write(f"{taxid}\t{rank}\t{taxpath}\t{taxpathsn}\t{abund}\n")

print(f"[DONE] OPAL species-only taxonomy profile written to: {output_file}")
