#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

if len(sys.argv) != 5:
    sys.exit(
        "Usage: python sylph_add_ncbi_taxid.py "
        "<input.sylphmpa> <assembly_summary_genbank.txt> "
        "<assembly_summary_refseq.txt> <output.tsv>"
    )

sylph_file = sys.argv[1]
genbank_file = sys.argv[2]
refseq_file = sys.argv[3]
output_file = sys.argv[4]

accession2taxid = {}

def load_assembly_summary(path):
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            accession = cols[0]
            taxid = cols[6]
            species_taxid = cols[7]
            accession2taxid[accession] = (taxid, species_taxid)

load_assembly_summary(genbank_file)
load_assembly_summary(refseq_file)

print(f"[INFO] Loaded {len(accession2taxid)} assemblies")

with open(sylph_file) as fin, open(output_file, "w") as fout:
    fout.write(
        "clade_name\trelative_abundance\tsequence_abundance\t"
        "ANI\tCoverage\tassembly_accession\ttaxid\tspecies_taxid\n"
    )

    for line in fin:
        if line.startswith("#") or line.startswith("clade_name"):
            continue

        cols = line.rstrip("\n").split("\t")
        clade = cols[0]

        accession = "NA"
        taxid = "NA"
        species_taxid = "NA"


        if "|t__" in clade:
            accession = clade.split("|t__")[-1]
            if accession in accession2taxid:
                taxid, species_taxid = accession2taxid[accession]

        fout.write(
            line.rstrip("\n")
            + f"\t{accession}\t{taxid}\t{species_taxid}\n"
        )

print(f"[DONE] Output written to: {output_file}")
