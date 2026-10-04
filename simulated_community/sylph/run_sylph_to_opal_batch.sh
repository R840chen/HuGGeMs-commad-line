#!/usr/bin/env bash
set -euo pipefail

# ========= CONFIG =========

# Sylph output files
INPUT_PATTERN="3-synthetic-community2017.12.04_18.45.54_sample_*-anonymous_reads.fq.sylphmpa"

# NCBI assembly summary files
ASSEMBLY_GENBANK="/ME4012-Vol01/chenc/MetaphIan/gut_bac-catalogue/NCBI/assembly_summary_genbank.txt"
ASSEMBLY_REFSEQ="/ME4012-Vol01/chenc/MetaphIan/gut_bac-catalogue/NCBI/assembly_summary_refseq.txt"

# Scripts
ADD_TAXID_SCRIPT="1-sylph_add_ncbi_taxid.py"
TO_OPAL_SCRIPT="2-sylph_ncbi_to_opal_profile.py"

# Output directories
NCBI_TSV_DIR="01.sylph_ncbi_tsv"
OPAL_DIR="02.opal_species"

mkdir -p "${NCBI_TSV_DIR}" "${OPAL_DIR}"

# ========= MAIN LOOP =========

for sylph_file in ${INPUT_PATTERN}; do
    echo "=============================================="
    echo "[INFO] Processing: ${sylph_file}"

    # basename without suffix
    base=$(basename "${sylph_file}" .sylphmpa)

    ncbi_tsv="${NCBI_TSV_DIR}/${base}.sylph.ncbi.tsv"
    opal_tsv="${OPAL_DIR}/${base}.species.opal.tsv"

    # Sample ID for OPAL (recommended: same as basename)
    sample_id="${base}"

    echo "[STEP 1] Add NCBI taxid"
    python "${ADD_TAXID_SCRIPT}" \
        "${sylph_file}" \
        "${ASSEMBLY_GENBANK}" \
        "${ASSEMBLY_REFSEQ}" \
        "${ncbi_tsv}"

    echo "[STEP 2] Convert to OPAL (species-only)"
    python "${TO_OPAL_SCRIPT}" \
        "${ncbi_tsv}" \
        "${sample_id}" \
        "${opal_tsv}"

    echo "[DONE] ${opal_tsv}"
done

echo "=============================================="
echo "[ALL DONE] Sylph ¡ú NCBI ¡ú OPAL pipeline finished."
