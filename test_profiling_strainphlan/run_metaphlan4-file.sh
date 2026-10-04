#!/bin/bash
set -euo pipefail

########################################
# Input / output paths
########################################

input_root="/mnt/ME5012-Vol01/lisc/cohort/PRJDB4176-colon_cancer/1"

output_base="/mnt/ME5012-Vol01/chenc/MetaphIan/paper-data/7.test_profiling_strainphlan/PRJDB4176-colon_cancer-strainphlan"

mkdir -p "$output_base"

########################################
# Resource parameters
########################################

threads=60
parallel_jobs=3

########################################
# Command list
########################################

cmd_file="${output_base}/run_PRJDB4176_metaphlan.txt"
: > "$cmd_file"

########################################
# MetaPhlAn / Bowtie2 DB
########################################

bowtie2db="/mnt/ME5012-Vol01/chenc/MetaphIan/7-seventh-all-maker-gene-db/20251219-filtered-gene"
index_name="GMB-longer-first-20260123_NEW"

########################################
# Find paired-end FASTQ files
########################################

echo "Searching FASTQ files under:"
echo "$input_root"

while IFS= read -r -d '' r1_file; do

    # Infer R2 filename
    r2_file="${r1_file/_1.fastq.gz/_2.fastq.gz}"

    # Check whether R2 exists
    if [ ! -f "$r2_file" ]; then
        echo "[WARN] R2 file not found, skipping:"
        echo "       $r1_file"
        continue
    fi

    ####################################
    # Sample name
    ####################################

    base_name="$(basename "$r1_file")"
    base_name="${base_name%%_1.fastq.gz}"

    ####################################
    # Output files
    ####################################

    output_file="${output_base}/mpa4-${base_name}.txt"

    bowtie_file="${output_base}/${base_name}.bowtie2out.txt"

    # SAM output for StrainPhlAn
    sam_file="${output_base}/${base_name}.sam.bz2"

    ####################################
    # Generate MetaPhlAn command
    ####################################

    printf '%s\n' \
"metaphlan \
\"$r1_file\",\"$r2_file\" \
--input_type fastq \
-o \"$output_file\" \
--bowtie2out \"$bowtie_file\" \
-s \"$sam_file\" \
--nproc $threads \
--bowtie2db \"$bowtie2db\" \
--force \
--index $index_name" >> "$cmd_file"

done < <(
    find "$input_root" \
    -type f \
    -name "*_1.fastq.gz" \
    -print0
)

########################################
# Execute
########################################

echo "Command list generated:"
echo "$cmd_file"

task_count=$(wc -l < "$cmd_file")

echo "Total number of tasks: $task_count"

if [ "$task_count" -eq 0 ]; then
    echo "[ERROR] No FASTQ files found."
    exit 1
fi

echo "Running MetaPhlAn with parallel -j $parallel_jobs"

parallel -j "$parallel_jobs" < "$cmd_file"

echo "All MetaPhlAn analyses completed."