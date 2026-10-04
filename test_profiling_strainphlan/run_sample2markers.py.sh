#!/bin/bash
set -euo pipefail

input_dir="/mnt/ME5012-Vol01/chenc/MetaphIan/paper-data/7.test_profiling_strainphlan/PRJDB4176-colon_cancer-strainphlan"

marker_dir="${input_dir}/markers"

mkdir -p "$marker_dir"

for sam_file in "$input_dir"/*.sam.bz2; do

    sample_name="$(basename "$sam_file" .sam.bz2)"

    echo "Processing: $sample_name"

    sample2markers.py -i "$sam_file" -o "$marker_dir"  -d /mnt/ME5012-Vol01/chenc/MetaphIan/7-seventh-all-maker-gene-db/20251219-filtered-gene/GMB-longer-first-20260123_NEW.pkl -n 100



done

echo "All sample2markers analyses completed."