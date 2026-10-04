#!/usr/bin/env bash
# Download the three genomes analysed in this folder from NCBI into ./data
#
#   bash fetch_genomes.sh
#
# The URL pattern follows the NCBI FTP layout
#   https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/xxx/yyy/zzz/<acc>_<asm>/<acc>_<asm>_genomic.fna.gz
# If a link ever stops working, download the same accessions through
#   https://www.ncbi.nlm.nih.gov/datasets/genome/   (or `datasets download genome accession <GCA_...>`)
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data

FTP=https://ftp.ncbi.nlm.nih.gov/genomes/all

fetch () {                    # $1 = path under genomes/all, $2 = output file name in ./data
  local rel="$1" out="data/$2"
  if [ -s "$out" ]; then
    echo "already present: $out"
    return
  fi
  echo "downloading $2 ..."
  curl -fsSL "$FTP/$rel.gz" | gunzip -c > "$out"
}

fetch GCA/025/567/015/GCA_025567015.1_ASM2556701v1/GCA_025567015.1_ASM2556701v1_genomic.fna \
      GCA_025567015.1_ASM2556701v1_genomic.fna
fetch GCA/020/687/025/GCA_020687025.1_ASM2068702v1/GCA_020687025.1_ASM2068702v1_genomic.fna \
      GCA_020687025.1_ASM2068702v1_genomic.fna
fetch GCA/014/297/375/GCA_014297375.1_ASM1429737v1/GCA_014297375.1_ASM1429737v1_genomic.fna \
      GCA_014297375.1_ASM1429737v1_genomic.fna

echo
echo "Done. The two marker sets"
echo "  Brotolimicola_acetigignens_species_2981769.ffn"
echo "  Dysosmobacter_hominis_species_2763041.ffn"
echo "are part of the HuGGeMs database / the Supplementary data of the manuscript;"
echo "copy them into ./data as well, then run the pipeline (see README.md)."
