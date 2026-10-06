# Positional distribution of marker genes within genomes

Analysis behind the Discussion paragraph that answers the reviewer's question of
**where the marker genes sit in a genome** and whether the **replication-associated
copy-number gradient** (more DNA copies near the origin than near the terminus in a
replicating cell) could bias coverage-based abundance estimates.

The marker set of each of two species (200 markers each, from the HuGGeMs database)
was mapped onto three publicly available complete genomes: two strains of *Waltera
acetigignens* and one *Dysosmobacter hominis*.

---

## 1. Requirements

Python 3.9+ with `numpy` and `openpyxl`:

```bash
pip install numpy openpyxl
```

No read mapper, no aligner package and no C compiler are needed (see §5).

## 2. Input files

Put the five FASTA files in `data/`, or point `MARKER_DATA_DIR` at wherever they are.
They are **not** distributed in this repository.

| File | Source |
|---|---|
| `Brotolimicola_acetigignens_species_2981769.ffn` | marker set of *Waltera acetigignens*, HuGGeMs database (Supplementary data) |
| `Dysosmobacter_hominis_species_2763041.ffn` | marker set of *Dysosmobacter hominis*, HuGGeMs database (Supplementary data) |
| `GCA_025567015.1_ASM2556701v1_genomic.fna` | NCBI, *W. acetigignens* strain 1 |
| `GCA_020687025.1_ASM2068702v1_genomic.fna` | NCBI, *W. acetigignens* strain 2 |
| `GCA_014297375.1_ASM1429737v1_genomic.fna` | NCBI, *D. hominis* |

The three genomes are public and are downloaded by accession:

```bash
bash fetch_genomes.sh          # downloads the three genomes into ./data
```

The first marker set keeps its file name as exported from the HuGGeMs database; the
species was reclassified from *Brotolimicola acetigignens* to *Waltera acetigignens*
(PMID 38722771). The scripts refer to the files by name, so the names are kept as is.

## 3. Scripts (run in this order)

| Script | Purpose |
|---|---|
| `mapper.py` | in-house seed-and-extend aligner (module imported by the others) |
| `02_map_markers.py` | maps the markers to the three genomes → `marker_records.json` |
| `03_position_stats.py` | loci count, spacing, window density, uniformity tests |
| `04_oric_axis_bias.py` | reconstructs an oriC→ter axis and estimates the replication bias |
| `05_gc_and_skew.py` | marker vs genome GC content, local GC skew (diagnostic) |
| `06_final_checks.py` | aligner self-check and controls (diagnostic) |

## 4. How to run

```bash
cd marker_position_analysis

# optional: where the .ffn / .fna files are, and where results are written
export MARKER_DATA_DIR=/path/to/input_fastas     # default: ./data
export MARKER_OUT_DIR=/path/to/output            # default: this directory

python 02_map_markers.py      # ~2 min for 400 markers x 3 genomes
python 03_position_stats.py
python 04_oric_axis_bias.py
python 05_gc_and_skew.py
python 06_final_checks.py
```

## 5. Notes

- **Why a custom aligner.** Neither `blastn` nor `minimap2` was available on the machine
  used, and no C compiler was available either (`pip install mappy|edlib|parasail` all
  failed to build). `mapper.py` therefore implements a minimal nucleotide aligner in pure
  numpy: a k = 15 k-mer index of the genome (both strands), seed clustering into candidate
  start positions, gapless extension with a ±8 bp shift search, and reporting of the best
  hit. Accuracy self-check: the markers whose own source genome is among the three analysed
  are all recovered from that same genome with an identity of exactly 1.0000.
- **Limitations.** All three genomes are drafts (36–125 contigs, N50 70–300 kb). The
  oriC→ter axis is reconstructed from cumulative GC skew (covering 89–97% of each genome)
  and is an approximation; the dispersion results (loci count, spacing, window density) do
  not depend on it. The copy-number-profile figures indicate the *order of magnitude* of
  the bias, not a measured value. Only two species / three genomes were examined, so the
  results should not be extrapolated to the whole database. The argument in the manuscript
  rests on the marker set being dispersed and, above all, on the fact that the **same marker
  set is used for every sample**, so any residual offset is a species-specific constant that
  cancels in case–control contrasts.
