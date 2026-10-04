# Positional distribution of HuGGeMs marker genes and replication-associated bias

This folder contains the analysis behind the paragraph in the Discussion that
addresses the reviewer's comment on **where marker genes sit in a genome** and
whether the **replication-associated copy-number gradient** (more DNA copies
near the origin than near the terminus in a replicating cell) could bias
coverage-based abundance estimates.

Two SGBs were examined:

| SGB | Species | Marker set | Genomes analysed |
|---|---|---|---|
| SGB_2981769 | *Waltera acetigignens* (formerly *Brotolimicola acetigignens*, PMID 38722771) | 200 markers | GCA_025567015.1, GCA_020687025.1 |
| SGB_2763041 | *Dysosmobacter hominis* | 200 markers | GCA_014297375.1 |

---

## 1. Requirements

- Python 3.9+ with `numpy` and `openpyxl`
- No aligner installation is needed (see §5)

```bash
pip install numpy openpyxl
```

## 2. Input files

Put the five input files in `data/` (or point `MARKER_DATA_DIR` at wherever
they are; see §4). They are **not** distributed in this repository.

| File | Source |
|---|---|
| `Brotolimicola_acetigignens_species_2981769.ffn` | marker set of SGB_2981769, HuGGeMs database (Supplementary data of the manuscript) |
| `Dysosmobacter_hominis_species_2763041.ffn` | marker set of SGB_2763041, HuGGeMs database (Supplementary data of the manuscript) |
| `GCA_025567015.1_ASM2556701v1_genomic.fna` | NCBI, *W. acetigignens* strain 1 |
| `GCA_020687025.1_ASM2068702v1_genomic.fna` | NCBI, *W. acetigignens* strain 2 |
| `GCA_014297375.1_ASM1429737v1_genomic.fna` | NCBI, *D. hominis* |

The three genomes are public and can be downloaded by accession, e.g.

```bash
bash fetch_genomes.sh          # downloads the three genomes into ./data
```

The marker `.ffn` files keep the original file name of SGB_2981769, which
still uses the pre-2024 name; the species was reclassified as
*Waltera acetigignens*.

## 3. Scripts (run in this order)

| Script | Purpose |
|---|---|
| `mapper.py` | in-house seed-and-extend aligner (module imported by the others) |
| `02_map_markers.py` | maps all 400 markers to all 3 genomes → `marker_records.json` |
| `03_position_stats.py` | loci count, spacing, window density, KS uniformity test, GC skew |
| `04_oric_axis_bias.py` | rebuilds an oriC→ter axis from GC skew and estimates the abundance bias |
| `05_gc_and_skew.py` | marker vs genome GC content, local GC skew relative to contig background |
| `06_final_checks.py` | aligner accuracy self-check, strain- vs species-level markers, multi-copy markers, cross-species control, annotation composition |
| `07_make_tableS21.py` | writes Supplementary Table S21 |

`05_gc_and_skew.py` and `06_final_checks.py` are diagnostic steps; they produce
the numbers used in the Supplementary Table S21 summary and in the
"correctness self-check" below.

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
python 07_make_tableS21.py
```

## 5. Why a custom aligner

The analysis started from a marker set and complete genome sequences, i.e. the
task is "where in this genome does this ~1.6–4.6 kb marker sit?". Neither
`blastn` nor `minimap2` was available on the machine used, and no C compiler
was available either (`pip install mappy|edlib|parasail` all failed to build).
`mapper.py` therefore implements a minimal nucleotide aligner in pure numpy:

1. a k = 15 k-mer index of the genome (both strands are searched);
2. seed hits are clustered into candidate start positions;
3. each candidate is verified by gapless extension with a ±8 bp shift search;
4. the best hit is reported with its identity, coverage, contig and position.

**Accuracy self-check** (`06_final_checks.py`, section A): the 7 markers whose
own source genome is one of the three analysed genomes are all recovered from
that same genome with an identity of exactly **1.0000**, i.e. the aligner
recovers a sequence from the very genome it came from without error.

## 6. Key results

| Quantity | Strain 1 | Strain 2 | *D. hominis* |
|---|---|---|---|
| Markers assigned to a genome position (≥90% id, ≥90% cov) | 179 / 200 | 178 / 200 | 195 / 200 |
| Independent loci (gap > 5 kb) | 130 | 130 | 132 |
| Median spacing between adjacent markers | 9.95 kb | 8.36 kb | 8.15 kb |
| Maximum markers in any 10-kb window | 4 (2.2%) | 4 (2.2%) | 4 (2.1%) |
| Uniformity within contigs (KS *P*) | 0.135 | 0.513 | 0.280 |
| Uniformity along the oriC→ter axis (KS *P*) | 0.003 | 0.058 | 0.661 |
| Estimated bias, oriC:ter = 2 | +4.31% | −2.24% | −1.06% |
| Estimated bias, oriC:ter = 4 | +8.03% | −4.84% | −2.19% |
| Estimated bias, oriC:ter = 8 | +11.04% | −7.61% | −3.32% |

Control: if the same 200 markers were confined to the 5% of the genome closest
to the origin, the same model would give **+36.3% / +78.6% / +125.7%** at
oriC:ter = 2 / 4 / 8.

Species vs strain: 176 / 200 markers (88.0%) were recovered in **both**
*W. acetigignens* strains, 5 (2.5%) were present in one strain only, and 19
(9.5%) were absent from both (mostly markers contributed by genomes at the edge
of the SGB). Median identity between the two strains was 0.988.

Cross-species control: **0 / 200** markers of either SGB matched the other
species' genome at ≥90% identity (median identity 0.272–0.274, i.e. the level
expected by chance between unrelated sequences).

## 7. Caveats

- All three genomes are **drafts** (36–125 contigs, N50 70–300 kb).
- The oriC→ter axis is **reconstructed** from cumulative GC skew and covers
  89–97% of each genome; it is an approximation, not a finished circular
  chromosome. The dispersion results (loci count, spacing, window density) do
  not depend on this reconstruction.
- The copy-number profile `C(x) = R^(1−x)` is a **model with an assumed
  parameter** R; the resulting percentages indicate the *order of magnitude* of
  the bias, not a measured value. Strain 1 deviates from uniformity along the
  reconstructed axis (KS *P* = 0.003), so a small offset may exist; this is why
  the argument rests on the marker set being dispersed and, above all, on the
  fact that the **same marker set is used for every sample**, so any residual
  offset is a species-specific constant that cancels in case–control contrasts.
- Only 2 SGBs / 3 genomes were examined; the results should not be extrapolated
  to the whole database without further sampling.
