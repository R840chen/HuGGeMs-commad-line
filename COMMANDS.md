# Complete command lines and parameters

Extracted verbatim from the scripts in this repository. Placeholders written
as `<...>` correspond to the script variables; the values actually used are
listed under each command.

---

## A. Synthetic / mock community benchmarks

### A1. Read simulation — CAMISIM

Scripts: `simulated_community/CAMISIM/1-generate-25-ini-file.py`,
`2-modified-ini-file.py`, `3-run-ini.sh`

Execution:

```bash
parallel -j 10 "python <CAMISIM>/metagenomesimulation.py default_config-simulated-{}.ini" ::: {1..16}
```

Key `.ini` parameters (identical across the 25 generated config files):

| Section | Parameter | Value |
|---|---|---|
| Main | `max_processors` | 15 |
| Main | `phase` | 0 (community design + read simulator) |
| Main | `gsa` | True (gold-standard assembly) |
| Main | `pooled_gsa` | True |
| Main | `anonymous` | True |
| Main | `compress` | 1 |
| Main | `dataset_id` | RL |
| ReadSimulator | `type` | art |
| ReadSimulator | `profile` | mbarc (MBARC-26, 150 bp) |
| ReadSimulator | `fragments_size_mean` | 270 |
| ReadSimulator | `fragment_size_standard_deviation` | 27 |
| ReadSimulator | `samtools` | tools/samtools-1.3/samtools |
| ReadSimulator | `readsim` | tools/art_illumina-2.3.6/art_illumina |
| CommunityDesign | `size` | 5 (Gbp) |
| CommunityDesign | `number_of_samples` | 1 |
| CommunityDesign | `num_communities` | 1 |
| community0 | `genomes_total` | 46 |
| community0 | `max_strains_per_otu` | 1 |
| community0 | `ratio` | 1 |
| community0 | `log_sigma` | 2 |
| community0 | `log_mu` | 1 |

### A2. Profiling of the simulated communities

**HuGGeMs** — `simulated_community/HuGGeMs/run-contrast-metagenomes.sh`

```bash
metaphlan "<reads.fastq.gz>" --input_type fastq \
  -o "<out>/mpa4-3-synthetic-sample-result-<sample>.txt" \
  --bowtie2out "<out>/bowtie2out/<sample>.bowtie2out.txt" \
  --nproc 90 \
  --bowtie2db "<HuGGeMs_INDEX_DIR>" \
  --force --index GMB-longer-first-20260123_NEW
```
Jobs dispatched with `parallel -j 3`.

**MetaPhlAn4 (default database)** — `simulated_community/MetaPhlAn4/run-contrast-metagenomes.sh`

```bash
metaphlan "<reads.fastq.gz>" --input_type fastq \
  -o "<out>/mpa4-3-synthetic-sample-result-<sample>.txt" \
  --bowtie2out "<out>/bowtie2out/<sample>.bowtie2out.txt" \
  --nproc 90 \
  --bowtie2db "/mnt/ME4012-Vol01/database/mpa4.2-202503" \
  --force --index mpa_vJan25_CHOCOPhlAnSGB_202503
```
Jobs dispatched with `parallel -j 3`.

**Profiling output → CAMI2 format** (for OPAL):
`simulated_community/{HuGGeMs,MetaPhlAn4}/1-metaphlan_to_cami2-test-parallel.py`

**Kraken2 + Bracken** — `simulated_community/Kraken2/`

Single-end (`run_kraken_bracken_batch-new.sh`):

```bash
kraken2 classify "<reads.fq>" \
  --db /ME4012-Vol01/database/kraken2 \
  --threads 20 --confidence 0 \
  --output "<out>/kraken2-<sample>.out" \
  --report "<out>/kraken2-<sample>.report"

bracken -d /ME4012-Vol01/database/kraken2 \
  -i "<out>/kraken2-<sample>.report" \
  -o "<out>/bracken-<sample>.out" \
  -w "<out>/bracken-<sample>.report" \
  -r 150 -t 20
```

Paired-end (`run_kraken_bracken_batch-new-paired.sh`):

```bash
kraken2 --db /ME4012-Vol01/database/kraken2 \
  --threads 20 --confidence 0 --paired \
  --output "<out>/kraken2-<sample>.out" \
  --report "<out>/kraken2-<sample>.report" \
  "<R1>" "<R2>"

bracken -d /ME4012-Vol01/database/kraken2 \
  -i "<out>/kraken2-<sample>.report" \
  -o "<out>/bracken-<sample>.out" \
  -w "<out>/bracken-<sample>.report" \
  -r 150 -t 20
```

Post-processing: `1-process_bracken_results.py` (keeps
`fraction_total_reads > 0.0001`, reformats to CAMI2 columns).

**mOTUs v3** — `simulated_community/mOTUs3/`

| Script | Command | `-g` | output dir |
|---|---|---|---|
| `run_motus_batch-precision-test.sh` | `motus profile -g 6 -t 50 -s <reads.fq> -o <out>` | 6 | `4.mOTUs-result/precision` |
| `run_motus_batch-recall-test.sh` | `motus profile -g 1 -t 50 -s <reads.fq> -o <out>` | 1 | `4.mOTUs-result/recall` |

**mOTUs v4 (4.1.0)** — `simulated_community/mOTUs4/`

| Script | Command | `-g` | output dir |
|---|---|---|---|
| `run_motus_batch-precision-test.sh` | `motus profile -s <reads.fq.gz> -g 6 -t 50 -o <out>` | 6 | `7.motus4-result` |
| `run_motus_batch-recall-test.sh` | `motus profile -s <reads.fq.gz> -g 1 -t 50 -o <out>` | 1 | `7.motus4-result/recall` |

Post-processing to CAMI format: `motu2cami.py`, `motu_batch2cami.py`.

**sylph + sylph-tax** — `simulated_community/sylph/`

```bash
sylph profile /mnt/ME4012-Vol01/database/sylph/gtdb-r220-c200-dbv1.syldb \
  -r "<reads.fq>" -t 150 -m 95 > profiling.tsv

sylph-tax taxprof -o 3-synthetic-community profiling.tsv \
  -t GTDB_r226 IMGVR_4.1 GTDB_r220
```

Conversion to OPAL profile (species only):
`run_sylph_to_opal_batch.sh` → `1-sylph_add_ncbi_taxid.py` →
`2-sylph_ncbi_to_opal_profile.py`

```bash
python 1-sylph_add_ncbi_taxid.py <in.sylphmpa> <assembly_summary_genbank.txt> \
       <assembly_summary_refseq.txt> <out.ncbi.tsv>
python 2-sylph_ncbi_to_opal_profile.py <in.ncbi.tsv> <sample_id> <out.opal.tsv>
```

### A3. Benchmark evaluation — OPAL

Scripts: `simulated_community/*/3-run_CAMI-opal*.sh`

```bash
opal.py "<prediction_profile>" -g "<gold_standard>" -o "<output_dir>"
```

- gold standard: `community-<i>.cami2.tsv` (synthetic) / `taxonomic_profile_<i>.txt` (CAMI2)
- prediction: the tool's CAMI-format profile
- samples iterated `i = 0..4` (synthetic) or `0..13` / `0..16` (CAMI2), depending on the tool

> Note: in `mOTUs3/3-run_CAMI-opal.sh` and `mOTUs4/3-run_CAMI-opal.sh` the two
> positional arguments are supplied in the reverse order compared with the other
> scripts. Please verify before reuse.

---

## B. Real metagenome profiling and association analysis

### B1. Profiling command (real cohorts)

Same MetaPhlAn invocation as in A2 (single-end), with the HuGGeMs index:

```bash
metaphlan "<reads.fastq.gz>" --input_type fastq -o "<out>.txt" \
  --bowtie2out "<out>.bowtie2out.txt" --nproc 90 \
  --bowtie2db "<HuGGeMs_INDEX_DIR>" --force --index GMB-longer-first-20260123_NEW
```

### B2. Species matrix construction

`real-metagenome-profiling/HuGGeMs/生成matrices文件/`

| Step | Script | Parameter |
|---|---|---|
| 1 | `1-metaphlan_txt_to_species_matrix.py` | MetaPhlAn txt → species × sample matrix |
| 1' | `1-1-batch_process_matrices.py` | batch wrapper |
| 2 | `2-filter_by_prevalence.py` | `prevalence_threshold = 0.05` (≥5% of samples) |
| 2' | `2-2-filter_by_prevalence-0.1.py` | `PREVALENCE_THRESHOLD = 0.1` (sensitivity check) |
| 3 | `3-clr_transform.py` | `pseudocount = 1e-6`, CLR = log(x+pc) − row/column geometric mean |
| 3' | `3-1-modified-matrices_file-name.py` | rename to `*_species_abundance_matrix.prev5pct.maaslin_ready.tsv` |
| 4 | `4-filter_and_rename_metaphlan_matrix.py` | final filtering / renaming |

### B3. MaAsLin3

**Input preparation** — `real-metagenome-profiling/HuGGeMs/run_maaslin3/`

```bash
python 00_unify_metadata.py       # unify metadata format across cohorts
python 01_make_maaslin_inputs.py  # explicit complete-case filtering per level
```

`01_make_maaslin_inputs.py` parameters:

| Parameter | Value | Meaning |
|---|---|---|
| `MIN_GROUP_N` | 10 | minimum samples per disease/control group |
| `MIN_COV_COVERAGE` | 0.70 | minimum non-missing fraction per covariate |
| `MIN_LEVEL_N` | 5 | minimum samples per categorical level |
| `MIN_NUM_UNIQUE` | 3 | minimum distinct values for a numeric covariate |
| `REFERENCE_GROUP` | Health | reference level |
| `GEO_MIN_N` | 10 | minimum samples per country (for the LG level) |

Model levels (formula):

| Level | Formula | Covariates added |
|---|---|---|
| L1 | `~ group` | none (baseline) |
| L2 | `~ group + age` | age |
| L3 | `~ group + sex` | sex |
| L4 | `~ group + BMI` | BMI |
| L5 | `~ group + age + sex + BMI` | all three |
| LG | `~ group + country` | country (multi-country cohorts only) |

**Model fitting**

- `run_maaslin3.R` — GMB profiles, `formula = "~ group"`
- `metaphlan4/RUN.R` — MetaPhlAn4 profiles, `formula = "~ group"`
- `run_maaslin3_by_level.R` — all levels L1–L5 + LG

`maaslin3()` call and parameters:

```r
maaslin3(
  input_data       = "<PID>__<LEVEL>__data.tsv",
  input_metadata   = meta_df,
  output           = "<out>/<PID>/<LEVEL>/",
  formula          = "~ group + age + sex + BMI",   # level-dependent
  normalization    = "TSS",
  transform        = "LOG",
  standardize      = TRUE,
  augment          = TRUE,
  min_prevalence   = 0.05,
  max_significance = 0.1,
  cores            = 4,
  save_models      = TRUE,
  reference        = "group,Health"
)
```

Additional safeguards in `run_maaslin3_by_level.R`: numeric coercion of
`age`/`BMI`, factorisation of `group`/`sex`/`country`, `relevel(group, ref =
"Health")`, and a rank check on the model matrix (`qr(mm)$rank < ncol(mm)` →
skip) to avoid collinear covariates.

**Result collection** — `02_collect_results.py`

| Parameter | Value |
|---|---|
| `Q_THRESHOLD` | 0.05 (on `qval_joint`) |
| `MODELS` | `abundance`, `prevalence` |
| `REF_GROUP` | Health |

Model filter: `metadata == "group"` and `value != "Health"`; direction from the
sign of `coef` (>0 enriched, <0 depleted).

Outputs: `_all_significant.tsv`, `_cross_level_summary.tsv`,
`_robust_species.tsv`, `_count_by_cohort.tsv`.

**Cross-cohort / covariate robustness** — `03_analyze_robustness.py`

| Parameter | Value |
|---|---|
| `QLVL` | 0.05 |
| `MAIN_IBD` | PRJEB1220, PRJNA389280, PRJNA398089, SRP057027 |
| `COV_COHORT` | PRJEB1220 |
| `COV_LEVEL` | L5 |
| `LEVEL_ORDER` | L1, L2, L3, L4, L5, LG |

Outputs: `_core_candidate_species.tsv`, `_cross_cohort_consistency.tsv`.

---

## C. Strain-level test (StrainPhlAn)

`test_profiling_strainphlan/`

**1. MetaPhlAn with SAM output** — `run_metaphlan4-file.sh`
(generates `run_PRJDB4176_metaphlan.txt`)

```bash
metaphlan "<R1.fastq.gz>,<R2.fastq.gz>" --input_type fastq \
  -o "<out>/mpa4-<sample>.txt" \
  --bowtie2out "<out>/<sample>.bowtie2out.txt" \
  -s "<out>/<sample>.sam.bz2" \
  --nproc 60 \
  --bowtie2db "<HuGGeMs_INDEX_DIR>" \
  --force --index GMB-longer-first-20260123_NEW
```
Jobs dispatched with `parallel -j 3`.

**2. Marker extraction** — `run_sample2markers.py.sh`

```bash
sample2markers.py -i "<out>/<sample>.sam.bz2" -o "<out>/markers" \
  -d <HuGGeMs_INDEX_DIR>/GMB-longer-first-20260123_NEW.pkl -n 100
```

**3. StrainPhlAn** — example output in `test-output/`
(`Phocaeicola dorei`, RAxML tree). The final `strainphlan` call was run
interactively and is not scripted in this folder.

Included reference/marker files:
- `reference/GCA_013009555.1_ASM1300955v1_genomic.fna`
- `clade_markers/Phocaeicola_dorei_species_357276.ffn`
