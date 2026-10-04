# HuGGeMs — command lines and parameters

This repository contains the scripts and the exact command lines used for
(i) the CAMI2 and synthetic-community benchmarks, (ii) the profiling of real
metagenomes, (iii) the strain-level test, and (iv) the MaAsLin3 differential
abundance / prevalence analysis reported in the HuGGeMs manuscript.

It is provided to satisfy the request for a *complete command-line parameter
list* and to allow full reproduction of the analyses.

A condensed overview of every tool and its parameters is given below; the
per-step detail is in [`COMMANDS.md`](COMMANDS.md).

---

## 1. Directory map

```
simulated_community/                 # CAMI2 + synthetic (mock) community benchmarks
├── CAMISIM/                         # simulation of 16 synthetic communities (ART reads)
├── HuGGeMs/                         # profiling with the HuGGeMs database
├── MetaPhlAn4/                      # profiling with the default MetaPhlAn4 database
├── Kraken2/                         # Kraken2 + Bracken
├── mOTUs3/                          # mOTUs v3
├── mOTUs4/                          # mOTUs v4 (4.1.0)
└── sylph/                           # sylph + sylph-tax

real-metagenome-profiling/           # real metagenomes
├── HuGGeMs/
│   ├── 生成matrices文件/            # MetaPhlAn txt -> species matrix -> prevalence filter -> CLR
│   └── run_maaslin3/                # MaAsLin3 input preparation, model fitting, result collection
└── metaphlan4/
    └── RUN.R                        # MaAsLin3 batch run on MetaPhlAn4 profiles

test_profiling_strainphlan/          # strain-level test (PRJDB4176, colon cancer)
```

## 2. Software used

| Software | Version | Role |
|---|---|---|
| CAMISIM (`metagenomesimulation.py`) | — | simulation of synthetic communities |
| ART (`art_illumina`) | 2.3.6 | Illumina read simulation |
| samtools | 1.3 | SAM/BAM handling during simulation |
| MetaPhlAn | 4.x | species-level profiling (HuGGeMs and default database) |
| Bowtie2 | bundled with MetaPhlAn | read alignment against the marker index |
| Kraken2 | — | k-mer based classification |
| Bracken | — | abundance re-estimation from Kraken2 reports |
| mOTUs | v3 and v4.1.0 | marker-gene based profiling |
| sylph / sylph-tax | — | k-mer containment based profiling |
| OPAL (`opal.py`) | — | benchmark evaluation against CAMI gold standards |
| StrainPhlAn (`sample2markers.py`) | 4.x | strain-level marker extraction |
| MaAsLin3 (R package) | 1.2.0 | differential abundance / prevalence models |

> Versions marked “—” were not recorded in the scripts; please see
> [§5 Notes](#5-notes-on-reproducibility) before citing them.

### Databases / indices

| Resource | Path used |
|---|---|
| HuGGeMs Bowtie2 index | `.../7-seventh-all-maker-gene-db/20251219-filtered-gene`, index name `GMB-longer-first-20260123_NEW` |
| Default MetaPhlAn4 database | `mpa4.2-202503`, index name `mpa_vJan25_CHOCOPhlAnSGB_202503` |
| Kraken2 database | `/ME4012-Vol01/database/kraken2` |
| sylph database | `gtdb-r220-c200-dbv1.syldb` (taxonomy: GTDB_r226, IMGVR_4.1, GTDB_r220) |

## 3. Quick reference — the commands that matter

**MetaPhlAn (single-end)**
```bash
metaphlan <reads.fq> --input_type fastq -o <out.txt> \
  --bowtie2out <out.bowtie2out.txt> --nproc 90 \
  --bowtie2db <INDEX_DIR> --force --index <INDEX_NAME>
```

**MetaPhlAn (paired-end, with SAM output for StrainPhlAn)**
```bash
metaphlan "<R1.fastq.gz>,<R2.fastq.gz>" --input_type fastq -o <out.txt> \
  --bowtie2out <out.bowtie2out.txt> -s <out.sam.bz2> --nproc 60 \
  --bowtie2db <INDEX_DIR> --force --index <INDEX_NAME>
```

**Kraken2 + Bracken**
```bash
kraken2 --db <KRAKEN_DB> --threads 20 --confidence 0 --paired \
  --output <out> --report <report> <R1> <R2>
bracken -d <KRAKEN_DB> -i <report> -o <out> -w <report.w> -r 150 -t 20
```

**mOTUs**
```bash
motus profile -s <reads.fq> -g <1|6> -t 50 -o <out.txt>
```

**sylph + sylph-tax**
```bash
sylph profile <db.syldb> -r <reads.fq> -t 150 -m 95 > profiling.tsv
sylph-tax taxprof -o <out_prefix> profiling.tsv -t GTDB_r226 IMGVR_4.1 GTDB_r220
```

**OPAL**
```bash
opal.py <prediction_profile> -g <gold_standard> -o <output_dir>
```

**StrainPhlAn marker extraction**
```bash
sample2markers.py -i <sample.sam.bz2> -o <marker_dir> -d <INDEX>.pkl -n 100
```

**MaAsLin3 (R)**
```r
maaslin3(input_data = <species_matrix.tsv>,
         input_metadata = <meta_df>,
         output = <out_dir>,
         formula = "~ group",
         normalization = "TSS", transform = "LOG",
         standardize = TRUE, augment = TRUE,
         min_prevalence = 0.05, max_significance = 0.1,
         cores = 4, save_models = TRUE,
         reference = "group,Health")
```

## 4. How to run

Each directory is self-contained and the scripts are numbered in execution
order. Absolute paths inside the scripts (e.g. `/ME5012-Vol01/...`,
`D:/Desktop/...`) must be replaced with paths valid on your own system.

Typical CAMI2 / synthetic-community workflow:

```bash
# 1) simulate communities
cd simulated_community/CAMISIM
python 1-generate-25-ini-file.py
python 2-modified-ini-file.py
bash 3-run-ini.sh

# 2) profile the simulated reads, one directory per tool, e.g.
cd ../HuGGeMs
bash run-contrast-metagenomes.sh

# 3) evaluate against the CAMI gold standard
bash 3-run_CAMI-opal.sh
```

Typical real-metagenome workflow:

```bash
cd real-metagenome-profiling/HuGGeMs/生成matrices文件
python 1-metaphlan_txt_to_species_matrix.py
python 2-filter_by_prevalence.py          # prevalence >= 5%
python 3-clr_transform.py                 # pseudocount 1e-6

cd ../run_maaslin3
python 00_unify_metadata.py
python 01_make_maaslin_inputs.py          # complete-case filtering, L1-L5 + LG
Rscript run_maaslin3_by_level.R
python 02_collect_results.py              # qval_joint < 0.05
python 03_analyze_robustness.py
```

## 5. Notes on reproducibility

- **Data files are not included.** Large intermediate files
  (`*.pkl` StrainPhlAn markers, `*.fna` reference genomes) are excluded by
  `.gitignore`; they are either publicly available (genomes) or reproducible
  from the commands here. The three test genomes are available from NCBI under
  the accessions listed in the manuscript.
- **Hard-coded absolute paths.** All scripts were written for a specific
  server layout and contain absolute paths. They are reproduced here verbatim
  so that the exact parameters are unambiguous; adjust the paths to run them.
- **Versions to be confirmed.** A few versions were not recorded in the
  scripts (CAMISIM, Kraken2, Bracken, sylph, OPAL, MetaPhlAn, StrainPhlAn).
  See the table in §2 and fill in from `--version` if required.
- **OPAL argument order.** In `simulated_community/mOTUs3/3-run_CAMI-opal.sh`
  and `simulated_community/mOTUs4/3-run_CAMI-opal.sh` the two arguments are
  passed in the reverse order relative to the other tools' scripts; please
  double-check before reuse.
- **StrainPhlAn step.** This directory contains the MetaPhlAn and
  `sample2markers.py` steps plus an example output; the final `strainphlan`
  call itself was run interactively and is not scripted here.
