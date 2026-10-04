# Input data

Put the five input files here, or set the environment variable
`MARKER_DATA_DIR` to the directory that contains them.

| File name | Where it comes from |
|---|---|
| `Brotolimicola_acetigignens_species_2981769.ffn` | marker set (200 sequences) of SGB_2981769, extracted from the HuGGeMs database; also provided as Supplementary data |
| `Dysosmobacter_hominis_species_2763041.ffn` | marker set (200 sequences) of SGB_2763041, extracted from the HuGGeMs database; also provided as Supplementary data |
| `GCA_025567015.1_ASM2556701v1_genomic.fna` | NCBI GenBank `GCA_025567015.1` (*Waltera acetigignens*, strain 1) |
| `GCA_020687025.1_ASM2068702v1_genomic.fna` | NCBI GenBank `GCA_020687025.1` (*Waltera acetigignens*, strain 2) |
| `GCA_014297375.1_ASM1429737v1_genomic.fna` | NCBI GenBank `GCA_014297375.1` (*Dysosmobacter hominis*) |

The three genomes are downloaded automatically by

```bash
bash ../fetch_genomes.sh
```

The file name of the first marker set still uses the pre-2024 species name
*Brotolimicola acetigignens*; the species was reclassified as
*Waltera acetigignens* (PMID 38722771). Scripts refer to it by file name, so
the name is kept as is.
