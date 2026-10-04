# -*- coding: utf-8 -*-
import os

# 原始固定配置文件内容（可存放在脚本中，也可以从文件读取）
base_config = """[Main]
# maximum number of processes
max_processors=15

# 0: community design + read simulator,
# 1: read simulator only
phase=0

# ouput directory, where the output will be stored (will be overwritten if set in from_profile)
output_directory=/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/result-simulate-reads

# temporary directory
temp_directory=/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/tmp

# gold standard assembly
gsa=True

# gold standard for all samples combined
pooled_gsa=True

# anonymize sequences?
anonymous=True

# compress data (levels 0-9, recommended is 1 the gain of higher levels is not too high)
compress=1

# id of dataset, used in foldernames and is prefix in anonymous sequences
dataset_id=RL

# Read Simulation settings, relevant also for from_profile
[ReadSimulator]
# which readsimulator to use:
#           Choice of 'art', 'wgsim', 'nanosim', 'pbsim'
type=art

# Samtools (http://www.htslib.org/) takes care of sam/bam files. Version 1.0 or higher required!
# file path to executable
samtools=tools/samtools-1.3/samtools

# file path to read simulation executable
readsim=tools/art_illumina-2.3.6/art_illumina

#error profiles:
#for ART:
#HiSeq 150bp: hi150
#MBARC-26 150bp: mbarc
#custom profile (see below): own
#for wgsim:
#error rate as <float> (e.g. 0.05 for 5% error rate)
#blank for nanosim and wgsim
profile=mbarc

# Directory containing error profiles (can be blank for wgsim)
error_profiles=tools/art_illumina-2.3.6/profiles/

# For supplying custom error profiles with "own" option:
# path to error profile files (without "[1/2].txt")
base_profile_name=
# read length for custom error profile
profile_read_length=

#paired end read, insert size (not applicable for nanosim)
fragments_size_mean=270
fragment_size_standard_deviation=27

# Only relevant if not from_profile is run:
[CommunityDesign]
# specify the samples size in Giga base pairs
size=5

distribution_file_paths=/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/abundance0.tsv

# how many different samples?
number_of_samples=1

# how many communities
num_communities=1

# directory containing the taxdump of ncbi, version from 22.02.2017 is shipped
# "nodes.dmp"
# "merged.dmp"
# "names.dmp"
ncbi_taxdump=tools/ncbi-taxonomy_20170222.tar.gz

# the strain simulator for de novo strain creation
strain_simulation_template=scripts/StrainSimulationWrapper/sgEvolver/simulation_dir/

# define communities: [community<integer>]
[community0]
# information about all included genomes:
# can be used for multiple samples
metadata=/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/metadata.tsv
id_to_genome_file=/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/genome_to_id.tsv

# how many genomes do you want to sample over all?
genomes_total=46
num_real_genomes=

# how many genomes per species taxon
#   (species taxon will be replaced by OTU-cluster later on)
max_strains_per_otu=1
ratio=1

# which kind of different samples do you need?
#   replicates / timeseries_lognormal / timeseries_normal / differential
mode=

# Part: community design
# Set parameters of log-normal and normal distribution, number of samples
# sigma > 0; influences shape (higher sigma -> smaller peak and longer tail),
log_sigma=2

# mu (real number) is a parameter for the log-scale
log_mu=1

# do you want to see a distribution before you decide to use it? yes/no
view=no
"""

# 创建文件
for i in range(1, 26):
    output_dir = f"/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/simulated-{i}-result-reads"
    abundance_file = f"/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/simulated-abundance{i}.tsv"
    metadata_file = f"/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/simulated-{i}-metadata.tsv"
    genome_id_file = f"/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/simulated-{i}-genome_to_id.tsv"

    # 替换内容
    config = base_config.replace(
        "/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/result-simulate-reads", output_dir
    ).replace(
        "/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/abundance0.tsv", abundance_file
    ).replace(
        "/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/metadata.tsv", metadata_file
    ).replace(
        "/ME4012-Vol01/chenc/software/CAMISIM/simulated-reads/20250726-ssimulated-16-community/genome_to_id.tsv", genome_id_file
    )

    # 写入文件
    filename = f"default_config-simulated-{i}.ini"
    with open(filename, 'w') as f:
        f.write(config)

print("All config files generated successfully.")
