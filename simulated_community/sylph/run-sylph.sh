sylph profile /mnt/ME4012-Vol01/database/sylph/gtdb-r220-c200-dbv1.syldb -r /mnt/ME5012-Vol01/chenc/CAMII-simulated-reads/all-simulated-reads/2017.12.04_18.45.54_sample_*-anonymous_reads.fq -t 150 -m 95 > profiling.tsv

sylph-tax taxprof -o 3-synthetic-community profiling.tsv -t GTDB_r226 IMGVR_4.1 GTDB_r220