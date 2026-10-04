#!/bin/bash

# 设置前缀时间戳
prefix="3-synthetic-community2017.12.04_18.45.54_sample"

# 循环处理 sample_0 到 sample_12
for i in {0..16}
do
    echo "Processing sample_$i..."

    # 设置输入文件
    gold_standard="taxonomic_profile_${i}.txt"
    prediction="${prefix}_${i}-anonymous_reads.fq.species.opal.tsv"
    output_dir="simulated_sample_CC${i}"

    # 检查文件是否存在
    if [[ -f "$gold_standard" && -f "$prediction" ]]; then
        opal.py "$prediction" -g "$gold_standard" -o "$output_dir"
        echo "Finished sample_$i."
    else
        echo "Skipping sample_$i: Missing input file(s)."
    fi
done

echo "All processing done."
