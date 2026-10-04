#!/bin/bash


# 循环处理 sample_0 到 sample_12
for i in {0..13}
do
    echo "Processing sample_$i..."

    # 设置输入文件
    gold_standard="2017.12.04_18.45.54_sample_${i}.cami"
    prediction="taxonomic_profile_${i}.txt"
    output_dir="cami-motus-2017.12.04_18.45.54_sample_${i}"

    # 检查文件是否存在
    if [[ -f "$gold_standard" && -f "$prediction" ]]; then
        opal.py "$gold_standard" -g "$prediction" -o "$output_dir"
        echo "Finished sample_$i."
    else
        echo "Skipping sample_$i: Missing input file(s)."
    fi
done

echo "All processing done."
