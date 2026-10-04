#!/bin/bash

# 设置前缀时间戳
prefix="standard_profile"

# 循环处理 sample_0 到 sample_12
for i in {0..4}
do
    echo "Processing sample_$i..."

    # 设置输入文件
    gold_standard="/ME5012-Vol01/chenc/real-metagenomes-validate/3-synthetic-mock-community/community-species/opal-file/community-${i}.cami2.tsv"
    prediction="cami2_gold_${prefix}-CC${i}_1.txt"
    output_dir="3-synthetic_sample_CC${i}"

    # 检查文件是否存在
    if [[ -f "$gold_standard" && -f "$prediction" ]]; then
        opal.py "$prediction" -g "$gold_standard" -o "$output_dir"
        echo "Finished sample_$i."
    else
        echo "Skipping sample_$i: Missing input file(s)."
    fi
done

echo "All processing done."
