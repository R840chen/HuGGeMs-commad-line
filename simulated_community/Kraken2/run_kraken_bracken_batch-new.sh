#!/bin/bash
set -euo pipefail

########################################
# 输入与输出路径
########################################

# 输入目录
input_dir="/ME5012-Vol01/chenc/CAMII-simulated-reads/all-simulated-reads"

# 统一输出目录
output_dir="/ME5012-Vol01/chenc/MetaphIan/paper-data/2.CAMI2-profiling/3.Kraken2-result/confidence-0"

# Kraken2 数据库
kraken_db="/ME4012-Vol01/database/kraken2"

########################################
# 参数设置
########################################

threads=20
read_length=150
confidence=0

########################################
# 准备环境
########################################

mkdir -p "$output_dir"

cd "$input_dir" || {
    echo "? 无法进入输入目录: $input_dir"
    exit 1
}

echo "?? 输入目录: $input_dir"
echo "?? 输出目录: $output_dir"
echo "?? 线程数: $threads"
echo "?? Kraken2 confidence: $confidence"
echo "----------------------------------"

########################################
# 批量运行 Kraken2 + Bracken
########################################

shopt -s nullglob
for fq_file in *.fq *.fq.gz; do

    # 去掉 .fq 或 .fq.gz 后缀
    base_name=$(basename "$fq_file")
    base_name=${base_name%.fq}
    base_name=${base_name%.fq.gz}

    echo "?? 处理文件: $fq_file"

    # Kraken2 输出
    kraken_output="${output_dir}/kraken2-${base_name}.out"
    kraken_report="${output_dir}/kraken2-${base_name}.report"

    echo "  ? Kraken2 运行中..."
    kraken2 classify "$fq_file" \
        --db "$kraken_db" \
        --threads "$threads" \
        --confidence "$confidence" \
        --output "$kraken_output" \
        --report "$kraken_report"

    # Bracken 输出
    bracken_output="${output_dir}/bracken-${base_name}.out"
    bracken_report="${output_dir}/bracken-${base_name}.report"

    echo "  ? Bracken 运行中..."
    bracken \
        -d "$kraken_db" \
        -i "$kraken_report" \
        -o "$bracken_output" \
        -w "$bracken_report" \
        -r "$read_length" \
        -t "$threads"

    echo "  ? $fq_file 分析完成"
    echo "----------------------------------"
done

echo "?? 所有 Kraken2 + Bracken 批量分析完成！"
