#!/bin/bash
set -euo pipefail

########################################
# 输入与输出路径
########################################

# 双端 reads 输入目录
input_dir="/ME5012-Vol01/chenc/real-metagenomes-validate/3-synthetic-mock-community/CP2025120400151/H101SC25122830/RSMD00104/X101SC25122830-Z01/X101SC25122830-Z01-J001/00.CleanData/all-fq"

# 统一输出目录
output_dir="/ME5012-Vol01/chenc/MetaphIan/paper-data/4.all-3-synthetic-community/3.Kraken2-result/confidence-0"

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
# 批量运行 Kraken2 + Bracken（双端）
########################################

shopt -s nullglob
for r1 in *_1.clean.fq *_1.clean.fq.gz; do

    # 样本名（去掉 _1.clean.fq(.gz)）
    sample=$(basename "$r1")
    sample=${sample%_1.clean.fq}
    sample=${sample%_1.clean.fq.gz}

    # 对应的 R2
    r2="${sample}_2.clean.fq"
    [ -f "$r2" ] || r2="${r2}.gz"

    # 检查 R2 是否存在
    if [[ ! -f "$r2" ]]; then
        echo "??  缺失 R2，跳过样本: $sample"
        continue
    fi

    echo "?? 处理样本: $sample"
    echo "   R1: $r1"
    echo "   R2: $r2"

    # Kraken2 输出
    kraken_output="${output_dir}/kraken2-${sample}.out"
    kraken_report="${output_dir}/kraken2-${sample}.report"

    echo "  ? Kraken2 运行中 (paired-end)..."
    kraken2 \
        --db "$kraken_db" \
        --threads "$threads" \
        --confidence "$confidence" \
        --paired \
        --output "$kraken_output" \
        --report "$kraken_report" \
        "$r1" "$r2"

    # Bracken 输出（不区分单/双端）
    bracken_output="${output_dir}/bracken-${sample}.out"
    bracken_report="${output_dir}/bracken-${sample}.report"

    echo "  ? Bracken 运行中..."
    bracken \
        -d "$kraken_db" \
        -i "$kraken_report" \
        -o "$bracken_output" \
        -w "$bracken_report" \
        -r "$read_length" \
        -t "$threads"

    echo "? 样本完成: $sample"
    echo "----------------------------------"
done

echo "?? 所有 Kraken2 + Bracken 双端分析完成！"
