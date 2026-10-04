#!/bin/bash
set -euo pipefail

# 输入根目录（递归搜索该目录及其所有子目录）
input_root="/mnt/ME5012-Vol01/chenc/real-metagenomes-validate/3-synthetic-mock-community/CP2025120400151/H101SC25122830/RSMD00104/X101SC25122830-Z01/X101SC25122830-Z01-J001/00.CleanData/all-fq/1"

# 输出基础目录（每个输入文件的子目录会在此下创建对应子目录）
output_base="./"
mkdir -p "$output_base"

# 输出文件前缀
output_prefix="3-synthetic-sample-result"

# 线程数（每个任务内部的线程数）
threads=90

# 每次同时运行的任务数（parallel -j）
parallel_jobs=3

# 生成命令列表文件（放在指定位置）
cmd_file="./run_metaphlan_commands.txt"
# 清空或创建文件
: > "$cmd_file"

# Bowtie2 DB（保持原样，如需更改请修改此变量）
bowtie2db="/mnt/ME4012-Vol01/database/mpa4.2-202503"
index_name="mpa_vJan25_CHOCOPhlAnSGB_202503"

# 找到所有 fastq.gz / fq.gz 文件并为每个生成一条 metaphlan 命令
# 使用 find 递归查找（不区分大小写）
while IFS= read -r -d '' fq_file; do
    # 计算输入文件所在的目录（绝对路径）
    fq_dir="$(dirname "$fq_file")"

    # 计算相对于 input_root 的子目录路径（去掉前缀 input_root/）
    rel_subdir="${fq_dir#$input_root/}"

    # 如果文件直接在 input_root 下，rel_subdir 会变为空字符串 -> 直接使用 output_base
    if [ -z "$rel_subdir" ] || [ "$rel_subdir" = "$fq_dir" ]; then
        out_subdir="$output_base"
    else
        out_subdir="$output_base/$rel_subdir"
    fi

    # 确保输出目录存在（包括 bowtie2out 子目录）
    mkdir -p "$out_subdir"
    mkdir -p "$out_subdir/bowtie2out"

    # 计算基础文件名（去掉 .fastq.gz 或 .fq.gz 或 .fastq 或 .fq 后缀）
    base_name="$(basename "$fq_file")"
    # 依次去掉可能的后缀（最长匹配先移除 .fastq.gz/.fq.gz）
    base_name="${base_name%%.fastq.gz}"
    base_name="${base_name%%.fq.gz}"
    base_name="${base_name%%.fastq}"
    base_name="${base_name%%.fq}"

    output_file="${out_subdir}/mpa4-${output_prefix}-${base_name}.txt"
    bowtie_file="${out_subdir}/bowtie2out/${base_name}.bowtie2out.txt"

    # 注意：这里将命令写成一行并适当引用路径，方便 parallel 执行
    printf '%s\n' "metaphlan \"$fq_file\" --input_type fastq -o \"$output_file\" --bowtie2out \"$bowtie_file\" --nproc $threads --bowtie2db \"$bowtie2db\" --force --index $index_name" >> "$cmd_file"

done < <(find "$input_root" -type f \( -iname "*.fastq.gz" -o -iname "*.fq.gz" -o -iname "*.fastq" -o -iname "*.fq" \) -print0)

echo "命令列表已生成： $cmd_file"
echo "将使用 parallel -j $parallel_jobs 执行每条命令（请确保系统已安装 GNU parallel）"

# 使用 parallel 执行命令列表
parallel -j "$parallel_jobs" < "$cmd_file"

echo "所有 MetaPhlAn 分析命令执行完毕（或已并发启动）。"
