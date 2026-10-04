#!/bin/bash
set -euo pipefail

########################
# 参数设置
########################

# 输入文件目录
input_dir="/ME5012-Vol01/chenc/CAMII-simulated-reads/CAMI2-reads-from-gut/all-simulated-reads"

# 输出目录
output_dir="/ME5012-Vol01/chenc/MetaphIan/paper-data/2.CAMI2-profiling/7.motus4-result/recall"

# 输出文件前缀
output_prefix="motus-recall-result"

# 每个 mOTUs 任务使用的线程数
threads=50

# 同时运行的最大任务数
max_jobs=6

########################
# 准备工作
########################

mkdir -p "$output_dir"

echo "===================================="
echo "mOTUs version:"
motus profile -h | head -5
echo "===================================="

echo "Input dir       : $input_dir"
echo "Output dir      : $output_dir"
echo "Threads per job : $threads"
echo "Max parallel    : $max_jobs"
echo "===================================="

########################
# 并行执行函数
########################

run_motus () {

    local fq_file="$1"

    # 去掉路径和 .fq / .fq.gz 后缀
    local base_name
    base_name="$(basename "$fq_file")"
    base_name="${base_name%.fq.gz}"
    base_name="${base_name%.fq}"

    local output_file="${output_dir}/${output_prefix}-${base_name}.txt"

    echo "[START] $(date '+%F %T')"
    echo "        Input : $fq_file"
    echo "        Output: $output_file"

    motus profile \
        -s "$fq_file" \
        -g 1 \
        -t "$threads" \
        -o "$output_file"

    echo "[DONE ] $(date '+%F %T')  $fq_file"
}

########################
# 主循环
########################

job_count=0

cd "$input_dir"

for fq_file in *.fq.gz; do

    [ -e "$fq_file" ] || continue

    run_motus "$fq_file" &

    ((job_count+=1))

    # 达到最大并发数后等待
    if (( job_count % max_jobs == 0 )); then
        wait
    fi

done

# 等待最后一批任务
wait

echo "===================================="
echo "所有 mOTUs 4.1.0 分析完成！"
echo "===================================="