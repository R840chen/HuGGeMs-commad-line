#!/bin/bash
set -euo pipefail

########################
# 参数设置
########################

# 输入文件目录
input_dir="/mnt/ME5012-Vol01/chenc/CAMII-simulated-reads/all-simulated-reads/20250726-simulated-reads"

# 输出目录
output_dir="/mnt/ME5012-Vol01/chenc/MetaphIan/paper-data/2.CAMI2-profiling/4.mOTUs-result/precision"

# 输出文件前缀
output_prefix="motus-recall-result"

# 每个 motus 任务使用的线程数
threads=50

# 同时运行的最大任务数（?? 很重要）
max_jobs=6

########################
# 准备工作
########################

cd "$input_dir" || exit 1
mkdir -p "$output_dir"

echo "Input dir : $input_dir"
echo "Output dir: $output_dir"
echo "Threads per job: $threads"
echo "Max parallel jobs: $max_jobs"
echo "------------------------------------"

########################
# 并行执行函数
########################

run_motus () {
    local fq_file="$1"
    local base_name
    base_name="$(basename "${fq_file%%.fq*}")"
    local output_file="${output_dir}/${output_prefix}-${base_name}.txt"

    echo "[START] $(date '+%F %T')  $fq_file"

    motus profile \
        -g 6 \
        -t "$threads" \
        -s "$fq_file" \
        -o "$output_file"

    echo "[DONE ] $(date '+%F %T')  $fq_file"
}

########################
# 主循环（并发控制）
########################

job_count=0

for fq_file in *.fq *.fq.gz; do
    [ -e "$fq_file" ] || continue

    run_motus "$fq_file" &

    ((job_count++))

    # 达到最大并发数就等待
    if (( job_count % max_jobs == 0 )); then
        wait
    fi
done

# 等待剩余任务
wait

echo "===================================="
echo "所有 mOTUs recall 分析完成！"
