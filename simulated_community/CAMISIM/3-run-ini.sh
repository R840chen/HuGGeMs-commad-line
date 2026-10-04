#!/bin/bash

# 设置 Python 脚本路径
SCRIPT_PATH="/ME4012-Vol01/chenc/software/CAMISIM/metagenomesimulation.py"

# 生成 1 到 15 的序列并并行执行
parallel -j 10 "python $SCRIPT_PATH default_config-simulated-{}.ini" ::: {1..16}
