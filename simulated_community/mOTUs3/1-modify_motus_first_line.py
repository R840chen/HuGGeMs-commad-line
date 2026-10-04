# -*- coding: utf-8 -*-
import os

# 要处理的文件名范围
for n in range(26):
    filename = f"motus-recall-result-2017.12.04_18.45.54_sample_{n}-anonymous_reads.txt"
    if os.path.exists(filename):
        # 构建替换的第一行内容
        new_first_line = f"@SampleID:2017.12.04_18.45.54_sample_{n}-result-reads\n"
        # 读取原文件内容
        with open(filename, 'r') as f:
            lines = f.readlines()
        # 替换第一行
        lines[4] = new_first_line
        # 写回原文件
        with open(filename, 'w') as f:
            f.writelines(lines)
        print(f"Updated: {filename}")
    else:
        print(f"File not found: {filename}")
