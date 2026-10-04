# -*- coding: utf-8 -*-
import os
import pandas as pd
import re

# 输入文件夹路径（包含所有 bracken-result-*.txt 文件）
input_dir = './'
# 输出文件夹路径
output_dir = './'
os.makedirs(output_dir, exist_ok=True)

# 获取所有符合命名规则的 Bracken 文件
for filename in os.listdir(input_dir):
    if filename.startswith("bracken-") and ".out" in filename:
        file_path = os.path.join(input_dir, filename)
        
        # 从文件名中提取 sample_id，例如 sample_1
        sample_id_match = re.search(r'CC(\d+)', filename)
        if not sample_id_match:
            continue
        sample_number = sample_id_match.group(1)
        sample_id = f"CC{sample_number}"

        # 读取 Bracken 文件
        df = pd.read_csv(file_path, sep='\t')

        # 筛选 fraction_total_reads > 0.0001
        df_filtered = df[df['fraction_total_reads'] > 0.0001].copy()

        # 构建新 DataFrame，按顺序为：TAXID, RANK, TAXPATH, TAXPATHSN, PERCENTAGE
        df_output = pd.DataFrame()
        df_output['TAXID'] = df_filtered['taxonomy_id']
        df_output['RANK'] = 'species'
        df_output['TAXPATH'] = df_filtered['taxonomy_id']
        df_output['TAXPATHSN'] = df_filtered['name']
        df_output['PERCENTAGE'] = df_filtered['fraction_total_reads'] * 100

        # 构建输出文件名
        output_filename = f"braken-{filename}-gold-file"
        output_path = os.path.join(output_dir, output_filename)

        # 写入文件
        with open(output_path, 'w') as f:
            f.write(f"@SampleID:{sample_id}\n")
            f.write("@Version: 0.9.1\n")
            f.write("@Ranks:superkingdom|phylum|class|order|family|genus|species\n")
            f.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
            df_output.to_csv(f, sep='\t', index=False, header=False)
        
        print(f"Processed: {filename}  {output_filename}")
