# -*- coding: utf-8 -*-
import pandas as pd
import re
import os
from glob import glob

# 只保留这七种rank
valid_ranks = [
    "superkingdom", "phylum", "class", "order", "family", "genus", "species"
]

# rank mapping
rank_map = {
    0: "superkingdom",
    1: "phylum",
    2: "class",
    3: "order",
    4: "family",
    5: "genus",
    6: "species"
}

# 遍历所有符合条件的输入文件
input_files = glob("mpa4-3-synthetic-sample-mpa4-result-*.clean.txt")

for input_file in input_files:
    # 提取中间部分
    match = re.search(r'mpa4-3-synthetic-sample-mpa4-result-(.*?).clean\.txt', input_file)
    if not match:
        continue
    mid_str = match.group(1)
    sample_id = f"{mid_str}"
    output_file = f"cami2_gold_standard_profile-{mid_str}.txt"

    # 读取MetaPhlAn结果
    df = pd.read_csv(input_file, sep="\t", comment="#", header=None)
    df.columns = ["clade_name", "NCBI_tax_id", "relative_abundance", "additional_species"]

    output_rows = []

    for _, row in df.iterrows():
        clade = row["clade_name"]
        ncbi_ids = str(row["NCBI_tax_id"])
        abundance = float(row["relative_abundance"])

        # 跳过菌株级别（带 t__ 或 多于7级）
        if clade.count("|") > 6 or "|t__" in clade:
            continue

        taxpath = ncbi_ids
        taxid = taxpath.strip().split("|")[-1]

        depth = clade.count("|")
        if depth > 6:
            continue  # 只保留到 species
        rank = rank_map[depth]

        name_parts = []
        for part in clade.split("|"):
            clean = re.sub(r'^[a-z]__', '', part)
            clean = clean.replace("_", " ")
            name_parts.append(clean)
        taxpathsn = "|".join(name_parts)

        output_rows.append([
            taxid, rank, taxpath, taxpathsn, f"{abundance:.4f}", "", ""
        ])

    # 写入输出文件
    with open(output_file, "w", encoding="utf-8") as out:
        out.write(f"@SampleID:{sample_id}\n")
        out.write(f"@Version:0.9.1\n")
        out.write(f"@Ranks:superkingdom|phylum|class|order|family|genus|species\n")
        out.write("\n")
        out.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n")
        for row in output_rows:
            out.write("\t".join(row) + "\n")

    print(f" Converted: {input_file} --> {output_file}")
