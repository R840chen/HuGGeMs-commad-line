# -*- coding: utf-8 -*-
import os

# 设置 .ini 文件所在目录
ini_directory = "./"  # 可以替换为你的实际路径

# 替换规则：键 -> 新路径
replacements = {
    "samtools=": "samtools=/ME4012-Vol01/chenc/software/CAMISIM/tools/samtools-1.3/samtools",
    "readsim=": "readsim=/ME4012-Vol01/chenc/software/CAMISIM/tools/art_illumina-2.3.6/art_illumina",
    "error_profiles=": "error_profiles=/ME4012-Vol01/chenc/software/CAMISIM/tools/art_illumina-2.3.6/profiles/",
    "strain_simulation_template=": "strain_simulation_template=/ME4012-Vol01/chenc/software/CAMISIM/scripts/StrainSimulationWrapper/sgEvolver/",
    "ncbi_taxdump=": "ncbi_taxdump=/ME4012-Vol01/chenc/software/CAMISIM/tools/ncbi-taxonomy_20170222.tar.gz"
}

# 遍历文件夹中所有以 .ini 结尾的文件
for filename in os.listdir(ini_directory):
    if filename.endswith(".ini"):
        filepath = os.path.join(ini_directory, filename)

        # 读取文件内容
        with open(filepath, 'r') as file:
            lines = file.readlines()

        # 替换路径
        new_lines = []
        for line in lines:
            modified = False
            for key, new_value in replacements.items():
                if line.startswith(key):
                    new_lines.append(new_value + "\n")
                    modified = True
                    break
            if not modified:
                new_lines.append(line)

        # 写回修改后的内容
        with open(filepath, 'w') as file:
            file.writelines(new_lines)

        print(f"Updated: {filename}")

print("All INI files updated successfully.")
