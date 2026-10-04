# -*- coding: utf-8 -*-
"""诊断：按 marker 来源基因组分组，看一致性与存在率"""
import sys, os, re, time, json, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mapper import *

# Input FASTA files live in $MARKER_DATA_DIR (default: ./data); outputs go to $MARKER_OUT_DIR (default: this directory)
BASE = os.environ.get('MARKER_DATA_DIR', os.path.join(HERE, 'data'))
OUTDIR = os.environ.get('MARKER_OUT_DIR', HERE)
GENOMES = {
    'W1_GCA_025567015': 'GCA_025567015.1_ASM2556701v1_genomic.fna',
    'W2_GCA_020687025': 'GCA_020687025.1_ASM2068702v1_genomic.fna',
    'D1_GCA_014297375': 'GCA_014297375.1_ASM1429737v1_genomic.fna',
}
MARKERS = {
    'Waltera_2981769': ('Brotolimicola_acetigignens_species_2981769.ffn',
                        ['W1_GCA_025567015', 'W2_GCA_020687025']),
    'Dysosmobacter_2763041': ('Dysosmobacter_hominis_species_2763041.ffn',
                              ['D1_GCA_014297375']),
}

idxs = {}
for gk, gf in GENOMES.items():
    contigs = [(n.split()[0], encode(s)) for n, s in read_fasta(os.path.join(BASE, gf))]
    idxs[gk] = GenomeIndex(contigs)

REC = []
for mk, (mf, own) in MARKERS.items():
    mks = read_fasta(os.path.join(BASE, mf))
    for i, (nm, sq) in enumerate(mks):
        # 解析来源基因组
        src = nm.split()[0]
        m = re.match(r'^(GCA_\d+\.\d+)', src)
        source = m.group(1) if m else src.split('_')[0]
        hypo = 'hypothetical protein' in nm
        rec = {'species': mk, 'index': i, 'header': nm, 'source': source,
               'len': len(sq), 'hypothetical': hypo}
        for gk in GENOMES:
            res = map_marker(idxs[gk], sq, nm)
            if res:
                b = res[0]
                rec[f'{gk}_id'] = round(b[3], 4)
                rec[f'{gk}_cov'] = round(b[4], 3)
                rec[f'{gk}_strand'] = b[0]
                cn, cs, cl = contig_of(idxs[gk], b[1])
                rec[f'{gk}_contig'] = cn
                rec[f'{gk}_pos'] = cs
                rec[f'{gk}_len_c'] = cl
                rec[f'{gk}_nloc'] = len([r for r in res if r[3] >= 0.90 and r[4] >= 0.90])
            else:
                rec[f'{gk}_id'] = None
        REC.append(rec)

json.dump(REC, open(os.path.join(OUTDIR, 'marker_records.json'), 'w'), indent=1)
print(f"记录数: {len(REC)}")

for mk, (mf, own) in MARKERS.items():
    sub = [r for r in REC if r['species'] == mk]
    print(f"\n===== {mk} =====")
    print(f"  marker 数: {len(sub)}  长度: min={min(r['len'] for r in sub)} "
          f"中位={int(np.median([r['len'] for r in sub]))} max={max(r['len'] for r in sub)}")
    print(f"  hypothetical protein 占比: {sum(r['hypothetical'] for r in sub)}/{len(sub)} "
          f"= {sum(r['hypothetical'] for r in sub)/len(sub)*100:.1f}%")
    # 来源分布
    print("  来源基因组分布(前6):")
    for s, c in collections.Counter(r['source'] for r in sub).most_common(6):
        print(f"     {s}: {c}")
    for gk in own:
        ids = np.array([r[f'{gk}_id'] if r[f'{gk}_id'] is not None else -1 for r in sub])
        pres = [r for r in sub if r[f'{gk}_id'] is not None and r[f'{gk}_id'] >= 0.90
                and r[f'{gk}_cov'] >= 0.90]
        print(f"  vs {gk}: 存在 {len(pres)}/{len(sub)}   中位一致性={np.median(ids[ids>0]):.3f}")
        # 按来源分组看低一致性
        low = [r for r in sub if r[f'{gk}_id'] is None or r[f'{gk}_id'] < 0.90]
        print(f"     低一致性/缺失 {len(low)} 个，其来源分布: "
              f"{dict(collections.Counter(r['source'] for r in low).most_common(5))}")

# 两株 Waltera 的存在重叠
W = [r for r in REC if r['species'] == 'Waltera_2981769']
both = sum(1 for r in W if r['W1_GCA_025567015_id'] and r['W1_GCA_025567015_id'] >= 0.9
           and r['W2_GCA_020687025_id'] and r['W2_GCA_020687025_id'] >= 0.9)
only1 = sum(1 for r in W if r['W1_GCA_025567015_id'] and r['W1_GCA_025567015_id'] >= 0.9
            and not (r['W2_GCA_020687025_id'] and r['W2_GCA_020687025_id'] >= 0.9))
only2 = sum(1 for r in W if r['W2_GCA_020687025_id'] and r['W2_GCA_020687025_id'] >= 0.9
            and not (r['W1_GCA_025567015_id'] and r['W1_GCA_025567015_id'] >= 0.9))
print(f"\n===== 两株 Waltera 之间的 marker 存在性 =====")
print(f"  两株都有(=种水平保守): {both}")
print(f"  仅 W1(GCA_025567015) 有: {only1}")
print(f"  仅 W2(GCA_020687025) 有: {only2}")
print(f"  两株都无: {len(W)-both-only1-only2}")
