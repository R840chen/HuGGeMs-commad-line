# -*- coding: utf-8 -*-
"""GC 组成 & 局部 GC skew（相对 contig 背景）——位置偏好性检验"""
import sys, os, json, collections, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mapper import read_fasta

# Input FASTA files live in $MARKER_DATA_DIR (default: ./data); outputs go to $MARKER_OUT_DIR (default: this directory)
BASE = os.environ.get('MARKER_DATA_DIR', os.path.join(HERE, 'data'))
OUTDIR = os.environ.get('MARKER_OUT_DIR', HERE)
GENOMES = {
    'W1': ('GCA_025567015.1_ASM2556701v1_genomic.fna', 'Waltera_2981769', 'W1_GCA_025567015'),
    'W2': ('GCA_020687025.1_ASM2068702v1_genomic.fna', 'Waltera_2981769', 'W2_GCA_020687025'),
    'D1': ('GCA_014297375.1_ASM1429737v1_genomic.fna', 'Dysosmobacter_2763041', 'D1_GCA_014297375'),
}
MARKERS = {'Waltera_2981769': 'Brotolimicola_acetigignens_species_2981769.ffn',
           'Dysosmobacter_2763041': 'Dysosmobacter_hominis_species_2763041.ffn'}
recs = json.load(open(os.path.join(OUTDIR, 'marker_records.json')))


def gc(s):
    s = s.upper()
    g = s.count('G'); c = s.count('C')
    return (g + c) / len(s) * 100 if s else 0


print("=== 基因组 GC% vs marker GC% ===")
for mk, mf in MARKERS.items():
    mks = read_fasta(os.path.join(BASE, mf))
    mks = [s for _, s in mks]
    mg = np.array([gc(s) for s in mks])
    print(f"  {mk}: marker GC% 中位={np.median(mg):.1f}  "
          f"(P10={np.percentile(mg,10):.1f}, P90={np.percentile(mg,90):.1f})")
for gk, (gf, mk, key) in GENOMES.items():
    seqs = read_fasta(os.path.join(BASE, gf))
    tot = sum(gc(s) * len(s) for _, s in seqs) / sum(len(s) for _, s in seqs)
    print(f"  {gk} ({mk}) 基因组 GC% = {tot:.1f}")

print()
print("=== 局部 GC skew（marker ±5kb）相对 contig 背景 ===")
for gk, (gf, mk, key) in GENOMES.items():
    cmap = {n.split()[0]: s for n, s in read_fasta(os.path.join(BASE, gf))}
    contig_sk = {}
    for c, s in cmap.items():
        s = s.upper(); n = len(s)
        if n < 10000: continue
        w = 5000
        G = sum(s[i:i+w].count('G') for i in range(0, n - w, w))
        C = sum(s[i:i+w].count('C') for i in range(0, n - w, w))
        if G + C: contig_sk[c] = (G - C) / (G + C)
    diffs = []
    for r in recs:
        if r['species'] != mk: continue
        if r.get(f'{key}_id') is None or r[f'{key}_id'] < 0.90: continue
        c = r[f'{key}_contig']; p = r[f'{key}_pos']
        if c not in contig_sk: continue
        G = cmap[c][max(0, p-5000):p+5000].upper().count('G')
        C = cmap[c][max(0, p-5000):p+5000].upper().count('C')
        if G + C == 0: continue
        diffs.append((G - C) / (G + C) - contig_sk[c])
    diffs = np.array(diffs)
    # Wilcoxon 符号秩的简易 z 检验
    pos = (diffs > 0).sum(); neg = (diffs < 0).sum(); nn = pos + neg
    z = (pos - nn / 2) / np.sqrt(nn / 4) if nn else 0
    print(f"  {gk}: n={len(diffs)}  局部skew−contig背景 中位={np.median(diffs):+.4f}  "
          f"正={pos} 负={neg}  (符号检验 z={z:+.2f})")
