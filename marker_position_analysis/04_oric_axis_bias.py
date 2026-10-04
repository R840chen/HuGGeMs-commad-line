# -*- coding: utf-8 -*-
"""
Part 2: 用 GC skew 给 contig 定向，重建 oriC->ter 轴，
        计算 marker 沿复制轴的分布 + 复制不均导致的丰度偏差估计 (R4-C4)
"""
import sys, os, json, collections, itertools, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mapper import read_fasta

# Input FASTA files live in $MARKER_DATA_DIR (default: ./data); outputs go to $MARKER_OUT_DIR (default: this directory)
BASE = os.environ.get('MARKER_DATA_DIR', os.path.join(HERE, 'data'))
OUTDIR = os.environ.get('MARKER_OUT_DIR', HERE)
GENOMES = {
    'W1': ('GCA_025567015.1_ASM2556701v1_genomic.fna', 'Waltera_2981769',
           'W1_GCA_025567015'),
    'W2': ('GCA_020687025.1_ASM2068702v1_genomic.fna', 'Waltera_2981769',
           'W2_GCA_020687025'),
    'D1': ('GCA_014297375.1_ASM1429737v1_genomic.fna', 'Dysosmobacter_2763041',
           'D1_GCA_014297375'),
}
recs = json.load(open(os.path.join(OUTDIR, 'marker_records.json')))
MIN_CONTIG = 20000


def skew_profile(seq, nwin=30):
    w = max(2000, len(seq) // nwin)
    n = len(seq) // w
    s = seq[:n * w].upper()
    G = np.array([s[i * w:(i + 1) * w].count('G') for i in range(n)], dtype=float)
    C = np.array([s[i * w:(i + 1) * w].count('C') for i in range(n)], dtype=float)
    return (G - C) / np.maximum(G + C, 1), w


out = {}
for gk, (gf, mk, key) in GENOMES.items():
    contigs = {n.split()[0]: s for n, s in read_fasta(os.path.join(BASE, gf))}
    # 收集 marker
    hits = []
    for r in recs:
        if r['species'] != mk: continue
        if r.get(f'{key}_id') is None: continue
        if r[f'{key}_id'] < 0.90 or r[f'{key}_cov'] < 0.90: continue
        hits.append({'contig': r[f'{key}_contig'], 'pos': r[f'{key}_pos']})
    print(f"\n{'='*72}\n{gk} {mk}")

    # 定向
    info = []
    used_len = 0
    for c, s in contigs.items():
        L = len(s)
        if L < MIN_CONTIG:
            continue
        sk, w = skew_profile(s)
        cum = np.cumsum(sk)
        if len(cum) < 4:
            continue
        x = np.arange(len(cum))
        slope = np.polyfit(x, cum, 1)[0]
        flip = slope < 0
        if flip:
            sk2 = -sk[::-1]
            cum2 = np.cumsum(sk2)
        else:
            cum2 = cum
        info.append({'contig': c, 'len': L, 'flip': bool(flip),
                     'level': float(np.mean(cum2)), 'slope': float(abs(slope)),
                     'cum_start': float(cum2[0]), 'cum_end': float(cum2[-1])})
        used_len += L
    # 按 level 排序 => oriC(低) -> ter(高)
    info.sort(key=lambda d: d['level'])
    offset = 0
    coord = {}
    for d in info:
        d['offset'] = offset
        coord[d['contig']] = d
        offset += d['len']
    total = offset
    print(f"  定向成功 contig: {len(info)}/{len(contigs)}  覆盖 {used_len}/{sum(len(s) for s in contigs.values())} bp "
          f"({used_len/sum(len(s) for s in contigs.values())*100:.0f}%)")

    # marker 沿轴的位置
    xs = []
    for h in hits:
        d = coord.get(h['contig'])
        if not d: continue
        p = (d['len'] - h['pos']) if d['flip'] else h['pos']
        xs.append((d['offset'] + p) / total)
    xs = np.array(sorted(xs))
    n = len(xs)
    print(f"  落在重建轴上的 marker: {n}/{len(hits)}")
    if n == 0: continue
    m = len(xs)
    dplus = (np.arange(1, m + 1) / m - xs).max()
    dminus = (xs - np.arange(0, m) / m).max()
    D = max(dplus, dminus)
    en = math.sqrt(m)
    lam = (en + 0.12 + 0.11 / en) * D
    pks = min(max(2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam)
                          for k in range(1, 101)), 0.0), 1.0)
    print(f"  沿 oriC->ter 轴归一化位置: KS D={D:.3f}, p={pks:.3f}  (n={m})")
    hist, _ = np.histogram(xs, bins=10, range=(0, 1))
    print(f"  十分位(oriC->ter): {list(hist.astype(int))}")
    print(f"  均值位置={xs.mean():.3f} (均匀期望 0.500)   中位={np.median(xs):.3f}")
    print(f"  oriC 侧一半 (x<0.5): {int((xs<0.5).sum())}  ter 侧一半: {int((xs>=0.5).sum())}")

    # 复制梯度模型: C(x) = R^(1-x)，x=0(oriC) -> 1(ter)
    print("  复制梯度模型 C(x)=R^(1-x) 下的相对丰度偏差:")
    # 基因组平均（按连续近似）
    res = {}
    for R in (1.0, 2.0, 4.0, 8.0):
        if R == 1.0:
            bias = 1.0
        else:
            Cm = np.mean(R ** (1 - xs))
            mean_genome = (R - 1) / math.log(R)      # ∫0^1 R^(1-x) dx
            bias = Cm / mean_genome
        res[R] = bias
        print(f"     oriC:ter = {R:.0f}:1  ->  偏差 = {(bias-1)*100:+.2f}%")
    # 如果所有 marker 都挤在最靠近 oriC 的 5% 区域（极端反例，用于对照）
    print(f"     [对照] 若 200 markers 全部集中在最靠近 oriC 的 5% 区域:")
    for R in (2.0, 4.0, 8.0):
        xs_ext = np.linspace(0.0, 0.05, 200)
        bias = np.mean(R ** (1 - xs_ext)) / ((R - 1) / math.log(R))
        print(f"        oriC:ter = {R:.0f}:1  ->  偏差 = {(bias-1)*100:+.1f}%")

    out[gk] = {'n_on_axis': n, 'ks_D': float(D), 'ks_p': float(pks),
               'hist': [int(v) for v in hist], 'mean_pos': float(xs.mean()),
               'bias': {str(k): float(v) for k, v in res.items()},
               'contigs_used': len(info), 'cov_frac': float(used_len / total)}

json.dump(out, open(os.path.join(OUTDIR, 'oric_bias.json'), 'w'), indent=1)
print("\n已保存 oric_bias.json")
