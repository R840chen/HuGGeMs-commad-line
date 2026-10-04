# -*- coding: utf-8 -*-
"""Part 1: marker 位置分布 / 聚集度 / GC skew 复制方向"""
import sys, os, json, collections, itertools, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mapper import read_fasta, encode

# Input FASTA files live in $MARKER_DATA_DIR (default: ./data); outputs go to $MARKER_OUT_DIR (default: this directory)
BASE = os.environ.get('MARKER_DATA_DIR', os.path.join(HERE, 'data'))
OUTDIR = os.environ.get('MARKER_OUT_DIR', HERE)
GENOMES = {
    'W1': ('GCA_025567015.1_ASM2556701v1_genomic.fna', 'Waltera_2981769'),
    'W2': ('GCA_020687025.1_ASM2068702v1_genomic.fna', 'Waltera_2981769'),
    'D1': ('GCA_014297375.1_ASM1429737v1_genomic.fna', 'Dysosmobacter_2763041'),
}
recs = json.load(open(os.path.join(OUTDIR, 'marker_records.json')))


def ks_pvalue(d, n):
    """Kolmogorov-Smirnov 渐近 p 值 (单样本 vs 均匀分布 D 统计量)"""
    if n < 1: return None, None
    en = math.sqrt(n)
    lam = (en + 0.12 + 0.11 / en) * d
    p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101))
    return d, min(max(p, 0.0), 1.0)


def chi2_pvalue(x2, df):
    """Wilson-Hilferty 近似"""
    if df <= 0: return None
    z = ((x2 / df) ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
    # 正态尾概率
    return 0.5 * math.erfc(z / math.sqrt(2))


summary = {}
for gk, (gf, mk) in GENOMES.items():
    contigs = read_fasta(os.path.join(BASE, gf))
    cmap = {n.split()[0]: s for n, s in contigs}
    clen = {n: len(s) for n, s in cmap.items()}
    sys.stderr.write(f"\n===== {gk} ({mk}) =====\n")

    # 收集该 genome 下所有"存在"的 marker 位置
    hits = []
    for r in recs:
        if r['species'] != mk: continue
        key = {'W1': 'W1_GCA_025567015', 'W2': 'W2_GCA_020687025',
               'D1': 'D1_GCA_014297375'}[gk]
        if r.get(f'{key}_id') is None: continue
        if r[f'{key}_id'] < 0.90 or r[f'{key}_cov'] < 0.90: continue
        hits.append({'marker': r['index'], 'contig': r[f'{key}_contig'],
                     'start': r[f'{key}_pos'], 'len_c': r[f'{key}_len_c'],
                     'strand': r[f'{key}_strand'], 'ident': r[f'{key}_id'],
                     'nloc': r[f'{key}_nloc']})
    n = len(hits)
    hits = sorted(hits, key=lambda x: (x['contig'], x['start']))
    ncontig_with = len(set(h['contig'] for h in hits))
    print(f"\n{'='*70}\n{gk}  ({mk})  基因组 {len(clen)} contigs / {sum(clen.values())} bp")
    print(f"  定位成功的 marker: {n}  分布在 {ncontig_with} 个 contig 上")

    # --- 每个 contig 的 marker 数 vs 长度（均匀性检验）---
    obs = collections.Counter(h['contig'] for h in hits)
    # 只对长度 >= 20kb 的 contig 做检验（排除碎片造成的统计噪声）
    big = [c for c, L in clen.items() if L >= 20000]
    exp_tot = sum(obs.get(c, 0) for c in big)
    Ltot = sum(clen[c] for c in big)
    x2 = 0.0
    for c in big:
        e = exp_tot * clen[c] / Ltot
        o = obs.get(c, 0)
        if e > 0:
            x2 += (o - e) ** 2 / e
    df = len(big) - 1
    print(f"  均匀性检验(contig 层面, 仅 >20kb 的 {len(big)} 条): "
          f"chi2={x2:.1f}, df={df}, p={chi2_pvalue(x2, df):.3f}")

    # --- contig 内归一化位置：KS 检验 & 直方图 ---
    norm = []
    for h in hits:
        if h['len_c'] and h['len_c'] > 0:
            norm.append(h['start'] / h['len_c'])
    norm = np.array(sorted(norm))
    # KS 统计量 vs 均匀
    m = len(norm)
    dplus = (np.arange(1, m + 1) / m - norm).max()
    dminus = (norm - np.arange(0, m) / m).max()
    D = max(dplus, dminus)
    D, p = ks_pvalue(D, m)
    print(f"  contig 内归一化位置 KS 检验: D={D:.3f}, p={p:.3f}  (n={m})")
    hist, _ = np.histogram(norm, bins=10, range=(0, 1))
    print(f"  十分位分布 (0-100%): {list(hist)}")

    # --- 聚集度 / 位点数 ---
    loci = 0
    for c, grp in itertools.groupby(sorted(hits, key=lambda x: (x['contig'], x['start'])),
                                   key=lambda x: x['contig']):
        pos = sorted(h['start'] for h in grp)
        loci += 1
        for a, b in zip(pos, pos[1:]):
            if b - a > 5000:
                loci += 1
    print(f"  位点数(间隔>5kb 计为新位点): {loci}")
    # 每 10kb / 20kb 窗口最大 marker 数
    for win in (10000, 20000, 50000):
        mx = 0
        for c, grp in itertools.groupby(hits, key=lambda x: x['contig']):
            pos = sorted(h['start'] for h in grp)
            for i, p0 in enumerate(pos):
                cnt = sum(1 for p in pos if p0 <= p < p0 + win)
                mx = max(mx, cnt)
        print(f"  任一个 {win//1000}kb 窗口内的最大 marker 数: {mx} "
              f"({mx/n*100:.1f}% of {n})")
    # 相邻 marker 间距
    gaps = []
    for c, grp in itertools.groupby(hits, key=lambda x: x['contig']):
        pos = sorted(h['start'] for h in grp)
        gaps += [b - a for a, b in zip(pos, pos[1:])]
    gaps = np.array(gaps)
    if len(gaps):
        print(f"  同 contig 相邻 marker 间距(bp): 中位={int(np.median(gaps))} "
              f"P25={int(np.percentile(gaps,25))} P75={int(np.percentile(gaps,75))}")
    # 覆盖跨度
    spans = []
    for c, grp in itertools.groupby(hits, key=lambda x: x['contig']):
        grp = list(grp)
        L = next(h['len_c'] for h in grp)
        spans.append((max(h['start'] for h in grp) - min(h['start'] for h in grp)) / L)
    print(f"  每条 contig 内 marker 覆盖率(跨度/contig长): 中位={np.median(spans):.2f}")

    # --- strand 分布 ---
    st = collections.Counter(h['strand'] for h in hits)
    print(f"  链分布: {dict(st)}")

    # --- GC skew：marker 所在局部区域（±5kb）的 skew 符号 ---
    def gcskew(seq, a, b):
        s = seq[max(0, a):min(len(seq), b)].upper()
        g = s.count('G'); c = s.count('C')
        return (g - c) / (g + c) if (g + c) > 0 else 0.0
    sk = [gcskew(cmap[h['contig']], h['start'] - 5000, h['start'] + 5000) for h in hits]
    sk = np.array(sk)
    print(f"  局部 GC skew: 正={int((sk>0).sum())} 负={int((sk<0).sum())} "
          f"≈0={int((sk==0).sum())}  中位={np.median(sk):+.3f}")

    # --- 多拷贝 ---
    multi = [h for h in hits if (h['nloc'] or 0) > 1]
    print(f"  多命中(>=2 个 >=90% 位点) marker: {len(multi)}")

    # --- 每条 contig 明细（仅 >20kb）---
    print("  各 contig 明细 (>20kb):")
    for c, grp in itertools.groupby(hits, key=lambda x: x['contig']):
        grp = list(grp)
        L = grp[0]['len_c']
        if L < 20000:
            continue
        pos = sorted(h['start'] for h in grp)
        print(f"     {c}  len={L:>7}  markers={len(pos):>3}  "
              f"密度={len(pos)/L*1e5:.1f}/100kb  跨度={pos[0]}-{pos[-1]}")

    summary[gk] = {'n': n, 'ncontig': ncontig_with, 'loci': loci,
                   'ks_D': D, 'ks_p': p, 'chi2': x2, 'chi2_p': chi2_pvalue(x2, df),
                   'hist': list(hist), 'skew_pos': int((sk > 0).sum()),
                   'skew_neg': int((sk < 0).sum()), 'multi': len(multi),
                   'gaps_median': int(np.median(gaps)) if len(gaps) else None}

json.dump(summary, open(os.path.join(OUTDIR, 'pos_summary.json'), 'w'), indent=1,
          default=lambda o: int(o) if hasattr(o, 'item') else str(o))
print("\n已保存 pos_summary.json")
