# -*- coding: utf-8 -*-
"""最终核对：比对器准确性 / 菌株特异 marker / 多拷贝 / 注释构成"""
import json, os, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.environ.get('MARKER_OUT_DIR', HERE)
recs = json.load(open(os.path.join(OUTDIR, 'marker_records.json')))

print("="*74)
print("A. 比对器准确性自检（marker 来源基因组 == 目标基因组 时的一致性）")
print("="*74)
SRC = {'W1_GCA_025567015': 'GCA_025567015.1', 'W2_GCA_020687025': 'GCA_020687025.1',
       'D1_GCA_014297375': 'GCA_014297375.1'}
for gk, acc in SRC.items():
    sel = [r for r in recs if r['source'] == acc and r.get(f'{gk}_id') is not None]
    if sel:
        ids = [r[f'{gk}_id'] for r in sel]
        print(f"  {gk}: {len(sel)} 个 marker 来自该基因组自身, 一致性="
              f"{[round(x,4) for x in ids]}")
    else:
        print(f"  {gk}: 无同名来源 marker")

print()
print("="*74)
print("B. 菌株水平 vs 种水平 marker（两株 Waltera 对比）")
print("="*74)
W = [r for r in recs if r['species'] == 'Waltera_2981769']
k1, k2 = 'W1_GCA_025567015', 'W2_GCA_020687025'


def pres(r, k):
    return r.get(f'{k}_id') is not None and r[f'{k}_id'] >= 0.90 and r[f'{k}_cov'] >= 0.90


both = [r for r in W if pres(r, k1) and pres(r, k2)]
o1 = [r for r in W if pres(r, k1) and not pres(r, k2)]
o2 = [r for r in W if pres(r, k2) and not pres(r, k1)]
none = [r for r in W if not pres(r, k1) and not pres(r, k2)]
print(f"  两株共有（种水平保守）: {len(both)}/{len(W)} = {len(both)/len(W)*100:.1f}%")
print(f"  仅株1(GCA_025567015): {len(o1)}   仅株2(GCA_020687025): {len(o2)}"
      f"   → 菌株可变合计 {len(o1)+len(o2)} = {(len(o1)+len(o2))/len(W)*100:.1f}%")
print(f"  两株均缺失: {len(none)} = {len(none)/len(W)*100:.1f}%")
print(f"\n  菌株可变 marker（株1独有）来源: "
      f"{[r['source'] for r in o1]}")
print(f"  菌株可变 marker（株2独有）来源: "
      f"{[r['source'] for r in o2]}")
print(f"  两株均缺失 marker 来源分布: "
      f"{dict(collections.Counter(r['source'] for r in none).most_common())}")

print()
print("  两株共有 marker 的一致性分布:")
for k, nm in ((k1, '株1'), (k2, '株2')):
    ids = np.array([r[f'{k}_id'] for r in both])
    print(f"    {nm}: 中位={np.median(ids):.4f}  P25={np.percentile(ids,25):.4f}  "
          f"min={ids.min():.4f}  =1.0000 的比例={np.mean(ids>=0.9999)*100:.1f}%")

print()
print("="*74)
print("C. 多拷贝 marker（可能造成丰度高估）")
print("="*74)
for gk in ('W1_GCA_025567015', 'W2_GCA_020687025', 'D1_GCA_014297375'):
    sel = [r for r in recs if r.get(f'{gk}_id') is not None
           and r[f'{gk}_id'] >= 0.90 and r[f'{gk}_cov'] >= 0.90]
    multi = [r for r in sel if (r.get(f'{gk}_nloc') or 0) > 1]
    print(f"  {gk}: 存在的 marker {len(sel)}, 其中 >=2 个 >90% 位点的: {len(multi)}"
          f"  ({[r['index'] for r in multi]})")

print()
print("="*74)
print("D. 跨物种特异性（对照）")
print("="*74)
for mk, gks in (('Waltera_2981769', ('D1_GCA_014297375',)),
                ('Dysosmobacter_2763041', ('W1_GCA_025567015', 'W2_GCA_020687025'))):
    sel = [r for r in recs if r['species'] == mk]
    for gk in gks:
        ids = np.array([r[f'{gk}_id'] for r in sel if r.get(f'{gk}_id') is not None])
        n90 = sum(1 for r in sel if r.get(f'{gk}_id') is not None and r[f'{gk}_id'] >= 0.90)
        print(f"  {mk} markers vs {gk}: n={len(ids)}  中位一致性={np.median(ids):.3f}  "
              f"max={ids.max():.3f}  >=90% 的个数={n90}")

print()
print("="*74)
print("E. 注释构成（hypothetical protein 比例）")
print("="*74)
for mk in ('Waltera_2981769', 'Dysosmobacter_2763041'):
    sel = [r for r in recs if r['species'] == mk]
    hy = sum(r['hypothetical'] for r in sel)
    L = np.array([r['len'] for r in sel])
    print(f"  {mk}: {len(sel)} markers, hypothetical protein = {hy} "
          f"({hy/len(sel)*100:.1f}%), 长度 中位={int(np.median(L))} "
          f"(P10={int(np.percentile(L,10))}, P90={int(np.percentile(L,90))})")
