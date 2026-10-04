# -*- coding: utf-8 -*-
"""
HuGGeMs marker gene 在基因组上的定位（纯 numpy seed-and-extend 比对器）
用于回答 R4-C4: marker 在染色体上的位置分布 / 复制不均偏差
"""
import numpy as np, os, sys, re, json, collections

K = 15          # seed 长度
MIN_IDENT = 0.80   # 最低一致性（用于判定"存在"）
MIN_COV = 0.90     # 最低覆盖度


def read_fasta(fp):
    seqs = []
    name = None
    buf = []
    with open(fp, encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs.append((name, ''.join(buf)))
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs.append((name, ''.join(buf)))
    return seqs


_BASE = {65: 0, 67: 1, 71: 2, 84: 3}


def encode(seq):
    a = np.frombuffer(seq.upper().encode('ascii', 'ignore'), dtype=np.uint8)
    c = np.full(a.shape, 4, dtype=np.uint8)
    for ch, v in _BASE.items():
        c[a == ch] = v
    return c


def revcomp_code(c):
    comp = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
    return comp[c[::-1]]


class GenomeIndex:
    def __init__(self, contigs):
        """contigs: list of (name, encoded_array)"""
        self.names = []
        self.codes_list = contigs
        pieces = []
        self.bounds = []   # (start, end) per contig
        off = 0
        sep = np.full(50, 4, dtype=np.uint8)
        for nm, c in contigs:
            self.names.append(nm)
            self.bounds.append((off, off + len(c)))
            pieces.append(c)
            pieces.append(sep)
            off += len(c) + 50
        codes = np.concatenate(pieces)
        self.codes = codes
        self.n = len(codes)
        # contig id per position
        cid = np.full(self.n, -1, dtype=np.int32)
        for i, (s, e) in enumerate(self.bounds):
            cid[s:e] = i
        self.cid = cid
        # k-mer index
        m = self.n - K + 1
        val = np.zeros(m, dtype=np.int64)
        valid = np.ones(m, dtype=bool)
        for i in range(K):
            seg = codes[i:i + m]
            val *= 4
            val += seg.astype(np.int64)
            valid &= (seg < 4)
        val[~valid] = -1
        order = np.argsort(val, kind='stable')
        self.sv = val[order]
        self.order = order
        self.val = val

    def kmer_hits(self, q):
        """返回 query k-mer 值在基因组中的起始位置"""
        if q < 0:
            return np.empty(0, dtype=np.int64)
        lo = np.searchsorted(self.sv, q, 'left')
        hi = np.searchsorted(self.sv, q, 'right')
        return self.order[lo:hi]

    def verify(self, start, qcode, L, shifts=range(-8, 9)):
        best_id = -1.0
        best = None
        for s in shifts:
            p = start + s
            if p < 0 or p + L > self.n:
                continue
            seg = self.codes[p:p + L]
            ok = (seg < 4) & (qcode < 4)
            cov = int(ok.sum())
            if cov < L * 0.5:
                continue
            match = int(((seg == qcode) & ok).sum())
            ident = match / cov
            if ident > best_id:
                best_id = ident
                best = (p, s, match, cov)
        return best_id, best


def encode_query(seq):
    return encode(seq)


def map_marker(idx, qseq, qname):
    """把一个 marker 比对到基因组；返回 (strand, contig, start, end, identity, coverage, n_loci)"""
    results = []
    for strand, qc in (('+', encode_query(qseq)), ('-', revcomp_code(encode_query(qseq)))):
        L = len(qc)
        if L < K + 5:
            continue
        step = max(30, L // 15)
        seeds = list(range(0, L - K + 1, step))
        if seeds[-1] != L - K:
            seeds.append(L - K)
        cand = collections.Counter()
        for so in seeds:
            qk = qc[so:so + K]
            if (qk >= 4).any():
                continue
            qv = 0
            for b in qk:
                qv = qv * 4 + int(b)
            hits = idx.kmer_hits(qv)
            for h in hits:
                cand[int(h) - so] += 1
        if not cand:
            continue
        # 合并相邻候选位点
        starts = sorted(cand.items(), key=lambda x: -x[1])[:8]
        for c, votes in starts:
            ident, best = idx.verify(c, qc, L)
            if best is None:
                continue
            p, s, mt, cov = best
            results.append((strand, p, p + L, ident, cov / L, votes))
    if not results:
        return None
    results.sort(key=lambda x: -x[3])
    return results


def contig_of(idx, pos):
    cid = idx.cid[pos]
    if cid < 0:
        return None, None, None
    s, e = idx.bounds[cid]
    return idx.names[cid], pos - s, e - s


def gc_skew_oric(codes):
    """用累积 GC skew 估计 oriC（在最长 contig 上）"""
    c = codes
    ok = c < 4
    g = ((c == 2) & ok).astype(np.float64)
    cc = ((c == 1) & ok).astype(np.float64)
    # 累积 skew，用窗口平滑
    win = 10000
    n = len(c)
    if n < win * 3:
        return None
    ng = g[:(n // win) * win].reshape(-1, win).sum(axis=1)
    nc = cc[:(n // win) * win].reshape(-1, win).sum(axis=1)
    skew = (ng - nc) / np.maximum(ng + nc, 1)
    cum = np.cumsum(skew)
    # 切换点：cum 的最大值到最小值的下降段的中点
    i_max = int(np.argmax(cum))
    i_min = int(np.argmin(cum))
    mid = ((i_max + i_min) // 2) * win
    return mid, skew, cum


if __name__ == '__main__':
    print('mapper module ready')
