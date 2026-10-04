#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================================
 批量把 mOTUs4 .relab 结果转成 CAMI 格式（给 OPAL 比对）
============================================================================
 与 motu2cami.py 的区别:
   - motu2cami.py        : 输入已是 motu_ncbi_taxid.tsv（单文件）
   - motu_batch2cami.py  : 输入是原始 .relab（整个目录批量），
                           自己完成 mOTU → 代表基因组 → NCBI taxid → CAMI 全流程

 完整流程:
   1. 读 .relab，提取 (mOTU, Taxonomy, Abundance)，剔除 unassigned
   2. 读 mOTUs 代表基因组表 (mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz)
        → mOTU → 代表基因组 GENOME
   3. 从代表基因组名解析 NCBI 线索:
        *_GCF-024329905-V1_GENO_*   → GCF_024329905.1 → assembly_summary 查 taxid
        *_SAMN15871286_MAG_*        → SAMN15871286    → biosample 列查 accession
        其它                          → 无 NCBI 线索 → taxid 留空
   4. 用 taxdump 补全完整谱系，写出 CAMI 文件
   5. taxid 找不到的行 → 保留但 TAXID 留空（--missing-mode keep）
        （用户明确要求：找不到就空着，OPAL 会算作假阳性，符合预期）

 用法:
   python3 motu_batch2cami.py \
       --input-dir "D:/Desktop/.../2026.9.25之后添加motus的分析" \
       --output-dir "D:/Desktop/.../cami_output" \
       --rep-file /ME4012-Vol01/database/motus4/db_mOTU/mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz \
       --ncbi-dir /ME4012-Vol01/database/NCBI_tax \
       --biosample-map /ME4012-Vol01/database/NCBI_tax/assembly_summary_genbank.txt \
       --biosample-map /ME4012-Vol01/database/NCBI_tax/assembly_summary_refseq.txt \
       --missing-mode keep \
       --recursive

 输出:
   <output-dir>/<组名>/<子目录>/<样本名>.cami      每个样本一个 CAMI 文件
   <output-dir>/_batch_summary.tsv                  总汇总表
   <output-dir>/_unmapped_all.tsv                   所有未映射 mOTU 汇总
============================================================================
"""
import os
import re
import sys
import gzip
import glob
import argparse
from collections import OrderedDict, Counter

# 复用 motu2cami.py 里的核心类与函数
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motu2cami import (  # noqa: E402
    CAMI_RANKS, NCBI_RANK_MAP, Taxdump, AssemblySummary,
    genome_id_to_lookup, cami_rank_of,
)


# ============================== .relab 读取 ==============================

def read_relab(path):
    """读 mOTUs .relab 文件。
    格式:
        第1行: #tool_version=... database_version=...
        第2行: mOTU\tTaxonomy\t<sample name>
        第3行起: <mOTU id>\t<GTDB taxonomy>\t<abundance>

    返回 rows: [{motu, tax_str, gtdb_species, gtdb_genus, abundance}]
    """
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()          # #tool_version...
        header = fh.readline()         # mOTU / Taxonomy / sample
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            motu = cols[0].strip()
            tax_str = cols[1].strip()
            try:
                abd = float(cols[2])
            except ValueError:
                abd = 0.0

            # 剔除 unassigned 汇总行
            if "unassigned" in motu.lower():
                continue

            # 解析 GTDB 分类串
            gtdb_species = ""
            gtdb_genus = ""
            for part in tax_str.split(";"):
                if part.startswith("s__"):
                    gtdb_species = part[3:].strip()
                elif part.startswith("g__"):
                    gtdb_genus = part[3:].strip()

            rows.append({
                "motu": motu,
                "tax_str": tax_str,
                "gtdb_species": gtdb_species,
                "gtdb_genus": gtdb_genus,
                "abundance": abd,
            })
    return rows


# ============================== 代表基因组表 ==============================

def load_rep_table(rep_file):
    """读 mOTUs 代表基因组表，返回 mOTU -> rep_genome。
    文件: mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz
    列: MOTU | GENOME | GTDB  (可能还有表头)
    """
    rep_map = {}
    if not os.path.exists(rep_file):
        print(f"[warn] 代表基因组表不存在: {rep_file}", file=sys.stderr)
        return rep_map

    opener = gzip.open if rep_file.endswith(".gz") else open
    with opener(rep_file, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            motu = cols[0].strip()
            genome = cols[1].strip()
            # 跳过表头行
            if motu.lower() == "motu":
                continue
            rep_map[motu] = genome

    print(f"[info] 代表基因组表已加载: {len(rep_map)} 个 mOTU", file=sys.stderr)
    return rep_map


# ============================== 解析 mOTU 的 NCBI taxid ==============================

def resolve_taxid(motu, gtdb_species, rep_map, asm, taxdump, cache):
    """给定 mOTU，返回 (taxid, name, rep_genome, how)。
    taxid 为空字符串表示找不到。
    how 说明来源：'rep_accession' / 'rep_biosample' / ''(失败)

    结果缓存，避免重复查询。
    """
    if motu in cache:
        return cache[motu]

    rep_genome = rep_map.get(motu, "")
    taxid = ""
    name = ""
    how = ""

    if rep_genome and asm is not None and asm.loaded:
        kind, val = genome_id_to_lookup(rep_genome)
        acc = ""
        if kind == "accession":
            acc = val
        elif kind == "biosample":
            acc = asm.lookup_biosample(val)

        if acc:
            rec = asm.lookup_accession(acc)
            if rec:
                t_taxid, t_species, t_name = rec
                use = t_species or t_taxid
                if use and use != "na":
                    taxid = use
                    name = t_name or ""
                    how = "rep_" + kind

    result = (taxid, name, rep_genome, how)
    cache[motu] = result
    return result


# ============================== 组装 CAMI ==============================

def build_cami_from_relab(rows, rep_map, asm, taxdump, sample_id,
                          missing_mode="keep", version="0.9.1"):
    """从 .relab 行直接生成 CAMI 文本行。"""

    cache = {}
    mapped = []      # 成功映射
    unmapped = []    # 未映射（taxid 空）

    for r in rows:
        taxid, name, rep_genome, how = resolve_taxid(
            r["motu"], r["gtdb_species"], rep_map, asm, taxdump, cache)
        if taxid:
            mapped.append({
                "motu": r["motu"],
                "taxid": str(taxid),
                "name": name or r["gtdb_species"] or f"taxid_{taxid}",
                "abundance": r["abundance"],
                "rep_genome": rep_genome,
                "how": how,
            })
        else:
            unmapped.append({
                "motu": r["motu"],
                "gtdb_species": r["gtdb_species"],
                "gtdb_genus": r["gtdb_genus"],
                "rep_genome": rep_genome,
                "abundance": r["abundance"],
            })

    # ---- 同一 taxid 合并 ----
    merged = OrderedDict()
    for e in mapped:
        t = e["taxid"]
        if t in merged:
            merged[t]["abundance"] += e["abundance"]
        else:
            merged[t] = dict(e)

    # ---- 建节点 + 沿路径滚丰度（与 motu2cami.py 同一算法）----
    nodes = OrderedDict()

    def ensure_node(tid, nm, rank):
        if tid not in nodes:
            nodes[tid] = {"taxid": tid, "name": nm, "rank": rank,
                          "abundance": 0.0, "taxpath": [], "taxpathsn": []}
        return nodes[tid]

    leaf_total = 0.0
    for t, e in merged.items():
        tid_list, name_list = taxdump.lineage_cami(t) if taxdump.loaded else ([], [])
        if not tid_list:
            tid_list, name_list = [t], [e["name"]]
        for i, (tt, nn) in enumerate(zip(tid_list, name_list)):
            cr = None
            if taxdump.loaded:
                cr = NCBI_RANK_MAP.get(taxdump.ranks.get(tt, ""), None)
            if cr is None:
                cr = "species" if i == len(tid_list) - 1 else "no rank"
            nd = ensure_node(tt, nn, cr)
            if not nd["taxpath"]:
                nd["taxpath"] = tid_list[: i + 1]
                nd["taxpathsn"] = name_list[: i + 1]
        for tt in tid_list:
            nodes[tt]["abundance"] += e["abundance"]
        leaf_total += e["abundance"]

    # ---- 分母：所有行（含未映射）的丰度和，保证百分比语义一致 ----
    all_total = leaf_total + sum(x["abundance"] for x in unmapped)
    scale = 100.0 / all_total if all_total > 0 else 1.0

    out = []
    out.append(f"@SampleID:{sample_id}")
    out.append(f"@Version:{version}")
    out.append("@Ranks:" + "|".join(CAMI_RANKS))
    out.append("")
    out.append("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t"
               "_CAMI_genomeID\t_CAMI_OTU")

    rank_bucket = OrderedDict((r, []) for r in CAMI_RANKS)
    for t, nd in nodes.items():
        if nd["rank"] in CAMI_RANKS and nd["abundance"] > 0:
            rank_bucket[nd["rank"]].append(nd)

    n_out = 0
    for rank in CAMI_RANKS:
        for nd in sorted(rank_bucket[rank], key=lambda x: -x["abundance"]):
            out.append("\t".join([
                nd["taxid"], rank,
                "|".join(nd["taxpath"]) if nd["taxpath"] else nd["taxid"],
                "|".join(nd["taxpathsn"]) if nd["taxpathsn"] else nd["name"],
                f"{nd['abundance'] * scale:.4f}", "", "",
            ]))
            n_out += 1

    # ---- 未映射行：按用户要求「空着」----
    if missing_mode == "keep" and unmapped:
        for x in sorted(unmapped, key=lambda y: -y["abundance"]):
            out.append("\t".join([
                "",            # TAXID 空
                "species",
                "",            # TAXPATH 空
                x["gtdb_species"] or "unknown",
                f"{x['abundance'] * scale:.4f}",
                "", "",
            ]))
            n_out += 1

    stats = {
        "n_input": len(rows),
        "n_mapped": len(mapped),
        "n_unmapped": len(unmapped),
        "n_nodes_out": n_out,
        "mapped_abd_pct": leaf_total * scale,
        "unmapped_abd_pct": sum(x["abundance"] for x in unmapped) * scale,
        "unmapped_list": unmapped,
        "all_total_raw": all_total,
    }
    return out, stats


# ============================== 样本名提取 ==============================

# 命名的三大类：
#   motus-recall-result-2017.12.04_18.45.54_sample_5.txt.relab
#   motus-precision-result-simulated-3-result-reads.txt.relab
#   motus-recall-result-CC1_1.clean.txt.relab
def sample_id_from_filename(fname):
    """从 .relab 文件名提取样本名。
    处理三种实际命名:
        motus-recall-result-2017.12.04_18.45.54_sample_5.txt.relab  -> 2017.12.04_18.45.54_sample_5
        motus-precision-result-simulated-3-result-reads.txt.relab    -> simulated-3
        motus-recall-result-CC1_1.clean.txt.relab                    -> CC1_1
    """
    base = os.path.basename(fname)
    base = re.sub(r'\.relab$', '', base)          # 去 .relab
    base = re.sub(r'\.txt$', '', base)            # 去 .txt
    base = re.sub(r'\.clean$', '', base)          # 去 .clean
    base = re.sub(r'^motus-(recall|precision)-result-', '', base)  # 去 motus-xxx-result-
    base = re.sub(r'-result-reads$', '', base)    # 去 -result-reads
    return base


# ============================== 主流程 ==============================

def find_relabs(input_dir, recursive=True):
    if recursive:
        hits = glob.glob(os.path.join(input_dir, "**", "*.relab"), recursive=True)
    else:
        hits = glob.glob(os.path.join(input_dir, "*.relab"))
    return sorted(hits)


def main():
    ap = argparse.ArgumentParser(
        description="批量 .relab -> CAMI（OPAL 可直接用）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("--input-dir", required=True,
                    help="含 .relab 的根目录（会递归找）")
    ap.add_argument("--output-dir", required=True,
                    help="输出 CAMI 的目录")
    ap.add_argument("--rep-file", required=True,
                    help="mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz 路径")
    ap.add_argument("--ncbi-dir", default=None,
                    help="taxdump 目录（补全谱系）")
    ap.add_argument("--biosample-map", action="append", default=[],
                    help="assembly_summary_genbank.txt / _refseq.txt；可给多次")
    ap.add_argument("--missing-mode", choices=["keep", "drop"], default="keep",
                    help="无法映射的行：keep=保留但TAXID留空(默认) / drop=删掉")
    ap.add_argument("--no-recursive", action="store_true",
                    help="只处理 input-dir 顶层，不递归")
    ap.add_argument("--flat", action="store_true",
                    help="所有 CAMI 平铺到 output-dir，不保留组名层级")
    args = ap.parse_args()

    # ---- 加载参考数据 ----
    taxdump = Taxdump(args.ncbi_dir)
    rep_map = load_rep_table(args.rep_file)

    asm = None
    if args.biosample_map:
        asm = AssemblySummary()
        for p in args.biosample_map:
            asm.load(p)
        if not asm.loaded:
            asm = None

    relabs = find_relabs(args.input_dir, recursive=not args.no_recursive)
    if not relabs:
        print(f"[error] 在 {args.input_dir} 没找到 .relab 文件", file=sys.stderr)
        sys.exit(1)
    print(f"[info] 找到 {len(relabs)} 个 .relab 文件", file=sys.stderr)

    os.makedirs(args.output_dir, exist_ok=True)
    summary_rows = []
    all_unmapped = []

    for i, relab in enumerate(relabs, 1):
        sample_id = sample_id_from_filename(relab)
        rel_path = os.path.relpath(relab, args.input_dir)
        group_dir = os.path.dirname(rel_path).replace(os.sep, "__") or "."

        # 输出路径
        if args.flat:
            out_sub = args.output_dir
        else:
            out_sub = os.path.join(args.output_dir, os.path.dirname(rel_path))
            if not os.path.dirname(rel_path):
                out_sub = args.output_dir
        os.makedirs(out_sub, exist_ok=True)
        out_path = os.path.join(out_sub, sample_id + ".cami")

        rows = read_relab(relab)
        lines, stats = build_cami_from_relab(
            rows, rep_map, asm, taxdump, sample_id,
            missing_mode=args.missing_mode)

        with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(lines) + "\n")

        for x in stats["unmapped_list"]:
            y = dict(x)
            y["sample"] = sample_id
            y["group"] = group_dir
            all_unmapped.append(y)

        summary_rows.append({
            "group": group_dir,
            "sample": sample_id,
            "n_input": stats["n_input"],
            "n_mapped": stats["n_mapped"],
            "n_unmapped": stats["n_unmapped"],
            "n_nodes_out": stats["n_nodes_out"],
            "mapped_abd_pct": f"{stats['mapped_abd_pct']:.2f}",
            "unmapped_abd_pct": f"{stats['unmapped_abd_pct']:.2f}",
            "cami_out": os.path.relpath(out_path, args.output_dir).replace(os.sep, "/"),
        })

        if i % 10 == 0 or i == len(relabs):
            print(f"[{i}/{len(relabs)}] {sample_id}: "
                  f"映射 {stats['n_mapped']}, 未映射 {stats['n_unmapped']}",
                  file=sys.stderr)

    # ---- 汇总表 ----
    if summary_rows:
        sum_path = os.path.join(args.output_dir, "_batch_summary.tsv")
        cols = ["group", "sample", "n_input", "n_mapped", "n_unmapped",
                "n_nodes_out", "mapped_abd_pct", "unmapped_abd_pct", "cami_out"]
        with open(sum_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in summary_rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"[done] 汇总表: {sum_path}", file=sys.stderr)

    # ---- 未映射汇总 ----
    if all_unmapped:
        un_path = os.path.join(args.output_dir, "_unmapped_all.tsv")
        cols = ["group", "sample", "motu", "gtdb_species",
                "rep_genome", "abundance"]
        with open(un_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(cols) + "\n")
            for x in sorted(all_unmapped,
                            key=lambda y: (y["group"], y["sample"], -y["abundance"])):
                fh.write("\t".join(str(x.get(c, "")) for c in cols) + "\n")
        print(f"[done] 未映射汇总: {un_path}", file=sys.stderr)

    # ---- 控制台总结 ----
    tot_in = sum(r["n_input"] for r in summary_rows)
    tot_map = sum(r["n_mapped"] for r in summary_rows)
    tot_un = sum(r["n_unmapped"] for r in summary_rows)
    print("", file=sys.stderr)
    print(f"===== 完成: {len(summary_rows)} 个样本 =====", file=sys.stderr)
    print(f"  总 mOTU 行数: {tot_in}", file=sys.stderr)
    print(f"  成功映射:     {tot_map} ({tot_map/max(tot_in,1)*100:.1f}%)",
          file=sys.stderr)
    print(f"  未映射(空着): {tot_un} ({tot_un/max(tot_in,1)*100:.1f}%)",
          file=sys.stderr)
    print(f"  输出目录:     {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
