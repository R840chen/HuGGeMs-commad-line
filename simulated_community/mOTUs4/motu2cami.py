#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================================
 motu_ncbi_taxid.tsv  ->  CAMI profiling format (Bioboxes 0.9.1)
============================================================================
 用途:
   把 motus_ncbi_pipeline.py 产出的 motu_ncbi_taxid.tsv
   转成 OPAL 能直接比对的 CAMI 格式（cami2_gold_standard_profile-*.txt 同款）。

 为什么必须转换:
   OPAL 只认 CAMI Bioboxes 格式，要求：
     1. TAXID 是 NCBI 数字 taxid（不是 mOTU 编号，不是物种名字符串）
     2. RANK 是 superkingdom/phylum/class/order/family/genus/species
     3. TAXPATH / TAXPATHSN 是从根到当前节点的完整路径（| 分隔）
     4. PERCENTAGE 是百分比（0-100），同一 rank 内部累积

 核心难点（也是本脚本存在的理由）:
   mOTUs 用 GTDB 分类，很多种的代表基因组是 MAG / 非 NCBI，
   在 NCBI 里根本没有 taxid。这些行怎么办？ —— 用 BioSample 二次补救：

   第 1 步：代表基因组 ID 内嵌 accession
       RSGB23-1_GCF-024329905-V1_GENO_10000001  -> GCF_024329905.1  ✅
       RSGB23-1_GCA-000009285-V2_GENO_10000001  -> GCA_000009285.2  ✅
   第 2 步（本版新增）：MAG 的代表基因组只有 BioSample，用
       --biosample-map <assembly_summary_genbank.txt> 把 BioSample -> GCA/GCF
       （assembly_summary_genbank.txt 的第 3 列就是 biosample，格式：
         assembly_accession  bioproject  biosample  ...  taxid  species_taxid ...）

   --missing-mode drop   (推荐)  直接删掉，并把占比归一到 100%
   --missing-mode keep           保留行但 TAXID 留空（OPAL 会当 root 悬挂）

   → drop 是最稳的：OPAL 内部会把这些未匹配物种算作 "假阴性"，
     这是真实反映 mOTUs 无法回收到 NCBI 物种的现状，正是你要评估的东西。

 用法:
   python motu2cami.py \
       --input  motu_ncbi_taxid.tsv \
       --output sample_5.cami \
       --ncbi-dir /ME4012-Vol01/database/NCBI_tax \
       --biosample-map /ME4012-Vol01/database/NCBI_tax/assembly_summary_genbank.txt \
       --sample-id 2017.12.04_18.45.54_sample_5 \
       --missing-mode drop

 依赖: 无（纯标准库）。taxdump 有则用，没有也能跑（自动跳过谱系补全）。
============================================================================
"""
import os
import re
import sys
import gzip
import argparse
from collections import OrderedDict, defaultdict

# ============================== 常量 ==============================

# CAMI 标准 rank 顺序（必须与 @Ranks 头一致）
CAMI_RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

# NCBI rank 名 -> CAMI rank 名（CAMI 只认 7 级；亚种/strain 归到 species）
NCBI_RANK_MAP = {
    "superkingdom": "superkingdom",
    "domain": "superkingdom",          # 有些 taxdump 用 domain
    "kingdom": None,                   # CAMI 无 kingdom，跳过
    "phylum": "phylum",
    "class": "class",
    "order": "order",
    "family": "family",
    "genus": "genus",
    "species": "species",
    "subspecies": "species",
    "strain": "species",
    "no rank": None,                   # 中间无 rank 节点，跳过
    "clade": None,
}

# 起始 taxid（Bacteria=2, Archaea=2157, Eukaryota=2759 的统一父节点是 1 root）
ROOT_TAXID = "1"
BACTERIA_TAXID = "2"


# ============================== taxdump 读取 ==============================

class Taxdump:
    """读取 NCBI taxdump，提供 taxid -> (name, rank, parent) 与完整谱系。"""

    def __init__(self, ncbi_dir=None):
        self.nodes = {}     # taxid -> parent_taxid
        self.ranks = {}     # taxid -> rank
        self.names = {}     # taxid -> scientific name
        self.merged = {}    # old taxid -> new taxid
        self.deleted = set()
        self.loaded = False
        if ncbi_dir:
            self.load(ncbi_dir)

    # -------- 加载 --------
    def _open(self, path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return open(path, "r", encoding="utf-8", errors="replace")

    def load(self, ncbi_dir):
        def find(name):
            p = os.path.join(ncbi_dir, name)
            if os.path.exists(p):
                return p
            gz = p + ".gz"
            if os.path.exists(gz):
                return gz
            return None

        nodes_p, names_p = find("nodes.dmp"), find("names.dmp")
        if not nodes_p or not names_p:
            print(f"[warn] 未在 {ncbi_dir} 找到 nodes.dmp / names.dmp，跳过谱系补全",
                  file=sys.stderr)
            return

        # 注意: NCBI dmp 文件每行的结构是
        #     "值1\t|\t值2\t|\t...\t|\r\n"
        # 即字段间用 "\t|\t" 分隔, 行尾还有一个多余的 "\t|"。
        # Python 文本模式会把 \r\n 统一成 \n, 所以只需:
        #   1) rstrip("\n")
        #   2) 去掉行尾的 "\t|"
        #   3) 用 "\t|\t" 切分
        def fields(line):
            s = line.rstrip("\n")
            if s.endswith("\t|"):
                s = s[:-2]
            return [x.strip() for x in s.split("\t|\t")]

        # nodes.dmp: taxid | parent | rank | ...
        with self._open(nodes_p) as fh:
            for line in fh:
                parts = fields(line)
                if len(parts) < 3:
                    continue
                t = parts[0]
                if t:
                    self.nodes[t] = parts[1]
                    self.ranks[t] = parts[2]

        # names.dmp: taxid | name | unique name | name class |
        with self._open(names_p) as fh:
            for line in fh:
                parts = fields(line)
                if len(parts) < 4:
                    continue
                if parts[3] == "scientific name":
                    t = parts[0]
                    if t:
                        self.names[t] = parts[1]

        # merged.dmp: old | new |
        mp = find("merged.dmp")
        if mp:
            with self._open(mp) as fh:
                for line in fh:
                    parts = fields(line)
                    if len(parts) >= 2 and parts[0] and parts[1]:
                        self.merged[parts[0]] = parts[1]

        # delnodes.dmp
        dp = find("delnodes.dmp")
        if dp:
            with self._open(dp) as fh:
                for line in fh:
                    parts = fields(line)
                    if parts and parts[0]:
                        self.deleted.add(parts[0])

        self.loaded = True
        print(f"[info] taxdump 已加载: {len(self.nodes)} nodes, "
              f"{len(self.names)} names, {len(self.merged)} merged",
              file=sys.stderr)

    # -------- 查询 --------
    def resolve(self, taxid):
        """把 merged/已删除的 taxid 归一到现存 taxid。"""
        taxid = str(taxid).strip()
        seen = set()
        while taxid in self.merged and taxid not in seen:
            seen.add(taxid)
            taxid = self.merged[taxid]
        if taxid in self.deleted:
            return None
        return taxid

    def lineage(self, taxid, max_depth=50):
        """返回从根(含)到 taxid(含) 的 [(taxid, name, rank), ...]。
        若 taxid 不在 rank 表里（缺 taxdump），返回单节点。"""
        taxid = self.resolve(taxid)
        if taxid is None:
            return []
        if not self.loaded or taxid not in self.nodes:
            return [(taxid, self.names.get(taxid, taxid), "species")]

        chain = []
        cur = taxid
        depth = 0
        while cur and depth < max_depth:
            chain.append((cur,
                          self.names.get(cur, cur),
                          self.ranks.get(cur, "no rank")))
            parent = self.nodes.get(cur)
            if parent is None or parent == cur:
                break
            cur = parent
            depth += 1
        chain.reverse()  # 根在前
        return chain

    def lineage_cami(self, taxid):
        """只保留 CAMI 7 级 rank，返回 (taxids[], names[]) 从 root 到该节点。"""
        chain = self.lineage(taxid)
        if not chain:
            return [], []
        tid_list, name_list = [], []
        rank_seen = {}  # rank -> index，用于同名 rank 去重（保留最深）
        for t, n, r in chain:
            cr = NCBI_RANK_MAP.get(r, None)
            if cr is None:
                continue
            rank_seen[cr] = len(tid_list)
            tid_list.append(t)
            name_list.append(n)
        return tid_list, name_list


# ============================== assembly_summary 读取 ==============================

class AssemblySummary:
    """读 assembly_summary_genbank.txt（或 _refseq.txt），提供：
         accession -> (taxid, species_taxid, organism_name)
         biosample -> accession        ← 本类新增，用于救 MAG
       列结构（真实表头，注意首行是 '##' 说明行，第 2 行才是 '#列名'）：
         assembly_accession  bioproject  biosample  wgs_master  refseq_category
         taxid  species_taxid  organism_name  ...
    """

    def __init__(self, path=None):
        self.acc2taxid = {}      # accession(不含版本) -> (taxid, species_taxid, name)
        self.acc2taxid_ver = {}  # accession(含版本)   -> (taxid, species_taxid, name)
        self.biosample2acc = {}  # biosample -> accession
        self.biosample_ambig = {}  # biosample -> 多个 accession 时记数，供报告用
        self.loaded = False
        if path:
            self.load(path)

    def _open(self, path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return open(path, "r", encoding="utf-8", errors="replace")

    def load(self, path):
        if not os.path.exists(path):
            gz = path + ".gz"
            if os.path.exists(gz):
                path = gz
            else:
                print(f"[warn] 未找到 {path}，跳过 accession/biosample 映射",
                      file=sys.stderr)
                return

        with self._open(path) as fh:
            header = None
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                # 首行 '##' 说明行；第二行 '#assembly_accession ...' 才是表头
                if line.startswith("#"):
                    if not line.startswith("##") and header is None:
                        header = line.lstrip("#").split("\t")
                    continue
                if header is None:
                    continue
                cols = line.split("\t")
                if len(cols) < len(header):
                    continue

                def col(name):
                    i = header.index(name) if name in header else None
                    return cols[i].strip() if i is not None and i < len(cols) else ""

                acc = col("assembly_accession")
                taxid = col("taxid")
                sp_taxid = col("species_taxid")
                org = col("organism_name")
                bs = col("biosample")

                if acc:
                    rec = (taxid, sp_taxid, org)
                    self.acc2taxid_ver[acc] = rec
                    base = acc.rsplit(".", 1)[0] if "." in acc else acc
                    self.acc2taxid[base] = rec

                if bs and bs != "na" and acc:
                    if bs in self.biosample2acc:
                        # 同一 BioSample 对多个 assembly：优先 GCF（RefSeq），其余记数
                        self.biosample_ambig[bs] = self.biosample_ambig.get(bs, 1) + 1
                        old = self.biosample2acc[bs]
                        if acc.startswith("GCF") and not old.startswith("GCF"):
                            self.biosample2acc[bs] = acc
                    else:
                        self.biosample2acc[bs] = acc

        self.loaded = True
        print(f"[info] assembly_summary 已加载: {len(self.acc2taxid_ver)} accession, "
              f"{len(self.biosample2acc)} biosample "
              f"(其中 {len(self.biosample_ambig)} 个对应多个 assembly)",
              file=sys.stderr)

    def lookup_accession(self, acc):
        """accession -> (taxid, species_taxid, organism_name)；找不到返回 None"""
        if not self.loaded or not acc:
            return None
        if acc in self.acc2taxid_ver:
            return self.acc2taxid_ver[acc]
        base = acc.rsplit(".", 1)[0] if "." in acc else acc
        return self.acc2taxid.get(base)

    def lookup_biosample(self, bs):
        """biosample -> accession；找不到返回 ''"""
        if not self.loaded or not bs:
            return ""
        return self.biosample2acc.get(bs, "")


# ============================== MAG / 基因组 ID 解析 ==============================

# mOTUs 基因组 ID 命名规则（官方 naming.html）
RE_REF_GENOME = re.compile(r'^(.+?)_(GC[AF])-(\d+)-V(\d+)_GENO_\d+$')
RE_MAG_GENOME = re.compile(r'^(.+?)_(SAM[NDE][A-Z]?\d+(?:-S\d+)?)_MAG_\d+$')

# 纯 BioSample 编号（可去掉 -S### 后缀）
RE_BIOSAMPLE = re.compile(r'^(SAM[NDE][A-Z]?\d+)')


def genome_id_to_lookup(gid):
    """从 mOTUs 基因组 ID 里提取可用于查 NCBI 的键。
    返回 (kind, value)：
       ('accession', 'GCF_024329905.1')  -> RefSeq/GenBank 分离株
       ('biosample', 'SAMN15871286')     -> MAG
       (None, None)                       -> 非 NCBI（JGI 等），无法查
    """
    if not gid:
        return (None, None)

    # 1) 参考基因组：内嵌 GCF/GCA
    m = RE_REF_GENOME.match(gid)
    if m:
        prefix, digits, ver = m.group(2), m.group(3), m.group(4)
        # 还原成标准 accession：GCF_024329905.1（digits 补零 9 位）
        num = digits.zfill(9)
        return ("accession", f"{prefix}_{num}.{ver}")

    # 2) MAG：内嵌 BioSample（可能带 -S002 后缀）
    m = RE_MAG_GENOME.match(gid)
    if m:
        bs_raw = m.group(2)
        # 去掉 -S### 后缀（那是 SRA 样本序号，不是 BioSample 本体）
        bs = bs_raw.split("-S")[0]
        # 校验是不是合法 BioSample
        mb = RE_BIOSAMPLE.match(bs)
        if mb:
            return ("biosample", mb.group(1))
        return (None, None)

    # 3) 其它（JGI GA####、SRA run 名等）无法查
    return (None, None)


# ============================== 读输入 TSV ==============================

def read_input(path, missing_mode):
    """读 motu_ncbi_taxid.tsv，返回 [row_dict, ...]（已过滤 unassigned）。"""
    rows = []
    with open(path, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}

        def get(cols, name, default=""):
            i = idx.get(name)
            if i is None or i >= len(cols):
                return default
            return cols[i].strip()

        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")

            motu = get(cols, "mOTU")
            is_un = get(cols, "is_unassigned", "no")

            # 跳过 unassigned 汇总行
            if "unassigned" in motu.lower() or is_un.lower() == "yes":
                continue

            taxid = get(cols, "NCBI_taxid")
            rank = get(cols, "NCBI_rank") or "species"
            name = get(cols, "NCBI_name")
            gtdb_s = get(cols, "GTDB_species")
            gtdb_g = get(cols, "GTDB_genus")
            rep = get(cols, "rep_genome")
            rep_type = get(cols, "rep_type")
            try:
                abd = float(get(cols, "Abundance") or 0.0)
            except ValueError:
                abd = 0.0

            rows.append({
                "motu": motu,
                "taxid": taxid,
                "rank": rank,
                "name": name,
                "gtdb_species": gtdb_s,
                "gtdb_genus": gtdb_g,
                "rep_genome": rep,
                "rep_type": rep_type,
                "abundance": abd,
            })
    return rows


# ============================== 分配 CAMI rank ==============================

def cami_rank_of(taxid, gtdb_rank_hint, taxdump):
    """决定这一行在 CAMI 里算什么 rank。"""
    if taxid and taxdump.loaded:
        r = taxdump.ranks.get(taxdump.resolve(taxid) or "", None)
        cr = NCBI_RANK_MAP.get(r, None)
        if cr:
            return cr
    return "species"   # 默认物种级


# ============================== 主体转换 ==============================

def build_cami(rows, taxdump, missing_mode, sample_id, version="0.9.1",
               asm=None):
    """返回 CAMI 文件的文本行列表。

    asm: AssemblySummary 实例（可选）。给了就能对 taxid 为空的行做二次补救：
         用 rep_genome 里的 BioSample / accession 反查 NCBI taxid。
    """

    # ---- 1. 把每一行的 abundance 归到某个 (taxid, rank) ----
    # taxid 可能为空 -> 先尝试用 assembly_summary 二次补救 -> 再走 missing 策略
    entries = []       # (taxid, rank, abundance, name)
    dropped = []       # 被丢弃的行（用于报告）
    unresolved = []    # taxid 空但保留的行
    rescued = []       # 被 BioSample/accession 二次补救成功的行

    for r in rows:
        taxid = r["taxid"]
        name = r["name"]
        rank_hint = r["rank"]

        # ---- 二次补救：taxid 为空时，从代表基因组反查 ----
        if not taxid and asm is not None and asm.loaded and r["rep_genome"]:
            kind, val = genome_id_to_lookup(r["rep_genome"])
            acc = ""
            if kind == "accession":
                acc = val
            elif kind == "biosample":
                acc = asm.lookup_biosample(val)

            if acc:
                rec = asm.lookup_accession(acc)
                if rec:
                    t_taxid, t_species, t_name = rec
                    # 优先 species_taxid（物种级评估更稳），退化到 taxid
                    use = t_species or t_taxid
                    if use and use != "na":
                        taxid = use
                        name = t_name or name
                        rank_hint = "species"
                        rescued.append({
                            "motu": r["motu"],
                            "rep_type": r["rep_type"],
                            "rep_genome": r["rep_genome"],
                            "key_kind": kind,
                            "key_val": val,
                            "assembly": acc,
                            "taxid": use,
                            "organism": t_name,
                        })

        if not taxid:
            if missing_mode == "drop":
                dropped.append(r)
                continue
            elif missing_mode == "keep":
                unresolved.append(r)
                continue

        taxid = str(taxid)
        rank = cami_rank_of(taxid, rank_hint, taxdump)
        if rank not in CAMI_RANKS:
            rank = "species"
        entries.append({
            "taxid": taxid,
            "rank": rank,
            "abundance": r["abundance"],
            "name": name or r["gtdb_species"] or f"taxid_{taxid}",
        })

    # ---- 2. 同一 taxid 出现多次 -> 合并丰度 ----
    merged = OrderedDict()
    for e in entries:
        t = e["taxid"]
        if t in merged:
            merged[t]["abundance"] += e["abundance"]
        else:
            merged[t] = dict(e)

    # ---- 3. 为每个叶子行解析完整谱系，并把丰度沿路径"向上累加一次" ----
    #
    # 关键设计（避免 double counting）:
    #   把每一行看成一个"叶子"，它的完整路径是 path = [root, ..., leaf]。
    #   叶子的原始丰度 leaf_abd 沿这条 path 的**每个**节点累加一次：
    #       nodes[root] += leaf_abd, nodes[path[1]] += leaf_abd, ..., nodes[leaf] += leaf_abd
    #   这样每个节点天然得到"自身+所有后代"的累积值，且不会有重复。
    #
    #   如果一行缺谱系（taxdump 不全），它的 path 只有 [self]，
    #   那么只累加到自己，不会被错误地加到 Bacteria 上。
    nodes = OrderedDict()

    def ensure_node(tid, name, rank):
        if tid not in nodes:
            nodes[tid] = {
                "taxid": tid, "name": name, "rank": rank,
                "abundance": 0.0,
                "taxpath": [], "taxpathsn": [],
            }
        return nodes[tid]

    leaf_abd_total = 0.0   # 所有成功映射行的原始丰度和（分母）

    for t, e in merged.items():
        tid_list, name_list = taxdump.lineage_cami(t) if taxdump.loaded else ([], [])
        if not tid_list:
            tid_list, name_list = [t], [e["name"]]

        # 记录每个节点的 path（同 taxid 只写一次，取第一次的路径）
        for i, (tt, nn) in enumerate(zip(tid_list, name_list)):
            cr = NCBI_RANK_MAP.get(taxdump.ranks.get(tt, ""), None) if taxdump.loaded else None
            if cr is None:
                cr = e["rank"] if i == len(tid_list) - 1 else "no rank"
            nd = ensure_node(tt, nn, cr)
            if not nd["taxpath"]:
                nd["taxpath"] = tid_list[: i + 1]
                nd["taxpathsn"] = name_list[: i + 1]

        # 这一行的丰度，沿路径每一级都累加一次
        abd = e["abundance"]
        for tt in tid_list:
            nodes[tt]["abundance"] += abd
        leaf_abd_total += abd

    # ---- 4. 归一化 ----
    # 分母 = 所有 leaf 原始丰度之和（= 成功映射行的 Abundance 总和）
    # 这样 species 级各行的和恰好 = 100，父节点 = 子节点累积。
    scale = 100.0 / leaf_abd_total if leaf_abd_total > 0 else 1.0

    # ---- 5. 写出 CAMI ----
    # 头格式严格对齐官方 gold standard（bioboxes profiling 0.9.1）:
    #   @SampleID / @Version / @Ranks / 空行 / @@列头
    out = []
    out.append(f"@SampleID:{sample_id}")
    out.append(f"@Version:{version}")
    out.append("@Ranks:" + "|".join(CAMI_RANKS))
    out.append("")
    out.append("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t"
               "_CAMI_genomeID\t_CAMI_OTU")

    # 按 rank 分组，rank 内按丰度降序（CAMI 惯例，便于阅读）
    rank_bucket = OrderedDict((r, []) for r in CAMI_RANKS)
    for t, nd in nodes.items():
        rank = nd["rank"]
        if rank not in CAMI_RANKS:
            continue
        if nd["abundance"] <= 0:
            continue
        rank_bucket[rank].append(nd)

    total_written = 0
    for rank in CAMI_RANKS:
        bucket = sorted(rank_bucket[rank], key=lambda x: -x["abundance"])
        for nd in bucket:
            pct = nd["abundance"] * scale
            out.append("\t".join([
                nd["taxid"],
                rank,
                "|".join(nd["taxpath"]) if nd["taxpath"] else nd["taxid"],
                "|".join(nd["taxpathsn"]) if nd["taxpathsn"] else nd["name"],
                f"{pct:.4f}",
                "",   # _CAMI_genomeID 留空
                "",   # _CAMI_OTU 留空
            ]))
            total_written += 1

    # ---- 6. 保留空 taxid 行（--missing-mode keep）----
    if unresolved:
        for r in unresolved:
            out.append("\t".join([
                "",                              # TAXID 空
                "species",
                "",                              # TAXPATH 空
                r["gtdb_species"] or "unknown",  # 名字放这里便于人工核对
                f"{r['abundance'] * scale:.4f}",
                "",
                "",
            ]))
            total_written += 1

    return out, {
        "n_rows_input": len(rows),
        "n_nodes_written": total_written,
        "n_dropped": len(dropped),
        "n_unresolved_kept": len(unresolved),
        "n_rescued": len(rescued),
        "dropped_abundance_pct": sum(x["abundance"] for x in dropped) * scale,
        "leaf_total_raw": leaf_abd_total,
        "dropped_list": dropped,
        "unresolved_list": unresolved,
        "rescued_list": rescued,
    }


# ============================== CLI ==============================

def main():
    ap = argparse.ArgumentParser(
        description="motu_ncbi_taxid.tsv -> CAMI profiling format (OPAL-ready)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("--input", required=True, help="motu_ncbi_taxid.tsv")
    ap.add_argument("--output", required=True, help="输出的 CAMI 文件")
    ap.add_argument("--ncbi-dir", default=None,
                    help="taxdump 目录（含 nodes.dmp/names.dmp/merged.dmp），"
                         "用于补全完整谱系；缺省则只用单节点路径")
    ap.add_argument("--biosample-map", default=None,
                    help="assembly_summary_genbank.txt（或 _refseq.txt）路径。"
                         "给了就能用 BioSample / accession 反查 taxid，"
                         "把 MAG 二次补救回来。强烈建议提供。")
    ap.add_argument("--sample-id", default="sample",
                    help="@SampleID 头，建议用样本名")
    ap.add_argument("--missing-mode", choices=["drop", "keep"], default="drop",
                    help="无 NCBI taxid 的行如何处理（默认 drop）")
    ap.add_argument("--report", default=None,
                    help="额外写出一个统计报告（可选）")
    args = ap.parse_args()

    taxdump = Taxdump(args.ncbi_dir)
    asm = AssemblySummary(args.biosample_map) if args.biosample_map else None
    rows = read_input(args.input, args.missing_mode)
    print(f"[info] 读入 {len(rows)} 行（已剔除 unassigned）", file=sys.stderr)

    lines, stats = build_cami(rows, taxdump, args.missing_mode, args.sample_id,
                              asm=asm)

    # newline="\n" 强制 Unix 换行（LF），OPAL 在 Linux 上跑，绝不能用 CRLF
    with open(args.output, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"[done] 写出 {args.output}", file=sys.stderr)
    print(f"       行数: 输入 {stats['n_rows_input']} -> 输出 {stats['n_nodes_written']}",
          file=sys.stderr)
    if stats["n_rescued"]:
        print(f"       BioSample 二次补救成功: {stats['n_rescued']} 行",
              file=sys.stderr)
    print(f"       丢弃(仍无 taxid): {stats['n_dropped']} 行, "
          f"占总丰度 {stats['dropped_abundance_pct']:.2f}%", file=sys.stderr)

    if args.report:
        with open(args.report, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# motu2cami 转换报告\n\n")
            fh.write(f"输入: {args.input}\n")
            fh.write(f"输出: {args.output}\n")
            fh.write(f"missing-mode: {args.missing_mode}\n")
            fh.write(f"assembly_summary: {args.biosample_map or '(未提供)'}\n\n")
            fh.write(f"输入行数: {stats['n_rows_input']}\n")
            fh.write(f"输出行数: {stats['n_nodes_written']}\n")
            fh.write(f"BioSample/accession 二次补救: {stats['n_rescued']} 行\n")
            fh.write(f"仍无 taxid 被丢弃: {stats['n_dropped']} 行 "
                     f"({stats['dropped_abundance_pct']:.2f}% 丰度)\n")
            fh.write(f"无 taxid 保留: {stats['n_unresolved_kept']} 行\n\n")

            if stats.get("rescued_list"):
                fh.write("## 通过 BioSample / accession 补救成功的 mOTU\n\n")
                fh.write("mOTU\trep_type\t查询键\t键值\t命中assembly\t"
                         "NCBI_taxid\torganism\n")
                for x in sorted(stats["rescued_list"], key=lambda y: y["motu"]):
                    fh.write(f"{x['motu']}\t{x['rep_type']}\t{x['key_kind']}\t"
                             f"{x['key_val']}\t{x['assembly']}\t{x['taxid']}\t"
                             f"{x['organism']}\n")
                fh.write("\n")

            if stats["dropped_list"]:
                fh.write("## 仍无法映射到 NCBI 的 mOTU（已丢弃）\n\n")
                fh.write("mOTU\tGTDB_species\trep_type\trep_genome\tAbundance\n")
                for r in sorted(stats["dropped_list"],
                                key=lambda x: -x["abundance"]):
                    fh.write(f"{r['motu']}\t{r['gtdb_species']}\t"
                             f"{r['rep_type']}\t{r['rep_genome']}\t"
                             f"{r['abundance']:.8f}\n")
        print(f"[done] 报告写出 {args.report}", file=sys.stderr)


if __name__ == "__main__":
    main()
