# -*- coding: utf-8 -*-
"""
Step 01  生成 MaAsLin3 的输入文件（按协变量层 L1–L5 + 地域层 LG）
================================================================
为什么需要这一步：
  MaAsLin3 内部用 lm/glm 拟合，R 默认 na.action = na.omit，
  → **只要样本在公式的任一变量上缺失，该样本会被静默丢弃**。
  → 不报错，但样本量会悄悄缩水，而且 L1 和 L5 的样本集不一致，无法直接比较。
  所以这里**显式做 complete-case 筛选**，并为每一层生成一对已对齐的
  data/meta 文件，让 R 端只管跑、不用操心缺失值。

输出：GMB/maaslin3_v2/maaslin_inputs/
        <PID>__L1__data.tsv / __meta.tsv
        <PID>__L2__data.tsv / __meta.tsv      ( + age )
        <PID>__L3__data.tsv / __meta.tsv      ( + sex )
        <PID>__L4__data.tsv / __meta.tsv      ( + BMI )
        <PID>__L5__data.tsv / __meta.tsv      ( + age + sex + BMI )
        <PID>__LG__data.tsv / __meta.tsv      ( + country，仅多国家队列 )
        _design.tsv                           所有任务的设计表

阈值都集中在下面的“参数区”，想改标准只改这里。
"""
import csv, os, collections

# ============================================================================
# 参数区（标准都在这，改这里即可）
# ============================================================================
MIN_GROUP_N = 10        # 每个疾病/对照组最少样本数，低于此的组被剔除
MIN_COV_COVERAGE = 0.70  # 协变量非缺失比例下限（0.7 = 70%）
MIN_LEVEL_N = 5         # 分类协变量每个水平最少样本数
MIN_NUM_UNIQUE = 3      # 数值协变量最少不同取值个数
REFERENCE_GROUP = "Health"   # 对照组名称
GEO_MIN_N = 10          # 地域层：每个国家最少样本数

# ============================================================================
# 路径
# ============================================================================
BASE = r"D:/Desktop/刘双江-实验室/硕士发表文章/各种材料/4.真实宏基因组样本信息/真实宏基因组样本分析/4.疾病健康队列分析/GMB"
META_UNIFIED = os.path.join(BASE, "maaslin3_v2", "metadata_unified")
TAX_DIR = os.path.join(BASE, "相对丰度matrix")
OUT_DIR = os.path.join(BASE, "maaslin3_v2", "maaslin_inputs")
SUF = "_species_abundance_matrix.prev5pct.maaslin_ready.tsv"

# 层的定义：层名 -> (公式, 需要的协变量列表)
LEVELS = [
    ("L1", "~ group",                              []),
    ("L2", "~ group + age",                        ["age"]),
    ("L3", "~ group + sex",                        ["sex"]),
    ("L4", "~ group + BMI",                        ["BMI"]),
    ("L5", "~ group + age + sex + BMI",            ["age", "sex", "BMI"]),
]

os.makedirs(OUT_DIR, exist_ok=True)


def read_table(path):
    with open(path, encoding="utf-8", errors="ignore", newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def write_table(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


design = []
log = []


def record(pid, level, formula, status, n=0, n_ref=0, n_oth=0,
           cov_used="", note=""):
    design.append({
        "project_id": pid, "level": level, "formula": formula, "status": status,
        "n_samples": n, "n_Health": n_ref, "n_other": n_oth,
        "covariates_used": cov_used, "note": note,
        "file_prefix": "%s__%s" % (pid, level) if status == "OK" else "",
    })


# 必须先有统一 metadata
if not os.path.isdir(META_UNIFIED):
    raise SystemExit("找不到 %s，请先运行 00_unify_metadata.py" % META_UNIFIED)

meta_files = sorted(f for f in os.listdir(META_UNIFIED)
                    if f.endswith(".tsv") and not f.startswith("_"))
print("待处理队列: %d 个" % len(meta_files))
print("=" * 110)

for mf in meta_files:
    pid = mf[:-4]
    mx_path = os.path.join(TAX_DIR, pid + SUF)
    if not os.path.exists(mx_path):
        log.append("[跳过] %s 找不到丰度矩阵" % pid)
        continue

    # ---- 读统一 metadata ----
    mrows = read_table(os.path.join(META_UNIFIED, mf))
    mhdr = mrows[0]
    mi = {c: i for i, c in enumerate(mhdr)}
    meta = {}
    for r in mrows[1:]:
        if not r or not r[0].strip():
            continue
        meta[r[0].strip()] = {
            "group": r[mi["group"]].strip(),
            "age": r[mi["age"]].strip(),
            "sex": r[mi["sex"]].strip(),
            "BMI": r[mi["BMI"]].strip(),
            "country": r[mi["country"]].strip(),
        }

    # ---- 读丰度矩阵 ----
    trows = read_table(mx_path)
    thdr = trows[0]
    tax_samples = [r[0].strip() for r in trows[1:] if r and r[0].strip()]

    # ---- 取交集 ----
    common = [s for s in tax_samples if s in meta]
    n_only_mx = len(set(tax_samples) - set(meta))
    n_only_mt = len(set(meta) - set(tax_samples))

    # ---- 按 MIN_GROUP_N 过滤小组 ----
    gcnt = collections.Counter(meta[s]["group"] for s in common if meta[s]["group"])
    keep_groups = {g for g, c in gcnt.items() if c >= MIN_GROUP_N}
    dropped_groups = {g: c for g, c in gcnt.items() if c < MIN_GROUP_N}

    base_samples = [s for s in common if meta[s]["group"] in keep_groups]
    gcnt2 = collections.Counter(meta[s]["group"] for s in base_samples)

    if REFERENCE_GROUP not in gcnt2:
        log.append("[跳过] %s 无对照(%s) 或样本不足 | 组分布=%s"
                   % (pid, REFERENCE_GROUP, dict(gcnt2)))
        record(pid, "-", "-", "SKIP_NO_CONTROL", note="组分布=%s" % dict(gcnt2))
        continue
    if len([g for g in gcnt2 if g != REFERENCE_GROUP]) == 0:
        log.append("[跳过] %s 只有对照组" % pid)
        record(pid, "-", "-", "SKIP_NO_DISEASE", note="组分布=%s" % dict(gcnt2))
        continue

    # ---- 判断各协变量可用性 ----
    n_base = len(base_samples)

    def cov_ok(col, kind):
        vals = [meta[s][col] for s in base_samples]
        nonmiss = [v for v in vals if v != ""]
        if len(nonmiss) / max(n_base, 1) < MIN_COV_COVERAGE:
            return False, "覆盖率 %.0f%% < %.0f%%" % (
                100 * len(nonmiss) / max(n_base, 1), 100 * MIN_COV_COVERAGE)
        if kind == "cat":
            lv = collections.Counter(nonmiss)
            if len(lv) < 2:
                return False, "只有1个水平"
            if min(lv.values()) < MIN_LEVEL_N:
                return False, "最小水平 n=%d < %d" % (min(lv.values()), MIN_LEVEL_N)
        else:
            if len(set(nonmiss)) < MIN_NUM_UNIQUE:
                return False, "不同取值 < %d" % MIN_NUM_UNIQUE
        return True, ""

    cov_status = {}
    for col, kind in (("age", "num"), ("sex", "cat"), ("BMI", "num")):
        ok, why = cov_ok(col, kind)
        cov_status[col] = (ok, why)

    # ---- 逐层生成 ----
    made = []
    for lvl, formula, need in LEVELS:
        missing = [c for c in need if not cov_status[c][0]]
        if missing:
            why = "; ".join("%s(%s)" % (c, cov_status[c][1]) for c in missing)
            record(pid, lvl, formula, "SKIP_COV", note=why)
            continue

        # complete-case 筛选
        sel = [s for s in base_samples
               if all(meta[s][c] != "" for c in need)]
        g = collections.Counter(meta[s]["group"] for s in sel)
        if g.get(REFERENCE_GROUP, 0) == 0 or len([k for k in g if k != REFERENCE_GROUP]) == 0:
            record(pid, lvl, formula, "SKIP_EMPTY_GROUP",
                   note="筛选后组分布=%s" % dict(g))
            continue
        # 筛选后小组再次检查
        bad = {k: v for k, v in g.items() if v < MIN_LEVEL_N}
        if bad:
            record(pid, lvl, formula, "SKIP_TINY_GROUP", note="组过小=%s" % bad)
            continue

        # 写文件
        pre = "%s__%s" % (pid, lvl)
        # data：保留 ID 列 + 全部物种列，只留 sel 样本
        sel_set = set(sel)
        d_rows = [r for r in trows[1:] if r and r[0].strip() in sel_set]
        # 统一按 sel 顺序
        d_map = {r[0].strip(): r for r in d_rows}
        d_rows = [d_map[s] for s in sel if s in d_map]
        write_table(os.path.join(OUT_DIR, pre + "__data.tsv"), thdr, d_rows)
        # meta
        m_out_hdr = ["ID", "group", "age", "sex", "BMI", "country"]
        m_out = [[s, meta[s]["group"], meta[s]["age"], meta[s]["sex"],
                  meta[s]["BMI"], meta[s]["country"]] for s in sel]
        write_table(os.path.join(OUT_DIR, pre + "__meta.tsv"), m_out_hdr, m_out)

        record(pid, lvl, formula, "OK", n=len(sel),
               n_ref=g.get(REFERENCE_GROUP, 0),
               n_oth=len(sel) - g.get(REFERENCE_GROUP, 0),
               cov_used=",".join(need) if need else "-")
        made.append((lvl, len(sel)))

    # ---- 地域层 LG：仅多国家队列 ----
    ctry = collections.Counter(meta[s]["country"] for s in base_samples)
    ctry_ok = {k: v for k, v in ctry.items() if v >= GEO_MIN_N}
    if len(ctry_ok) >= 2:
        sel = [s for s in base_samples if meta[s]["country"] in ctry_ok]
        g = collections.Counter(meta[s]["group"] for s in sel)
        if g.get(REFERENCE_GROUP, 0) >= MIN_LEVEL_N and \
           len([k for k in g if k != REFERENCE_GROUP]) > 0:
            pre = "%s__LG" % pid
            sel_set = set(sel)
            d_map = {r[0].strip(): r for r in trows[1:] if r and r[0].strip() in sel_set}
            d_rows = [d_map[s] for s in sel if s in d_map]
            write_table(os.path.join(OUT_DIR, pre + "__data.tsv"), thdr, d_rows)
            m_out_hdr = ["ID", "group", "age", "sex", "BMI", "country"]
            m_out = [[s, meta[s]["group"], meta[s]["age"], meta[s]["sex"],
                      meta[s]["BMI"], meta[s]["country"]] for s in sel]
            write_table(os.path.join(OUT_DIR, pre + "__meta.tsv"), m_out_hdr, m_out)
            record(pid, "LG", "~ group + country", "OK", n=len(sel),
                   n_ref=g.get(REFERENCE_GROUP, 0),
                   n_oth=len(sel) - g.get(REFERENCE_GROUP, 0),
                   cov_used="country", note="国家=%s" % dict(ctry_ok))
    else:
        record(pid, "LG", "~ group + country", "SKIP_COV",
               note="多国家队列不足（%s）" % dict(ctry))

    log.append("[完成] %-12s 交集=%-5d 剔除小组=%-18s 生成层=%s"
               % (pid, len(common), str(dropped_groups)[:18],
                  ",".join("%s(%d)" % (l, n) for l, n in made)))

# ============================================================================
# 输出设计表与日志
# ============================================================================
design_path = os.path.join(OUT_DIR, "_design.tsv")
cols = ["project_id", "level", "formula", "status", "n_samples", "n_Health",
        "n_other", "covariates_used", "note", "file_prefix"]
with open(design_path, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in design:
        w.writerow(r)

print("\n".join(log))
print("=" * 110)
ok = [r for r in design if r["status"] == "OK"]
print("生成任务数: %d" % len(ok))
by_lvl = collections.Counter(r["level"] for r in ok)
print("按层统计: %s" % dict(sorted(by_lvl.items())))
print("输入文件目录: %s" % OUT_DIR)
print("设计表: %s" % design_path)
print("\n提示：SKIP 原因可通过 _design.tsv 的 status/note 列查看")
