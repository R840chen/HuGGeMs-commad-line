# -*- coding: utf-8 -*-
# =============================================================================
#  03_analyze_robustness.py
#  在 02_collect_results.py 的汇总结果之上，做「跨队列一致性 × 协变量校正」交叉分析
# =============================================================================
#  输入：maaslin3_v2/summary/_all_significant.tsv
#        maaslin3_v2/summary/_robust_species.tsv
#  输出：maaslin3_v2/summary/_core_candidate_species.tsv   （核心候选物种）
#        maaslin3_v2/summary/_cross_cohort_consistency.tsv （跨队列一致性矩阵）
#
#  用法：python 03_analyze_robustness.py
# =============================================================================
import csv, os, collections, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUM = os.path.join(HERE, "summary")

# ---- 参数区 ----------------------------------------------------------------
QLVL = 0.05                       # 显著性阈值（与已发表 261+237 口径一致）
MAIN_IBD = ["PRJEB1220", "PRJNA389280", "PRJNA398089", "SRP057027"]  # IBD 4 主队列
COV_COHORT = "PRJEB1220"          # 协变量最全的队列（L1-L5+LG）
COV_LEVEL = "L5"                  # 全协变量层（~ group + age + sex + BMI）
LEVEL_ORDER = ["L1", "L2", "L3", "L4", "L5", "LG"]
# ---------------------------------------------------------------------------

def load(name):
    fp = os.path.join(SUM, name)
    if not os.path.exists(fp):
        sys.exit(f"[错误] 缺少 {fp}，请先运行 02_collect_results.py")
    with open(fp, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

sig = load("_all_significant.tsv")
rob = load("_robust_species.tsv")
print(f"读入显著记录 {len(sig)} 条，稳健物种 {len(rob)} 行")

# ---------- 1) 跨队列一致性（L1，主队列） ----------
cross = collections.defaultdict(lambda: collections.defaultdict(dict))  # model -> feature -> {pid: dir}
for r in sig:
    if r["level"] == "L1" and r["project_id"] in MAIN_IBD:
        cross[r["model"]][r["feature"]][r["project_id"]] = r["direction"]

consistent = {}
for model, feats in cross.items():
    for f, v in feats.items():
        if len(v) == len(MAIN_IBD) and len(set(v.values())) == 1:
            consistent[(model, f)] = list(v.values())[0]
print(f"\n[1] 跨 {len(MAIN_IBD)} 个 IBD 主队列方向一致的 (模型,物种): {len(consistent)}")

# ---------- 2) 协变量校正表现 ----------
p = collections.defaultdict(dict)     # (model,feature) -> {level: dir}
for r in sig:
    if r["project_id"] == COV_COHORT:
        p[(r["model"], r["feature"])][r["level"]] = r["direction"]

rows, n_keep = [], 0
for (m, f), d in sorted(consistent.items()):
    cells = []
    for l in LEVEL_ORDER:
        dd = p.get((m, f), {}).get(l)
        cells.append("" if dd is None else ("+" if dd == d else "x"))
    keep = cells[LEVEL_ORDER.index(COV_LEVEL)] == "+"
    if keep:
        n_keep += 1
    rows.append((m, f, d, keep, cells))

print(f"[2] 其中在 {COV_COHORT} 的 {COV_LEVEL} 层仍显著且方向一致: {n_keep}")
print(f"\n{'模型':<11}{'物种':<42}{'方向':<7}" + "".join(f"{l:>5}" for l in LEVEL_ORDER))
for m, f, d, keep, cells in rows:
    mark = " ★" if keep else ""
    print(f"{m:<11}{f:<42}{d:<7}" + "".join(f"{c:>5}" for c in cells) + mark)

# ---------- 3) 输出跨队列一致性矩阵 ----------
out1 = os.path.join(SUM, "_cross_cohort_consistency.tsv")
with open(out1, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["model", "feature", "direction_in_" + COV_COHORT] +
               [f"IBD4_{pid}" for pid in MAIN_IBD] +
               [f"{COV_COHORT}_{l}" for l in LEVEL_ORDER] + ["robust_after_covariate"])
    for m, f, d, keep, cells in rows:
        v = cross[m][f]
        w.writerow([m, f, d] + [v.get(pid, "") for pid in MAIN_IBD] + cells + ["Y" if keep else "N"])
print(f"\n已写出: {out1}")

# ---------- 4) 核心候选：跨队列一致 且 协变量校正后稳健 ----------
core = [(m, f, d) for m, f, d, keep, _ in rows if keep]
out2 = os.path.join(SUM, "_core_candidate_species.tsv")
with open(out2, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["model", "feature", "direction", f"n_IBD_cohorts_consistent", f"{COV_COHORT}_{COV_LEVEL}_significant"])
    for m, f, d in core:
        w.writerow([m, f, d, len(MAIN_IBD), "Y"])
print(f"已写出: {out2}   （{len(core)} 个核心候选物种）")

# ---------- 5) 附：PRJEB1220 L1 与 L5 双稳健 ----------
l1 = {(r["model"], r["feature"]): r["direction"] for r in sig
      if r["project_id"] == COV_COHORT and r["level"] == "L1"}
l5 = {(r["model"], r["feature"]): r["direction"] for r in sig
      if r["project_id"] == COV_COHORT and r["level"] == COV_LEVEL}
both = {k for k in set(l1) & set(l5) if l1[k] == l5[k]}
na = len([1 for m, f in both if m == "abundance"])
npr = len([1 for m, f in both if m == "prevalence"])
print(f"\n[3] {COV_COHORT}: L1 与 {COV_LEVEL} 均显著且方向一致 —— 丰度 {na} 个、患病率 {npr} 个（共 {len(both)}）")

print("\n=== 分析完成 ===")
