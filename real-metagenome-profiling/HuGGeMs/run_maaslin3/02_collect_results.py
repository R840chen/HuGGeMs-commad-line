# -*- coding: utf-8 -*-
"""
Step 02  汇总 MaAsLin3 分层结果，找出「跨层稳健」的差异物种
=========================================================
输入：GMB/maaslin3_v2/Maaslin3-result-v2/<PID>/<LEVEL>/all_results.tsv
输出：GMB/maaslin3_v2/summary/
        _all_significant.tsv    所有显著结果明细
        _cross_level_summary.tsv 物种 × 队列 × 各层显著性（宽表）
        _robust_species.tsv     在「所有可用层」都显著的稳健物种（核心产出）
        _count_by_cohort.tsv    每个 队列×层 的显著物种数

判定口径：
  - 只看疾病主效应：metadata == "group" 且 value != "Health"
  - 显著性：qval_joint < Q_THRESHOLD（默认 0.05，与已发表口径一致）
  - 方向一致性：coef 符号（>0 富集 / <0 缺失）
"""
import csv, os, collections

# ============================================================================
Q_THRESHOLD = 0.05     # 显著性阈值（与已发表 261+237 的口径一致）
MODELS = ["abundance", "prevalence"]   # 两个模型分别看
# ============================================================================

BASE = r"D:/Desktop/刘双江-实验室/硕士发表文章/各种材料/4.真实宏基因组样本信息/真实宏基因组样本分析/4.疾病健康队列分析/GMB"
RES_DIR = os.path.join(BASE, "maaslin3_v2", "Maaslin3-result-v2")
OUT_DIR = os.path.join(BASE, "maaslin3_v2", "summary")
REF_GROUP = "Health"

os.makedirs(OUT_DIR, exist_ok=True)


def read_tsv(path):
    with open(path, encoding="utf-8", errors="ignore", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return [], []
    return rows[0], rows[1:]


# ---------------------------------------------------------------------------
# 1. 扫描所有 all_results.tsv
# ---------------------------------------------------------------------------
found = []
for pid in sorted(os.listdir(RES_DIR)) if os.path.isdir(RES_DIR) else []:
    pdir = os.path.join(RES_DIR, pid)
    if not os.path.isdir(pdir):
        continue
    for lvl in sorted(os.listdir(pdir)):
        f = os.path.join(pdir, lvl, "all_results.tsv")
        if os.path.exists(f):
            found.append((pid, lvl, f))

print("找到 all_results.tsv: %d 个" % len(found))
if not found:
    raise SystemExit("没有结果文件，请先运行 run_maaslin3_by_level.R")

# ---------------------------------------------------------------------------
# 2. 逐个读取，筛选疾病主效应
# ---------------------------------------------------------------------------
sig_rows = []       # 明细
by_level = collections.defaultdict(set)   # (pid, lvl) -> {(feature, model, sign)}
level_counts = []

for pid, lvl, fp in found:
    hdr, rows = read_tsv(fp)
    ix = {c: i for i, c in enumerate(hdr)}

    def g(r, name):
        i = ix.get(name)
        return r[i].strip() if (i is not None and i < len(r)) else ""

    n_ab = n_pv = 0
    for r in rows:
        if not r or len(r) < 2:
            continue
        meta_var = g(r, "metadata")
        val = g(r, "value")
        model = g(r, "model")
        err = g(r, "error")

        if meta_var != "group":
            continue
        if val == REF_GROUP:
            continue
        if err not in ("", "NA", "nan"):
            continue
        try:
            qj = float(g(r, "qval_joint"))
        except ValueError:
            continue
        if not (qj < Q_THRESHOLD):
            continue

        feat = g(r, "feature")
        try:
            coef = float(g(r, "coef"))
        except ValueError:
            continue
        sign = "UP" if coef > 0 else "DOWN"

        sig_rows.append({
            "project_id": pid, "level": lvl, "feature": feat, "model": model,
            "coef": coef, "direction": sign, "qval_joint": qj,
            "qval_individual": g(r, "qval_individual"),
        })
        by_level[(pid, lvl)].add((feat, model, sign))
        if model == "abundance":
            n_ab += 1
        elif model == "prevalence":
            n_pv += 1

    level_counts.append({"project_id": pid, "level": lvl,
                         "n_abundance_sig": n_ab, "n_prevalence_sig": n_pv,
                         "n_features": len({f for f, m, s in by_level[(pid, lvl)]})})
    print("  %-12s %-3s  abundance=%-5d prevalence=%-5d" % (pid, lvl, n_ab, n_pv))

# ---------------------------------------------------------------------------
# 3. 写明细
# ---------------------------------------------------------------------------
det_path = os.path.join(OUT_DIR, "_all_significant.tsv")
with open(det_path, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
        "project_id", "level", "feature", "model", "coef", "direction",
        "qval_joint", "qval_individual"])
    w.writeheader()
    for r in sig_rows:
        w.writerow(r)
print("\n明细: %s (%d 行)" % (det_path, len(sig_rows)))

# ---------------------------------------------------------------------------
# 4. 汇总每个队列×层的计数
# ---------------------------------------------------------------------------
cnt_path = os.path.join(OUT_DIR, "_count_by_cohort.tsv")
with open(cnt_path, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
        "project_id", "level", "n_abundance_sig", "n_prevalence_sig", "n_features"])
    w.writeheader()
    for r in sorted(level_counts, key=lambda x: (x["project_id"], x["level"])):
        w.writerow(r)
print("计数: %s" % cnt_path)

# ---------------------------------------------------------------------------
# 5. 跨层稳健性：对每个队列每个模型，找出在「所有可用层」都显著的物种
#    注意：L1 是所有队列都有的层，所以以「有 L1 的队列」为基准
# ---------------------------------------------------------------------------
# 该队列在该模型下有哪些层产生了结果
cohort_levels = collections.defaultdict(set)
for (pid, lvl) in by_level:
    cohort_levels[pid].add(lvl)

robust = []
cross = []

for pid in sorted(cohort_levels):
    lvls = sorted(cohort_levels[pid])
    for model in MODELS:
        # 该队列该模型下，每个层显著物种集合
        sets = {}
        for lvl in lvls:
            sets[lvl] = {f for (f, m, s) in by_level[(pid, lvl)] if m == model}
        if not any(sets.values()):
            continue
        # 取所有层都有结果的层
        usable = [l for l in lvls if l in sets]
        if not usable:
            continue
        inter = set.intersection(*[sets[l] for l in usable]) if usable else set()
        for feat in inter:
            # 各层方向是否一致
            dirs = []
            for lvl in usable:
                d = [s for (f, m, s) in by_level[(pid, lvl)]
                     if f == feat and m == model]
                dirs.append(d[0] if d else "NA")
            consistent = len(set(dirs)) == 1
            robust.append({
                "project_id": pid, "model": model, "feature": feat,
                "n_levels_significant": len(usable),
                "levels": ";".join(usable),
                "direction": dirs[0] if dirs else "NA",
                "direction_consistent": "Y" if consistent else "N",
            })
        # 每层的出现情况（宽表）
        for feat in set().union(*sets.values()) if sets else set():
            rec = {"project_id": pid, "model": model, "feature": feat}
            for lvl in lvls:
                rec[lvl] = ""
                d = [s for (f, m, s) in by_level[(pid, lvl)]
                     if f == feat and m == model]
                if d:
                    rec[lvl] = d[0]
            rec["n_levels"] = sum(1 for lvl in lvls if rec.get(lvl))
            cross.append(rec)

cross_lvls = sorted({r for rec in cross for r in rec
                     if r.startswith("L") and r not in ("L",)})
cross_path = os.path.join(OUT_DIR, "_cross_level_summary.tsv")
with open(cross_path, "w", encoding="utf-8", newline="") as fh:
    cols = ["project_id", "model", "feature"] + cross_lvls + ["n_levels"]
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for r in sorted(cross, key=lambda x: (x["project_id"], x["model"], -x["n_levels"])):
        w.writerow(r)
print("跨层宽表: %s" % cross_path)

rob_path = os.path.join(OUT_DIR, "_robust_species.tsv")
with open(rob_path, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
        "project_id", "model", "feature", "n_levels_significant",
        "levels", "direction", "direction_consistent"])
    w.writeheader()
    for r in sorted(robust, key=lambda x: (x["project_id"], x["model"], x["feature"])):
        w.writerow(r)
print("稳健物种: %s (%d 行)" % (rob_path, len(robust)))

# ---------------------------------------------------------------------------
# 6. 控制台摘要
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("稳健物种（在所有可用层都显著）Top 20:")
robust_sorted = sorted(robust, key=lambda x: -x["n_levels_significant"])
for r in robust_sorted[:20]:
    print("  %-12s %-11s %-42s 层=%s 方向=%s%s"
          % (r["project_id"], r["model"], r["feature"][:42],
             r["levels"], r["direction"],
             "" if r["direction_consistent"] == "Y" else " (方向不一致!)"))
