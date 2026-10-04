# -*- coding: utf-8 -*-
"""
Step 00  统一 GMB 队列 metadata 格式
======================================================
输入：GMB/metadata的信息/ 下的三批 metadata tsv（两套表头格式）
输出：GMB/maaslin3_v2/metadata_unified/
        <PROJECT>.tsv        每个队列统一格式的 metadata
        _cohort_summary.tsv  队列级覆盖度/可用性报告
        _unmapped_values.tsv 未能映射的原始取值（供人工核查）

统一后的列：
  ID  group  group_raw  age  sex  BMI  country  continent  disease_icd  project_id  source_folder
"""
import csv, os, collections

# ----------------------------------------------------------------------------
# 路径配置（按需修改）
# ----------------------------------------------------------------------------
BASE = r"D:/Desktop/刘双江-实验室/硕士发表文章/各种材料/4.真实宏基因组样本信息/真实宏基因组样本分析/4.疾病健康队列分析/GMB"
META_ROOT = os.path.join(BASE, "metadata的信息")
OUT_DIR = os.path.join(BASE, "maaslin3_v2", "metadata_unified")

META_DIRS = [
    os.path.join(META_ROOT, "完成分析了的"),
    os.path.join(META_ROOT, "group仅有一列信息"),
    os.path.join(META_ROOT, "group仅有一列信息", "sample-less-50"),
]

# ----------------------------------------------------------------------------
# 取值标准化规则
# ----------------------------------------------------------------------------
# 缺失值标记
NA_TOKENS = {"", "na", "n/a", "nan", "none", "null", "-", "unknown", "not applicable"}


def is_na(v):
    return v is None or str(v).strip().lower() in NA_TOKENS


# 疾病/分组标准化
GROUP_MAP = {
    # --- 对照 ---
    "health": "Health", "healthy": "Health", "control": "Health",
    "controls": "Health", "control group": "Health", "control groups": "Health",
    "normal": "Health", "h": "Health", "healthy control": "Health",
    "healthy controls": "Health", "non-ibd control": "Health",
    # --- IBD ---
    "ibd": "IBD", "inflammatory bowel diseases": "IBD",
    "inflammatory bowel disease": "IBD",
    "crohn's disease": "Crohn's disease", "crohns disease": "Crohn's disease",
    "crohn disease": "Crohn's disease",
    "ulcerative colitis": "Ulcerative colitis", "uc": "Ulcerative colitis",
    # --- 结直肠 ---
    "colorectal neoplasms": "CRC", "colorectal cancer": "CRC", "crc": "CRC",
    "colon cancer": "CRC", "colorectal carcinoma": "CRC",
    "adenoma": "Adenoma", "advanced adenoma": "Adenoma",
    "colorectal adenoma": "Adenoma",
    # --- 代谢 ---
    "obesity": "Obesity", "obese": "Obesity",
    "type 2 diabetes": "T2D", "type 2 diabetes mellitus": "T2D",
    "t2d": "T2D", "t2dm": "T2D", "diabetes": "T2D",
    "thinness": "Thinness", "overweight": "Overweight",
    # --- 其他 ---
    "arthritis, rheumatoid": "RA", "rheumatoid arthritis": "RA",
    "spondylitis, ankylosing": "AS", "ankylosing spondylitis": "AS",
}

# 性别标准化
SEX_MAP = {
    "male": "Male", "m": "Male", "man": "Male",
    "female": "Female", "f": "Female", "woman": "Female",
}

# 国家标准化（合并同一国家的不同写法）
COUNTRY_MAP = {
    "united states of america": "USA",
    "the united states of america": "USA",
    "united states": "USA", "usa": "USA", "us": "USA",
    "united kingdom of great britain and northern ireland": "UK",
    "united kingdom of great britain and northern irela": "UK",
    "united kingdom": "UK", "uk": "UK", "england": "UK",
    "china": "China", "japan": "Japan", "denmark": "Denmark",
    "netherlands": "Netherlands", "germany": "Germany", "spain": "Spain",
    "austria": "Austria", "italy": "Italy", "sweden": "Sweden",
    "luxembourg": "Luxembourg", "australia": "Australia",
    "france": "France", "finland": "Finland", "norway": "Norway",
    "iceland": "Iceland", "estonia": "Estonia", "hungary": "Hungary",
    "slovakia": "Slovakia", "poland": "Poland", "ireland": "Ireland",
    "belgium": "Belgium", "canada": "Canada", "india": "India",
    "korea": "South Korea", "republic of korea": "South Korea",
    "south korea": "South Korea", "singapore": "Singapore",
    "brazil": "Brazil", "mexico": "Mexico", "israel": "Israel",
    "switzerland": "Switzerland", "portugal": "Portugal",
    "greece": "Greece", "czech republic": "Czechia", "czechia": "Czechia",
    "russia": "Russia", "new zealand": "New Zealand",
}

# 国家 → 大洲
COUNTRY2CONT = {
    "China": "Asia", "Japan": "Asia", "India": "Asia", "South Korea": "Asia",
    "Singapore": "Asia", "Israel": "Asia",
    "Denmark": "Europe", "Netherlands": "Europe", "Germany": "Europe",
    "Spain": "Europe", "Austria": "Europe", "Italy": "Europe",
    "Sweden": "Europe", "Luxembourg": "Europe", "UK": "Europe",
    "Slovakia": "Europe", "Hungary": "Europe", "Estonia": "Europe",
    "Finland": "Europe", "Iceland": "Europe", "Norway": "Europe",
    "France": "Europe", "Poland": "Europe", "Ireland": "Europe",
    "Belgium": "Europe", "Switzerland": "Europe", "Portugal": "Europe",
    "Greece": "Europe", "Czechia": "Europe", "Russia": "Europe",
    "USA": "North America", "Canada": "North America", "Mexico": "North America",
    "Brazil": "South America",
    "Australia": "Oceania", "New Zealand": "Oceania",
}

# ----------------------------------------------------------------------------
# 列名识别（兼容两套表头）
# ----------------------------------------------------------------------------
COL_ALIASES = {
    "ID": ["id", "sample", "runid", "run_id", "sampleid", "sample_id"],
    "group": ["group"],
    "age": ["age", "host_age"],
    "sex": ["sex", "gender"],
    "BMI": ["bmi"],
    "country": ["country", "site", "country_or_site"],
    "continent": ["continent"],
    "disease_icd": ["disease"],
    "project_id": ["project_id", "project"],
    "known_species": ["knownspecies"],
    "total_species": ["totalspecies"],
    "unknown_species": ["unknownspecies"],
}

unmapped = []      # (project, field, raw_value)
cohort_rows = []


def find_col(header, key):
    """在表头里找到 key 对应的列索引，找不到返回 None"""
    lower = [h.strip().lower() for h in header]
    for alias in COL_ALIASES[key]:
        # 精确优先
        if alias in lower:
            return lower.index(alias)
    return None


# ----------------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------------
os.makedirs(OUT_DIR, exist_ok=True)

# 建 project -> 文件 映射（后出现的目录优先级更高：sample-less-50 > group仅有一列信息 > 完成分析了的）
proj_files = {}
for d in META_DIRS:
    if not os.path.isdir(d):
        print("[警告] 目录不存在，跳过: %s" % d)
        continue
    for f in sorted(os.listdir(d)):
        if f.endswith(".tsv"):
            proj_files[f[:-4]] = (os.path.join(d, f), os.path.basename(d))

print("发现 metadata 文件: %d 个" % len(proj_files))
print("-" * 100)

OUT_COLS = ["ID", "group", "group_raw", "age", "sex", "BMI",
            "country", "continent", "disease_icd", "project_id", "source_folder"]

for pid in sorted(proj_files):
    fp, src = proj_files[pid]
    with open(fp, encoding="utf-8", errors="ignore", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if len(rows) < 2:
        print("[跳过] %s 行数不足" % pid)
        continue

    header = rows[0]
    ci = {k: find_col(header, k) for k in COL_ALIASES}

    # 逐行转换
    out_rows = []
    grp_c = collections.Counter()
    geo_c = collections.Counter()
    sex_c = collections.Counter()
    n_age = n_sex = n_bmi = n_all3 = 0
    for r in rows[1:]:
        if not r or not str(r[0]).strip():
            continue

        def get(key):
            i = ci.get(key)
            if i is None or i >= len(r):
                return ""
            return str(r[i]).strip()

        sid = get("ID")
        if not sid:
            continue

        # group 标准化
        g_raw = get("group")
        if is_na(g_raw):
            g_std = ""
        else:
            g_std = GROUP_MAP.get(g_raw.lower(), g_raw)   # 未知值原样保留
            if g_raw.lower() not in GROUP_MAP:
                unmapped.append((pid, "group", g_raw))

        # 性别
        s_raw = get("sex")
        if is_na(s_raw):
            sex = ""
        else:
            sex = SEX_MAP.get(s_raw.strip().lower(), "")
            if not sex:
                sex = s_raw
                unmapped.append((pid, "sex", s_raw))

        # 国家
        c_raw = get("country")
        if is_na(c_raw):
            country = ""
        else:
            country = COUNTRY_MAP.get(c_raw.strip().lower(), c_raw)
            if c_raw.strip().lower() not in COUNTRY_MAP:
                unmapped.append((pid, "country", c_raw))

        # 大洲：优先用原表的 continent 列，否则按国家映射
        cont_raw = get("continent")
        if not is_na(cont_raw):
            cont = cont_raw.strip().title()
        else:
            cont = COUNTRY2CONT.get(country, "")

        # 年龄 / BMI（只保留数值）
        def num(v):
            v = v.strip()
            if is_na(v):
                return ""
            try:
                f = float(v)
                return ("%g" % f)
            except ValueError:
                unmapped.append((pid, "numeric", v))
                return ""

        age = num(get("age"))
        bmi = num(get("BMI"))

        icd = get("disease_icd")
        if is_na(icd):
            icd = ""

        if g_std:
            grp_c[g_std] += 1
        if country:
            geo_c[country] += 1
        if sex:
            sex_c[sex] += 1
        if age:
            n_age += 1
        if sex:
            n_sex += 1
        if bmi:
            n_bmi += 1
        if age and sex and bmi:
            n_all3 += 1

        out_rows.append([sid, g_std, g_raw, age, sex, bmi,
                         country, cont, icd, pid, src])

    if not out_rows:
        print("[跳过] %s 无有效样本" % pid)
        continue

    outp = os.path.join(OUT_DIR, pid + ".tsv")
    with open(outp, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(OUT_COLS)
        w.writerows(out_rows)

    n = len(out_rows)
    n_ctrl = grp_c.get("Health", 0)
    n_disease = n - n_ctrl - sum(v for k, v in grp_c.items() if k == "")
    has_ctrl = "Y" if n_ctrl > 0 else "N"
    n_disease_types = len([k for k in grp_c if k and k != "Health"])
    usable = "Y" if (n_ctrl > 0 and n_disease_types > 0) else "N"

    cohort_rows.append({
        "project_id": pid, "source_folder": src, "n_total": n,
        "n_Health": n_ctrl, "n_other": n - n_ctrl,
        "n_age": n_age, "n_sex": n_sex, "n_BMI": n_bmi, "n_age_sex_BMI": n_all3,
        "has_control": has_ctrl, "n_group_levels": len(grp_c),
        "n_geo": len(geo_c), "geo_list": ";".join(sorted(geo_c))[:60],
        "groups": ";".join("%s(%d)" % (k, v) for k, v in grp_c.most_common()),
        "usable_for_differential": usable,
    })

    print("%-14s n=%-5d Health=%-5d age=%-5d sex=%-5d BMI=%-5d 全=%-5d %s | %s"
          % (pid, n, n_ctrl, n_age, n_sex, n_bmi, n_all3,
             "可用" if usable == "Y" else "不可用", src))

# ----------------------------------------------------------------------------
# 报告
# ----------------------------------------------------------------------------
sum_path = os.path.join(OUT_DIR, "_cohort_summary.tsv")
cols = ["project_id", "source_folder", "n_total", "n_Health", "n_other",
        "n_age", "n_sex", "n_BMI", "n_age_sex_BMI", "has_control",
        "n_group_levels", "n_geo", "geo_list", "groups", "usable_for_differential"]
with open(sum_path, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in sorted(cohort_rows, key=lambda x: x["project_id"]):
        w.writerow(r)

if unmapped:
    um_path = os.path.join(OUT_DIR, "_unmapped_values.tsv")
    with open(um_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["project_id", "field", "raw_value"])
        seen = set()
        for rec in unmapped:
            if rec in seen:
                continue
            seen.add(rec)
            w.writerow(rec)

print("-" * 100)
n_usable = sum(1 for r in cohort_rows if r["usable_for_differential"] == "Y")
print("完成: %d 个队列写入 %s" % (len(cohort_rows), OUT_DIR))
print("其中可用于差异分析(有对照+有疾病): %d 个" % n_usable)
print("汇总报告: %s" % sum_path)
if unmapped:
    print("未映射取值: %s （请人工核查）" % um_path)
