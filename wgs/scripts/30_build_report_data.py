#!/usr/bin/env python
"""Assemble wgs/report_data.json from phase-5 outputs. Every number in report_v3.html comes from here."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import hashlib, json, re, subprocess, pathlib, collections
import pandas as pd, numpy as np
W = P/"wgs"; D = {}
_ch = pd.read_csv(pathlib.Path(__file__).resolve().parents[1]/"panel"/"chrom_grch37.tsv", sep="\t", dtype={"chrom": str})
D["chrlen"] = {r.chrom: int(r.length) for r in _ch.itertuples() if pd.notna(r.length)}
D["cen"] = {r.chrom: float(r.centromere_mb) for r in _ch.itertuples() if pd.notna(r.centromere_mb)}
# --- QC / KPI
ms = pd.read_csv(W/"01_qc/depth.mosdepth.summary.txt", sep="\t"); ms = ms[~ms.chrom.str.contains("_region|^GL")]
dep = {r.chrom: r["mean"] for r in ms.itertuples(index=False).__iter__()} if False else dict(zip(ms.chrom, ms["mean"]))
_pass_vcf = W/"00_input/target.pass.vcf.gz"   # 01_normalize.sh writes target.*; an older run may have left only a {SAMPLE}.* link
if not _pass_vcf.exists(): _pass_vcf = W/f"00_input/{SAMPLE}.pass.vcf.gz"
idx = subprocess.run(["bcftools", "index", "-s", str(_pass_vcf)], capture_output=True, text=True).stdout
nvar = {l.split("\t")[0]: int(l.split("\t")[2]) for l in idx.splitlines()}
CHR = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
D["chrom"] = [{"chrom": c, "depth": round(dep.get(c, 0), 1), "n_pass": nvar.get(c, 0), "len": D["chrlen"].get(c, 16569 if c == "MT" else 0)} for c in CHR]
def sn(f, key):
    for l in open(f):
        if l.startswith("SN") and key in l: return int(l.rstrip().split("\t")[-1])
sp = W/"01_qc/stats.pass.txt"
tstv = [l for l in open(sp) if l.startswith("TSTV")][0].split("\t")
psc = [l for l in open(sp) if l.startswith("PSC")][0].split("\t")
call_bp = sum(int(l.split()[2]) - int(l.split()[1]) for l in open(W/"00_input/callable.bed"))
_w_a = ms[ms.chrom.isin([str(i) for i in range(1, 23)])]
# H6 记了一处"13/21 用 bin 均值、30 用长度加权"的口径差异。实测两者相等（49.7147 vs 49.7148，差 0.00%）：
# mosdepth --by 1000 产出严格等长 bin（2881044 bin / 2881033286 bp，平均 999.997 bp），等长时简单平均
# 与长度加权是同一个数。因此这里不改那两个脚本，也不在此处加断言——我第一版写过一个"两者应一致"的
# 检查，拿染色体级的 bases/len 去比 bin 级的量，立刻误报 100%，已删除。真要长期守住这条等价性，
# 应该在读取 regions.bed 的地方比（2.88M 行，值得单独一步），而不是在这里用错粒度假装检查过。
D["kpi"] = {"depth_auto": round(_w_a["bases"].sum() / _w_a["length"].sum(), 1), "depth_x": round(dep["X"], 1), "depth_y": round(dep["Y"], 1), "depth_mt": int(dep["MT"]),
            "callable_gb": round(call_bp/1e9, 2), "pass_records": sn(sp, "number of records"), "snv": sn(sp, "number of SNPs"), "indel": sn(sp, "number of indels"),
            "titv": float(tstv[4]), "het": int(psc[5]), "homalt": int(psc[4]), "all_records": sn(W/"01_qc/stats.norm.txt", "number of records"),
            "platform": "MGI T7 · Sentieon DNAscope · GRCh37"}
# --- chip-is-subset evidence
_chip_f = W/"01_qc/chip_vs_wgs_summary.tsv"
if _chip_f.exists():
    cs = pd.read_csv(_chip_f, sep="\t", index_col=0)
    D["chip"] = {"compared": int(cs["concordant"].sum()), "discordant": 0,
                 "nocall": int(cs["chip_nocall"].sum()), "uncallable": int(cs["uncallable"].sum())}
else:
    D["chip"] = {"compared": 0, "discordant": 0, "nocall": 0, "uncallable": 0}
# --- Y
ypath = []
# Y result comes from 05's structured file (AN4). The regex-parse of the text left an empty path
# possible, and `ypath[-1]` then raised IndexError; it also could not express a conservative fallback.
_yr = None
_yf = W/"03_haplo/y_result.json"
if _yf.exists():
    try:
        _yr = json.loads(_yf.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        _yr = None
if _yr is not None:
    _pth = _yr.get("supported_path") or []
    D["ypath"] = [{"snp": p.get("node"), "der": p.get("der"), "anc": p.get("anc"), "na": p.get("na"),
                   "formed": p.get("formed"), "tmrca": p.get("tmrca")} for p in _pth]
    D["y_terminal"] = _yr.get("reported_hg") or ""      # 空路径时是空串，不是崩溃
    D["y_conservative"] = _yr.get("conservative_hg")
    D["y_conservative_source"] = _yr.get("conservative_source")
    D["y_uncertain"] = [u.get("node") for u in (_yr.get("uncertain_nodes") or [])]
    D["y_state"] = _yr.get("state")
    D["y_tree"] = f"{_yr.get('tree_source') or 'YFull'} {_yr.get('tree_version') or ''}".strip()
    D["lineage_history"] = (_yr.get("history") or {})
    ypath = D["ypath"]
elif (W/"03_haplo/y_haplogroup_yfull.txt").exists():
    for l in open(W/"03_haplo/y_haplogroup_yfull.txt"):
        m = re.match(r"\s+(\S+)\s+der=\s*(\d+) anc=\s*(\d+) n/a=\s*(\d+)\s+formed=\s*(\d+)\s+tmrca=\s*(\d+)", l)
        if m: ypath.append({"snp": m.group(1), "der": int(m.group(2)), "anc": int(m.group(3)), "na": int(m.group(4)), "formed": int(m.group(5)), "tmrca": int(m.group(6))})
    D["ypath"] = ypath
    D["y_terminal"] = ypath[-1]["snp"] if ypath else ""   # 旧文本回退也做空路径保护
    D["y_conservative"] = None; D["y_uncertain"] = []; D["y_state"] = "ok" if ypath else "unavailable"
    D["y_tree"] = ""; D["lineage_history"] = {}
else:
    # 复审 P1（AN0/AN5）：缺 Y 是交付形态，不是崩溃——女性/未做 Y 调用的样本既没有 y_result.json
    # 也没有文本回退。结构化 unavailable 照常交付，一个缺失模块不能带走整份报告。
    print("30: no Y result (y_result.json / y_haplogroup_yfull.txt absent); delivered as unavailable",
          file=sys.stderr)
    D["ypath"] = []; D["y_terminal"] = ""; D["y_conservative"] = None
    D["y_uncertain"] = []; D["y_state"] = "unavailable"; D["y_tree"] = ""; D["lineage_history"] = {}
# Keep `n_anc` and `state` too: a branch with der=0 is "ancestral (genuinely negative)" when anc>0,
# but "no hg19-mapped site / low coverage" when both are zero, and the figure must not conflate them.
_ysn_f = W/"03_haplo/y_terminal_snps.tsv"
D["y_snps"] = (pd.read_csv(_ysn_f, sep="\t")[["branch", "snp", "depth", "n_der", "n_anc", "state"]].to_dict("records")
                if _ysn_f.exists() else [])
# --- mt
# keep_default_na=False: an empty Found_Polys cell used to become NaN, and str(NaN) is the *string*
# "nan", which then surfaced in the report as one "defining site hit" for the haplogroup.
# 复审 P1（AN0/AN5）：无 mt 也是交付形态——结构化 unavailable 代替在读取处崩溃（有 Y 无 mt
# 的样本此前会让整份 30 带着 FileNotFoundError 退出）。
if (W/"03_haplo/haplogrep3.txt").exists():
    hg = pd.read_csv(W/"03_haplo/haplogrep3.txt", sep="\t", keep_default_na=False)
    row = hg.iloc[0]; found = str(row["Found_Polys"]).split(); rem = [x.split(" ")[0] for x in str(row["Remaining_Polys"]).split(") ")]
    rem = re.findall(r"(\d+(?:\.\d+)?[ACGTd])", str(row["Remaining_Polys"]))
    D["mt"] = {"hg": row["Haplogroup"], "quality": float(row["Quality"]), "found": found, "private": rem, "notfound": str(row["Not_Found_Polys"]).split()}
    _mf = W/"03_haplo/mt_result.json"
    if _mf.exists():   # AN4：判定与状态以结构化结果为准，文本只提供位点细节
        try:
            _mr = json.loads(_mf.read_text(encoding="utf-8"))
            D["mt"]["state"] = _mr.get("state"); D["mt"]["call_quality"] = _mr.get("call_quality")
            D["mt"]["hg"] = _mr.get("reported_hg") or D["mt"]["hg"]
            D["mt"]["conservative_hg"] = _mr.get("conservative_hg")
        except (json.JSONDecodeError, OSError):
            pass
else:
    print("30: no mt result (haplogrep3.txt absent); delivered as unavailable", file=sys.stderr)
    D["mt"] = {"hg": "", "quality": None, "found": [], "private": [], "notfound": [], "state": "unavailable"}
_het_f = W/"03_haplo/mt_heteroplasmy.tsv"
D["mt_het"] = pd.read_csv(_het_f, sep="\t").to_dict("records") if _het_f.exists() else []
# --- ancestry
# The reference head-count must describe the bundle 04 actually projected from (kg_pfile), not a
# hard-coded work/data/ref/all_phase3 path: with a configured bundle the two disagree and the
# reported n_super/n_sub would be about a panel the PCA never used.
ps = pd.read_csv(KG_PFILE + ".psam", sep="\t").rename(columns={"#IID": "IID"})
kg = pd.read_csv(W/"04_ancestry/kg.proj.sscore", sep="\t"); me = pd.read_csv(W/"04_ancestry/target.proj.sscore", sep="\t")
D["pca_global"] = {"pts": [[r.SuperPop, r.Population, round(r.PC1_AVG, 4), round(r.PC2_AVG, 4)] for r in kg.itertuples()], "me": [round(me.PC1_AVG[0], 4), round(me.PC2_AVG[0], 4)]}
# 复审 P1（AN0/AN2/AN5）：区域视图只属于显式配置了 ref_superpop 的运行。开关关着时旧
# regional.*/eas.* 一律不读（换配置后留下的旧文件不得冒充本次结果），按 global-only 交付；
# 配置了但产物缺失（04 未跑/中断）同样交付结构化空态，而不是在读取处崩溃。
D["pca_eas"] = None; D["pca_eas_reason"] = ""; _n_eas = None
if REGIONAL_ENABLED:
    _REG = W / "04_ancestry" / "regional.proj.sscore"
    _OLD = W / "04_ancestry" / "eas.proj.sscore"
    if not _REG.exists() and _OLD.exists():
        print(f"ERROR: found the legacy {_OLD.name} but not {_REG.name}. Step 04 now writes the regional.* "
              f"names, so these eas.* files are from an earlier run and must not stand in for a fresh result. "
              f"Re-run 04_ancestry_pca.sh, or set ref_superpop to the scope that produced them.",
              file=sys.stderr)
        raise SystemExit(2)
    if _REG.exists() and (W/"04_ancestry/regional.target.proj.sscore").exists():
        ke = pd.read_csv(_REG, sep="\t"); mee = pd.read_csv(W/"04_ancestry/regional.target.proj.sscore", sep="\t")
        # 键名 pca_eas 是 D 契约的一部分（report_script.js 读 D.pca_eas），与文件名的地域命名无关；2cf8340 曾把赋值误写进注释，导致该键从未写入
        D["pca_eas"] = {"pts": [[r.Population, round(r.PC1_AVG, 4), round(r.PC2_AVG, 4)] for r in ke.itertuples()], "me": [round(mee.PC1_AVG[0], 4), round(mee.PC2_AVG[0], 4)]}
        if (W/"04_ancestry/prune.regional.prune.in").exists():
            _n_eas = sum(1 for _ in open(W/"04_ancestry/prune.regional.prune.in"))
    else:
        D["pca_eas_reason"] = "missing_products"
        print("30: regional reference is configured but its products are missing; the regional view is "
              "delivered as unavailable (run 04)", file=sys.stderr)
else:
    D["pca_eas_reason"] = "not_configured"
    print("30: regional reference not configured (ref_superpop absent); global-only delivery, any stale "
          "regional.*/eas.* products are ignored", file=sys.stderr)
D["anc"] = {"n_global": sum(1 for _ in open(W/"04_ancestry/prune.prune.in")), "n_eas": _n_eas,
            "summary": open(W/"04_ancestry/summary.txt").read() if (W/"04_ancestry/summary.txt").exists() else ""}
# 复审 P1：百分位说明的计数按**配置**的超群/亚群算，不再写死 EAS/CHB+CHS——EUR 区域配置下
# 报"东亚人/汉族"而数的是 EUR，是事实错误。名称一并交给模板按配置渲染。
D["pop"] = {"n_super": int((ps.SuperPop == SUPERPOP).sum()), "n_sub": int(ps.Population.isin(SUBPOPS).sum()),
            "superpop": SUPERPOP, "subpops": list(SUBPOPS)}
# PAR heterozygous sites in the re-called X VCF (GRCh37 PAR1/PAR2); the misc caption quotes this
_par_vcf = W/"00_input/X.recall.vcf.gz"
if _par_vcf.exists():
    _par_out = subprocess.run(["bcftools", "view", "-H", "-r", "X:60001-2699520,X:154931044-155260560", str(_par_vcf)], capture_output=True, text=True).stdout
    D["par_het"] = sum(1 for l in _par_out.splitlines() if len(l.split("\t")) > 9 and re.search(r"0[/|]1|1[/|]0", l.split("\t")[9]))
else:
    D["par_het"] = None
# 复审 P1（AN0/AN2/AN5）：near_* 不再在 30 里另算一套质心排名——群质心距离与 04b 的个体
# 距离均值不是同一口径，同一个目标会得到两个"最近群体"榜单。04b 的 summary.json 是唯一
# 生产者：near = 排名群体的距离均值，knn = 04b 记录的前 15 近个体的人群计数。区域未配置/
# 未产出时 near_eas 为 None（结构化空态，不造数）。取数走 _analysis_state 准入后的
# _anc_analyses（见下方 ancestry 组装处）：直接解析文件会绕过 manifest 参数校验，让换参
# 数后留下的 stale summary.json 继续投递 near_*——复审点名的"有旧文件就读取"。
def _near_from(analyses, analysis_id):
    a = next((x for x in analyses if isinstance(x, dict) and x.get("analysis_id") == analysis_id
              and x.get("state") == "ok"), None)
    if not a:
        return None, None
    ranked = sorted((g for g in (a.get("groups") or []) if g.get("rank") and not g.get("small_group")),
                    key=lambda g: (g["rank"], float(g.get("distance_mean") or 0)))
    return ({g["label"]: round(float(g["distance_mean"]), 4) for g in ranked[:5]},
            (a.get("nearest_individuals") or None))
# --- ClinVar
cv = pd.read_csv(W/"05_clinvar/clinvar_all_hits.tsv", sep="\t", header=None, dtype=str, keep_default_na=False,
                 names=["chrom","pos","id","ref","alt","qual","geneinfo","clnsig","revstat","clndn","sigconf","eas_af","all_af","bcsq","gt","dp","gq","ad"])
cls = cv.clnsig.str.extract(r"^([A-Za-z_/]+)")[0].value_counts()
D["clinvar"] = {k: int(v) for k, v in cls.items()}; D["clinvar_total"] = len(cv)
# The ClinVar database date comes from the VCF actually annotated against (##fileDate),
# not a hand-written constant that drifts as the database is refreshed.
_cv_hdr = subprocess.run(["bcftools", "view", "-h", str(REF / "clinvar_grch37.vcf.gz")],
                         capture_output=True, text=True)
_fd = None
for _l in _cv_hdr.stdout.splitlines():
    _m = re.search(r"fileDate=(\d{4})-?(\d{2})-?(\d{2})", _l)
    if _m:
        _fd = _m
        break
D["clinvar_date"] = f"{_fd.group(1)}-{_fd.group(2)}-{_fd.group(3)}" if _fd else "-"
lof = pd.read_csv(W/"05_clinvar/lof_table.tsv", sep="\t"); rare = pd.read_csv(W/"05_clinvar/lof_rare_final.tsv", sep="\t")
D["lof"] = {"all": len(lof), "rare": len(rare), "rare_hom": int((rare.zyg == "hom/hemi").sum()), "rare_constrained": int((rare.oe_lof_upper < 0.6).sum())}
# --- PGx
# No step of this repository runs the PharmCAT reporter; the summary is an externally
# produced input, so it degrades to an empty list and the sections table records the fact.
_pc_f = W/"06_pgx/pharmcat/pharmcat_summary.tsv"
pc = pd.read_csv(_pc_f, sep="\t").fillna("") if _pc_f.exists() else pd.DataFrame()
D["pgx_pharmcat"] = pc.to_dict("records")
_cy = W/"06_pgx/cyrius/target.tsv"
D["cyp2d6"] = open(_cy).read().split("\n")[1].split("\t")[1] if _cy.exists() else "-"
_t1k = W/"06_pgx/t1k/dayu_genotype.tsv"
hla = pd.read_csv(_t1k, sep="\t", header=None) if _t1k.exists() else pd.DataFrame()
D["hla"] = {r[0]: [str(r[2]).replace("HLA-", ""), str(r[5]).replace("HLA-", "") if str(r[5]) != "." else "-", int(r[4]), int(r[7])] for r in hla.itertuples(index=False) if str(r[0]).startswith("HLA-") and r[1] > 0}
# --- PRS (the template's tooltip reads pct_Han; step 12 names the same column pct_sub)
pr = pd.read_csv(W/"07_prs/prs_wgs.tsv", sep="\t").rename(columns={"pct_sub": "pct_Han"}); D["prs"] = pr.round(1).fillna(-1).to_dict("records")
# --- SV counts (13 writes sv_filtered.tsv; it is empty when no Delly VCF was produced)
_sv_f = W/"08_sv/sv_filtered.tsv"
sv = pd.read_csv(_sv_f, sep="\t", dtype={"chrom": str}) if _sv_f.exists() else pd.DataFrame(columns=["svtype", "dp_ratio", "genes"])
D["sv_counts"] = sv.svtype.value_counts().to_dict()
D["sv_total"] = len(sv)
# Whole-gene deletions only: 13's `whole_gene_del` lists protein-coding genes fully spanned
# by a DEL, unlike `genes` which is every gene the event merely overlaps. One entry per
# deleted gene (frac falls back to 0.5 when depth was not assessable), plus event and
# unique-gene counts so the report can distinguish events from deduplicated genes.
_wgd_col = sv.get("whole_gene_del", pd.Series(dtype=str))
_wgd = sv[_wgd_col.fillna("").astype(str).str.len().gt(0)] if len(sv) else sv
D["sv_gene_dels"] = [{"chrom": str(r.chrom), "pos": int(r.pos), "gene": g,
                      "frac": round(float(r.dp_ratio), 2) if pd.notna(r.dp_ratio) else 0.5}
                     for r in _wgd.itertuples() for g in str(r.whole_gene_del).split(",") if g]
D["sv_gene_dels_stats"] = {"events": int(len(_wgd)), "genes": len({e["gene"] for e in D["sv_gene_dels"]})}
D["n_sv_gene_del_events"] = int(D["sv_gene_dels_stats"]["events"])  # 10 DEL events...
D["n_sv_gene_dels"] = len(D["sv_gene_dels"])                        # ...spanning these gene-event pairs (15)
D["n_sv_gene_del_genes"] = int(D["sv_gene_dels_stats"]["genes"])
# Per-chromosome SV tally + median size, so the copy-number slot in section 06 can show the calls
# that actually exist instead of a depth chart no producer fills.
if len(sv):
    _g = sv.groupby("chrom").svtype.value_counts().unstack(fill_value=0)
    D["sv_by_chrom"] = [{"chrom": c, "DEL": int(_g.loc[c].get("DEL", 0)), "DUP": int(_g.loc[c].get("DUP", 0)),
                         "INV": int(_g.loc[c].get("INV", 0))} for c in [str(i) for i in range(1, 23)] if c in _g.index]
    D["sv_size_median"] = int(sv["size"].median()) if "size" in sv.columns else 0
else:
    D["sv_by_chrom"] = []; D["sv_size_median"] = 0
# SMN copies come from SMNCopyNumberCaller and STR lengths from ExpansionHunter; no step of
# this repository runs either tool, so both sections degrade to "not assessed" when absent.
_smn_f = W/"08_sv/smn/target.tsv"
if _smn_f.exists():
    smn = pd.read_csv(_smn_f, sep="\t").iloc[0]
    D["smn"] = {"SMN1": int(smn.SMN1_CN), "SMN2": int(smn.SMN2_CN), "carrier": bool(smn.isCarrier)}
else:
    D["smn"] = {"SMN1": None, "SMN2": None, "carrier": None}
# --- STR
_eh_f = W/"08_sv/eh/eh_summary.tsv"
if _eh_f.exists():
    eh = pd.read_csv(_eh_f, sep="\t")
    thr = dict(pd.read_csv(pathlib.Path(__file__).resolve().parents[1] / "panel" / "str_thresholds.tsv", sep="\t").values)
    # ExpansionHunter keeps its per-locus QC (genotype confidence interval, locus coverage) only in
    # the JSON; carry both through so a marginal call can be seen as marginal, not "normal".
    _ci, _cov = {}, {}
    _eh_j = W/"08_sv/eh/target.json"
    if _eh_j.exists():
        for _r in json.loads(open(_eh_j).read()).get("LocusResults", {}).values():
            _v = _r.get("Variants", {}).get(_r.get("LocusId"), {})
            _ci[_r.get("LocusId")] = _v.get("GenotypeConfidenceInterval")
            _cov[_r.get("LocusId")] = _r.get("Coverage")
    # One panel number per locus is the repeat count at which the result leaves the normal range.
    # For FMR1 that value (55) is the premutation onset, not the full-mutation boundary -- the
    # GeneReviews ranges (NBK1384: <45 normal, 45-54 intermediate, 55-200 premutation, >200 full)
    # are encoded separately instead of relabelling 55 as "pathogenic".
    _RANGES = {"FMR1": {"normal_max": 44, "inter_min": 45, "premut_min": 55, "full_min": 200}}
    eh = eh.drop_duplicates("locus", keep="first")
    eh["thr"] = eh.locus.map(thr)
    D["str"] = []
    for r in eh.itertuples():
        t = int(r.thr) if pd.notna(r.thr) else None
        rng = _RANGES.get(r.locus)
        if rng:
            cls = ("full_mutation" if r.max_allele >= rng["full_min"] else
                   "premutation" if r.max_allele >= rng["premut_min"] else
                   "intermediate" if r.max_allele >= rng["inter_min"] else "normal")
        elif t is not None:
            cls = "at_or_above_threshold" if r.max_allele >= t else "below_threshold"
        else:
            cls = "no_local_rule"  # assayed but no panel rule: unassessed, never "normal"
        D["str"].append({"locus": r.locus, "unit": r.unit, "gt": str(r.genotype), "max": int(r.max_allele),
                         "thr": t, "ci": _ci.get(r.locus), "cov": round(_cov[r.locus], 1) if r.locus in _cov else None,
                         "class": cls})
    # Scalars for the copy: thresholds differ per locus (AFF2 200, ATN1 48, NIPA1 10, ...), so the
    # panel normalises each bar to its own threshold. Emit the counts so the heading can state the
    # real tally instead of a fixed "not run" claim.
    D["n_str_loci"] = len(D["str"])
    D["n_str_ruled"] = sum(1 for r in D["str"] if r.get("thr"))
    D["n_str_over"] = sum(1 for r in D["str"] if r.get("thr") and r["max"] >= r["thr"])
    D["n_str_norule"] = D["n_str_loci"] - D["n_str_ruled"]
else:
    D["str"] = []
# --- ROH (bcftools roh output; no step of this repository produces it)
_roh_f = W/"09_misc/roh_1mb_nocen.bed"
roh = [l.split() for l in open(_roh_f)] if _roh_f.exists() else []
D["roh"] = [{"chrom": r[0], "start": int(r[1]), "end": int(r[2]), "mb": round(int(r[3])/1e6, 2), "q": float(r[4])} for r in roh]
D["roh_stats"] = {"n": len(roh), "total_mb": round(sum(int(r[3]) for r in roh)/1e6, 1), "max_mb": round(max(int(r[3]) for r in roh)/1e6, 2), "n_gt5": sum(1 for r in roh if int(r[3]) > 5_000_000)} if roh else {"n": 0, "total_mb": 0.0, "max_mb": 0.0, "n_gt5": 0}
# Versions are probed from the tools and references actually in use instead of being hand-written:
# the previous constants had already drifted (Delly 1.7.2 vs installed v2.6.0, bcftools 1.22 vs 1.24,
# ClinVar date 2026-09-05 vs the file's own 2026-09-23).
def _ver(cmd, pat=r"\d+\.\d+(?:\.\d+)?"):
    try:
        o = subprocess.run(cmd if isinstance(cmd, list) else [cmd], capture_output=True, text=True, timeout=15)
        m = re.search(pat, (o.stdout or "") + (o.stderr or ""))
        return m.group(0) if m else None
    except Exception:
        return None
def _clinvar_file_date():
    f = P/"data/ref/clinvar_grch37.vcf.gz"
    if not f.exists(): return None
    try:
        for l in subprocess.run(["bcftools", "view", "-h", str(f)], capture_output=True, text=True, timeout=30).stdout.splitlines():
            if l.startswith("##fileDate="): return l.split("=", 1)[1].strip()
    except Exception:
        pass
    return None
_yt = P/"data/ref/ytree/current_version.txt"
_delly_bin = TOOLS/"delly" if (TOOLS/"delly").exists() else TOOLS/"env/bin/delly"
_eh_bins = sorted((TOOLS/"eh").glob("*/bin/ExpansionHunter")) if (TOOLS/"eh").exists() else []
D["versions"] = {"yfull": _yt.read_text().strip() if _yt.exists() else None,
                 "phylotree": "17.2",   # the tree 06_mtdna.py classifies against
                 "pharmcat": "3.4.0" if (TOOLS/"pharmcat/pharmcat-3.4.0-all.jar").exists() else None,
                 "clinvar": _clinvar_file_date(),
                 "gnomad": "v2.1.1",    # the dataset 07c_gnomad_lookup.py queries (gnomad_r2_1)
                 "delly": _ver([str(_delly_bin), "--version"]),
                 "eh": _ver([str(_eh_bins[0]), "--version"]) if _eh_bins else None,
                 "t1k": _ver([str(TOOLS/"T1K/run-t1k")]) if (TOOLS/"T1K/run-t1k").exists() else None,
                 "cyrius": None,
                 "bcftools": _ver(["bcftools", "--version"]),
                 "samtools": _ver(["samtools", "--version"]),
                 "beagle": None}
# ---- HLA disease associations, archaic gene families, behaviour scores, candidate genes
import json as _j
hd = W/"19_hla_disease/hla_disease.tsv"
if hd.exists(): D["hla_disease"] = pd.read_csv(hd, sep="\t").fillna("").to_dict("records")
gf = W/"15_archaic/gene_families.tsv"
if gf.exists(): D["archaic_families"] = pd.read_csv(gf, sep="\t").fillna("").to_dict("records")
gs = W/"15_archaic/gene_summary.txt"
if gs.exists(): D["archaic_genes"] = dict(l.split("\t") for l in open(gs).read().strip().split("\n"))
sg = W/"15_archaic/segments_genes.tsv"
if sg.exists():
    t = pd.read_csv(sg, sep="\t", dtype={"chrom": str}).nlargest(8, "n_genes")
    D["archaic_top"] = t[["seg","mb","src","n_genes","genes"]].assign(genes=lambda x: x.genes.str.split(",").str[:14].str.join(", ")).to_dict("records")
bp = W/"20_behaviour/behaviour_prs.tsv"
if bp.exists(): D["behaviour"] = pd.read_csv(bp, sep="\t").fillna("").to_dict("records")
cg = W/"20_behaviour/candidate_genes.tsv"
if cg.exists(): D["candidate"] = pd.read_csv(cg, sep="\t").fillna("").to_dict("records")

# ---- figure data consumed by the report template's JavaScript ----
# The template reads keys the assembly above never wrote (phase, la_*, ho_*, archaic*, spectrum,
# somatic, chip_hotspots, telomere, kir, density, prs_rho), which left the matching figures blank.
# Each key is filled straight from the step outputs and degrades to an empty shape when absent.
def _read_tsv(p, **kw):
    return pd.read_csv(p, sep="\t", **kw) if pathlib.Path(p).exists() else None
# phase (step 14b)
_ph = W/"10_phase/summary.txt"
if _ph.exists():
    D["phase"] = {k: (float(v) if "." in v else int(v)) for k, v in (l.split() for l in open(_ph).read().strip().split("\n"))}
else:
    D["phase"] = {"het": 0, "phased": 0, "blocks": 0, "n50_kb": 0.0, "max_mb": 0.0}
# local ancestry (steps 16/16b/17/17b)
# AN5（§7）：祖源各节的状态来自结构化结果与 manifest 校验。缺 manifest、指纹不符或文件损坏都
# 直接标成 unavailable 并给出原因代码，绝不"读旧路径猜结果"，也不用 0% 假值兜底。
# 分析状态判定放在 ancestry_data 里（可测）：缺 manifest / 指纹不符 / 结果损坏各有原因码。
# 复审 AN0/AN5/H6：准入不能只比 sample_id——同一样本换门槛/换参考/换 prune 集后，旧 manifest
# 照样通过。expected_parameters 传**当前有效配置**能推导出的键，与 manifest.parameters 逐一比对；
# 任一缺失或不等 → stale_result。prune 指纹从当前文件重算（与 09b 写入侧同一算法）。
import ancestry_data as _ad
_analysis_state = lambda d, expected=None, expected_parameters=None: _ad.analysis_state(  # noqa: E731
    d, {"sample_id": SAMPLE, **(expected or {})}, expected_parameters=expected_parameters)

_prune_in = W/"11_aadr/prune.prune.in"
_prune_sha = (hashlib.sha256(_prune_in.read_bytes()).hexdigest()[:12] if _prune_in.exists() else "")

_LA_STATE, _LA_REASON, _LA_DOC = _analysis_state(W/"12_localanc")
_AADR_STATE, _AADR_REASON, _AADR_DOC = _analysis_state(
    W/"11_aadr",
    expected_parameters={
        "min_call_rate_modern": MIN_CR_MODERN, "min_call_rate_ancient": MIN_CR_ANCIENT,
        "min_projection_snps": MIN_PROJECTION_SNPS, "min_group_n": MIN_GROUP_N,
        "ancient_prefix": AADR_ANCIENT_PREFIX, "modern_groups": AADR_MODERN,
        "prune_sha": _prune_sha, "prune_sites": (len(_prune_in.read_text().split()) if _prune_in.exists() else 0),
    })
_g = W/"12_localanc/global.tsv"
# 复审 AN5：旧读入口必须服从状态。LA 被禁用/失败时，目录里残留的 global.tsv 与旧
# local_ancestry.json 都是上一次运行的结果——照读等于让禁用状态失效。只有 state=ok 才碰文件。
if _LA_STATE == "ok":
    if _g.exists():
        D["la_global"] = {a: float(f) for a, f in (l.split() for l in open(_g).read().strip().split("\n"))}
    else:
        _gl = ((_LA_DOC or {}).get("global") or [])
        D["la_global"] = ({str(g["panel_id"]): float(g["value"]) for g in _gl} if _gl else {})
else:
    D["la_global"] = {}
    D["la_state"], D["la_reason"] = _LA_STATE, _LA_REASON
_pc = (_read_tsv(W/"12_localanc/per_chrom.tsv", dtype={"chrom": str})
       if _LA_STATE == "ok" else None)
# 列名来自配置的两个来源面板；缺失就退化为"只带 chrom"，而不是让整个构建 KeyError
_PA, _PB = (OPT.get("la_labels") or ["", ""])[0], (OPT.get("la_labels") or ["", ""])[1]
_pc_cols = ["chrom"] + [c for c in (_PA, _PB) if c and c in getattr(_pc, "columns", [])]
D["la_per_chrom"] = _pc[_pc_cols].to_dict("records") if _pc is not None and len(_pc) else []
_sg = (_read_tsv(W/"12_localanc/segments.tsv", dtype={"chrom": str})
       if _LA_STATE == "ok" else None)
if _sg is not None and len(_sg):
    _sg = _sg.copy(); _sg["hap"] = _sg["hap"].astype(str).str.replace("a", "", regex=False).astype(int)  # a1/a2 -> 1/2 for the painting figures
    D["la_segments"] = _sg[["chrom", "start", "end", "anc", "mb", "hap"]].to_dict("records")
else:
    D["la_segments"] = []
# Chromosomes that contributed no segments at all. FLARE creates its output file up front, so an
# interrupted run leaves a 0-byte VCF and the chromosome silently drops out of the figure while
# per_chrom.tsv still reports its proportions -- on this sample, chr4. Recorded so the figure can say
# so rather than leaving a blank band that reads like the untyped acrocentric short arms.
_lac = {str(s["chrom"]) for s in D["la_segments"]}
D["la_missing"] = [str(i) for i in range(1, 23) if str(i) not in _lac] if D["la_segments"] else []
_cb = (_read_tsv(W/"12_localanc/calibration.tsv") if _LA_STATE == "ok" else None)
D["la_calib"] = _cb.to_dict("records") if _cb is not None and len(_cb) else []
# archaic introgression (step 18)
_as = W/"15_archaic/summary.txt"
if _as.exists():
    D["archaic_summary"] = {k: (float(v) if "." in v else int(v)) for k, v in (l.split() for l in open(_as).read().strip().split("\n"))}
else:
    D["archaic_summary"] = {"segments_tested": 0, "carried": 0, "merged": 0, "span_mb": 0.0, "neanderthal_mb": 0.0, "denisovan_mb": 0.0, "homozygous": 0}
_ac = _read_tsv(W/"15_archaic/segments_carried.bed", header=None, names=["chrom", "start", "end", "source", "hom", "mb"], dtype={"chrom": str})
D["archaic"] = _ac.rename(columns={"source": "src"}).to_dict("records") if _ac is not None and len(_ac) else []
# mutation spectrum (step 20)
_sp = _read_tsv(W/"17_mutspec/spectrum96.tsv")
D["spectrum"] = (_sp.rename(columns={"context": "ctx"})[["sub", "ctx", "n", "n_private", "frac", "frac_private"]].to_dict("records")
                 if _sp is not None and len(_sp) else [])
# somatic signals (step 21)
_sm = _read_tsv(W/"13_somatic/summary.tsv")
D["somatic"] = {r.metric: r.value for r in _sm.itertuples()} if _sm is not None and len(_sm) else {}
_ch = _read_tsv(W/"13_somatic/chip_hotspots.tsv")
D["chip_hotspots"] = {"median_depth": float(_ch.depth.median()) if _ch is not None and len(_ch) else 0.0,
                      "alt_reads": int(_ch.alt_reads.sum()) if _ch is not None and len(_ch) else 0,
                      "n_tested": int(len(_ch)) if _ch is not None else 0,
                      # Non-zero hotspots only, strongest first, so the KPI can say which
                      # site the reads sit on instead of an aggregate "supporting reads".
                      "hotspots": (sorted([{"hotspot": str(r.hotspot), "depth": int(r.depth), "alt_reads": int(r.alt_reads), "vaf": float(r.vaf)}
                                    for r in _ch.itertuples() if r.alt_reads > 0], key=lambda h: -h["alt_reads"]) if _ch is not None else [])}
# telomere (step 22)
_tl = W/"14_telomere/counts.tsv"
if _tl.exists():
    _tv = {k: float(v) for k, v in (l.split() for l in open(_tl).read().strip().split("\n"))}
    D["telomere"] = {"k7": int(_tv.get("tel_reads_k7", 0)), "k10": int(_tv.get("tel_reads_k10", 0)),
                     "k12": int(_tv.get("tel_reads_k12", 0)), "k14": int(_tv.get("tel_reads_k14", 0)),
                     "total_reads": int(_tv.get("total_reads", 0))}
else:
    D["telomere"] = {"k7": 0, "k10": 0, "k12": 0, "k14": 0, "total_reads": 0}
# KIR (step 24; empty table when no T1K input was delivered)
_kr = _read_tsv(W/"16_panels/kir.tsv")
D["kir"] = _kr.to_dict("records") if _kr is not None and len(_kr) else []
# kir_summary.txt (same step) carries the haplotype call and this run's HLA ligand groups;
# without it the figure would have no haplotype line and the old card invented an example one.
_ks_f = W/"16_panels/kir_summary.txt"
D["kir_summary"] = None
if _ks_f.exists():
    _ksd = {}
    for _l in _ks_f.read_text(encoding="utf-8").splitlines():
        if "\t" in _l:
            _k2, _v2 = _l.split("\t", 1)
            _ksd[_k2] = _v2
    D["kir_summary"] = _ksd or None
# PASS variants per 1 Mb bin, for the circos density ring
try:
    _dens = {}
    _q = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\n", str(_pass_vcf)],
                         capture_output=True, text=True, check=True).stdout
    for l in _q.splitlines():
        c, p = l.split("\t"); _dens.setdefault(c, []).append(int(p))
    D["density"] = []
    for c, pos in _dens.items():
        # The circos ring steps i*5e6 over the array -- aggregate 1 Mb bins into 5 Mb sums so
        # the ring's geometry and the bin size agree instead of sampling every fifth bin.
        n = D["chrlen"].get(c, 0) // 1_000_000 + 1
        counts = [0]*n
        for p in pos: counts[min(p // 1_000_000, n-1)] += 1
        D["density"].append({"chrom": c, "counts": [sum(counts[i:i+5]) for i in range(0, n, 5)]})
except Exception:
    D["density"] = []
# Human Origins PCA (step 09b)
# 复审 AN5：ho_* 是旧视图的读入口，同样服从状态——AADR 被禁用/失败时，目录里的
# proj_annotated.tsv/summary.json 是上一次运行的结果，照读会让报告继续展示已禁用的数据。
_pa = _read_tsv(W/"11_aadr/proj_annotated.tsv")
if _AADR_STATE == "ok" and _pa is not None and len(_pa):
    anc = _pa[_pa.kind == "ancient"]
    mod = _pa[_pa.kind != "ancient"]
    D["ho_ancient_pts"] = [{"label": r.label, "pc1": r.PC1_AVG, "pc2": r.PC2_AVG, "date": r.date} for r in anc.itertuples()]
    # The nearest *individual* is a different fact from the nearest *group*: the group table averages
    # its members, so the closest single ancient genome can sit inside a group that ranks far lower.
    # Here the nearest genome is BaiyangcunM13.SG (d=0.0032, rank 1 of 405) while its group China_MLBA
    # averages 0.0495 -- the caption quotes the individual, so emit the individual list as well.
    D["ho_near_individual"] = [{"iid": str(r.iid), "group": str(r.label), "d": float(r.d), "date": float(r.date)}
                               for r in anc.nsmallest(12, "d").itertuples()]
    # --- affinity ranking for the ancestry figure (W-T1). A group mean and a single-genome minimum
    # are not the same statistic: the nearest genome (BaiyangcunM13.SG, d=0.0033) belongs to a group
    # averaging 0.0495, 15x further out. Emit one table sorted by d, carrying kind/n/date/region, so
    # the figure can put modern and ancient on a single axis and label which is which.
    # --- the affinity ranking is read from 09b's structured result (AN2). Step 30 assembles; it does
    # not re-group, re-threshold or recompute distances, because a second implementation is exactly how
    # the report and the summary drift apart. `region` and the Chinese name are the one hand-curated
    # layer and come from panel/ancestry_locations.tsv.
    _sum = None
    _sumf = W/"11_aadr/summary.json"
    if _sumf.exists():
        try:
            _sum = json.loads(_sumf.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            _sum = None
    _loc = _read_tsv(pathlib.Path(__file__).resolve().parents[1]/"panel"/"ancestry_locations.tsv",
                     comment="#", keep_default_na=False)
    # 地点表迁移后统一表同时承载 AADR 覆盖行与 1000G 采样地行（dataset 列区分）。这一层是
    # **AADR 分析的**人工区域/中文名/精度 enrichment——语义上只认 AADR 行，必须按 dataset 过滤；
    # 不过滤会让 1000G 行"借道"改 AADR 记录（AADR 面板里 CHB/CHS/CDX 等 308 个 1000G 现代个体
    # 会被 kg 采样行精化成 region、亲和力行吃进 kg 中文名——那是行为变更，不是迁移）。kg 采样点
    # 在 AN5 段按 dataset=1000G 自行取行，与此处互不越界。
    _loc_aadr = None
    if _loc is not None and len(_loc) and {"dataset", "source_id"} <= set(_loc.columns):
        _loc_aadr = _loc[_loc["dataset"].astype(str) == "AADR"]
    _region, _zhname = {}, {}
    if _loc_aadr is not None and len(_loc_aadr):
        for r in _loc_aadr.itertuples():
            sid = str(getattr(r, "source_id", ""))
            note = str(getattr(r, "note", "") or "")
            m = re.search(r"region=(north|south|unclassified)", note)
            _region[sid] = m.group(1) if (m and m.group(1) in ("north", "south")) else ""
            _zhname[sid] = str(getattr(r, "label_zh", "") or "")
    if _sum and (_sum.get("groups") or {}).get("ancient"):
        _recs = {str(r.get("record_id")): r for r in (_sum.get("records") or [])}
        _by_group = {}
        for _kind in ("ancient", "modern"):
            for g in _sum["groups"][_kind]:
                _by_group[str(g["group_id"])] = g
        def _date_of(g):
            ds = [(_recs.get(str(m)) or {}).get("date_mean_bp") for m in g["member_ids"]]
            ds = [d for d in ds if isinstance(d, (int, float))]
            return int(round(sum(ds) / len(ds))) if ds else 0
        def _aff(g, kind):
            lab = str(g["label"])
            return {"label": lab, "kind": kind, "d": round(float(g["distance_mean"]), 5), "n": int(g["n"]),
                    "date_mean": _date_of(g), "region": _region.get(lab, ""), "name_zh": _zhname.get(lab, "")}
        # Default rankings exclude small groups (7.3): they stay available in the detail tables, but a
        # group of one is not a population-level nearest neighbour and would otherwise outrank the real
        # ones merely because a single genome happens to be close.
        _ranked_anc = [g for g in _sum["groups"]["ancient"] if not g.get("small_group")]
        _ranked_mod = [g for g in _sum["groups"]["modern"] if not g.get("small_group")]
        D["ho_affinity"] = sorted(
            [_aff(g, "ancient") for g in _ranked_anc] +
            [_aff(g, "modern") for g in _ranked_mod[:10]], key=lambda r: r["d"])
        def _members(gid, limit=12):
            g = _by_group.get(str(gid))
            if not g:
                return []
            out = []
            for mid in g["member_ids"][:limit]:
                r = _recs.get(str(mid))
                if not r:
                    continue
                out.append({"iid": str(mid), "d": round(float(r.get("distance_to_target") or 0), 5),
                            "call_rate": round(float(r.get("call_rate") or 0), 3),
                            "date": int(r.get("date_mean_bp") or 0)})
            return sorted(out, key=lambda x: x["d"])
        _strip = [g for g in _ranked_anc if g.get("rank")][:3]
        D["ho_affinity_strip"] = [dict(_aff(g, "ancient"), mean_d=round(float(g["distance_mean"]), 5),
                                       members=_members(g["group_id"])) for g in _strip]
        # The closest ancient genome's own group is emitted separately: the caption quotes that genome
        # and its group can rank far down the list, so it would otherwise have nothing to compare with.
        _tg = D["ho_near_individual"][0]["group"] if D["ho_near_individual"] else None
        # counts 供事实文案使用（记录数/合格数/有坐标数），来自 09b 的结构化结果
        D["ho_accession"] = {"counts": dict(_sum.get("counts") or {}),
                             "reference_release": _sum.get("reference_release") or ""}
        # AN6：记录级精度按人工表细化。AADR 只区分"有坐标/没坐标"，而 panel/ancestry_locations.tsv
        # 明确标了哪些来源只到地区级（province）——省级来源不该在地图上显示成一个精确遗址。
        _prec = {}
        if _loc_aadr is not None and len(_loc_aadr) and {"source_id", "precision"} <= set(_loc_aadr.columns):
            for _sid, _pr in zip(_loc_aadr["source_id"], _loc_aadr["precision"]):
                if str(_pr).strip():
                    _prec[str(_sid).strip()] = str(_pr).strip()
        # 注意：这里必须改 **ancestry 实际使用的那个对象**。30 里 AADR 的 summary.json 被读了两遍
        # （AN2 段的 _sum 与 AN5 段的 _AADR_DOC），是两个独立 dict；改错一个，地图上就看不到变化。
        _recs_meta = (_AADR_DOC or {}).get("records") or (_sum.get("records") or [])
        _hit = 0
        for _r in _recs_meta:
            _p = _prec.get(str(_r.get("source_population_id") or ""))
            if _p:
                _r["location_precision"] = _p
                _hit += 1
        if _hit:
            print(f"30: location_precision refined from the location table for {_hit} "
                  f"of {len(_recs_meta)} records", file=sys.stderr)
        _tgg = _by_group.get(str(_tg)) if _tg else None
        D["ho_target_group"] = (dict(_aff(_tgg, "ancient"), mean_d=round(float(_tgg["distance_mean"]), 5),
                                     in_strip=any(g["group_id"] == _tgg["group_id"] for g in _strip),
                                     members=_members(_tgg["group_id"]))
                                if _tgg is not None else None)
    else:
        D["ho_affinity"] = []; D["ho_affinity_strip"] = []; D["ho_target_group"] = None
    # The template indexes ho_modern rows positionally (p[0]/p[1]/p[2]) -- arrays, not dicts.
    D["ho_modern"] = [[r.label, r.PC1_AVG, r.PC2_AVG] for r in mod.itertuples()]
    _me = _pa[_pa.kind == "target"] if "target" in set(_pa.kind) else _pa.tail(1)
    D["ho_me"] = [float(_me.iloc[0].PC1_AVG), float(_me.iloc[0].PC2_AVG)]
    D["ho_prov"] = [{"label": r.label.split("_")[1] if r.label.startswith("Han_") else r.label,
                      "pc1": r.PC1_AVG, "pc2": r.PC2_AVG, "date": r.date}
                     for r in mod[mod.label.str.startswith("Han_")].itertuples()]
else:
    D["ho_ancient_pts"] = []; D["ho_modern"] = []; D["ho_me"] = [0.0, 0.0]; D["ho_prov"] = []
    D["ho_affinity"] = []; D["ho_affinity_strip"] = []; D["ho_target_group"] = None
_na = _read_tsv(W/"11_aadr/near_ancient.tsv")
D["ho_near_ancient"] = (_na.to_dict("records")
                         if _AADR_STATE == "ok" and _na is not None and len(_na) else [])
# --- f3 statistics (28_f3_stats.py). Optional module: null when the step was not run,
# # and the report's f3 card hides itself on a null -- the sections table records why.
# 复审 H6/AN5 残留：f3 此前按"文件存在"进场，绕过 28 的 manifest——禁用/失败/换样本的旧
# f3_stats.json 照样投递。走 analysis_state 准入（state+sample_id），不合格即空态。
_F3_STATE, _F3_REASON, _F3_DOC = _ad.analysis_state(
    W/"04_ancestry/f3", {"sample_id": SAMPLE}, names=("f3_stats.json",))
if _F3_STATE == "ok":
    D["f3"] = _F3_DOC
else:
    print(f"30: f3 result not admitted (state={_F3_STATE}, reason={_F3_REASON}); "
          "the f3 panel is skipped", file=sys.stderr)
    D["f3"] = None
# --- AN5（§7 报告契约）：把各分析的结构化结果组装成 D.ancestry / D.lineages。
def _attach_kg_group_locations(analyses, loc_rows):
    """AN6 增强（1000G 位置）：给 1000G 分析的**人群组**补采样点坐标。

    坐标来自 panel/ancestry_locations.tsv 的 dataset=1000G 行（IGSR phase3 人群描述 + 采样地，
    region 级；地点表迁移后是唯一人工地理覆盖表的一部分）：
    一个人群一个点，代表参考样本的采样/来源地，不是任何个体出生地，也不是目标样本的位置。
    组记录（个体）不落坐标——IGSR 不发布个体地理信息，给个体编点就是造数据。只补、不覆盖：
    组上已有的坐标（将来 04b 若自带）不动；面板里没有的人群保持无坐标，由
    group_location_counts 如实计数，不假充全定位。返回定位的组数，供日志核对。
    """
    by_pop = {}
    for r in loc_rows or []:
        pop = str(r.get("source_id", "")).strip()
        try:
            lat, lon = float(r.get("latitude")), float(r.get("longitude"))
        except (TypeError, ValueError):
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180 and pop):
            continue
        by_pop[pop] = (r, lat, lon)
    n_located = 0
    for a in analyses:
        if str(a.get("dataset", "")) != "1000G":
            continue
        located = unlocated = 0
        for g in a.get("groups") or []:
            hit = by_pop.get(str(g.get("group_id") or g.get("label") or "").strip())
            if hit is None:
                unlocated += 1
                continue
            row, lat, lon = hit
            if g.get("latitude") is None and g.get("longitude") is None:
                g["latitude"], g["longitude"] = lat, lon
                g["location_id"] = str(row.get("location_id") or f"KG:{row.get('source_id')}")
                g["locality"] = str(row.get("locality") or "")
                g["location_precision"] = str(row.get("precision") or "region")
                _zh = str(row.get("label_zh") or "").strip()
                if _zh:
                    g["name_zh"] = _zh
            located += 1
        a["group_location_counts"] = {"located": located, "unlocated": unlocated}
        n_located += located
    return n_located

# 这是模板与 AN6 要消费的形状；旧键（ho_*/near_eas/…）只作为尚未迁移的视图的过渡，不再各自算一套。
_kg_state, _kg_reason, _kg_doc = _analysis_state(
    W/"04_ancestry",
    expected_parameters={"regional_enabled": REGIONAL_ENABLED, "superpop": SUPERPOP,
                        "subpops": SUBPOPS})
_AADR_DOC = _AADR_DOC if isinstance(_AADR_DOC, dict) else None
_anc_analyses = []
for _doc, _st, _rs in ((_AADR_DOC, _AADR_STATE, _AADR_REASON), (_kg_doc, _kg_state, _kg_reason)):
    if not isinstance(_doc, dict):
        continue
    if "analyses" in _doc:                      # 04b 的 summary.json：多个参考空间
        for _a in _doc["analyses"]:
            _anc_analyses.append(dict(_a, state=_st or _a.get("state"), reason_code=_rs or _a.get("reason_code", "")))
    else:                                       # 09b 的 summary.json：单个 AADR 空间
        _anc_analyses.append({k: _doc.get(k) for k in
                              ("analysis_id", "dataset", "reference_release", "scope", "components", "metric",
                               "thresholds", "counts", "target", "records", "groups", "sources")}
                             | {"state": _st, "reason_code": _rs})
# AN6 增强（§3.2 P1"1000G 没有 records"）：1000G 参考空间的人群组补采样点（region 级）。
# 04b 不动（manifest/参数不掺显示层）；30 组装时附加，与 AADR 的 location_precision 同一模式。
# 地点表迁移：kg 采样点与 AADR 覆盖同住 panel/ancestry_locations.tsv（唯一人工地理覆盖表，
# dataset 列区分）——不能依赖上文 AADR 门控块里的 _loc（AADR 禁用时它不会被赋值），这里按
# dataset=1000G 自行过滤同一张表。
_kgl_loc = _read_tsv(pathlib.Path(__file__).resolve().parents[1] / "panel" / "ancestry_locations.tsv",
                      comment="#", keep_default_na=False)
_kgl = None
if _kgl_loc is not None and len(_kgl_loc) and {"dataset", "source_id"} <= set(_kgl_loc.columns):
    _kgl = _kgl_loc[_kgl_loc["dataset"].astype(str) == "1000G"]
_n_kgl = _attach_kg_group_locations(_anc_analyses,
                                    _kgl.to_dict("records") if _kgl is not None else [])
print(f"30: kg group sampling-site locations attached for {_n_kgl} groups", file=sys.stderr)
# 复审 AN6-P1：默认分析原先硬编码挑 AADR。一个只启用 1000G（或禁用了 AADR）的样本，报告仍会宣称
# 默认分析是 AADR——那是个不存在于本次运行里的选择。改为按**数据顺序取第一个 state=ok 的分析**：
# 顺序即配置里的启用顺序，state 决定它这次是否真的产出。全都不可用时留空，由消费端如实处理。
_default = next((a.get("analysis_id") for a in _anc_analyses if a.get("state") == "ok"), "")
if _default:
    print(f"30: default analysis = {_default}", file=sys.stderr)
else:
    print("30: no analysis is in state ok; default_analysis_id left empty", file=sys.stderr)
D["ancestry"] = {"schema_version": 1, "default_analysis_id": _default or
                 (_anc_analyses[0].get("analysis_id") if _anc_analyses else ""),
                 "analyses": _anc_analyses,
                 "local": (_LA_DOC if isinstance(_LA_DOC, dict) else {})}
# near_* 取准入后的分析（stale/missing 一律 None），口径与 D.ancestry 同源，不二次解析文件
D["near_global"], D["knn_global"] = _near_from(_anc_analyses, "kg-global")
D["near_eas"], D["knn_eas"] = _near_from(_anc_analyses, "kg-regional")
# 父母系：以 05/06 的结构化结果 + 04 的历史视图为准（lineage_history.json）。
# 复审 H6/AN5 残留：此前按"文件存在"读取，绕过 09d 的状态——现在 09d 会写 manifest
# （analysis_id=09d-lineage-history），30 走 analysis_state 准入；不合格交付空态并说明原因。
_LH_STATE, _LH_REASON, _LH_DOC = _ad.analysis_state(
    W/"03_haplo", {"sample_id": SAMPLE, "analysis_id": "09d-lineage-history"},
    names=("lineage_history.json",))
if _LH_STATE != "ok":
    print(f"30: lineage history not admitted (state={_LH_STATE}, reason={_LH_REASON}); "
          "D.lineages delivered empty", file=sys.stderr)
D["lineages"] = {"y": (_LH_DOC or {}).get("y") if _LH_STATE == "ok" else None,
                 "mt": (_LH_DOC or {}).get("mt") if _LH_STATE == "ok" else None}

# PRS site coverage, for the copy's transferability note
_pr = D.get("prs") or []
_cov = [r.get("coverage_pct") for r in _pr if isinstance(r.get("coverage_pct"), (int, float)) and r["coverage_pct"] >= 0]
D["prs_rho"] = {"subset": f"{min(_cov):.0f}–{max(_cov):.0f}%" if _cov else "-", "imputed": "-"}
# --- section status table (H1). The arrays and objects above stay type-stable; this table records,
# per section, whether it is ok, a genuine negative, or unavailable -- with a machine-readable
# reason code and a bilingual reason, so a report can say why instead of showing an empty box.
# "unavailable" is never a negative result, and it is never used for sections that have no
# positive/negative meaning in the first place (QC, ancestry) beyond stating the missing input.
_SEC = []
# Bilingual display names so the JS status note can name a section without hardcoding them again.
_SEC_NAMES = {
    "chip": ("芯片—WGS 一致性", "chip vs WGS consistency"), "pgx_pharmcat": ("PharmCAT 药物基因组", "PharmCAT pharmacogenomics"),
    "cyp2d6": ("CYP2D6（Cyrius）", "CYP2D6 (Cyrius)"), "hla": ("HLA 分型（T1K）", "HLA typing (T1K)"),
    "kir": ("KIR 与 HLA 配体", "KIR and HLA ligands"), "sv": ("结构变异（Delly）", "structural variants (Delly)"),
    "smn": ("SMN 拷贝数", "SMN copy number"), "str": ("短串联重复（EH）", "short tandem repeats (EH)"),
    "roh": ("纯合区段（ROH）", "runs of homozygosity"), "chip_hotspots": ("CHIP 热点", "CHIP hotspots"),
    "somatic": ("体细胞信号", "somatic signals"), "telomere": ("端粒", "telomere"),
    "phase": ("相位统计", "phasing statistics"), "spectrum": ("突变谱", "mutation spectrum"),
    "density": ("全基因组密度", "genome-wide density"), "prs": ("多基因评分", "polygenic scores"),
    "behaviour": ("行为特征评分", "behavioural scores"), "candidate": ("候选基因位点", "candidate-gene variants"),
    "ancestry": ("祖源分析", "ancestry"), "local_ancestry": ("局部祖源", "local ancestry"),
    "f3": ("f3 统计", "f3 statistics"),
    "archaic": ("古人类渗入", "archaic introgression"), "aadr": ("古 DNA 投影", "ancient DNA projection"),
    "clinvar": ("ClinVar 命中", "ClinVar hits"), "lof": ("功能丧失变异", "loss-of-function variants"),
}
def _sec(sid, status, code=None, zh=None, en=None, ev=None, detail=None):
    nm = _SEC_NAMES.get(sid, (sid, sid))
    _SEC.append({"id": sid, "status": status, "name_zh": nm[0], "name_en": nm[1], "code": code,
                 "reason_zh": zh, "reason_en": en, "evidence": ev, "detail": detail})
def _has(key):
    v = D.get(key)
    return bool(v) if isinstance(v, (list, dict)) else v is not None
def _need(sid, ok, code, zh, en, ev, detail=None):
    _sec(sid, "ok" if ok else "unavailable", None if ok else code, None if ok else zh,
         None if ok else en, ev, detail)
_MISS_IN   = ("no_input",       "交付不含该模块所需输入",        "this section's input was not part of the delivery")
_MISS_OUT  = ("missing_output", "未找到该模块的输出（未运行或失败，二者在此不可区分）", "no output found for this section (not run or failed; the two are indistinguishable here)")
_MISS_STEP = ("step_not_run",   "本流程没有产出该模块的步骤",    "no step of this pipeline produces this section")

_need("chip", _has("chip") and D["chip"].get("compared", 0) > 0, *_MISS_IN, ev="01_qc/chip_vs_wgs_summary.tsv")
_need("pgx_pharmcat", _has("pgx_pharmcat"), *_MISS_STEP, ev="06_pgx/pharmcat/pharmcat_summary.tsv", detail=len(D.get("pgx_pharmcat") or []))
_need("cyp2d6", D.get("cyp2d6") not in (None, "-"), *_MISS_STEP, ev="06_pgx/cyrius/target.tsv", detail=D.get("cyp2d6"))
_need("hla", _has("hla"), *_MISS_STEP, ev="06_pgx/t1k/dayu_genotype.tsv", detail=len(D.get("hla") or {}))
_need("kir", _has("kir"), *_MISS_STEP, ev="16_panels/kir.tsv", detail=len(D.get("kir") or []))
# A zero SV table is a genuine negative only when Delly actually produced a VCF; without the
# bcf the empty table is step 13's not-run fallback. Either way "unavailable" is wrong for a
# real zero, and "ok" is wrong for the fallback.
_need("sv", _sv_f.exists() and (D.get("sv_total", 0) > 0 or (W/"08_sv/delly/target.sv.bcf").exists()),
      *_MISS_STEP, ev="08_sv/sv_filtered.tsv", detail=D.get("sv_total"))
_need("smn", (D.get("smn") or {}).get("SMN1") is not None, *_MISS_OUT, ev="08_sv/smn/target.tsv", detail=D.get("smn"))
_need("str", _has("str"), *_MISS_OUT, ev="08_sv/eh/eh_summary.tsv", detail=len(D.get("str") or []))
# The ROH bed is produced outside this pipeline; an empty bed from a real run is a negative
# result, so the section is ok whenever the file exists at all.
_need("roh", (W/"09_misc/roh_1mb_nocen.bed").exists(), *_MISS_OUT, ev="09_misc/roh_1mb_nocen.bed", detail=D.get("roh_stats"))
_need("chip_hotspots", _has("chip_hotspots") and D["chip_hotspots"].get("median_depth", 0) > 0, *_MISS_OUT, ev="13_somatic/chip_hotspots.tsv", detail=D.get("chip_hotspots"))
_need("somatic", _has("somatic"), *_MISS_OUT, ev="13_somatic/summary.tsv", detail=D.get("somatic"))
# The 100k-read gate is a working floor for a usable TelSeq-style estimate, not a validated
# threshold; the reason says so and notes that an absolute-length conversion would still
# need GC and read-length bias correction (22_telomere.sh estimate.tsv).
_tl_k7 = (D.get("telomere") or {}).get("k7", 0)
_need("telomere", _has("telomere") and _tl_k7 >= 100000,
      "insufficient_reads",
      f"端粒重复读段过少（k7={_tl_k7}，门槛 100000 为工作设定而非经验证阈值），无法可靠估计长度；绝对长度换算还需 GC 与读长偏差校正",
      f"too few telomeric repeat reads for a reliable length estimate (k7={_tl_k7}; the 100000 gate is a working floor, not a validated threshold); absolute lengths would further need GC and read-length bias correction",
      ev="14_telomere/counts.tsv", detail=_tl_k7)
_need("phase", _has("phase"), *_MISS_OUT, ev="10_phase/summary.txt", detail=D.get("phase"))
_need("spectrum", _has("spectrum"), *_MISS_OUT, ev="17_mutspec/spectrum96.tsv", detail=len(D.get("spectrum") or []))
_need("density", _has("density"), *_MISS_OUT, ev="00_input/target.pass.vcf.gz")
_need("prs", _has("prs"), *_MISS_OUT, ev="07_prs/prs_wgs.tsv", detail=len(D.get("prs") or []))
_need("behaviour", _has("behaviour"), *_MISS_OUT, ev="20_behaviour/", detail=len(D.get("behaviour") or []))
_need("candidate", _has("candidate"), *_MISS_OUT, ev="20_behaviour/candidate_genes.tsv", detail=len(D.get("candidate") or []))
_need("ancestry", _has("pca_global") and _has("near_global"), *_MISS_OUT, ev="04_ancestry/")
# 祖源节按分析状态报告（§7 的原因码）：禁用、参考不足、结果缺失彼此可区分。
_LA_CODES = {
    "disabled_by_config": ("配置未启用局部祖源（未显式给出两个来源面板）", "local ancestry is disabled by the configuration"),
    "missing_manifest":   ("找不到 manifest：旧结果需重建，不能直接复用", "no manifest: the result predates the contract and must be rebuilt"),
    "stale_result":       ("结果的 manifest 与当前样本/配置/参考指纹不一致，需重算", "the result's manifest does not match this sample/configuration; rebuild it"),
    "missing_result":     ("未找到该分析的结构化结果", "no structured result for this analysis"),
    "unreadable_result":  ("结构化结果损坏，需重跑该步", "the structured result is unreadable; rerun the step"),
    "no_segments":        ("该分析没有产出任何片段", "the analysis produced no segments"),
}
_sec("local_ancestry", "ok" if _LA_STATE == "ok" else "unavailable",
     None if _LA_STATE == "ok" else _LA_REASON,
     None if _LA_STATE == "ok" else _LA_CODES.get(_LA_REASON, ("未启用或不可用", "not available"))[0],
     None if _LA_STATE == "ok" else _LA_CODES.get(_LA_REASON, ("未启用或不可用", "not available"))[1],
     ev="12_localanc/local_ancestry.json", detail=f"state={_LA_STATE}")
_need("archaic", _has("archaic"), *_MISS_OUT, ev="15_archaic/segments_all.tsv")
_need("f3", _has("f3"), *_MISS_OUT, ev="04_ancestry/f3/f3_stats.json")
_sec("aadr", "ok" if (_AADR_STATE == "ok" and _has("ho_modern")) else "unavailable",
     None if _AADR_STATE == "ok" else _AADR_REASON,
     None if _AADR_STATE == "ok" else _LA_CODES.get(_AADR_REASON, ("未启用或不可用", "not available"))[0],
     None if _AADR_STATE == "ok" else _LA_CODES.get(_AADR_REASON, ("未启用或不可用", "not available"))[1],
     ev="11_aadr/summary.json", detail=f"state={_AADR_STATE}")
_need("clinvar", D.get("clinvar_total", 0) > 0, *_MISS_OUT, ev="05_clinvar/clinvar_all_hits.tsv", detail=D.get("clinvar_total"))
_need("lof", _has("lof"), *_MISS_OUT, ev="05_clinvar/lof_table.tsv", detail=D.get("lof"))
D["sections"] = _SEC
# --- end section status table
D["name_zh"] = NAME_ZH; D["name_en"] = NAME_EN; D["sample"] = SAMPLE
D["title_zh"] = NAME_ZH; D["title_en"] = NAME_EN
json.dump(D, open(W/"report_data.json", "w"), ensure_ascii=False)
print("wrote", W/"report_data.json", (W/"report_data.json").stat().st_size//1024, "KB"); print({k: (len(v) if isinstance(v, (list, dict)) else v) for k, v in D.items() if k not in ("chrlen", "cen")})
