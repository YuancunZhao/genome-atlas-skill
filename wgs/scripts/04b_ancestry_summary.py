#!/usr/bin/env python
"""Nearest 1000 Genomes reference groups in the global and (if configured) the regional PCA (AN2).

Distance convention is the same as 09b: each reference individual's distance to the target is computed
first, and a group's value is the mean over its members. The previous version compared group centroids
with the target, so the 1000G and AADR summaries were not answering the same question; that also made
the ranking change for reasons a reader could not see. The regional analysis runs only when the
configuration names a super-population -- the module default is never treated as a choice, and the
output files are named after the actual scope instead of the hard-coded "eas".

Writes 04_ancestry/summary.json (the Analysis shape of section 7.3) and keeps summary.txt as the audit
log.
"""
import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...
from wgsconfig import (SAMPLE, NAME_EN, MIN_CR_TARGET, MIN_CR_MODERN, MIN_PROJECTION_SNPS,
                       MIN_GROUP_N, SUBPOPS, SUPERPOP, REGIONAL_ENABLED, OPT)

import pandas as pd, numpy as np
import ancestry_data as ad

W = f"{P}/wgs/04_ancestry"
# 失效先于计算（复审 §3.2 P0 失败生命周期）：中途崩溃/退出 1 时，上一轮的 ok manifest 与
# summary.json 不得继续充当本次结果；成功路径的最终 manifest 原子覆盖这份记录。
ad.begin_run_manifest(f"{W}/manifest.json", SAMPLE, "04b-ancestry-summary")
log = []


def say(line):
    print(line)
    log.append(line)


# "regional" is the scope name in the configuration; older runs wrote the same products as "eas.*".
# Reading the configured name first keeps a stale file from being mistaken for the current scope.
def _first_existing(*names):
    for n in names:
        if (pathlib.Path(W) / n).exists():
            return n
    return None


spaces = [("global", "global", "kg.proj.sscore", "target.proj.sscore", "prune.prune.in")]
if REGIONAL_ENABLED:
    reg = _first_existing("regional.proj.sscore", "eas.proj.sscore")
    reg_t = _first_existing("regional.target.proj.sscore", "eas.target.proj.sscore")
    # 旧 eas.* 是硬编码东亚范围的产物：只有当前配置的区域**就是**东亚（SUPERPOP 与 psam 的
    # SuperPop 列同名，即 "EAS"）时它们才是本配置的结果。别的区域配置读 eas 文件，等于把
    # 错面板当区域分析发布——宁可缺区域，不可错区域。
    _legacy_eas = (reg or "").startswith("eas") or (reg_t or "").startswith("eas")
    if _legacy_eas and SUPERPOP.strip().upper() != "EAS":
        say(f"legacy eas.* products do not belong to the configured regional scope ({SUPERPOP}); "
            f"run 04 with the current scope name (regional.*) before reading them")
        reg = reg_t = None
    if reg and reg_t:
        _tag = reg.split(".")[0]                      # regional 或 eas：与 prune.<tag>.prune.in 对应
        spaces.append(("regional", _tag, reg, reg_t, f"prune.{_tag}.prune.in"))
    else:
        say(f"regional PCA is configured (ref_superpop={SUPERPOP}) but its products are missing; "
            f"run 04 with the current scope name before reading them")


def _plink_missing(pfile_prefix, prune_in, out_prefix):
    import subprocess
    r = subprocess.run([PLINK2, "--pfile", str(pfile_prefix), "--extract", str(prune_in),
                        "--missing", "--out", str(out_prefix),
                        "--threads", str(THREADS), "--memory", str(int(float(MEM_GB) * 1000))],
                       capture_output=True, text=True)
    return r.returncode == 0 and pathlib.Path(str(out_prefix) + ".smiss").exists()


def _space_coverage(tag, prune_name):
    """本空间在**最终剪枝位点**上的每个体覆盖：参考面板与目标分开算（04 的 PCA 正是这样用的），
    每次运行重算并写 qc.<tag>.{ref,target}.smiss（不做存在性缓存，见调用处注释）。缺任一输入
    （联合 pgen、剪枝集、目标 1kg）或 plink 失败时返回 None——覆盖算不出来就不过门槛，不回退成
    "人人合格"。"""
    w = pathlib.Path(W)
    prune_in = w / prune_name
    refs_out, tgt_out = w / f"qc.{tag}.ref.smiss", w / f"qc.{tag}.target.smiss"
    tgt_prefix = pathlib.Path(P) / "wgs" / "02_complete" / f"{SAMPLE}.1kg"
    refs_prefix = w / "kg.common"
    if not prune_in.exists():
        say(f"[{tag}] skipped: {prune_name} not found -- coverage (and therefore the gate) "
            f"cannot be computed on the final sites")
        return None
    if not all(pathlib.Path(f"{refs_prefix}.{e}").exists() for e in ("pgen", "psam", "pvar")):
        say(f"[{tag}] skipped: kg.common.* not found -- reference coverage cannot be computed")
        return None
    if not all(pathlib.Path(f"{tgt_prefix}.{e}").exists() for e in ("pgen", "psam", "pvar")):
        say(f"[{tag}] skipped: 02_complete/{SAMPLE}.1kg.* not found -- target coverage cannot be computed")
        return None
    # 复审 AN2/H6：每次重算，不做存在性缓存——qc.<tag>.*.smiss 是哪个 prune 集/哪个目标算出来的，
    # 事后无法从文件得知；换目标或换位点集后命中旧缓存，门槛就成了摆设。--missing 数秒即完。
    if not _plink_missing(refs_prefix, prune_in, w / f"qc.{tag}.ref"):
        say(f"[{tag}] skipped: reference missingness failed (plink2)")
        return None
    if not _plink_missing(tgt_prefix, prune_in, w / f"qc.{tag}.target"):
        say(f"[{tag}] skipped: target missingness failed (plink2)")
        return None
    cov = pd.concat([
        ad.coverage_from_smiss(pd.read_csv(refs_out, sep=r"\s+")),
        ad.coverage_from_smiss(pd.read_csv(tgt_out, sep=r"\s+")),
    ], ignore_index=True)
    return cov[["iid", "n_called_snps", "call_rate"]]


analyses = []
skip_reason = "no_reference_space"          # 一个空间都没有时的具体原因，写进 summary 与 manifest
for scope, tag, kg_file, target_file, prune_name in spaces:
    if not (pathlib.Path(W) / kg_file).exists() or not (pathlib.Path(W) / target_file).exists():
        say(f"[{tag}] skipped: {kg_file} / {target_file} not found")
        continue
    k = pd.read_csv(f"{W}/{kg_file}", sep="\t").rename(columns={"#IID": "iid", "IID": "iid"})
    d = pd.read_csv(f"{W}/{target_file}", sep="\t")
    pcs = [c for c in k.columns if c.endswith("_AVG")][:4]
    if len(pcs) < 4:
        say(f"[{tag}] skipped: only {len(pcs)} components in {kg_file}")
        continue
    D = d[pcs].values[0]
    k["distance_to_target"] = np.sqrt(((k[pcs].values - D) ** 2).sum(1))
    k["kind"] = "modern"
    k["group_id"] = k["Population"]
    # 覆盖来自最终剪枝位点上的真实 smiss（不再是 NaN+关门）：1000G 没有例外，参考与目标同一口径。
    cov = _space_coverage(tag, prune_name)
    if cov is None:
        skip_reason = "missing_qc"
        continue
    k = k.drop(columns=[c for c in ("n_called_snps", "call_rate") if c in k.columns]) \
         .merge(cov, on="iid", how="left")
    k["call_rate"] = k["call_rate"].fillna(-1.0)      # 不在 smiss 里 = 覆盖未知，必不过门槛
    _t = cov[cov["iid"] == SAMPLE]
    _t_rate = float(_t["call_rate"].iloc[0]) if len(_t) else -1.0
    _t_n = float(_t["n_called_snps"].iloc[0]) if len(_t) else -1.0
    if _t_rate < float(MIN_CR_TARGET) or _t_n < float(MIN_PROJECTION_SNPS):
        # 目标不过门槛：本空间不出排名（HANDOFF AN2 与 09b 同一条验收）。
        say(f"[{tag}] skipped: target call_rate={_t_rate:.4f} (min {MIN_CR_TARGET}), "
            f"n_called={_t_n:.0f} (min {MIN_PROJECTION_SNPS}) on the final pruned sites")
        skip_reason = "target_below_coverage_gate"
        continue
    recs = k.to_dict("records")
    groups = ad.group_summaries(
        ad.eligible_records(recs, kind="modern", min_rate=MIN_CR_MODERN, min_snps=MIN_PROJECTION_SNPS),
        min_group_n=MIN_GROUP_N)
    _elig = sum(1 for r in recs if r.get("eligible") is not False)

    say(f"== {tag} (PC1-4, mean individual distance): {NAME_EN} {np.round(D, 4)}")
    top = groups[:6]
    for g in top:
        say(f"   {g['label']:<12} n={g['n']:<4} d={g['distance_mean']:.4f}")
    nearest15 = k.nsmallest(15, "distance_to_target").Population.value_counts().to_dict()
    say(f"   nearest 15 individuals: {nearest15}")
    # 复审 P1（AN0/AN2/AN5）：30 的 near_*/knn_* 不再另算一套质心排名，从这里转录同一口径——
    # 个体距离均值（groups）与前 15 近个体的人群计数（nearest_individuals）。
    # Percentile position of the target inside the configured sub-populations (not a fixed list).
    for p in [x for x in SUBPOPS if x in set(k.Population)]:
        sub = k[k.Population == p]
        say(f"   {p}: PC1 pct {(sub[pcs[0]] < D[0]).mean() * 100:.0f}, "
            f"PC2 pct {(sub[pcs[1]] < D[1]).mean() * 100:.0f}")

    analyses.append({
        "analysis_id": f"kg-{tag}", "dataset": "1000G", "reference_release": "phase3",
        "scope": scope, "state": "ok", "reason_code": "", "components": [1, 2, 3, 4],
        "metric": "mean_individual_pc_distance",
        "thresholds": {"min_group_n": MIN_GROUP_N, "subpops": SUBPOPS, "superpop": SUPERPOP if scope != "global" else "",
                       "min_call_rate_modern": MIN_CR_MODERN, "min_call_rate_target": MIN_CR_TARGET,
                       "min_projection_snps": MIN_PROJECTION_SNPS},
        "counts": {"selected": len(recs), "eligible": _elig, "excluded": len(recs) - _elig,
                   "mapped": 0, "unmapped": len(recs)},
        "target": {"sample_id": SAMPLE, "pcs": [float(x) for x in D],
                   "call_rate": _t_rate, "n_called_snps": _t_n},
        "groups": groups,
        "nearest_individuals": nearest15,
        "sources": [{"scores": kg_file, "target_scores": target_file,
                     "coverage": f"qc.{tag}.ref.smiss+qc.{tag}.target.smiss over {prune_name}"},
        ],
    })

pathlib.Path(f"{W}/summary.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
state = "ok" if analyses else "unavailable"
pathlib.Path(f"{W}/summary.json").write_text(json.dumps({
    "schema_version": 1, "default_analysis_id": analyses[0]["analysis_id"] if analyses else "",
    "analyses": analyses, "local": {},
    "state": state, "reason_code": "" if analyses else skip_reason,
}, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
ad.write_manifest(f"{W}/manifest.json", ad.build_manifest(
    SAMPLE, "04b-ancestry-summary", state=state,
    reason_code="" if analyses else skip_reason,
    parameters={"regional_enabled": REGIONAL_ENABLED, "superpop": SUPERPOP, "subpops": SUBPOPS},
    outputs=["04_ancestry/summary.json", "04_ancestry/summary.txt"]))
print(f"wrote {W}/summary.json ({len(analyses)} reference space(s)); manifest written")
