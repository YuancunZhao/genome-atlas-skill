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
from wgsconfig import (SAMPLE, NAME_EN, MIN_CR_TARGET, MIN_GROUP_N, SUBPOPS, SUPERPOP,
                       REGIONAL_ENABLED, OPT)

import pandas as pd, numpy as np
import ancestry_data as ad

W = f"{P}/wgs/04_ancestry"
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


spaces = [("global", "global", "kg.proj.sscore", "target.proj.sscore")]
if REGIONAL_ENABLED:
    reg = _first_existing("regional.proj.sscore", "eas.proj.sscore")
    reg_t = _first_existing(f"regional.target.proj.sscore", "eas.target.proj.sscore")
    if reg and reg_t:
        spaces.append(("regional", "regional", reg, reg_t))
    else:
        say(f"regional PCA is configured (ref_superpop={SUPERPOP}) but its products are missing; "
            f"run 04 with the current scope name before reading them")

analyses = []
for scope, tag, kg_file, target_file in spaces:
    if not (pathlib.Path(W) / kg_file).exists() or not (pathlib.Path(W) / target_file).exists():
        say(f"[{tag}] skipped: {kg_file} / {target_file} not found")
        continue
    k = pd.read_csv(f"{W}/{kg_file}", sep="\t").rename(columns={"#IID": "iid"})
    d = pd.read_csv(f"{W}/{target_file}", sep="\t")
    pcs = [c for c in k.columns if c.endswith("_AVG")][:4]
    if len(pcs) < 4:
        say(f"[{tag}] skipped: only {len(pcs)} components in {kg_file}")
        continue
    D = d[pcs].values[0]
    k["distance_to_target"] = np.sqrt(((k[pcs].values - D) ** 2).sum(1))
    k["kind"] = "modern"
    k["group_id"] = k["Population"]
    k["n_called_snps"] = np.nan          # filled in by the PCA step on the final pruned site set
    k["call_rate"] = np.nan
    recs = k.to_dict("records")
    groups = ad.group_summaries(
        ad.eligible_records(recs, kind="modern", min_rate=None, min_snps=None),
        min_group_n=MIN_GROUP_N)

    say(f"== {tag} (PC1-4, mean individual distance): {NAME_EN} {np.round(D, 4)}")
    top = groups[:6]
    for g in top:
        say(f"   {g['label']:<12} n={g['n']:<4} d={g['distance_mean']:.4f}")
    nearest15 = k.nsmallest(15, "distance_to_target").Population.value_counts().to_dict()
    say(f"   nearest 15 individuals: {nearest15}")
    # Percentile position of the target inside the configured sub-populations (not a fixed list).
    for p in [x for x in SUBPOPS if x in set(k.Population)]:
        sub = k[k.Population == p]
        say(f"   {p}: PC1 pct {(sub[pcs[0]] < D[0]).mean() * 100:.0f}, "
            f"PC2 pct {(sub[pcs[1]] < D[1]).mean() * 100:.0f}")

    analyses.append({
        "analysis_id": f"kg-{tag}", "dataset": "1000G", "reference_release": "phase3",
        "scope": scope, "state": "ok", "reason_code": "", "components": [1, 2, 3, 4],
        "metric": "mean_individual_pc_distance",
        "thresholds": {"min_group_n": MIN_GROUP_N, "subpops": SUBPOPS, "superpop": SUPERPOP if scope != "global" else ""},
        "counts": {"selected": len(recs), "eligible": len(recs), "excluded": 0,
                   "mapped": 0, "unmapped": len(recs)},
        "target": {"sample_id": SAMPLE, "pcs": [float(x) for x in D], "call_rate": None, "n_called_snps": None},
        "groups": groups,
        "sources": [{"scores": kg_file, "target_scores": target_file}],
    })

pathlib.Path(f"{W}/summary.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
state = "ok" if analyses else "unavailable"
pathlib.Path(f"{W}/summary.json").write_text(json.dumps({
    "schema_version": 1, "default_analysis_id": analyses[0]["analysis_id"] if analyses else "",
    "analyses": analyses, "local": {},
    "state": state, "reason_code": "" if analyses else "no_reference_space",
}, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
ad.write_manifest(f"{W}/manifest.json", ad.build_manifest(
    SAMPLE, "04b-ancestry-summary", state=state,
    reason_code="" if analyses else "no_reference_space",
    parameters={"regional_enabled": REGIONAL_ENABLED, "superpop": SUPERPOP, "subpops": SUBPOPS},
    outputs=["04_ancestry/summary.json", "04_ancestry/summary.txt"]))
print(f"wrote {W}/summary.json ({len(analyses)} reference space(s)); manifest written")
