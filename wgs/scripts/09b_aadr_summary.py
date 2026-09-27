#!/usr/bin/env python
"""Nearest present-day and ancient groups in the AADR Human Origins PCA space (AN2).

Distance convention: every individual's distance to the target is computed first, and a group's value is
the mean over its members. The previous version compared each group's *centroid* against the target,
which is a different statistic -- the same group can rank differently under the two conventions, so the
two summaries could disagree about who is nearest. Both 04b and this step now use mean-of-members.

The eligible set comes from ancestry_data.eligible_records(), so this summary, the report and any map
see exactly the same people. The target is found by its internal kind, never by display name: a
reference group may legitimately share the sample's display name.
"""
import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...
from wgsconfig import (SAMPLE, AADR, AADR_ANNOTATION, MIN_CR_ANCIENT, MIN_CR_MODERN,
                       MIN_PROJECTION_SNPS, MIN_GROUP_N, AADR_ANCIENT_PREFIX, AADR_MODERN)

import pandas as pd, numpy as np
import ancestry_data as ad

W = f"{P}/wgs/11_aadr"
ANALYSIS_ID = "aadr-human-origins"


def _state(state, reason, detail):
    ad.write_manifest(f"{W}/manifest.json", ad.build_manifest(
        SAMPLE, "09b-aadr-summary", state=state, reason_code=reason, detail=detail))
    sys.exit(f"09b {state}: {reason} ({detail})")


if not (pathlib.Path(W) / "proj.sscore").exists():
    _state("unavailable", "missing_projection", f"{W}/proj.sscore not found")
if not (pathlib.Path(W) / "samples.tsv").exists():
    _state("unavailable", "missing_samples", f"{W}/samples.tsv not found")

s = pd.read_csv(f"{W}/proj.sscore", sep="\t").rename(columns={"#FID": "label", "IID": "iid"})
meta = pd.read_csv(f"{W}/samples.tsv", sep="\t")
# Both tables carry a 'label' column with different meanings (panel group vs. .fam FID); keep them
# apart instead of letting pandas suffix them into label_x/label_y and hoping the right one is picked.
if "label" in meta.columns:
    meta = meta.rename(columns={"label": "sample_label"})
s = s.merge(meta, on="iid", how="left")
if "label" not in s.columns and "label_x" in s.columns:
    s = s.rename(columns={"label_x": "label"})

pcs = [f"PC{i}_AVG" for i in range(1, 5)]
present = [c for c in pcs if c in s.columns]
if len(present) < 4:
    _state("unavailable", "insufficient_components", f"proj.sscore carries {present}; four are needed")

me = s[s.kind == "target"]
if len(me) != 1:
    _state("unavailable", "target_not_unique", f"{len(me)} rows carry kind=target")
D = me[present].values[0]
s["distance_to_target"] = np.sqrt(((s[present].values - D) ** 2).sum(1))

# n_called_snps / call_rate must be measured on the final pruned site set, not on the extraction-stage
# prefilter: the old code compared a 30%-filter call rate against the 50% ancient gate, which are
# different denominators. plink2 --missing over aadr.bed gives the real per-individual counts.
_MISS = pathlib.Path(f"{W}/aadr.smiss")   # plink2 writes .smiss/.vmiss; plink1 would write .imiss
if not _MISS.exists():
    if not pathlib.Path(f"{W}/aadr.bed").exists():
        _state("unavailable", "missing_panel", f"{W}/aadr.bed not found; run 08 first")
    import subprocess
    r = subprocess.run([PLINK2, "--bfile", f"{W}/aadr", "--missing", "--out", f"{W}/aadr",
                        "--threads", str(THREADS), "--memory", str(int(float(MEM_GB) * 1000))],
                       capture_output=True, text=True)
    if r.returncode != 0 or not _MISS.exists():
        _state("unavailable", "missingness_failed", (r.stderr or r.stdout or "")[-400:])
imiss = pd.read_csv(_MISS, sep=r"\s+").rename(columns={"#FID": "label", "IID": "iid"})
# plink2 sample-missing report: OBS_CT = non-missing sites, MISSING_CT = missing ones, F_MISS = ratio.
imiss["n_called_snps"] = imiss["OBS_CT"].astype(int)
imiss["call_rate"] = 1.0 - imiss["F_MISS"].astype(float)
s = s.drop(columns=[c for c in ("n_called_snps", "call_rate") if c in s.columns]) \
     .merge(imiss[["iid", "n_called_snps", "call_rate"]], on="iid", how="left")
if s["call_rate"].isna().any():
    # A record in proj.sscore that is not in the panel: keep it, but it cannot pass a coverage gate.
    s["call_rate"] = s["call_rate"].fillna(-1.0)
say = print
if "source_population_id" not in s.columns:
    s["source_population_id"] = s.get("sample_label", s.get("label"))
if "group_id" not in s.columns:
    # Group key is (dataset, population, location) per 7.3; the location half keeps the same population
    # at different sites apart instead of averaging them together.
    s["group_id"] = s["source_population_id"]
if "record_id" not in s.columns:
    s["record_id"] = s["iid"]

# Geographic and dating fields come from the normalised .anno metadata that 08 wrote in full (7.3
# requires the history views to read it, and the earlier code never carried it into this summary at
# all -- which is why every record looked unmapped). Only columns this frame lacks are merged in.
META = pathlib.Path(f"{W}/reference_metadata.tsv")
if META.exists():
    m = pd.read_csv(META, sep="\t", low_memory=False)
    take = [c for c in ("record_id", "locality", "location_id", "latitude", "longitude",
                        "location_precision", "location_source", "date_min_bp", "date_max_bp",
                        "y_hg_raw", "mt_hg_raw", "publication") if c in m.columns]
    m = m[take].rename(columns={"record_id": "iid"})
    add = [c for c in m.columns if c != "iid" and c not in s.columns]
    s = s.merge(m[["iid"] + add], on="iid", how="left")
for c in ("latitude", "longitude", "n_called_snps", "call_rate"):
    if c not in s.columns:
        s[c] = np.nan
# samples.tsv written before AN1 carries "date"; the AN1 field is "date_mean_bp". Accept either, so a
# panel produced by an older run still reports real dates instead of zeros.
_DATE_COL = "date_mean_bp" if "date_mean_bp" in s.columns else ("date" if "date" in s.columns else None)
if _DATE_COL and "date_mean_bp" not in s.columns:
    s["date_mean_bp"] = pd.to_numeric(s[_DATE_COL], errors="coerce")

recs = s.to_dict("records")
groups_mod = ad.group_summaries(
    ad.eligible_records(recs, kind="modern", min_rate=MIN_CR_MODERN, min_snps=MIN_PROJECTION_SNPS),
    min_group_n=MIN_GROUP_N)
groups_anc = ad.group_summaries(
    ad.eligible_records(recs, kind="ancient", min_rate=MIN_CR_ANCIENT, min_snps=MIN_PROJECTION_SNPS),
    min_group_n=MIN_GROUP_N)

# --- audit tables keep their previous file names and columns. Step 30 still reads them, so the
# legacy column names (d, date) are written next to the explicit ones rather than replacing them;
# AN5 moves the consumers over and then these can go.
pd.DataFrame([{"label": g["label"], "n": g["n"], "d": g["distance_mean"]} for g in groups_mod]
             ).to_csv(f"{W}/near_modern.tsv", sep="\t", index=False)
# The table has always been the n>=2 view, and its "date" column is the mean date of the group (not the
# lower bound of a range that most records do not have). Keep both, so existing readers are unaffected.
# per-group mean date, taken from the same rows this summary is built on
_anc_rows = s[s.kind == "ancient"]
_pub = {g["label"]: g for g in groups_anc if g["n"] >= MIN_GROUP_N}
_dates = {}
for g in _pub.values():
    sub = pd.to_numeric(_anc_rows[_anc_rows.source_population_id == g["group_id"]].get(_DATE_COL),
                        errors="coerce").dropna() if _DATE_COL else pd.Series(dtype=float)
    _dates[g["label"]] = float(sub.mean()) if len(sub) else 0.0
pd.DataFrame([{"label": g["label"], "n": g["n"], "d": g["distance_mean"], "date": _dates[g["label"]]}
              for g in _pub.values()]).to_csv(f"{W}/near_ancient.tsv", sep="\t", index=False)
s["d"] = s["distance_to_target"]
if "date_mean_bp" in s.columns:
    s["date"] = s["date_mean_bp"]            # legacy alias; 30 reads 'date'
s.to_csv(f"{W}/proj_annotated.tsv", sep="\t", index=False)
if not s["date"].notna().any() if "date" in s.columns else True:
    print("warning: no dates available in samples.tsv (run 08 again to get date_mean_bp)")

def _clean(v):
    if isinstance(v, float) and (np.isnan(v) or np.isinf(v)):
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    return v

excluded = [r for r in recs if r.get("eligible") is False]
by_reason = pd.Series([r.get("exclusion_reason") for r in excluded]).value_counts().to_dict() if excluded else {}
loc_ok = sum(1 for r in recs if r.get("latitude") is not None)
prev = ad.read_manifest(f"{W}/manifest.json") or {}
release = prev.get("reference_release") or ""
summary = {
    "schema_version": 1, "analysis_id": ANALYSIS_ID, "dataset": "AADR",
    "reference_release": release, "scope": "global", "state": "ok", "reason_code": "",
    "components": [1, 2, 3, 4], "metric": "mean_individual_pc_distance",
    "thresholds": {"min_call_rate_modern": MIN_CR_MODERN, "min_call_rate_ancient": MIN_CR_ANCIENT,
                   "min_projection_snps": MIN_PROJECTION_SNPS, "min_group_n": MIN_GROUP_N,
                   "ancient_prefix": AADR_ANCIENT_PREFIX, "modern_groups": AADR_MODERN},
    "counts": {"selected": len(recs), "eligible": len(recs) - len(excluded), "excluded": len(excluded),
               "mapped": loc_ok, "unmapped": len(recs) - loc_ok},
    "excluded_reasons": by_reason,
    "target": {"sample_id": SAMPLE, "pcs": [float(x) for x in D],
               "call_rate": _clean(me.call_rate.iloc[0]) if "call_rate" in me.columns else None,
               "n_called_snps": _clean(me.n_called_snps.iloc[0]) if "n_called_snps" in me.columns else None},
    "records": [{k: _clean(v) for k, v in r.items()} for r in recs],
    "groups": {"modern": groups_mod, "ancient": groups_anc},
    "sources": [{"aadr_prefix": AADR, "annotation": AADR_ANNOTATION or (AADR + ".anno")}],
}
pathlib.Path(f"{W}/summary.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
ad.write_manifest(f"{W}/manifest.json", ad.build_manifest(
    SAMPLE, "09b-aadr-summary", state="ok", reference_release=release,
    parameters=summary["thresholds"],
    outputs=["11_aadr/summary.json", "11_aadr/proj_annotated.tsv", "11_aadr/near_modern.tsv",
             "11_aadr/near_ancient.tsv"]))

print(f"== nearest present-day groups (mean individual distance, n>={MIN_GROUP_N})")
print(pd.DataFrame(groups_mod[:12]).to_string(index=False) if groups_mod else "  (none)")
print(f"\n== nearest ancient groups (call rate >= {MIN_CR_ANCIENT}, n>={MIN_GROUP_N})")
print(pd.DataFrame(groups_anc[:20]).to_string(index=False) if groups_anc else "  (none)")
print(f"\nexcluded reasons: {by_reason}")
print(f"target {SAMPLE} PCs {np.round(D, 4)}")
print(f"wrote {W}/summary.json")
