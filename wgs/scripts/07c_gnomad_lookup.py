#!/usr/bin/env python
"""Query gnomAD v2.1.1 (GRCh37) exome+genome AF / EAS AF for candidate variants (rare LoF + ClinVar P/LP/conflicting)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import pandas as pd, numpy as np, requests, time, json, os
W = f"{P}/wgs/05_clinvar"; CACHE = f"{W}/gnomad_cache.json"
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
Q = "{ variant(variantId: \"%s\", dataset: gnomad_r2_1) { exome { ac an populations { id ac an } } genome { ac an populations { id ac an } } } }"
def look(vid):
    # a cached failure is retried on the next run; only a definitive answer is final
    if vid in cache and not cache[vid].get("error"): return cache[vid]
    js = None
    for attempt in range(6):
        try:
            r = requests.post("https://gnomad.broadinstitute.org/api", json={"query": Q % vid}, timeout=60)
            if r.status_code == 200: js = r.json(); break
        except Exception: pass
        time.sleep(10 * (attempt + 1))
    if js is None:
        cache[vid] = {"found": False, "error": True}
        return cache[vid]
    d = js.get("data", {}).get("variant")
    if d is None: res = {"found": False}
    else:
        ac = an = eac = ean = 0
        for part in ("exome", "genome"):
            x = d.get(part)
            if x:
                ac += x["ac"]; an += x["an"]
                e = [p for p in x["populations"] if p["id"] == "eas"]
                if e: eac += e[0]["ac"]; ean += e[0]["an"]
        res = {"found": True, "af": ac / an if an else None, "eas_af": eac / ean if ean else None, "ac": ac, "an": an, "eas_ac": eac, "eas_an": ean}
    cache[vid] = res; time.sleep(1.2); return res
lof = pd.read_csv(f"{W}/lof_table.tsv", sep="\t", dtype={"chrom": str})
cand = lof[(lof.eas_af.fillna(0) < 0.01) & (lof.all_af.fillna(0) < 0.01)]
plp = pd.read_csv(f"{W}/clinvar_PLP.tsv", sep="\t", dtype={"chrom": str}); conf = pd.read_csv(f"{W}/clinvar_conflicting_with_P.tsv", sep="\t", dtype={"chrom": str})
vids = sorted(set(f"{r.chrom}-{r.pos}-{r.ref}-{r.alt}" for df in (cand, plp, conf) for r in df.itertuples()))
print("querying", len(vids), "variants")
for i, v in enumerate(vids):
    look(v)
    if i % 25 == 0: json.dump(cache, open(CACHE, "w")); print(i, flush=True)
json.dump(cache, open(CACHE, "w"))
_unresolved = [v for v in vids if v not in cache or cache[v].get("error")]
if _unresolved:
    print(f"WARNING: {len(_unresolved)} of {len(vids)} queries failed; re-run this script to retry them", file=sys.stderr)
def add(df):
    g = [cache.get(f"{r.chrom}-{r.pos}-{r.ref}-{r.alt}") for r in df.itertuples()]
    df = df.copy()
    df["gnomad_found"] = [bool(x and x.get("found")) for x in g]
    df["gnomad_af"] = [x.get("af") if x else None for x in g]
    df["gnomad_eas_af"] = [x.get("eas_af") if x else None for x in g]
    df["gnomad_ac"] = [x.get("ac") if x else None for x in g]
    # never queried (or the query failed) must not look like "absent from gnomAD"
    df["gnomad_unresolved"] = [x is None or bool(x.get("error")) for x in g]
    return df
lof2 = add(lof); lof2.to_csv(f"{W}/lof_table.tsv", sep="\t", index=False)
add(plp).to_csv(f"{W}/clinvar_PLP.tsv", sep="\t", index=False); add(conf).to_csv(f"{W}/clinvar_conflicting_with_P.tsv", sep="\t", index=False)
rare = lof2[(lof2.eas_af.fillna(0) < 0.01) & (lof2.all_af.fillna(0) < 0.01) & (lof2.gnomad_af.fillna(0) < 0.01) & (lof2.gnomad_eas_af.fillna(0) < 0.01)]
pd.set_option("display.width", 250)
print(f"\nrare LoF after gnomAD (<1% all & EAS): {len(rare)}; not in gnomAD at all: {(~rare.gnomad_found & ~rare.gnomad_unresolved).sum()}; unresolved (re-run to retry): {rare.gnomad_unresolved.sum()}; hom: {(rare.zyg=='hom/hemi').sum()}")
print("in constrained genes (LOEUF<0.6):"); print(rare[rare.oe_lof_upper < 0.6][["chrom","pos","id","ref","alt","gene","csq","zyg","gnomad_af","gnomad_eas_af","gnomad_ac","gnomad_unresolved","pLI","oe_lof_upper","dp","ad"]].to_string(index=False))
rare.to_csv(f"{W}/lof_rare_final.tsv", sep="\t", index=False)
