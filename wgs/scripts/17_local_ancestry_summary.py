#!/usr/bin/env python
"""Genome-wide and per-segment local ancestry from FLARE output."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import glob, gzip, io, json, os, subprocess, pandas as pd, numpy as np
P=str(P); W=f"{P}/wgs/12_localanc"
# 失效先于计算（复审 §3.2 P0 失败生命周期）：读段/汇总中途崩溃时，旧 ok manifest 与旧
# local_ancestry.json 不得继续充当本次结果；成功路径的最终 manifest 原子覆盖这份记录。
import ancestry_data as _ad_begin
_ad_begin.begin_run_manifest(f"{W}/manifest.json", SAMPLE, "17-local-ancestry")


def an_labels(work):
    """AN 索引 → 面板标签，来自 FLARE 的 .model（VCF 头只写 "Ancestry of first haplotype"，没有标签）。

    这是唯一权威来源：写死 NorthEA/SouthEA/European/SouthAsian 的列表在换参考面板时会静默错位，
    让百分比挂到错误的面板上。
    """
    for m in sorted(glob.glob(f"{work}/la.*.model")) + sorted(glob.glob(f"{work}/calib.*.model")):
        try:
            lines = open(m, encoding="utf-8", errors="replace").read().split("\n")
        except OSError:
            continue
        for i, l in enumerate(lines):
            if l.strip().lower().startswith("# list of ancestries"):
                for j in range(i + 1, min(i + 4, len(lines))):
                    row = lines[j].strip()
                    if row and not row.startswith("#"):
                        return row.split("\t")
    return []


ANC = an_labels(W)
if not ANC:
    print("no FLARE .model with an ancestry list found; cannot label AN indices", file=sys.stderr)
    raise SystemExit(2)
segs=[]; tot=np.zeros(len(ANC)); n_mark=0; missing=[]
for f in sorted(glob.glob(f"{W}/la.*.anc.vcf.gz"),key=lambda x:int(x.split("la.")[1].split(".")[0])):
    c=f.split("la.")[1].split(".")[0]
    q=subprocess.run(["bcftools","query","-f","%POS\t[%AN1\t%AN2]\n",f],capture_output=True,text=True).stdout
    d=pd.read_csv(io.StringIO(q),sep="\t",header=None,names=["pos","a1","a2"])
    if not len(d):
        # An empty per-chromosome VCF means FLARE never finished for that chromosome -- it creates the
        # output up front, so an interrupted run leaves a 0-byte file. This used to be skipped in
        # silence: chr4 was missing from the figure while per_chrom.tsv still carried its southern
        # share, and nothing on screen compared the two. Warn loudly and name the chromosome.
        print(f"WARNING: chr{c} has no ancestry calls ({f} is empty); skipping -- rerun step 16 for it", file=sys.stderr)
        missing.append(c); continue
    n_mark+=len(d)
    for col in ("a1","a2"):
        v=d[col].values
        for i in range(len(ANC)): tot[i]+=(v==i).sum()
        # segments
        chg=np.flatnonzero(np.diff(v))+1; bounds=np.concatenate([[0],chg,[len(v)]])
        for s,e in zip(bounds[:-1],bounds[1:]):
            segs.append((c,col,int(d.pos.iloc[s]),int(d.pos.iloc[e-1]),ANC[v[s]],e-s))
    print(f"chr{c}: {len(d)} markers",flush=True)
sg=pd.DataFrame(segs,columns=["chrom","hap","start","end","anc","n_markers"])
sg["mb"]=(sg.end-sg.start)/1e6
sg.to_csv(f"{W}/segments.tsv",sep="\t",index=False)
if missing:
    print(f"\nWARNING: {len(missing)} chromosome(s) contributed no segments at all: {', '.join(missing)}."
          f" Their proportions are absent, not zero -- do not read the figure as complete coverage.",
          file=sys.stderr)
frac=tot/tot.sum()
print("\ngenome-wide ancestry fractions (both haplotypes):")
for a,f_ in zip(ANC,frac): print(f"  {a:12s} {f_*100:6.2f}%")
# 第二个来源面板的名字来自配置（la_labels），不是写死的 SouthEA
_P1=(OPT.get("la_labels") or ["", ""])[1]
big=sg[(sg.anc==_P1)&(sg.mb>=2)].sort_values("mb",ascending=False)
print(f"\n{_P1} segments >=2 Mb: {len(big)}; longest {big.mb.max():.1f} Mb" if len(big) else f"\nno {_P1} segment >=2 Mb")
print(big.head(10).to_string(index=False))
sw=sg.groupby("anc").agg(n=("mb","size"),mb=("mb","sum"),median_mb=("mb","median"))
print("\nsegment counts:"); print(sw.round(2).to_string())
pd.Series(dict(zip(ANC,frac))).to_csv(f"{W}/global.tsv",sep="\t",header=False)

# --- structured result (7.3 Local). Three different quantities are kept apart on purpose: marker
# counts, length-weighted posterior, and segment spans. Merging them into one percentage is how a
# figure ends up disagreeing with its own table.
def _panel_rows():
    try:
        rows = [l.split("\t") for l in open(f"{W}/ref.panel", encoding="utf-8").read().split("\n") if l.strip()]
    except OSError:
        return []
    n = {}
    for r in rows:
        if len(r) >= 2:
            n[r[1]] = n.get(r[1], 0) + 1
    # A panel is a source if it carries one of the configured source labels (the control panel keeps its
    # own population names and is not relabelled as one of the sources).
    _src = {x for x in (os.environ.get("LA_LAB_A", ""), os.environ.get("LA_LAB_B", "")) if x}
    return [{"id": lab, "label": lab, "reference_pops": [lab], "n_reference": int(n.get(lab, 0)),
             "role": "source" if lab in _src else "control"} for lab in ANC]


_chroms = sorted({str(c) for c in sg.chrom}, key=lambda x: int(x)) if len(sg) else []
json.dump({
    "state": "ok" if len(sg) else "unavailable",
    "reason_code": "" if len(sg) else "no_segments",
    "panels": _panel_rows(),
    "global": [{"panel_id": a, "value": float(v)} for a, v in zip(ANC, frac)],
    "per_chrom": [{"panel_id": r.anc, "chrom": str(r.chrom),
                   "value": float(r.mb / sg[sg.chrom == r.chrom].mb.sum())} for r in
                  sg.groupby(["chrom", "anc"]).mb.sum().reset_index().itertuples()],
    "calibration": [],                       # filled by 17b when the holdout calibration ran
    "missing_chroms": missing,
    "aggregation_method": "marker_counts_for_global; segment_spans_for_per_chrom",
    "n_markers": int(n_mark),
    "chroms_contributing": _chroms,
}, open(f"{W}/local_ancestry.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"wrote {W}/local_ancestry.json ({len(ANC)} panels, chroms {_chroms})")

# AN5：把这次汇总的状态写进 manifest，30 据此判断能不能用这份结果（而不是看文件在不在）。
import ancestry_data as _ad
_ad.write_manifest(f"{W}/manifest.json", _ad.build_manifest(
    SAMPLE, "17-local-ancestry", state=("ok" if len(sg) else "unavailable"),
    reason_code=("" if len(sg) else "no_segments"),
    parameters={"panels": list(ANC), "missing_chroms": list(missing),
                "aggregation": "marker_counts_for_global; segment_spans_for_per_chrom"},
    outputs=["12_localanc/local_ancestry.json", "12_localanc/segments.tsv", "12_localanc/global.tsv",
             "12_localanc/calibration.tsv"]))
print(f"wrote {W}/manifest.json (state={'ok' if len(sg) else 'unavailable'})")
