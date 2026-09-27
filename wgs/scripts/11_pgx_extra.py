#!/usr/bin/env python
"""PharmGKB-style extra drug-gene variants that PharmCAT does not call, read straight from the WGS."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, pandas as pd, os, collections, bisect
from gt_alleles import gt_alleles
P=str(P); V=f"{P}/wgs/00_input/{SAMPLE}.norm.vcf.gz"
os.makedirs(f"{P}/wgs/16_panels", exist_ok=True)
PV=KG_PFILE+".pvar"
PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel" / "pgx_extra.tsv"
M = [(r.rsid, r.gene, r.allele, r.drug) for r in pd.read_csv(PANEL, sep="\t").fillna("").itertuples()]
seen=set(); want=[]
for rs,g,a,d in M:
    if rs in seen: continue
    seen.add(rs); want.append((rs,g,a,d))
pos={}
with open(PV) as f:
    for line in f:
        if line.startswith("#"): continue
        c,p,i,r,al=line.rstrip("\n").split("\t")[:5]
        for tok in i.split(";"):
            if tok in seen and tok not in pos: pos[tok]=(c,int(p),r,al)
# Callable mask, same convention as 23_bloodgroups.py: a site with no VCF record is a reference
# homozygote only when it lies inside the callable bed. Outside the mask, or when the query
# itself fails, the record must stay a no-call instead of silently becoming a hom-ref.
_call = collections.defaultdict(list)
_call_bed = pathlib.Path(P)/"wgs/00_input/callable.bed"
if _call_bed.exists():
    for l in open(_call_bed):
        f_ = l.split()[:3]
        if len(f_) == 3: _call[f_[0]].append((int(f_[1]), int(f_[2])))
def callable_(c, p):
    v = _call.get(c, []); i = bisect.bisect_right([x[0] for x in v], p-1)-1
    return i >= 0 and v[i][0] <= p-1 < v[i][1]

rows=[]
for rs,g,a,d in want:
    if rs not in pos: rows.append((g,rs,a,d,"-","not in 1000G panel")); continue
    c,p,r,al=pos[rs]
    q=subprocess.run(["bcftools","query","-r",f"{c}:{p}-{p}","-f","%REF\t%ALT\t%FILTER\t[%GT\t%DP]\n",V],capture_output=True,text=True)
    if q.returncode != 0:
        rows.append((g,rs,a,d,"-",f"query failed (rc={q.returncode})")); continue
    t=q.stdout.strip().splitlines()
    if t:
        rec=[x.split("\t") for x in t]; hit=[x for x in rec if x[0]==r and x[1]==al] or rec
        ref,alt,flt,gt,dp=hit[0]; ga=gt_alleles(gt,ref,alt)
        if ga is None:
            rows.append((g,rs,a,d,"-","no call"))
        else:
            ev=f"{flt} DP={dp}"
            if any(x is None for x in ga): ev+="; partial no call"
            rows.append((g,rs,a,d,"/".join(x if x is not None else "." for x in ga),ev))
    elif not callable_(c,p):
        rows.append((g,rs,a,d,"-","no call (outside the callable mask)"))
    else:
        rows.append((g,rs,a,d,f"{r}/{r}","hom-ref (callable)"))
df=pd.DataFrame(rows,columns=["gene","rsid","allele","drug","genotype","evidence"]).sort_values(["gene","rsid"])
df.to_csv(f"{P}/wgs/16_panels/pgx_extra.tsv",sep="\t",index=False)
carr=df[df.evidence.str.startswith("PASS")
        & df.genotype.str.split("/").apply(lambda xs: len(xs)==2 and "." not in xs and xs[0]!=xs[1])]
print(f"{len(df)} markers checked; {len(carr)} heterozygous / non-reference")
pd.set_option("display.width",220); print(df[df.evidence.str.startswith("PASS")].to_string(index=False))
