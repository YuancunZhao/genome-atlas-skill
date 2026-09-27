#!/usr/bin/env python
"""The classic 'personality gene' variants, read straight from the WGS, with what the evidence actually supports."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, pandas as pd, os, collections, bisect
from gt_alleles import gt_alleles
P=str(P); PV=KG_PFILE+".pvar"; W=f"{P}/wgs/20_behaviour"
os.makedirs(W,exist_ok=True)
# Step 01 writes target.*; fall back to the {SAMPLE}.* name of older run layouts.
V=f"{P}/wgs/00_input/target.norm.vcf.gz"
if not os.path.exists(V): V=f"{P}/wgs/00_input/{SAMPLE}.norm.vcf.gz"
PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel" / "candidate_genes.tsv"
_cg = pd.read_csv(PANEL, sep="\t").fillna("")
M = [(r.gene, r.rsid, r.variant_zh, r.claim_zh, r.evidence_zh) for r in _cg.itertuples()]
EN = {r.rsid: (r.variant_en, r.claim_en, r.evidence_en) for r in _cg.itertuples()}
pos={}
want=set(x[1] for x in M)
with open(PV) as f:
    for line in f:
        if line.startswith("#"): continue
        c,p,i,r,a=line.rstrip("\n").split("\t")[:5]
        for tok in i.split(";"):
            if tok in want and tok not in pos: pos[tok]=(c,int(p),r,a)
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

# b37 pseudoautosomal regions (same bounds as 30_build_report_data). A male's non-PAR X is
# hemizygous, so a diploid-homozygous record there is the caller's representation of one allele.
PAR_X = ((60001, 2699520), (154931044, 155260560))
def hemizygous_x(c, p):
    return c == "X" and SEX == "male" and not any(s <= p <= e for s, e in PAR_X)

rows=[]
for gene,rs,var,claim,truth in M:
    if rs not in pos: rows.append((gene,rs,var,"-","未在 1000G 面板中",claim,truth)); continue
    c,p,r,a=pos[rs]
    q=subprocess.run(["bcftools","query","-r",f"{c}:{p}-{p}","-f","%REF\t%ALT\t%FILTER\t[%GT\t%DP]\n",V],capture_output=True,text=True)
    if q.returncode != 0:
        rows.append((gene,rs,var,f"{c}:{p}",f"query failed (rc={q.returncode})",claim,truth)); continue
    t=q.stdout.strip().splitlines()
    if t:
        rec=[x.split("\t") for x in t]; hit=[x for x in rec if x[0]==r and x[1]==a] or rec
        ref,alt,flt,gt,dp=hit[0]; ga=gt_alleles(gt,ref,alt)
        if ga is None:
            geno="no call"
        else:
            if hemizygous_x(c,p) and len(ga)==2 and ga[0] is not None and ga[0]==ga[1]:
                ga=ga[:1]  # diploid-homozygous record of a hemizygous non-PAR X call
            geno="/".join(x if x is not None else "." for x in ga)
        rows.append((gene,rs,var,f"{c}:{p}",geno,claim,truth))
    elif not callable_(c,p):
        rows.append((gene,rs,var,f"{c}:{p}","no call (outside the callable mask)",claim,truth))
    else:
        rows.append((gene,rs,var,f"{c}:{p}", r if hemizygous_x(c,p) else f"{r}/{r}",claim,truth))
df=pd.DataFrame(rows,columns=["gene","rsid","variant","locus","genotype","popular_claim","what_evidence_supports"])
df["variant_en"]=[EN.get(r,("","",""))[0] for r in df.rsid]
df["popular_claim_en"]=[EN.get(r,("","",""))[1] for r in df.rsid]
df["evidence_en"]=[EN.get(r,("","",""))[2] for r in df.rsid]
df.to_csv(f"{W}/candidate_genes.tsv",sep="\t",index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",34)
print(df[["gene","rsid","variant","genotype","what_evidence_supports"]].to_string(index=False))
