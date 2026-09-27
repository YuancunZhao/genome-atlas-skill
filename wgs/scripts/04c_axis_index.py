#!/usr/bin/env python
"""Per-chromosome position on the configured axis (first configured population vs the second).

This is a *relative position* between two named reference populations, not an ancestry proportion, and
the axis populations and scope come from the configuration (AN2).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, io, os, pandas as pd, numpy as np
P=str(P); W=f"{P}/wgs/04_ancestry"; O=f"{P}/wgs/18_ns"; os.makedirs(O,exist_ok=True)
# The regional axis is only meaningful when the configuration explicitly names two populations: the
# module's own default must not be mistaken for a choice (7.3). Products are named after the scope.
SCOPE="regional"
if not AXIS or len(AXIS) != 2 or not REGIONAL_ENABLED:
    print("axis not configured (axis_pops / ref_superpop); run 04 to build the regional space first")
    raise SystemExit(0)
A,B=AXIS[0],AXIS[1]
_wts=f"{W}/{SCOPE}.pca.eigenvec.allele"
if not os.path.exists(_wts):
    raise SystemExit(f"{_wts} not found: the regional PCA was built under a different scope name")
wts=pd.read_csv(_wts,sep="\t",dtype={"#CHROM":str})
wts.columns=["chrom","id","ref","alt","a1"]+[f"PC{i}" for i in range(1,11)]
af=pd.read_csv(f"{W}/{SCOPE}.pca.afreq",sep="\t",dtype={"#CHROM":str}); af.columns=["chrom","id","ref","alt","alt_freq","n"]
kg=pd.read_csv(f"{W}/{SCOPE}.proj.sscore",sep="\t").rename(columns={"#IID":"iid"})
me=pd.read_csv(f"{W}/{SCOPE}.target.proj.sscore",sep="\t")
# per-chromosome scores: rerun --score restricted to each chromosome
rows=[]
for c in [str(i) for i in range(1,23)]:
    ids=wts[wts.chrom==c].id.drop_duplicates()
    idf=f"{O}/ids.{c}.txt"; ids.to_csv(idf,index=False,header=False)
    for who,pf,out in [("kg",f"{W}/kg.common",f"{O}/kg.{c}"),("target",f"{P}/wgs/02_complete/{SAMPLE}.1kg",f"{O}/target.{c}")]:
        subprocess.run([PLINK2,"--pfile",pf,"--extract",idf,"--read-freq",f"{W}/{SCOPE}.pca.afreq",
                        "--score",f"{W}/{SCOPE}.pca.eigenvec.allele","2","5","header-read","no-mean-imputation","variance-standardize",
                        "--score-col-nums","6","--out",out,"--threads",str(THREADS),"--memory",str(int(float(MEM_GB) * 1000))],capture_output=True)
    k=pd.read_csv(f"{O}/kg.{c}.sscore",sep="\t"); d=pd.read_csv(f"{O}/target.{c}.sscore",sep="\t")
    k=k.rename(columns={"#IID":"iid"})
    chb=k[k.Population==A].PC1_AVG; chs=k[k.Population==B].PC1_AVG
    x=d.PC1_AVG.iloc[0]
    idx=(x-chs.mean())/(chb.mean()-chs.mean())
    rows.append(dict(chrom=c,n_snp=len(ids),pc1=x,index=round(idx,3),sd_from_chb=round((x-chb.mean())/chb.std(),2),sd_from_chs=round((x-chs.mean())/chs.std(),2),axis_a=A,axis_b=B))
    print(rows[-1],flush=True)
    for f in [f"{O}/kg.{c}.sscore",f"{O}/target.{c}.sscore",idf]: os.remove(f)
df=pd.DataFrame(rows); df.to_csv(f"{O}/ns_index.tsv",sep="\t",index=False)
# genome-wide
chb=kg[kg.Population==A].PC1_AVG; chs=kg[kg.Population==B].PC1_AVG; x=me.PC1_AVG.iloc[0]
g=dict(index=round((x-chs.mean())/(chb.mean()-chs.mean()),3),sd_from_chb=round((x-chb.mean())/chb.std(),2),sd_from_chs=round((x-chs.mean())/chs.std(),2))
print("genome-wide:",g); pd.Series(g).to_csv(f"{O}/ns_genome.tsv",sep="\t",header=False)
print(f"per-chromosome index mean {df['index'].mean():.2f}, sd {df['index'].std():.2f}, range {df['index'].min():.2f}-{df['index'].max():.2f}")
