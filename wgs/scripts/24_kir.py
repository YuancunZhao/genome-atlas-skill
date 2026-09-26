#!/usr/bin/env python
"""KIR gene content from T1K plus the HLA class I ligand groups (C1/C2, Bw4/Bw6) that pair with them."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import pandas as pd
P=str(P); W=f"{P}/wgs/16_panels"
import os as _os
_os.makedirs(W, exist_ok=True)
# No step of this repository runs T1K; without its kir_genotype.tsv there is nothing to type.
import os as _os
if not _os.path.exists(f"{P}/wgs/06_pgx/t1k/kir_genotype.tsv"):
    pd.DataFrame(columns=["gene","present","call","ab1","ab2"]).to_csv(f"{W}/kir.tsv",sep="\t",index=False)
    open(f"{W}/kir_summary.txt","w").write("haplotype\tunavailable (T1K was not run)\npresent\t\nC_ligands\t\nB_epitopes\t\n")
    print("WARNING: 06_pgx/t1k/kir_genotype.tsv not found (no step runs T1K) -- empty kir.tsv written", file=sys.stderr)
    sys.exit(0)
k=pd.read_csv(f"{P}/wgs/06_pgx/t1k/kir_genotype.tsv",sep="\t",header=None,
              names=["gene","n","a1","ab1","q1","a2","ab2","q2","other"])
k["present"]=(k.n>0)&(k.ab1>=5)
k["call"]=[(f"{r.a1}" + (f" / {r.a2}" if r.n>1 and r.ab2>=5 else "")) if r.present else "absent" for r in k.itertuples()]
k[["gene","present","call","ab1","ab2"]].to_csv(f"{W}/kir.tsv",sep="\t",index=False)
pres=set(k[k.present].gene)
Ahap={"KIR3DL3","KIR2DL3","KIR2DP1","KIR2DL1","KIR3DP1","KIR2DL4","KIR3DL1","KIR2DS4","KIR3DL2"}
Bonly={"KIR2DL2","KIR2DL5A","KIR2DL5B","KIR2DS1","KIR2DS2","KIR2DS3","KIR2DS5","KIR3DS1"}
hap="AA (two A haplotypes)" if not (pres & Bonly) else "Bx (at least one B haplotype)"
print("KIR genes present:", ", ".join(sorted(pres)))
print("B-haplotype-specific genes present:", ", ".join(sorted(pres & Bonly)) or "none")
print("KIR haplotype:", hap)
# HLA ligands come from this run's T1K HLA typing; the previous hard-coded genotype and the
# concluding sentence described the example sample, not this one.
import re as _re
_hla_f = f"{P}/wgs/06_pgx/t1k/dayu_genotype.tsv"
hla = {"A": [], "B": [], "C": []}
if _os.path.exists(_hla_f):
    for _line in open(_hla_f):
        _f = _line.rstrip("\n").split("\t")
        _g = _f[0].replace("HLA-", "")
        if _g not in hla: continue
        _vals = [_f[2]] + ([_f[5]] if len(_f) > 5 and _f[1].isdigit() and int(_f[1]) > 1 and _f[5] != "." else [])
        hla[_g] = [_re.sub(r"^HLA-[A-Z0-9]+\*", "", v).rsplit(":", 1)[0] for v in _vals]
C1={"01","03","07","08","12","14","16"}; C2={"02","04","05","06","15","17","18"}
cg=[("C1" if a.split(":")[0] in C1 else "C2" if a.split(":")[0] in C2 else "?") for a in hla["C"]]
BW4={"13:01","27:05","44:02","44:03","51:01","52:01","53:01","57:01","58:01","38:01","37:01","47:01","49:01","57:03","59:01","63:01","77:01"}
bw=["Bw4" if a in BW4 else "Bw6" for a in hla["B"]]
print("HLA-C ligands: " + (", ".join(f"C*{a} = {g}" for a, g in zip(hla["C"], cg)) or "no C allele typed") + f"  ->  {'/'.join(cg) or '?'}")
print("HLA-B epitopes: " + (", ".join(f"B*{a} = {e}" for a, e in zip(hla["B"], bw)) or "no B allele typed") + f"  ->  {'/'.join(bw) or '?'}")
_act = sorted(pres & {"KIR2DS1","KIR2DS2","KIR2DS3","KIR2DS5","KIR3DS1"})
print("KIR-ligand check: KIR2DL3 " + ("present with a C1 ligand" if ("KIR2DL3" in pres and "C1" in cg) else "present without a detected C1 ligand" if "KIR2DL3" in pres else "absent")
      + "; KIR3DL1 " + ("present with a Bw4 ligand" if ("KIR3DL1" in pres and "Bw4" in bw) else "present without a detected Bw4 ligand" if "KIR3DL1" in pres else "absent")
      + "; activating KIRs: " + (", ".join(_act) or "none"))
open(f"{W}/kir_summary.txt","w").write(f"haplotype\t{hap}\npresent\t{','.join(sorted(pres))}\nC_ligands\t{'/'.join(cg)}\nB_epitopes\t{'/'.join(bw)}\n")
