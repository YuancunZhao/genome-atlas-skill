#!/usr/bin/env python
"""MT: heteroplasmy from pileup allele depths (alt fraction >= 3%, depth >= 100), haplogrep3 from PASS haploid calls
(plus haplogrep3 on the VCF directly for cross-check)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import os, subprocess, pandas as pd, numpy as np
PROJ = str(P); W = f"{PROJ}/wgs/03_haplo"; os.makedirs(W, exist_ok=True)
# Step 01 writes target.*; fall back to the {SAMPLE}.* name of older run layouts.
PASS_VCF = f"{PROJ}/wgs/00_input/target.pass.vcf.gz"
if not os.path.exists(PASS_VCF): PASS_VCF = f"{PROJ}/wgs/00_input/{SAMPLE}.pass.vcf.gz"
rows = []
for line in open(f"{W}/mt_pileup_ad.tsv"):
    pos, ref, alt, ad = line.rstrip("\n").split("\t"); alts = alt.split(","); ads = list(map(int, ad.split(",")))
    dp = sum(ads)
    for a, n in zip(alts, ads[1:]):
        if a == "<*>": continue
        rows.append((int(pos), ref, a, dp, n, n/dp if dp else 0))
df = pd.DataFrame(rows, columns=["pos", "ref", "alt", "dp", "alt_reads", "af"])
het = df[(df.dp >= 100) & (df.af >= 0.03) & (df.af <= 0.97) & (df.alt_reads >= 10)]
hom = df[(df.af > 0.97) & (df.dp >= 100)]   # same coverage floor as the heteroplasmy call
het.to_csv(f"{W}/mt_heteroplasmy.tsv", sep="\t", index=False)
print("homoplasmic (AF>97%):", len(hom), " heteroplasmic (3-97%, DP>=100, >=10 reads):", len(het))
print(het.to_string(index=False))
# HSD for haplogrep3, built from the pileup's homoplasmic calls. The PASS VCF carries no MT records
# (the call set has no MT interval), so building it from those left an empty mutation list and
# haplogrep3 classified on no evidence at all -- only the low quality score hinted at it. Haplogrep
# expects the mutations relative to rCRS (which is the pileup reference) as a space-separated list.
muts = [f"{int(r.pos)}{r.alt}" for r in hom.itertuples() if len(r.ref) == 1 and len(r.alt) == 1]
print("HSD mutations from the pileup (homoplasmic, DP>=100):", len(muts), file=sys.stderr)
with open(f"{W}/target_mt.hsd", "w") as f:
    f.write(f"SampleId\tRange\tHaplogroup\tPolymorphisms\n{SAMPLE}\t1-16569\t?\t" + " ".join(muts) + "\n")
# haplogrep3 3.3.2 writes its report but then never exits: bound it with a timeout and
_hg = subprocess.run(["timeout","180", HAPLOGREP3, "classify", "--tree", "phylotree-rcrs@17.2",
                      "--in", f"{W}/target_mt.hsd", "--out", f"{W}/haplogrep3.txt",
                      "--extend-report"], capture_output=True)  # accept exit 124 if output was written
if _hg.returncode not in (0, 124):
    print(_hg.stderr.decode(errors="replace")[-500:], file=sys.stderr)
print(open(f"{W}/haplogrep3.txt").read()[:800])
