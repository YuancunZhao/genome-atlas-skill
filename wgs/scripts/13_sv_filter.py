#!/usr/bin/env python
"""Filter Delly SVs (PASS, non-ref, 50bp-5Mb, PRECISE or strong PE support), depth-filter DEL/DUP >=2kb by read depth,
annotate overlapping genes/exons (Ensembl 87 GFF3)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import _sex_evidence as _se
EFFECTIVE_SEX, _SEX_REC = _se.evidence_sex()   # 证据优先；回退声明时已在记录里说明原因

import subprocess, gzip, io, collections, bisect, os
import pandas as pd, numpy as np
PROJ = str(P); W = f"{PROJ}/wgs/08_sv"
# The autosomal mean depth is measured here from the step-01 mosdepth bed, the same way
# 21_somatic.py does it; a fixed constant would misstate dp_ratio and invert the DEL/DUP
# read-depth check for any sample whose coverage differs from the one it was measured on.
_auto = subprocess.run(f"zcat {PROJ}/wgs/01_qc/depth.regions.bed.gz | awk '$1 ~ /^[0-9]+$/ && $1<23 {{s+=$4; n++}} END{{print (n ? s/n : \"\")}}'",
                       shell=True, capture_output=True, text=True).stdout.strip()
MEAN_DP = float(_auto) if _auto else None
if MEAN_DP is None:
    print("WARNING: autosomal mean depth unreadable from 01_qc/depth.regions.bed.gz -- DEL/DUP read-depth check will be reported as n/a", file=sys.stderr)
# No step of this repository runs Delly; without an externally produced target.sv.bcf there is
# nothing to filter, so write an empty table and let the pipeline continue instead of dying.
import os as _os
if not _os.path.exists(f"{W}/delly/target.sv.bcf"):
    _os.makedirs(W, exist_ok=True)
    _cols = "chrom pos end svtype svlen precise pe sr mapq chr2 pos2 geno gq rc rcl rcr dr dv rr rv size dp_ratio depth_check geno_dp genes cds_overlap whole_gene_del".split()
    pd.DataFrame(columns=_cols).to_csv(f"{W}/sv_filtered.tsv", sep="\t", index=False)
    print("WARNING: 08_sv/delly/target.sv.bcf not found (no step runs Delly) -- empty sv_filtered.tsv written", file=sys.stderr)
    sys.exit(0)
q = subprocess.run(["bcftools", "query", "-i", 'FILTER="PASS" && GT="alt"', "-f",
    "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/PRECISE\t%INFO/PE\t%INFO/SR\t%INFO/MAPQ\t%INFO/CHR2\t%INFO/POS2\t[%GT\t%GQ\t%RC\t%RCL\t%RCR\t%DR\t%DV\t%RR\t%RV]\n",
    f"{W}/delly/target.sv.bcf"], capture_output=True, text=True, check=True).stdout
sv = pd.read_csv(io.StringIO(q), sep="\t", header=None, names="chrom pos end svtype svlen precise pe sr mapq chr2 pos2 geno gq rc rcl rcr dr dv rr rv".split(), dtype={"chrom": str, "chr2": str}, na_values=".")
sv["precise"] = sv.precise.fillna(0).astype(int)
sv["size"] = np.where(sv.svtype == "BND", np.nan, (sv.end - sv.pos).abs())
main = [str(i) for i in range(1, 23)] + ["X", "Y"]
sv = sv[sv.chrom.isin(main) & (sv.gq >= 20)]
nonbnd = sv[(sv.svtype != "BND") & (sv["size"] >= 50) & (sv["size"] <= 5e6) & ((sv.precise == 1) | ((sv.pe >= 5) & (sv.mapq >= 40)))].copy()
# read-depth validation for DEL/DUP >= 2kb via mosdepth 1kb bins
def region_depth(c, s, e):
    if MEAN_DP is None:
        return np.nan
    out = subprocess.run(["tabix", f"{PROJ}/wgs/01_qc/depth.regions.bed.gz", f"{c}:{max(1,s)}-{e}"], capture_output=True, text=True).stdout
    v = [float(l.split("\t")[3]) for l in out.splitlines()]
    # X depth baseline follows the sample's sex (b37 non-PAR X is hemizygous in males); a fixed 0.5
    # factor would double every female X ratio and call normal coverage a DUP. Which sex this is comes
    # from the depth evidence (_sex_evidence), not from the declaration alone: a declared-male sample
    # that was sequenced without a Y would otherwise have every Y deletion scored against a wrong baseline.
    ploidy = 0.5 if (c == "Y" or (c == "X" and EFFECTIVE_SEX.startswith("m"))) else 1.0
    base = MEAN_DP * ploidy
    return np.median(v) / base if v else np.nan
big = nonbnd[(nonbnd.svtype.isin(["DEL", "DUP"])) & (nonbnd["size"] >= 2000)]
nonbnd["dp_ratio"] = np.nan
nonbnd.loc[big.index, "dp_ratio"] = [region_depth(r.chrom, r.pos, r.end) for r in big.itertuples()]
def rd_ok(r):
    if np.isnan(r.dp_ratio): return "n/a"
    if r.svtype == "DEL": return "ok" if r.dp_ratio < 0.75 else "mismatch"
    return "ok" if r.dp_ratio > 1.25 else "mismatch"
def geno_dp(r):
    if np.isnan(r.dp_ratio): return r.geno
    if r.svtype == "DEL": return "1/1" if r.dp_ratio < 0.2 else "0/1"
    return "1/1" if r.dp_ratio > 1.8 else "0/1"
nonbnd["depth_check"] = nonbnd.apply(rd_ok, axis=1)
# Keep the caller's GT in `geno`; the read-depth estimate is a separate column so the
# crude depth threshold never silently replaces the original Delly call.
nonbnd["geno_dp"] = nonbnd.apply(geno_dp, axis=1)
# Depth consistency is a filter, not validation: count and export what it removes before
# dropping depth-inconsistent large CNVs (and large imprecise inversions), so pre-filter,
# dropped and not-assessable numbers stay auditable instead of vanishing before the stats.
n_pre = len(nonbnd)
n_mm = int((nonbnd.depth_check == "mismatch").sum())
n_na = int((nonbnd.depth_check == "n/a").sum())
_inv_drop = (nonbnd.svtype == "INV") & (nonbnd["size"] > 100000) & ~((nonbnd.precise == 1) & (nonbnd.sr >= 5))
n_inv = int(_inv_drop.sum())
nonbnd[nonbnd.depth_check == "mismatch"].to_csv(f"{W}/sv_depth_mismatch.tsv", sep="\t", index=False)
nonbnd = nonbnd[(nonbnd.depth_check != "mismatch") & ~_inv_drop]
# gene annotation
genes = collections.defaultdict(list); exons = collections.defaultdict(list)
with gzip.open(f"{PROJ}/data/ref/annot/Homo_sapiens.GRCh37.87.gff3.gz", "rt") as f:
    for l in f:
        if l.startswith("#"): continue
        p = l.split("\t")
        if p[2] == "gene" and "biotype=protein_coding" in p[8]:
            name = [x for x in p[8].split(";") if x.startswith("Name=")]; genes[p[0]].append((int(p[3]), int(p[4]), name[0][5:] if name else p[8][:20]))
        elif p[2] == "CDS":
            exons[p[0]].append((int(p[3]), int(p[4])))
for d in (genes, exons):
    for c in d: d[c].sort()
def overlap(d, c, s, e, with_name=True):
    hits = []
    for a, b, *n in d.get(c, []):
        if a > e: break
        if b >= s: hits.append(n[0] if n else True)
    return hits
nonbnd["genes"] = [",".join(sorted(set(overlap(genes, r.chrom, r.pos, r.end)))) for r in nonbnd.itertuples()]
nonbnd["cds_overlap"] = [len(overlap(exons, r.chrom, r.pos, r.end)) > 0 for r in nonbnd.itertuples()]
nonbnd["whole_gene_del"] = [",".join(g for a, b, g in genes.get(r.chrom, []) if r.svtype == "DEL" and r.pos <= a and r.end >= b) for r in nonbnd.itertuples()]
nonbnd = nonbnd.sort_values(["chrom", "pos"], key=lambda s: s.map(lambda x: main.index(x) if x in main else x) if s.name == "chrom" else s)
nonbnd.to_csv(f"{W}/sv_filtered.tsv", sep="\t", index=False)
print("filtered non-BND SVs:", len(nonbnd)); print(nonbnd.svtype.value_counts().to_dict())
print("size classes:", pd.cut(nonbnd["size"], [0, 100, 1000, 10000, 100000, 5e6]).value_counts().sort_index().to_dict())
print(f"depth check DEL/DUP>=2kb: {n_pre} pre-filter events -> kept {len(nonbnd)} (mismatch dropped {n_mm}, imprecise INV dropped {n_inv}, depth not assessable {n_na}); survivors:", nonbnd.depth_check.value_counts().to_dict())
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
cds = nonbnd[nonbnd.cds_overlap & (nonbnd.depth_check != "mismatch")]
print(f"\nSVs overlapping CDS (depth-consistent): {len(cds)}; hom by depth estimate: {(cds.geno_dp=='1/1').sum()}")
print(cds[cds.geno_dp == "1/1"][["chrom", "pos", "end", "svtype", "size", "geno", "geno_dp", "dp_ratio", "genes", "whole_gene_del"]].to_string(index=False))
print("\nSVs deleting/duplicating whole genes or >=10kb with CDS:")
print(cds[(cds.geno_dp != "1/1") & ((cds.whole_gene_del != "") | (cds["size"] >= 10000))][["chrom", "pos", "end", "svtype", "size", "geno", "geno_dp", "dp_ratio", "genes", "whole_gene_del"]].to_string(index=False))
print("\ndp_ratio distribution by type/depth-estimated genotype (DEL/DUP >= 2kb):")
b = nonbnd[nonbnd.dp_ratio.notna()]
print(b.groupby(["svtype", "geno_dp"]).dp_ratio.describe()[["count", "25%", "50%", "75%"]].round(2))
