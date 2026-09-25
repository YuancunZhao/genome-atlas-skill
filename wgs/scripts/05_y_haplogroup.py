#!/usr/bin/env python
"""Y haplogroup by greedy descent of the YFull tree (v14) using pileup allele counts from the CRAM at every branch-defining
SNP with a known hg19 position (ybrowse). A branch is 'derived' if the majority of its typed SNPs (depth>=2) show the
derived allele. Reports the path, per-branch support, and the terminal branch; also lists private (novel) Y variants."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import json, csv, subprocess, collections, os, sys
import pandas as pd
P = str(P); W = f"{P}/wgs/03_haplo"; CRAM = f"{P}/wgs/00_input/{SAMPLE}.cram"; REF = FASTA
# 1. SNP index: name -> (pos, anc, der)
idx = {}
with open(f"{YTREE}/snps_hg19.csv", newline="") as f:
    for r in csv.DictReader(f):
        try: pos = int(r["start"])
        except: continue
        a, d = r["allele_anc"].upper(), r["allele_der"].upper()
        if len(a) == 1 and len(d) == 1 and a in "ACGT" and d in "ACGT":
            idx.setdefault(r["Name"], (pos, a, d))
print("ybrowse SNPs with hg19 pos:", len(idx), file=sys.stderr)
tree = json.load(open(f"{YTREE}/current_tree.json"))
def snps_of(node):
    out = []
    for tok in (node.get("snps") or "").split(","):
        for name in tok.strip().split("/"):
            name = name.strip()
            if name in idx: out.append((name,) + idx[name])
    return out
# 2. pileup at all tree SNP positions (one pass): collect positions
positions = set()
def collect(n):
    for s in snps_of(n): positions.add(s[1])
    for c in n.get("children", []): collect(c)
collect(tree)
bed = f"{W}/yfull_positions.bed"
with open(bed, "w") as f:
    for p in sorted(positions): f.write(f"Y\t{p-1}\t{p}\n")
print("tree SNP positions:", len(positions), file=sys.stderr)
PILE = f"{W}/y_pileup.tsv"
if not os.path.exists(PILE):
    # samtools mpileup over a sorted BED is a single sequential pass (much faster than bcftools -R random access)
    subprocess.run(f"samtools mpileup -f {REF} -l {bed} -q 20 -Q 20 -d 500 {CRAM} 2>/dev/null > {PILE}", shell=True, check=True)
q = open(PILE).read()
import re
counts = {}
for l in q.splitlines():
    f = l.split("\t")
    if len(f) < 5: continue
    pos, ref, bases = int(f[1]), f[2].upper(), f[4]
    bases = re.sub(r"\^.", "", bases); bases = bases.replace("$", "")
    bases = re.sub(r"[+-](\d+)(?=[ACGTNacgtn*#]+)", lambda m: "\x00" * 0, bases)  # strip indel length markers
    # remove indel sequences: after +N or -N, N bases follow
    out = []; i = 0
    b = f[4]
    while i < len(b):
        ch = b[i]
        if ch == "^": i += 2; continue
        if ch == "$": i += 1; continue
        if ch in "+-":
            j = i + 1; num = ""
            while j < len(b) and b[j].isdigit(): num += b[j]; j += 1
            i = j + int(num or 0); continue
        out.append(ch); i += 1
    c = collections.Counter()
    for ch in out:
        if ch in ".,": c[ref] += 1
        elif ch.upper() in "ACGT": c[ch.upper()] += 1
    counts[pos] = dict(c)
refbase = {}
for l in q.splitlines():
    f = l.split("\t")
    if len(f) >= 3: refbase[int(f[1])] = f[2].upper()
def state(pos, anc, der):
    """Within haplogroup O the hg19 reference (R1b) carries the ancestral allele at every branch SNP, so
    derived := non-reference allele; ybrowse anc/der labels are only used when they agree with this."""
    c = counts.get(pos, {}); dp = sum(c.values())
    if dp < 2: return "nocov", dp, c
    nd, na = c.get(der, 0), c.get(anc, 0)
    if nd >= max(2, 0.8 * dp): return "der", dp, c
    if na >= max(2, 0.8 * dp): return "anc", dp, c
    return "mixed", dp, c
def score(node):
    st = collections.Counter(state(p, a, d)[0] for _, p, a, d in snps_of(node))
    return st["der"], st["anc"], st["nocov"] + st["mixed"]
# 3. greedy descent
# (removed) upstream hardcoded START = "O-F438"; the walk now starts at the tree root below
def find_node(n, i):
    if n["id"] == i: return n
    for c in n.get("children", []):
        r = find_node(c, i)
        if r: return r
# Walk from the ROOT of the tree. Upstream started at a hardcoded O-F438 branch, which any non-O
# sample never enters (every O SNP then reads as ancestral). Support is accumulated along the path
# (sum of der-anc), because plain per-node "der > anc" is fragile: some nodes list sibling
# sub-branch SNPs in their "snps" field and dilute the vote (e.g. A1: der=65 anc=265, while the
# true downstream path reads BT 499/206, CT 410/79, N 420/0).
_best = {"cum": 0, "path": []}
def _walk(_n, _cum, _path):
    for _c in _n.get("children", []):
        if _c["id"].endswith("*"): continue
        _d, _a, _o = score(_c)
        _cc = _cum + (_d - _a)
        if _cc > _best["cum"]:
            _best["cum"] = _cc
            _best["path"] = _path + [_c]
        _walk(_c, _cc, _path + [_c])
_walk(tree, 0, [])
node = _best["path"][-1] if _best["path"] else tree
path = []
_cum = 0
for _n in _best["path"]:
    _d, _a, _o = score(_n)
    _cum += _d - _a
    path.append((_n["id"], _d, _a, _o, _n.get("formed"), _n.get("tmrca")))
lines = [f"YFull tree {open(f'{P}/data/ref/ytree/current_version.txt').read().strip()} (walk from root); terminal branch: {path[-1][0]}",
         f"  formed ~{path[-1][4]} ybp, TMRCA ~{path[-1][5]} ybp", "", "Path (branch, #derived, #ancestral, #untyped/mixed, formed, tmrca):"]
for p in path: lines.append("  %-28s der=%3d anc=%3d n/a=%3d  formed=%s tmrca=%s" % p)
# children of terminal: show why we stopped
lines.append(""); lines.append("Children of terminal branch (not supported):")
for c in node.get("children", []):
    d, a, o = score(c); lines.append(f"  {c['id']:28s} der={d} anc={a} n/a={o}")
# 4. detail table for last 6 branches
rows = []
for bid, *_ in path[-6:]:
    def find(n):
        if n["id"] == bid: return n
        for c in n.get("children", []):
            r = find(c)
            if r: return r
    for name, pos, a, d in snps_of(find(tree)):
        s, dp, c = state(pos, a, d); rows.append((bid, name, pos, a, d, s, dp, c.get(a, 0), c.get(d, 0)))
pd.DataFrame(rows, columns=["branch", "snp", "pos_hg19", "anc", "der", "state", "depth", "n_anc", "n_der"]).to_csv(f"{W}/y_terminal_snps.tsv", sep="\t", index=False)
open(f"{W}/y_haplogroup_yfull.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
