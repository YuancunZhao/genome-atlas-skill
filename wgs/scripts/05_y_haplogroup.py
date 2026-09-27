#!/usr/bin/env python
"""Y haplogroup by greedy descent of the configured YFull tree (current_version.txt records the exact
version, 14.06.0 on this machine) using pileup allele counts from the CRAM at every branch-defining SNP with a known hg19
position (ybrowse). A branch is 'derived' if the majority of its typed SNPs (depth>=2) show the derived allele.
Reports the path, per-branch support and the terminal branch, and writes y_result.json for the report."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import json, csv, subprocess, collections, os, sys
import pandas as pd
P = str(P); W = f"{P}/wgs/03_haplo"; CRAM = READS if os.path.exists(READS) else f"{P}/wgs/00_input/{SAMPLE}.cram"; REF = FASTA
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
    # Tightened for platforms whose base qualities run optimistic (MGI/DNBSEQ):
    # NOTE -q/-Q only matter if this pileup is ever rebuilt; an existing pileup keeps its own
    # filtering, but base quality can still be re-applied below because column 6 carries it.
    subprocess.run(f"samtools mpileup -f {REF} -l {bed} -q 30 -Q 25 -d 500 {CRAM} 2>/dev/null > {PILE}", shell=True, check=True)
q = open(PILE).read()
import re
# --- quality gates (tightened; base qualities on MGI/DNBSEQ runs are optimistic) ---
# Column 6 of the pileup carries per-base qualities, so a higher BQ floor is re-applied here
# WITHOUT rebuilding the pileup. MAPQ is not recoverable from an existing pileup file.
MIN_BQ = 25      # base-quality floor applied while counting (pileup itself was built with -Q 20)
MIN_DP = 5       # minimum depth for a site to be callable at all
MIN_VOTES = 3    # minimum number of agreeing reads
CONC = 0.90      # minimum fraction of reads that must agree
print(f"quality gates: BQ>={MIN_BQ} DP>={MIN_DP} votes>={MIN_VOTES} conc>={CONC}", file=sys.stderr)
counts, refbase = {}, {}
with open(PILE) as _fh:
    for l in _fh:
        f = l.rstrip("\n").split("\t")
        if len(f) < 6: continue
        try: pos = int(f[1])
        except ValueError: continue
        ref, bases, quals = f[2].upper(), f[4], f[5]
        refbase[pos] = ref
        c = collections.Counter()
        bp = qp = 0; nb = len(bases)
        while bp < nb:
            ch = bases[bp]
            if ch == "^": bp += 2; continue          # read start + mapq char (no qual consumed)
            if ch == "$": bp += 1; continue          # read end (no qual consumed)
            if ch in "+-":                           # indel: digit(s) + sequence (no qual consumed)
                j = bp + 1; num = ""
                while j < nb and bases[j].isdigit(): num += bases[j]; j += 1
                bp = j + (int(num) if num else 0); continue
            _q = (ord(quals[qp]) - 33) if qp < len(quals) else 0
            bp += 1; qp += 1
            if _q < MIN_BQ: continue
            if ch in ".,": c[ref] += 1
            elif ch.upper() in "ACGT": c[ch.upper()] += 1
        counts[pos] = dict(c)
def state(pos, anc, der):
    """Trust the ybrowse anc/der labels as-is -- the old "hg19 reference == ancestral" shortcut only
    holds inside haplogroup O (the hg19 Y reference is an R1b individual, whose own branch SNPs carry
    the derived allele). A site counts only at depth >= MIN_DP with a >= CONC majority backed by at
    least MIN_VOTES reads."""
    c = counts.get(pos, {}); dp = sum(c.values())
    if dp < MIN_DP: return "nocov", dp, c
    nd, na = c.get(der, 0), c.get(anc, 0)
    need = max(MIN_VOTES, CONC * dp)
    if nd >= need: return "der", dp, c
    if na >= need: return "anc", dp, c
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

# --- structured result (7.3 Lineage). The text above stays as the audit trail; step 30 reads this
# file instead of regex-matching the text, which is how an empty path became an IndexError("'SAMPLE'
# is not in list") and a 1-5 site branch looked like a confirmed terminal.
import json as _json
import lineage_history as _lh
TREE_VERSION = open(f"{P}/data/ref/ytree/current_version.txt").read().strip()
# The conservative call is a rule, not a per-sample exception: walk the path from the end and fall back
# to the deepest node that still has solid support. Nodes that rest on a handful of sites are recorded
# as uncertain and never treated as a proven terminal on their own.
SOLID = 5
TAIL = 4          # 分辨率判据只看末端四级（与 HANDOFF 记录的"末端四级分辨率有限"一致）
# 复核结论优先（§7：conservative 来自证据或带理由的复核记录）；没有复核记录时退回通用规则
# conservative_from_path（连续弱链 → 上一级），并把来源写清楚，免得两种来源被当成一回事。
_HIST = _lh.load_history(_json.loads((pathlib.Path(__file__).resolve().parents[1] / "panel" /
                                      "lineage_history.json").read_text(encoding="utf-8")))
_CONS, _CONS_WHY = _lh.reviewed_call(_HIST, "y", path[-1][0] if path else None)
_if_reviewed = _CONS is not None
if _CONS:
    _CONS_SRC = "reviewed:" + str(((_HIST.get("reviewed_calls") or {}).get(
        "y:" + str(path[-1][0]), {}) or {}).get("tree_version") or _HIST.get("tree_source") or "manual")
else:
    _CONS = _lh.conservative_from_path(path, solid=SOLID, tail=TAIL)
    _CONS_SRC = f"rule:weak_chain_in_last_{TAIL}_levels(solid={SOLID})"
    _CONS_WHY = ""
    if _CONS and path and _CONS != path[-1][0]:
        _CONS_WHY = "the terminal sits at the end of a weakly supported chain"
if not path:
    _res = _lh.unavailable_lineage("y", "no_supported_path",
                                   "the tree walk found no supported branch (no usable Y pileup or no derived sites)")
else:
    _res = {"kind": "y", "state": "ok", "reason_code": "",
            "reported_hg": path[-1][0],
            "conservative_hg": _CONS, "conservative_source": _CONS_SRC, "conservative_reason": _CONS_WHY,
            "tree_source": "YFull", "tree_version": TREE_VERSION,
            "call_quality": "automatic",
            "supported_path": [{"node": b, "der": d, "anc": a, "na": o, "formed": f, "tmrca": t}
                               for b, d, a, o, f, t in path],
            # 只标末端四级：主干上支持位点少的节点（HIJK/K2 等）是因为那些 SNP 不属于本样本的
            # 谱系或未覆盖，不是"分辨率不确定"；HANDOFF 记录的分辨率极限正是末端四级。
            "uncertain_nodes": [{"node": b, "reason": f"only {d} supporting site(s); below the {SOLID}-site floor"}
                                for b, d, a, o, f, t in path[-TAIL:] if d < SOLID],
            "conflicts": [], "route_review": "automatic"}
_json.dump(_res, open(f"{W}/y_result.json", "w"), ensure_ascii=False, indent=1)
print(f"wrote {W}/y_result.json: reported={_res.get('reported_hg')} conservative={_res.get('conservative_hg')} "
      f"uncertain={len(_res.get('uncertain_nodes') or [])} state={_res['state']}")
