#!/usr/bin/env python
"""28_f3_stats.py -- f3 statistics on the pruned-SNP panels (outgroup-f3 ranking, adjacent
contrasts, pooled admixture-f3), for both the modern 1000G side and the ancient AADR side.

Inputs already produced by earlier steps:
  modern:  04_ancestry/kg.common.{psam,pvar} + 04_ancestry/prune.prune.in
           + 02_complete/{SAMPLE}.1kg / target.1kg (the target's 1kg-space genotypes)
  ancient: 11_aadr/aadr.{bed,bim,fam} + 11_aadr/prune.prune.in
           + panel/ancestry_locations.tsv (source_id -> region north/south via its note)

Group definitions live in panel/f3_groups.tsv (modern) and the region token of
ancestry_locations.tsv notes (ancient); nothing sample-specific is written here.

Writes:
  04_ancestry/f3/f3_results.tsv    modern tests (one row per test, f3/SE/Z)
  11_aadr/anc_f3_results.tsv       ancient tests
  04_ancestry/f3/f3_stats.json     combined payload 30_build_report_data.py reads

f3(A;B,C) per site = (pA-pB)(pA-pC) - pA(1-pA)/(nA-1), where nA is the number of called
alleles for the FIRST population at that site. The subtracted term is the finite-sample
(finite-sample-corrected) estimator convention qp3pop uses: E[(pAhat-pBhat)(pAhat-pChat)]
= f3 + pA(1-pA)/nA, so the sample-frequency product of a small panel (e.g. a single diploid
target, nA=2) is positively biased by its own sampling variance. Without the correction a
one-diploid target with zero true f3 yields a strictly positive statistic.

SE/Z come from a delete-one-block jackknife (5 Mb blocks) over the same corrected per-site
values; the correction term cancels inside same-site paired contrasts (same first population),
which is why contrast rows keep their exact difference. Estimator verified against analytic
zero/known-value conditions (tests), not against a bundled ADMIXTOOLS run. Missing policy:
a site needs every set present with >=2 called alleles; REF/ALT-flipped sites are dropped
and counted, never mixed. SE of exactly 0, a single block, or any non-finite value yields
null f3/se/z entries ("unavailable"), never NaN or Infinity in the JSON.

Conventions follow the ADMIXTOOLS qp3pop documentation:
  https://uqrmaie1.github.io/admixtools/reference/qp3pop.html
  https://uqrmaie1.github.io/admixtools/reference/f3blockdat_from_geno.html
  https://uqrmaie1.github.io/admixtools/articles/fstats.html
A significantly negative admixture-f3 is the signature of A descending from a mix of B and
C; outgroup-f3 rises with shared drift between the target and the profiled group.
"""
import json
import math
import pathlib
import re
import subprocess
import sys

import numpy as np
import pandas as pd

PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel"
BLOCK_MB = 5_000_000
MIN_GROUP_N = 20          # an ancient region pool below this is too small to pool quietly
# 复审 §3.2 W-T2：估计量串随加权公式的引入而改写——它同时是 28 manifest/30 准入的比对键，
# 改串让按旧等块公式算出的旧 f3_stats.json 判 stale，不复算前不认证旧 Z。
ESTIMATOR_ID = ("site-mean[(pA-pB)(pA-pC) - pA(1-pA)/(nA-1)] "
                "+ weighted delete-one-block(5Mb) jackknife SE (unequal blocks)")


def _num(v, nd=6):
    """Round for the payload; None/non-finite stay None so the JSON says 'unavailable'."""
    if v is None:
        return None
    v = float(v)
    return round(v, nd) if math.isfinite(v) else None


def _cfg():
    """Bind the run configuration lazily so the pure helpers stay importable in tests
    (importing wgsconfig creates the work directories -- a side effect tests must not cause)."""
    global W, PLINK2, SAMPLE, THREADS
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import wgsconfig
    W, PLINK2, SAMPLE, THREADS = wgsconfig.W, wgsconfig.PLINK2, wgsconfig.SAMPLE, wgsconfig.THREADS


# ---------------------------------------------------------------- pure helpers (tested)
def jackknife(x, bcode, nblk):
    """Delete-one-block mean and SE of the per-site statistic x (one value per site).

    Unequal blocks use the weighted block jackknife (admixtools R/resampling.R
    jack_vec_stats2): with h_b = n/n_b, est = weighted.mean(loo, 1 - 1/h_b) and
    var = mean((est - loo)^2 * (h_b - 1)). The plain equal-block formula this
    replaced inflated the SE whenever block site counts differed (5 Mb blocks at
    chromosome edges are never equal): on per-site [0 x100, 0.1, 0.2] with block
    sizes [100, 1, 1] both give mean 0.002941176, but the SE is 0.09901155 vs the
    official 0.01741909. With equal blocks the two formulas coincide, so the
    hand-checked equal-block cases are unchanged.

    SE == 0 (constant statistic / single block) or a non-finite theta is reported as
    (theta, se, None) / (None, None, None): a zero-SE point estimate has no usable Z, and
    Infinity must never reach the JSON.
    """
    x = np.asarray(x, dtype=float)
    if len(x) == 0 or not np.isfinite(x).all() or nblk < 2:
        return None, None, None
    S = x.sum()
    Sb = np.bincount(bcode, weights=x, minlength=nblk)
    Nb = np.bincount(bcode, minlength=nblk).astype(float)
    theta = S / len(x)
    loo = (S - Sb) / (len(x) - Nb)
    h = len(x) / Nb
    w = 1.0 - 1.0 / h                       # 0 only when a block holds every site (nblk==1)
    est = (w * loo).sum() / w.sum() if w.sum() > 0 else theta
    var = (((est - loo) ** 2) * (h - 1.0)).sum() / nblk
    se = math.sqrt(var)
    if not (math.isfinite(theta) and math.isfinite(se)):
        return None, None, None
    z = theta / se if se > 0 else None
    return float(theta), float(se), (float(z) if z is not None else None)


def block_codes(chroms, positions, block=BLOCK_MB):
    """Site -> 0-based jackknife block id (chrom + floor(pos/block))."""
    chroms = pd.Series(chroms).astype(str).str.replace("chr", "", regex=False)
    return pd.factorize((chroms + "_" + (pd.Series(positions) // block).astype(str)).values)[0]


def merge_freqs(freqs, pos, min_alleles=2):
    """Join per-set afreq tables on variant ID, validating REF/ALT across sources.

    `freqs` maps set name -> raw plink2 .afreq DataFrame; `pos` is ID -> (CHROM_N, POS).
    Besides the set's ALT_FREQS column, each set's OBS_CT (called allele count at the site)
    is carried through as an ``N:<set>`` column -- the finite-sample correction needs the
    per-site denominator, and dropping it is exactly the uncorrected-estimator bug.
    A flipped allele would invert the site's contribution, so mismatched sites are dropped
    and counted, never mixed. A site where any set has fewer than `min_alleles` called
    alleles is also dropped and counted: nA < 2 cannot estimate its sampling variance.
    Returns (value frame sorted by chrom/pos, n_allele_mismatch, n_low_allele_dropped).
    """
    keys = list(freqs)
    alt = {k: freqs[k].set_index("ID")["ALT_FREQS"] for k in keys}
    nct = {f"N:{k}": freqs[k].set_index("ID")["OBS_CT"] for k in keys}
    ra0 = freqs[keys[0]].set_index("ID")[["REF", "ALT"]]
    df = pd.DataFrame(alt).join(pd.DataFrame(nct)).join(ra0)
    bad = pd.Series(False, index=df.index)
    for k in keys[1:]:
        ra = freqs[k].set_index("ID")[["REF", "ALT"]]
        ok = ra.reindex(df.index)
        # A site absent from this source is a missing site (dropped by dropna below), not an
        # allele mismatch; only a present-but-different REF/ALT flips a site's contribution.
        miss = ok["REF"].isna() | ok["ALT"].isna()
        bad |= (~miss & ((df["REF"] != ok["REF"]) | (df["ALT"] != ok["ALT"]))).fillna(False)
    df = df.join(pos, how="inner")
    df = df[~bad.reindex(df.index, fill_value=False)].dropna(subset=keys)
    ncol = [c for c in df.columns if c.startswith("N:")]
    low = (df[ncol] < min_alleles).any(axis=1) if ncol else pd.Series(False, index=df.index)
    df = df[~low.fillna(True)]
    return (df.sort_values(["CHROM_N", "POS"]), int(bad.sum()), int(low.sum()))


def f3_matrix(df, outgroup, target, profiles, pools):
    """Jackknifed, finite-sample-corrected f3 rows for the panel.

    Per-site statistic (a is the f3 "first" population whose sampling variance is removed):
    (pa-pb)(pa-pc) - pa(1-pa)/(na-1), na = that population's called-allele count at the
    site (the ``N:<set>`` columns merge_freqs carries). E[(pahat-pbhat)(pahat-pchat)] =
    f3 + pa(1-pa)/na, so without the term a single-diploid first population is positively
    biased by its own sampling noise. outgroup may be None (ancient side: no outgroup group
    exists in the AADR panel -- pooled ancients are sources, not outgroups -- so only
    admixture rows are produced there). Ranking rows: f3(outgroup; target, set) per profiled
    set, ranked desc. Contrast rows: same-site difference between adjacent ranked pairs; the
    first population (and hence the correction) is common to both terms and cancels exactly.
    """
    skip = ("CHROM_N", "POS", "REF", "ALT")
    arr = {c: df[c].values for c in df.columns if c not in skip and not c.startswith("N:")}
    ncnt = {c[2:]: df[c].values.astype(float)
            for c in df.columns if c.startswith("N:") and c[2:] not in skip}
    bcode = block_codes(df.CHROM_N, df.POS)
    nblk = int(bcode.max()) + 1

    def f3(a, b, c):
        pa, na = arr[a], ncnt[a]
        stat = (pa - arr[b]) * (pa - arr[c]) - pa * (1 - pa) / (na - 1.0)
        return jackknife(stat, bcode, nblk)

    def row(t, **kw):
        return {"f3": t[0], "se": t[1], "z": t[2], **kw}

    ranked, contrasts, admixture = [], [], []
    if outgroup:
        vals = {p: f3(outgroup, target, p) for p in profiles}
        order = sorted(profiles, key=lambda p: -(vals[p][0] if vals[p][0] is not None else 0.0))
        ranked = [row(vals[p], set=p) for p in order]
        for p, q in zip(order, order[1:]):
            d = jackknife((arr[outgroup] - arr[target]) * (arr[outgroup] - arr[p])
                          - (arr[outgroup] - arr[target]) * (arr[outgroup] - arr[q]), bcode, nblk)
            contrasts.append(row(d, a=p, b=q))
    for i, a in enumerate(pools):
        for b in pools[i + 1:]:
            admixture.append(row(f3(target, a, b), a=a, b=b))
    return {"sites": len(df), "blocks": nblk, "ranked": ranked,
            "contrasts": contrasts, "admixture": admixture}


# ---------------------------------------------------------------- plink2 wrappers
def plink_freq(out_prefix, *, pfile=None, bfile=None, extract=None, keep=None):
    cmd = [PLINK2, "--freq", "--memory", "8000", "--threads", str(THREADS), "--out", str(out_prefix)]
    cmd += ["--pfile", str(pfile)] if pfile else ["--bfile", str(bfile)]
    if extract:
        cmd += ["--extract", str(extract)]
    if keep:
        cmd += ["--keep", str(keep)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    af = pathlib.Path(f"{out_prefix}.afreq")
    if r.returncode != 0 or not af.exists():
        sys.exit(f"28: plink2 --freq failed for {out_prefix.name}:\n{r.stderr[-2000:]}")
    return pd.read_csv(af, sep="\t")


def read_pos(pvar_path, extract_ids):
    # kg.common.pvar is VCF-styled (#CHROM POS ID REF ALT QUAL FILTER INFO): position is
    # column 2. Names match the selected columns exactly, so the mapping cannot drift if
    # the pvar gains extra trailing columns.
    pv = pd.read_csv(pvar_path, sep="\t", comment="#", header=None,
                     names=["CHROM", "POS", "ID"], usecols=[0, 1, 2], dtype={"CHROM": str})
    pv = pv[pv.ID.isin(extract_ids)].copy()
    pv["CHROM_N"] = pv.CHROM.str.replace("chr", "", regex=False)
    return pv.set_index("ID")[["CHROM_N", "POS"]]


def bim_pos(bim_path, extract_ids):
    bim = pd.read_csv(bim_path, sep="\t", header=None,
                      names=["CHROM", "ID", "POS"], usecols=[0, 1, 3], dtype={"CHROM": str})
    bim = bim[bim.ID.isin(extract_ids)].copy()
    bim["CHROM_N"] = bim.CHROM.str.replace("chr", "", regex=False)
    return bim.set_index("ID")[["CHROM_N", "POS"]]


# ---------------------------------------------------------------- modern block
def modern_block(groups):
    d04 = W / "04_ancestry"
    prune = d04 / "prune.prune.in"
    if not (prune.exists() and (d04 / "kg.common.psam").exists() and (d04 / "kg.common.pvar").exists()):
        print("28: modern inputs missing (04_ancestry kg.common/prune.prune.in); skipping modern f3",
              file=sys.stderr)
        return None
    tgt_prefix = next((p for p in (W / f"02_complete/{SAMPLE}.1kg", W / "02_complete/target.1kg")
                       if pathlib.Path(f"{p}.pgen").exists()), None)
    if tgt_prefix is None:
        print("28: no 02_complete target 1kg pfile; skipping modern f3", file=sys.stderr)
        return None

    outdir = d04 / "f3"
    outdir.mkdir(exist_ok=True)
    extract_ids = set(prune.read_text().split())
    ps = pd.read_csv(d04 / "kg.common.psam", sep="\t").rename(columns={"#FID": "FID", "#IID": "IID"})
    pops = set(ps.Population)
    missing = sorted({p for row in groups for p in row["members"]} - pops)
    if missing:
        sys.exit(f"28: panel/f3_groups.tsv names populations absent from kg.common.psam: {missing}")

    def afreq_for(row):
        members = ps[ps.Population.isin(row["members"])]
        (outdir / f"{row['set']}.ids").write_text(
            "".join(f"{r.FID}\t{r.IID}\n" for r in members.itertuples()))
        return plink_freq(outdir / row["set"], pfile=d04 / "kg.common",
                          extract=prune, keep=outdir / f"{row['set']}.ids"), len(members)

    freqs, meta = {}, {}
    by_kind = {k: [r for r in groups if r["kind"] == k] for k in ("outgroup", "profile", "pool")}
    if not by_kind["outgroup"] or not by_kind["profile"]:
        sys.exit("28: panel/f3_groups.tsv needs at least one outgroup and one profile row")
    for row in by_kind["outgroup"]:
        f, n = afreq_for(row)
        freqs[row["set"]] = f
        meta[row["set"]] = {**row, "n": n}
    for row in by_kind["profile"] + by_kind["pool"]:
        f, n = afreq_for(row)
        freqs[row["set"]] = f
        meta[row["set"]] = {**row, "n": n}
    freqs["target"] = plink_freq(outdir / "target", pfile=tgt_prefix, extract=prune)
    meta["target"] = {"kind": "target", "set": "target", "members": [],
                      "label_zh": "本样本", "label_en": "this sample", "n": 1}

    outgroup = by_kind["outgroup"][0]["set"]
    df, nbad, nlow = merge_freqs(freqs, read_pos(d04 / "kg.common.pvar", extract_ids))
    res = f3_matrix(df, outgroup, "target",
                    [r["set"] for r in by_kind["profile"]], [r["set"] for r in by_kind["pool"]])
    tsv = ([{"test": f"outgroup_f3({outgroup};target,{r['set']})", "f3": r["f3"], "SE": r["se"], "Z": r["z"]}
            for r in res["ranked"]]
           + [{"test": f"contrast_f3({c['a']}_minus_{c['b']})", "f3": c["f3"], "SE": c["se"], "Z": c["z"]}
              for c in res["contrasts"]]
           + [{"test": f"admixture_f3(target;{a['a']},{a['b']})", "f3": a["f3"], "SE": a["se"], "Z": a["z"]}
              for a in res["admixture"]])
    _write_tsv(outdir / "f3_results.tsv", tsv)
    print(f"28: modern f3 sites={res['sites']} blocks={res['blocks']} "
          f"allele_mismatch={nbad} low_allele_dropped={nlow}")
    return {
        "outgroup": outgroup, "sites": res["sites"], "blocks": res["blocks"],
        "allele_mismatch": nbad, "low_allele_dropped": nlow,
        "outgroup_f3": [{"set": r["set"],
                         "label_zh": meta[r["set"]]["label_zh"], "label_en": meta[r["set"]]["label_en"],
                         "kind": meta[r["set"]]["kind"], "n": meta[r["set"]]["n"],
                         "f3": _num(r["f3"]), "se": _num(r["se"]), "z": _num(r["z"], 2)}
                        for r in res["ranked"]],
        "contrasts": [{"a": c["a"], "b": c["b"], "diff": _num(c["f3"]), "z": _num(c["z"], 2)}
                      for c in res["contrasts"]],
        "admixture": [{"a": a["a"], "b": a["b"], "f3": _num(a["f3"]),
                       "se": _num(a["se"]), "z": _num(a["z"], 2)} for a in res["admixture"]],
        "groups": {s: {"kind": m["kind"], "label_zh": m["label_zh"], "label_en": m["label_en"], "n": m["n"]}
                   for s, m in meta.items()},
    }


# ---------------------------------------------------------------- ancient block
def regions_from_locations(loc):
    """source_id -> 'north'/'south' (else '') from the unified manual table's note token.

    地点表迁移（§3.2 P1 AN1/AN6）：旧的独立区域表已删除，南北区域判定的唯一出处是
    ancestry_locations.tsv 的 note `region=north|south`（秦岭—淮河界）。没有该 token 的行——
    非 China 前缀、跨区域聚合（region=unclassified）、1000G 采样地行——一律 ''，池分组本来就
    把 '' 当不可分类跳过，与旧表空 region 行为一致（迁移时 94 行逐值比对过）。
    """
    out = {}
    for sid, note in zip(loc["source_id"].astype(str), loc["note"].astype(str)):
        m = re.search(r"region=(north|south)\b", str(note))
        out[sid] = m.group(1) if m else ""
    return out


def ancient_block():
    d11 = W / "11_aadr"
    for f in ("aadr.bed", "aadr.bim", "aadr.fam", "prune.prune.in"):
        if not (d11 / f).exists():
            print(f"28: ancient input {f} missing (11_aadr); skipping ancient f3", file=sys.stderr)
            return None
    loc = pd.read_csv(PANEL / "ancestry_locations.tsv", sep="\t", comment="#", dtype=str)
    region_of = regions_from_locations(loc)
    fam = pd.read_csv(d11 / "aadr.fam", sep=r"\s+", header=None,
                      names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"], dtype=str)
    anc = fam[(fam.IID != SAMPLE) & (fam.FID != SAMPLE)]
    pools = {r: g for r, g in anc.assign(region=anc.FID.map(region_of).fillna(""))
             .groupby("region") if r and len(g) >= MIN_GROUP_N}
    if len(pools) < 2:
        print(f"28: fewer than two ancient region pools with n>={MIN_GROUP_N}; skipping ancient f3",
              file=sys.stderr)
        return None

    prune = d11 / "prune.prune.in"
    extract_ids = set(prune.read_text().split())
    freqs = {}
    for r, g in pools.items():
        (d11 / f"anc_{r}.ids").write_text("".join(f"{row.FID}\t{row.IID}\n" for row in g.itertuples()))
        freqs[r] = plink_freq(d11 / f"anc_{r}", bfile=d11 / "aadr", extract=prune,
                              keep=d11 / f"anc_{r}.ids")
    (d11 / "anc_target.ids").write_text(f"{SAMPLE}\t{SAMPLE}\n")
    freqs["target"] = plink_freq(d11 / "anc_target", bfile=d11 / "aadr", extract=prune,
                                 keep=d11 / "anc_target.ids")

    df, nbad, nlow = merge_freqs(freqs, bim_pos(d11 / "aadr.bim", extract_ids))
    # No outgroup exists in the AADR panel (pooled ancients are sources, not outgroups), so
    # only the pooled admixture rows are produced; a "ranking" against a pooled ancient would
    # not be an outgroup-f3 at all.
    res = f3_matrix(df, None, "target", [], sorted(pools))
    _write_tsv(d11 / "anc_f3_results.tsv",
               [{"test": f"admixture_f3(target;{a['a']},{a['b']})", "f3": a["f3"],
                 "SE": a["se"], "Z": a["z"]} for a in res["admixture"]])
    print(f"28: ancient f3 sites={res['sites']} blocks={res['blocks']} "
          f"allele_mismatch={nbad} low_allele_dropped={nlow}")
    zh = {"north": "北方古代池", "south": "南方古代池"}
    return {
        "sites": res["sites"], "blocks": res["blocks"],
        "allele_mismatch": nbad, "low_allele_dropped": nlow,
        "admixture": [{"a": a["a"], "b": a["b"],
                       "label_zh": f"{zh.get(a['a'], a['a'])}×{zh.get(a['b'], a['b'])}",
                       "label_en": f"{a['a']} x {a['b']}",
                       "n_a": int(len(pools[a["a"]])), "n_b": int(len(pools[a["b"]])),
                       "f3": _num(a["f3"]), "se": _num(a["se"]), "z": _num(a["z"], 2)}
                      for a in res["admixture"]],
    }


def _write_tsv(path, rows):
    # Unavailable rows (zero-SE etc.) write NA, not an empty cell that reads as 0 downstream.
    rows = [{k: ("NA" if v is None else v) for k, v in r.items()} for r in rows]
    pd.DataFrame(rows, columns=["test", "f3", "SE", "Z"]).to_csv(
        path, sep="\t", index=False, float_format="%.6f")


def main():
    _cfg()
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import ancestry_data as ad
    outdir = W / "04_ancestry" / "f3"
    # 失效先于计算（复审 §3.2 P0 失败生命周期）：真实复现里 28 因参考群缺失在下方
    # sys.exit(1)，旧 manifest=ok 与旧 f3_stats.json 原样保留、30 照旧准入。开工先写
    # 失效记录，成功后的最终 manifest 原子覆盖。
    ad.begin_run_manifest(outdir / "manifest.json", SAMPLE, "28-f3-stats")
    groups = pd.read_csv(PANEL / "f3_groups.tsv", sep="\t", comment="#", dtype=str)
    groups["members"] = groups.members.str.split()
    glist = groups.to_dict("records")
    out = {"sample_id": SAMPLE, "estimator": ESTIMATOR_ID, "block_mb": BLOCK_MB // 1_000_000,
           "modern": modern_block(glist), "ancient": ancient_block()}
    if out["modern"] is None and out["ancient"] is None:
        ad.write_manifest(outdir / "manifest.json", ad.build_manifest(
            SAMPLE, "28-f3-stats", state="unavailable", reason_code="no-inputs",
            detail="neither modern (04_ancestry) nor ancient (11_aadr) inputs are present"))
        sys.exit("28: neither modern nor ancient inputs are present; nothing to compute")
    dst = outdir / "f3_stats.json"
    # allow_nan=False: a leaked NaN/Infinity must fail the build here, not ship invalid JSON
    # that renders as "null" and silently reads as a number downstream.
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1, allow_nan=False), encoding="utf-8")
    outputs = [str(dst.relative_to(W))]
    for extra in (outdir / "f3_results.tsv", W / "11_aadr" / "anc_f3_results.tsv"):
        if extra.exists():
            outputs.append(str(extra.relative_to(W)))
    ad.write_manifest(outdir / "manifest.json", ad.build_manifest(
        SAMPLE, "28-f3-stats", state="ok",
        parameters={"estimator": ESTIMATOR_ID, "block_mb": BLOCK_MB // 1_000_000,
                    "min_group_n": MIN_GROUP_N,
                    # 复审 §3.2 P0 指纹绑定：分组面板与 modern 侧 prune 集的内容指纹——换
                    # f3_groups 或换 04 的修剪集后，旧 f3_stats.json 不得继续冒充本次结果。
                    "panel_sha": ad.file_sha(PANEL / "f3_groups.tsv"),
                    "prune_sha": ad.file_sha(W / "04_ancestry" / "prune.prune.in")},
        outputs=outputs,
        sides=[k for k in ("modern", "ancient") if out[k] is not None]))
    print("28: wrote", dst)


if __name__ == "__main__":
    main()
