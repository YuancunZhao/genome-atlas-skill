#!/usr/bin/env python
"""28_f3_stats.py -- f3 statistics on the pruned-SNP panels (outgroup-f3 ranking, adjacent
contrasts, pooled admixture-f3), for both the modern 1000G side and the ancient AADR side.

Inputs already produced by earlier steps:
  modern:  04_ancestry/kg.common.{psam,pvar} + 04_ancestry/prune.prune.in
           + 02_complete/{SAMPLE}.1kg / target.1kg (the target's 1kg-space genotypes)
  ancient: 11_aadr/aadr.{bed,bim,fam} + 11_aadr/prune.prune.in
           + panel/aadr_site_regions.tsv (fam FID -> region north/south)

Group definitions live in panel/f3_groups.tsv (modern) and the region column of
aadr_site_regions.tsv (ancient); nothing sample-specific is written here.

Writes:
  04_ancestry/f3/f3_results.tsv    modern tests (one row per test, f3/SE/Z)
  11_aadr/anc_f3_results.tsv       ancient tests
  04_ancestry/f3/f3_stats.json     combined payload 30_build_report_data.py reads

f3(A;B,C) = unweighted per-site mean of (pA-pB)(pA-pC); SE/Z from a delete-one-block
jackknife (5 Mb blocks), the same estimator AdmixTools uses. A significantly negative
admixture-f3 is the signature of A descending from a mix of B and C; outgroup-f3 rises
with shared drift between the target and the profiled group.
"""
import json
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd

PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel"
BLOCK_MB = 5_000_000
MIN_GROUP_N = 20          # an ancient region pool below this is too small to pool quietly


def _cfg():
    """Bind the run configuration lazily so the pure helpers stay importable in tests
    (importing wgsconfig creates the work directories -- a side effect tests must not cause)."""
    global W, PLINK2, SAMPLE, THREADS
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import wgsconfig
    W, PLINK2, SAMPLE, THREADS = wgsconfig.W, wgsconfig.PLINK2, wgsconfig.SAMPLE, wgsconfig.THREADS


# ---------------------------------------------------------------- pure helpers (tested)
def jackknife(x, bcode, nblk):
    """Delete-one-block mean and SE of the per-site statistic x (one value per site)."""
    S = x.sum()
    Sb = np.bincount(bcode, weights=x, minlength=nblk)
    Nb = np.bincount(bcode, minlength=nblk).astype(float)
    theta = S / len(x)
    loo = (S - Sb) / (len(x) - Nb)
    se = np.sqrt((nblk - 1) / nblk * ((loo - loo.mean()) ** 2).sum())
    return float(theta), float(se), float(theta / se)


def block_codes(chroms, positions, block=BLOCK_MB):
    """Site -> 0-based jackknife block id (chrom + floor(pos/block))."""
    chroms = pd.Series(chroms).astype(str).str.replace("chr", "", regex=False)
    return pd.factorize((chroms + "_" + (pd.Series(positions) // block).astype(str)).values)[0]


def merge_freqs(freqs, pos):
    """Join per-set afreq tables on variant ID, validating REF/ALT across sources.

    `freqs` maps set name -> raw plink2 .afreq DataFrame; `pos` is ID -> (CHROM_N, POS).
    Returns (value frame sorted by chrom/pos, n_allele_mismatch). A flipped allele would
    invert the site's contribution, so mismatched sites are dropped and counted, never mixed.
    """
    keys = list(freqs)
    alt = {k: freqs[k].set_index("ID")["ALT_FREQS"] for k in keys}
    ra0 = freqs[keys[0]].set_index("ID")[["REF", "ALT"]]
    df = pd.DataFrame(alt).join(ra0)
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
    return df.sort_values(["CHROM_N", "POS"]), int(bad.sum())


def f3_matrix(df, outgroup, target, profiles, pools):
    """Jackknifed f3 rows for the panel.

    outgroup may be None (ancient side: no outgroup group exists in the AADR panel -- pooled
    ancients are sources, not outgroups -- so only admixture rows are produced there).
    Ranking rows: f3(outgroup; target, set) per profiled set, ranked desc. Contrast rows:
    same-site difference between adjacent ranked pairs; the target's own drift error cancels
    in the difference, which is why the ranking quotes these Z values, not per-row SEs.
    """
    arr = {c: df[c].values for c in df.columns if c not in ("CHROM_N", "POS", "REF", "ALT")}
    bcode = block_codes(df.CHROM_N, df.POS)
    nblk = int(bcode.max()) + 1

    def f3(a, b, c):
        return jackknife((arr[a] - arr[b]) * (arr[a] - arr[c]), bcode, nblk)

    ranked, contrasts, admixture = [], [], []
    if outgroup:
        vals = {p: f3(outgroup, target, p) for p in profiles}
        order = sorted(profiles, key=lambda p: -vals[p][0])
        ranked = [{"set": p, "f3": vals[p][0], "se": vals[p][1], "z": vals[p][2]} for p in order]
        for p, q in zip(order, order[1:]):
            d = jackknife((arr[outgroup] - arr[target]) * (arr[outgroup] - arr[p])
                          - (arr[outgroup] - arr[target]) * (arr[outgroup] - arr[q]), bcode, nblk)
            contrasts.append({"a": p, "b": q, "f3": d[0], "se": d[1], "z": d[2]})
    for i, a in enumerate(pools):
        for b in pools[i + 1:]:
            v = f3(target, a, b)
            admixture.append({"a": a, "b": b, "f3": v[0], "se": v[1], "z": v[2]})
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
    df, nbad = merge_freqs(freqs, read_pos(d04 / "kg.common.pvar", extract_ids))
    res = f3_matrix(df, outgroup, "target",
                    [r["set"] for r in by_kind["profile"]], [r["set"] for r in by_kind["pool"]])
    tsv = ([{"test": f"outgroup_f3({outgroup};target,{r['set']})", "f3": r["f3"], "SE": r["se"], "Z": r["z"]}
            for r in res["ranked"]]
           + [{"test": f"contrast_f3({c['a']}_minus_{c['b']})", "f3": c["f3"], "SE": c["se"], "Z": c["z"]}
              for c in res["contrasts"]]
           + [{"test": f"admixture_f3(target;{a['a']},{a['b']})", "f3": a["f3"], "SE": a["se"], "Z": a["z"]}
              for a in res["admixture"]])
    _write_tsv(outdir / "f3_results.tsv", tsv)
    print(f"28: modern f3 sites={res['sites']} blocks={res['blocks']} allele_mismatch={nbad}")
    return {
        "outgroup": outgroup, "sites": res["sites"], "blocks": res["blocks"],
        "allele_mismatch": nbad,
        "outgroup_f3": [{"set": r["set"],
                         "label_zh": meta[r["set"]]["label_zh"], "label_en": meta[r["set"]]["label_en"],
                         "kind": meta[r["set"]]["kind"], "n": meta[r["set"]]["n"],
                         "f3": round(r["f3"], 6), "se": round(r["se"], 6), "z": round(r["z"], 2)}
                        for r in res["ranked"]],
        "contrasts": [{"a": c["a"], "b": c["b"], "diff": round(c["f3"], 6), "z": round(c["z"], 2)}
                      for c in res["contrasts"]],
        "admixture": [{"a": a["a"], "b": a["b"], "f3": round(a["f3"], 6),
                       "se": round(a["se"], 6), "z": round(a["z"], 2)} for a in res["admixture"]],
        "groups": {s: {"kind": m["kind"], "label_zh": m["label_zh"], "label_en": m["label_en"], "n": m["n"]}
                   for s, m in meta.items()},
    }


# ---------------------------------------------------------------- ancient block
def ancient_block():
    d11 = W / "11_aadr"
    for f in ("aadr.bed", "aadr.bim", "aadr.fam", "prune.prune.in"):
        if not (d11 / f).exists():
            print(f"28: ancient input {f} missing (11_aadr); skipping ancient f3", file=sys.stderr)
            return None
    reg = pd.read_csv(PANEL / "aadr_site_regions.tsv", sep="\t", comment="#", dtype=str)
    region_of = dict(zip(reg.label, reg.region.fillna("")))
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

    df, nbad = merge_freqs(freqs, bim_pos(d11 / "aadr.bim", extract_ids))
    # No outgroup exists in the AADR panel (pooled ancients are sources, not outgroups), so
    # only the pooled admixture rows are produced; a "ranking" against a pooled ancient would
    # not be an outgroup-f3 at all.
    res = f3_matrix(df, None, "target", [], sorted(pools))
    _write_tsv(d11 / "anc_f3_results.tsv",
               [{"test": f"admixture_f3(target;{a['a']},{a['b']})", "f3": a["f3"],
                 "SE": a["se"], "Z": a["z"]} for a in res["admixture"]])
    print(f"28: ancient f3 sites={res['sites']} blocks={res['blocks']} allele_mismatch={nbad}")
    zh = {"north": "北方古代池", "south": "南方古代池"}
    return {
        "sites": res["sites"], "blocks": res["blocks"], "allele_mismatch": nbad,
        "admixture": [{"a": a["a"], "b": a["b"],
                       "label_zh": f"{zh.get(a['a'], a['a'])}×{zh.get(a['b'], a['b'])}",
                       "label_en": f"{a['a']} x {a['b']}",
                       "n_a": int(len(pools[a["a"]])), "n_b": int(len(pools[a["b"]])),
                       "f3": round(a["f3"], 6), "se": round(a["se"], 6), "z": round(a["z"], 2)}
                      for a in res["admixture"]],
    }


def _write_tsv(path, rows):
    pd.DataFrame(rows, columns=["test", "f3", "SE", "Z"]).to_csv(
        path, sep="\t", index=False, float_format="%.6f")


def main():
    _cfg()
    groups = pd.read_csv(PANEL / "f3_groups.tsv", sep="\t", comment="#", dtype=str)
    groups["members"] = groups.members.str.split()
    glist = groups.to_dict("records")
    out = {"block_mb": BLOCK_MB // 1_000_000,
           "modern": modern_block(glist), "ancient": ancient_block()}
    if out["modern"] is None and out["ancient"] is None:
        sys.exit("28: neither modern nor ancient inputs are present; nothing to compute")
    (W / "04_ancestry" / "f3").mkdir(parents=True, exist_ok=True)
    dst = W / "04_ancestry" / "f3" / "f3_stats.json"
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print("28: wrote", dst)


if __name__ == "__main__":
    main()
