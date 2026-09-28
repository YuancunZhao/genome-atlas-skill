#!/usr/bin/env python
"""Genome-wide north/south ancestry proportions for {NAME_EN}, expressed against held-out CHB and CHS individuals
scored with the same panels (FLARE posterior global ancestry, weighted by chromosome length)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...
import ancestry_data as _ad  # chromosome_lengths: 优先读 .fai，回退内置表并告警

import glob, gzip, pandas as pd, numpy as np, os
P=str(P); W=f"{P}/wgs/12_localanc"
# 仅在读不到任何**同 build**索引时使用；它只对 hg19/GRCh37 正确。候选必须声明各自索引的 build：
# b37 缺席的机器上顺着候选读到 hg38 的 .fai 是跨 build 错误，chromosome_lengths 会整条跳过它。
_HG19=[249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]
LEN=_ad.chromosome_lengths([(CFG.get('reference_fai') or '', BUILD),
                            (REF/'b37/human_g1k_v37.fasta.fai', 'GRCh37'),
                            (REF/'hg38/Homo_sapiens_assembly38.fasta.fai', 'GRCh38')],
                     build=BUILD,
                     fallback={str(i): l for i, l in zip(range(1, 23), _HG19)})
rows=[]
for f in glob.glob(f"{W}/la.*.global.anc.gz"):
    c=os.path.basename(f).split(".")[1]
    d=pd.read_csv(f,sep="\t"); d["chrom"]=c; rows.append(d)
dy=pd.concat(rows); dy["w"]=dy.chrom.map(LEN)
# Labels come from the FLARE .model, in the same order as the data columns (AN3): a hard-coded list
# would attach percentages to the wrong panel as soon as the configured panels change.
import glob as _glob
def _an_labels(work):
    for m in sorted(_glob.glob(f"{work}/la.*.model")) + sorted(_glob.glob(f"{work}/calib.*.model")):
        try:
            lines = open(m, encoding="utf-8", errors="replace").read().split("\n")
        except OSError:
            continue
        for i, l in enumerate(lines):
            if l.strip().lower().startswith("# list of ancestries"):
                for j in range(i + 1, min(i + 4, len(lines))):
                    row = lines[j].strip()
                    if row and not row.startswith("#"):
                        return row.split("\t")
    raise SystemExit("no FLARE .model with an ancestry list; cannot label the columns")

anc = _an_labels(W)
# §7：比较的是**配置里的两个来源面板**，不是 FLARE 的列顺序。列顺序由参考面板首次出现的次序决定
# （这里是 European/SouthEA/SouthAsian/NorthEA），拿 anc[0]/anc[1] 当"北方/南方"会把欧洲面板
# 当成来源面板——数值看着合理，含义完全错了。
_PANELS = [str(x) for x in (OPT.get("la_labels") or [])]
if len(_PANELS) != 2:
    raise SystemExit("local_ancestry_labels must name the two source panels to compare")
A_COL, B_COL = _PANELS[0], _PANELS[1]
for _c in (A_COL, B_COL):
    if _c not in anc and _c not in dy.columns and False:
        raise SystemExit(f"{_c} is not among the calibration columns {anc}")
gw={a:np.average(dy[a],weights=dy.w) for a in anc}
print(f"{NAME_EN}, length-weighted over 22 autosomes:")
for a in anc: print(f"  {a:12s} {gw[a]*100:5.2f}%")
ps=pd.read_csv(KG_PFILE+".psam",sep="\t").rename(columns={"#IID":"SAMPLE","IID":"SAMPLE"})
cal=[]
for f in glob.glob(f"{W}/calib.*.global.anc.gz"):
    c=os.path.basename(f).split(".")[1]
    d=pd.read_csv(f,sep="\t"); d["chrom"]=c; cal.append(d)
if cal:
    cd=pd.concat(cal).merge(ps[["SAMPLE","Population"]],on="SAMPLE"); cd["w"]=cd.chrom.map(LEN)
    print(f"\ncalibration on chromosomes {sorted(set(cd.chrom))} (held-out reference individuals):")
    out=[]
    for pop,g in cd.groupby("Population"):
        per=g.groupby("SAMPLE").apply(lambda x: pd.Series({a:np.average(x[a],weights=x.w) for a in anc}),include_groups=False)
        out.append((pop, len(per), per[A_COL].mean(), per[A_COL].std(), per[B_COL].mean(), per[B_COL].std()))
        print(f"  {pop}: n={len(per)}  {A_COL} {per[anc[0]].mean()*100:.1f}% (sd {per[anc[0]].std()*100:.1f})  {B_COL} {per[B_COL].mean()*100:.1f}% (sd {per[B_COL].std()*100:.1f})")
    # {NAME_EN} restricted to the same chromosomes for a fair comparison
    same=dy[dy.chrom.isin(set(cd.chrom))]
    dn=np.average(same[A_COL],weights=same.w); ds=np.average(same[B_COL],weights=same.w)
    print(f"  {NAME_EN} (same chromosomes): {A_COL} {dn*100:.1f}%  {B_COL} {ds*100:.1f}%")
    for pop,n,m,s,ms,ss in out:
        print(f"    vs {pop}: {NAME_EN} is {(dn-m)/s:+.1f} sd on {A_COL}")
    pd.DataFrame(out,columns=["pop","n","north_mean","north_sd","south_mean","south_sd"]).to_csv(f"{W}/calibration.tsv",sep="\t",index=False)
    # §7/AN6：把校准结果写回结构化结果，报告才能用**真实的**群体、人数与染色体，而不是写死 CHB/CHS
    # 与 1/2/6/22。这里同时记下实际参与汇总的染色体集合。
    import json as _json
    _lj = pathlib.Path(f"{W}/local_ancestry.json")
    if _lj.exists():
        try:
            _d = _json.loads(_lj.read_text(encoding="utf-8"))
            _d["calibration"] = [{"population": p, "n": n, "north_mean": m, "north_sd": sd,
                                  "south_mean": ms, "south_sd": ss} for p, n, m, sd, ms, ss in out]
            _d["calibration_chroms"] = sorted(set(str(c) for c in cd.chrom), key=lambda x: int(x))
            _d["calibration_panels"] = [A_COL, B_COL]   # 参与比较的两个来源面板（来自配置）
            _lj.write_text(_json.dumps(_d, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
            print(f"wrote calibration for {len(out)} population(s) over chromosomes "
                  f"{_d['calibration_chroms']} into local_ancestry.json")
        except (OSError, _json.JSONDecodeError) as e:
            print(f"warning: could not attach calibration to local_ancestry.json ({e})", file=sys.stderr)
pd.Series(gw).to_csv(f"{W}/dayu_global.tsv",sep="\t",header=False)
dy.to_csv(f"{W}/per_chrom.tsv",sep="\t",index=False)
print(f"\nper-chromosome {B_COL} share ({NAME_EN}):")
print(dy.sort_values("chrom",key=lambda s:s.astype(int))[["chrom",A_COL,B_COL]].round(3).to_string(index=False))
