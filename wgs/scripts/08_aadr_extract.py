#!/usr/bin/env python
"""Extract a chosen set of AADR Human Origins individuals (TGENO packed format) plus {NAME_EN}'s genotypes
at the same SNPs, and write a plink1 .bed/.bim/.fam set for PCA."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import os, re, sys, io, json, subprocess, numpy as np, pandas as pd
import ancestry_data as ad
P=str(P); A=str(pathlib.Path(AADR).parent); W=f"{P}/wgs/11_aadr"; os.makedirs(W,exist_ok=True)
# The panel location comes from the configuration, not from a hard-coded v66 filename.
PREF=str(AADR)
ANNO = AADR_ANNOTATION or f"{PREF}.anno"
if not pathlib.Path(ANNO).exists():
    _alt = PREF.replace(".patch.PUB", ".PUB") + ".anno"
    ANNO = _alt if pathlib.Path(_alt).exists() else ANNO
if not pathlib.Path(ANNO).exists():
    sys.exit(f"no annotation file for aadr_prefix={PREF!r}: tried {ANNO} -- refusing to guess")
for _ext in (".geno", ".snp", ".ind"):
    if not pathlib.Path(PREF + _ext).exists():
        sys.exit(f"incomplete AADR panel: {PREF}{_ext} is missing")
RELEASE = (re.match(r"(v\d+(?:\.\d+)?)", pathlib.Path(PREF).name) or [None, "unknown"])[1]

def _fail(reason, detail):
    # 复审 AN1：提取失败不是一行日志。bcftools 可以先吐一条合法 stdout 再 exit 1，旧代码
    # 照常把 manifest 写成 ok，09b 就拿着半份基因型出结论。失败必须落 failed manifest
    # （原子覆盖上一次运行留下的 ok）并以非零码退出。
    ad.write_manifest(f"{W}/manifest.json", ad.build_manifest(
        SAMPLE, "08-aadr-extract", state="failed", reason_code=reason, reference_release=RELEASE,
        parameters={"aadr_prefix": PREF, "detail": str(detail)}, outputs=[]))
    sys.exit(f"08-aadr-extract failed ({reason}): {detail}")

def validate_tgeno_header(path, n_ind, n_snp):
    """返回 None 表示 .geno 的 TGENO 头与 .ind/.snp 维度一致；否则返回原因字符串。

    TGENO 头是 48 字节：'TGENO' 魔数后跟空格分隔的十进制 n_ind、n_snp（个体优先布局，
    每行 rlen=(n_snp+3)//4 字节）。旧代码只断言文件总尺寸——任何等长的损坏头都放行，
    而维度对不上的读取就是静默把两个人的基因型错位拼在一起。"""
    try:
        with open(path, "rb") as f:
            raw = f.read(48)
    except OSError as e:
        return f"unreadable .geno: {e}"
    if len(raw) < 48:
        return f"truncated header: {len(raw)} bytes"
    if not raw.startswith(b"TGENO"):
        return "bad magic: expected b'TGENO'"
    fields = raw.rstrip(b"\x00").decode("ascii", "replace").split()
    if len(fields) < 3:
        return f"malformed header fields: {fields!r}"
    try:
        h_ind, h_snp = int(fields[1]), int(fields[2])
    except ValueError:
        return f"non-numeric header dimensions: {fields[1:3]!r}"
    if h_ind != n_ind:
        return f"header n_ind {h_ind} != .ind rows {n_ind}"
    if h_snp != n_snp:
        return f"header n_snp {h_snp} != .snp rows {n_snp}"
    rlen = (n_snp + 3) // 4
    size = os.path.getsize(path)
    if 48 + rlen * n_ind != size:
        return f"size {size} != 48 + {rlen} * {n_ind}"
    return None

# Original geno row numbers, captured before any merge: reading packed genotypes by a post-merge
# DataFrame index is how a merge that reorders rows silently mixes individuals up.
ind=pd.read_csv(f"{PREF}.ind",sep=r"\s+",header=None,names=["iid","sex","pop"])
ind["geno_row"]=np.arange(len(ind),dtype=np.int64)
snp=pd.read_csv(f"{PREF}.snp",sep=r"\s+",header=None,names=["rsid","chrom","cm","pos","a1","a2"],dtype={"chrom":str})
n_ind,n_snp=len(ind),len(snp); rlen=(n_snp+3)//4
# 复审 AN1：断言不能代替输入校验——校验不过要带着原因落 failed manifest 退出，不是 AssertionError 裸栈。
_bad_tgeno = validate_tgeno_header(f"{PREF}.geno", n_ind, n_snp)
if _bad_tgeno:
    _fail("tgeno_header_invalid", _bad_tgeno)

# Panel membership comes from the configuration (aadr_modern / aadr_ancient_prefix); the age of a
# record never decides its kind, and no 500 BP cutoff is applied. Everything the panel does not name
# stays unknown and out of both PCA sets.
anno=pd.read_csv(ANNO,sep="\t",dtype=str,low_memory=False)
anno.columns=[c.strip() for c in anno.columns]
recs=ad.assign_kind(ad.normalize_metadata(anno.to_dict("records"),dataset="AADR",release=RELEASE),
                    ancient_prefixes=AADR_ANCIENT_PREFIX,modern_groups=AADR_MODERN)
# 复审 AN1：配置点名但面板里一条记录都没有的条目要点名道姓——一个拼错的群体名会静默掏空
# 那一侧的 PCA 集合，分析照常出结果但训练集合悄悄变了。计数行不区分"没配置"与"没命中"。
_um=ad.unmatched_panel_entries(recs,AADR_MODERN,AADR_ANCIENT_PREFIX)
if _um["modern"] or _um["ancient_prefix"]:
    print("warning: configured panel entries matched no records: "
          f"modern={_um['modern']} ancient_prefix={_um['ancient_prefix']}",file=sys.stderr)
if not any(r.get("kind")=="modern" for r in recs):
    print("warning: no modern references selected -- the modern PCA/projection side is empty",file=sys.stderr)
if not any(r.get("kind")=="ancient" for r in recs):
    print("warning: no ancient references selected -- the ancient projection side is empty",file=sys.stderr)
# 复审 AN1：同一个人可以有多种技术表示（不同 call 版本、重复记录）。之前只在测试里有去重，生产路径
# 不去重，于是同一个人可能以两条记录各自计入分组与计数。这里接上：保留可用的那一条，被丢弃的带原因
# 落盘供审计，而不是悄悄消失。
_dd = ad.dedupe_by_master_id(recs)
recs = _dd["kept"]
if _dd["dropped"]:
    pathlib.Path(f"{W}/dedup_dropped.tsv").write_text(
        "record_id\tmaster_id\treason_code\tdetail\n" + "".join(
            f"{d['record_id']}\t{d['master_id']}\t{d['reason_code']}\t{d['detail']}\n" for d in _dd["dropped"]),
        encoding="utf-8")
    print(f"dedupe: kept {len(recs)}, dropped {len(_dd['dropped'])} duplicate representation(s) "
          f"-> {W}/dedup_dropped.tsv", file=sys.stderr)

# .ind and .anno are both keyed by the Genetic ID (the .ind's own "pop" column is a patch artefact).
by_record={r["record_id"]: r for r in recs}
def _f(iid, field, default=None):
    r=by_record.get(str(iid))
    return r.get(field) if r else default
ind["source_population_id"]=[_f(i,"source_population_id") for i in ind.iid]
ind["kind"]=[_f(i,"kind","unknown") for i in ind.iid]
ind["locality"]=[_f(i,"locality") for i in ind.iid]
ind["date_mean_bp"]=[_f(i,"date_mean_bp") for i in ind.iid]
bad=ind["source_population_id"].fillna("").str.contains("Ignore|QCremove|DontUse|_dup|Discovery|_o$|-o$|outlier|contam|_lc$|_1d|_2d",case=False,na=False)
keep=ind[(ind.kind.isin(["modern","ancient"]))&~bad].copy()
keep["pop"]=keep["source_population_id"]
keep["label"]=keep["pop"]
# Every .anno record goes out for the history views -- not just the ones this PCA happens to use.
# Genotypes are only extracted for the names the panel selected.
pd.DataFrame(recs).to_csv(f"{W}/reference_metadata.tsv",sep="\t",index=False)
_selected=set(keep.iid.astype(str))
print(f"panel selection: {len(_selected)} of {len(recs)} annotated records "
      f"(modern {(keep.kind=='modern').sum()}, ancient {(keep.kind=='ancient').sum()}, "
      f"unknown {len(recs)-len(_selected)})",file=sys.stderr)
print(f"keeping {len(keep)} individuals: modern {(keep.kind=='modern').sum()}, ancient {(keep.kind=='ancient').sum()}",file=sys.stderr)
print(keep[keep.kind=="modern"].label.value_counts().head(20).to_dict(),file=sys.stderr)

# --- SNP selection: autosomal, biallelic ACGT, present in {NAME_EN}'s complete set
snp["idx"]=np.arange(n_snp)
sel=snp[(snp.chrom.isin([str(i) for i in range(1,23)]))&snp.a1.isin(list("ACGT"))&snp.a2.isin(list("ACGT"))].copy()
sel[["chrom","pos"]].to_csv(f"{W}/ho_sites.tsv",sep="\t",header=False,index=False)
q=subprocess.run(["bcftools","query","-R",f"{W}/ho_sites.tsv","-f","%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n",f"{P}/wgs/02_complete/{SAMPLE}.1kg_sites.vcf.gz"],capture_output=True,text=True)
# 复审 AN1：子进程失败立即终止并使本次状态失效。exit 1 就是失败——哪怕 stdout 里已有合法行
# （半份输出配上 ok manifest 会让 09b 拿截断的基因型照常出结论）；exit 0 但一行都没有同样
# 是输入坏了（错的 VCF/错的 build），不是"目标恰好全缺失"。
if q.returncode:
    _fail("bcftools_query_failed", f"bcftools exit {q.returncode}: {q.stderr.strip()[-500:]}")
if not q.stdout.strip():
    _fail("bcftools_query_empty", "bcftools exited 0 but returned no rows for ho_sites.tsv")
q=q.stdout
dg=pd.read_csv(io.StringIO(q),sep="\t",header=None,names=["chrom","pos","ref","alt","gt"],dtype={"chrom":str}).drop_duplicates(["chrom","pos"])
# LEFT join: a site the target cannot be called at keeps its reference individuals and is encoded as
# missing for the target. An inner join here would drop training SNPs because of the target's gaps,
# and the remaining sites would silently differ between runs and references.
sel=sel.merge(dg,on=["chrom","pos"],how="left")
_ok=((sel.a1==sel.ref)&(sel.a2==sel.alt))|((sel.a1==sel.alt)&(sel.a2==sel.ref))
_target_gt=sel["gt"]
sel["dayu_a1"]=[ad.a1_dosage(g,r,a,al) for g,r,a,al in zip(_target_gt.fillna("."),sel["ref"],sel["alt"],sel["a1"])]
sel.loc[~_ok,"dayu_a1"]=None
_usable=int(sel["dayu_a1"].notna().sum())
print(f"SNPs kept: {len(sel)} of {n_snp}; target called at {_usable} "
      f"({len(sel)-_usable} encoded missing)",file=sys.stderr)
if not len(sel):
    sys.exit("no reference SNPs survived filtering; refusing to write an empty panel")

# --- read packed rows for the kept individuals
G=np.zeros((len(keep),len(sel)),dtype=np.int8)
cols=sel.idx.values
with open(f"{PREF}.geno","rb") as f:
    for k,(row_i) in enumerate(keep.index.values):
        f.seek(48+rlen*row_i); buf=np.frombuffer(f.read(rlen),dtype=np.uint8)
        bits=np.unpackbits(buf).reshape(-1,2); vals=(bits[:,0]*2+bits[:,1])[:n_snp]
        v=vals[cols].astype(np.int8); v[v==3]=-1  # 3 = missing
        G[k]=v
        if k%200==0: print(f"  {k}/{len(keep)}",file=sys.stderr,flush=True)
miss=(G<0).mean(axis=1)
keep["call_rate"]=1-miss
keep2=keep[keep.call_rate>=0.3].copy(); G=G[keep.call_rate.values>=0.3]
print(f"after call-rate filter: {len(keep2)} individuals",file=sys.stderr)

# --- write plink1 bed (SNP-major): 00=hom a1a1? plink codes: 00 hom A1, 01 missing, 10 het, 11 hom A2
# We set A1 = eigen a1 (counted allele). geno value g = count of a1. plink 2-bit: g==2 -> 00, g==1 -> 10, g==0 -> 11, missing -> 01
# pandas 把缺失值读成 NaN，`v is None` 对 NaN 为 False，随后 int(nan) 抛 ValueError 带走整步。
# 用 pd.isna 同时覆盖 None 与 NaN —— 真实 AADR 的 GT 列本来就有缺失（女性和性别未知样本被标成 n/a）。
_tgt=np.array([(-1 if pd.isna(v) else int(v)) for v in sel["dayu_a1"]],dtype=np.int8)
allg=np.vstack([G,_tgt[None,:]])
# 复审 AN1：目标的 fid/iid 用**稳定的 sample_id**，不用显示名。显示名会随 name_en 配置改变，而下游按
# 名字匹配——改一次显示名就等于把目标换成了另一个人（身份随文案漂移）。sample_id 是样本的标识。
_TGT_ID = SAMPLE
labels=list(keep2.label)+[_TGT_ID]; kinds=list(keep2.kind)+["target"]; iids=list(keep2.iid)+[_TGT_ID]
code=np.select([allg==2,allg==1,allg==0],[0b00,0b10,0b11],default=0b01).astype(np.uint8)
n=len(labels); pad=(4-n%4)%4
with open(f"{W}/aadr.bed","wb") as f:
    f.write(bytes([0x6c,0x1b,0x01]))
    for j in range(allg.shape[1]):
        col=np.concatenate([code[:,j],np.zeros(pad,dtype=np.uint8)]).reshape(-1,4)
        f.write((col[:,0]|(col[:,1]<<2)|(col[:,2]<<4)|(col[:,3]<<6)).astype(np.uint8).tobytes())
pd.DataFrame({"chrom":sel.chrom,"rsid":sel.rsid,"cm":sel.cm,"pos":sel.pos,"a1":sel.a1,"a2":sel.a2}).to_csv(f"{W}/aadr.bim",sep="\t",header=False,index=False)
pd.DataFrame({"fid":labels,"iid":iids,"pat":0,"mat":0,"sex":0,"phe":-9}).to_csv(f"{W}/aadr.fam",sep=" ",header=False,index=False)
# The target's call rate is left empty on purpose: it is only meaningful once measured on the final
# pruned site set, and writing 1.0 here made an unchecked number look like an assayed one.
pd.DataFrame({"iid":iids,"label":labels,"kind":kinds,"date":list(keep2.date_mean_bp)+[None],
              "call_rate":list(keep2.call_rate)+[None],"geno_row":list(keep2.geno_row)+[-1]}).to_csv(f"{W}/samples.tsv",sep="\t",index=False)
ad.write_manifest(f"{W}/manifest.json", ad.build_manifest(
    SAMPLE,"08-aadr-extract",state="ok",reference_release=RELEASE,
    parameters={"aadr_prefix":PREF,"annotation":ANNO,"modern":AADR_MODERN,"ancient_prefix":AADR_ANCIENT_PREFIX,
                "unmatched_modern":_um["modern"],"unmatched_ancient_prefix":_um["ancient_prefix"],
                "min_call_rate_ancient":MIN_CR_ANCIENT,"min_projection_snps":MIN_PROJECTION_SNPS},
    input_fingerprints={"geno_size":os.path.getsize(f"{PREF}.geno"),"n_snp":int(n_snp),"n_ind":int(n_ind)},
    outputs=["11_aadr/reference_metadata.tsv","11_aadr/samples.tsv","11_aadr/aadr.bed","11_aadr/aadr.bim","11_aadr/aadr.fam"]))
print("wrote",W,"n_ind",n,"n_snp",allg.shape[1],"release",RELEASE)
