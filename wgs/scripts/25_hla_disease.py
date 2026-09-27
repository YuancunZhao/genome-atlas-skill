#!/usr/bin/env python
"""Disease and drug associations of the HLA types called by T1K, checked allele by allele."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import pandas as pd, re, json, os
P=str(P); W=f"{P}/wgs/19_hla_disease"; os.makedirs(W,exist_ok=True)
# No step of this repository runs T1K; without its dayu_genotype.tsv there are no HLA types to check.
if not os.path.exists(f"{P}/wgs/06_pgx/t1k/dayu_genotype.tsv"):
    json.dump({"genotype":{},"checks":[]},open(f"{W}/hla_disease.json","w"),ensure_ascii=False)
    print("WARNING: 06_pgx/t1k/dayu_genotype.tsv not found (no step runs T1K) -- empty hla_disease.json written", file=sys.stderr)
    sys.exit(0)
t=pd.read_csv(f"{P}/wgs/06_pgx/t1k/dayu_genotype.tsv",sep="\t",header=None,
              names=["gene","n","a1","ab1","q1","a2","ab2","q2","other"])
def two(a):  # keep two fields: 04:05
    m=re.match(r"(?:HLA-)?[A-Z0-9]+\*(\d+:\d+)",str(a)); return m.group(1) if m else None
G={}
for r in t.itertuples():
    if r.n<1: continue
    g=r.gene.replace("HLA-","")
    G[g]=[x for x in (two(r.a1),two(r.a2) if r.n>1 else None) if x]
print("typed:",{k:v for k,v in G.items() if k in ("A","B","C","DRB1","DQA1","DQB1","DPB1","DRB3","DRB4","DRB5")})

# A rule that queries a gene this run did not type must come out unknown, not negative:
# G.get(gene,[]) silently turned "no call" into "allele absent".
_UNTYPED=set()
def has(gene,allele):
    if gene not in G: _UNTYPED.add(gene); return False
    return allele in G[gene]
def hasgrp(gene,pref):
    if gene not in G: _UNTYPED.add(gene); return False
    return any(a.startswith(pref) for a in G[gene])

PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel" / "hla_disease.tsv"
CHECKS = []
for r in pd.read_csv(PANEL, sep="	").fillna("").itertuples():
    _UNTYPED.clear()
    pos = bool(eval(r.rule, {"has": has, "hasgrp": hasgrp}))
    CHECKS.append((r.category, r.condition_zh, r.condition_en, r.requirement_zh, r.requirement_en,
                   pos, ";".join(sorted(_UNTYPED)), r.note_zh, r.note_en))
EN={'乳糜泻 Coeliac disease': ('Coeliac disease', 'DQ2.5 (DQA1*05:01+DQB1*02:01) or DQ8 (DQA1*03+DQB1*03:02)', 'over 99% of patients carry DQ2 or DQ8; co-occurrence of the component alleles does not prove haplotype phase (no family or phase data here)'), '1 型糖尿病 Type 1 diabetes': ('Type 1 diabetes', 'DRB1*04:05-DQB1*04:01 (main East Asian risk haplotype)', 'reported odds ratios of 3 to 6 in East Asians; a risk shift, not a diagnosis'), '1 型糖尿病（另一风险单倍型）': ('Type 1 diabetes (other haplotype)', 'DRB1*09:01-DQB1*03:03', 'the second commonest East Asian risk haplotype'), '类风湿关节炎 Rheumatoid arthritis': ('Rheumatoid arthritis', 'shared-epitope alleles DRB1*04:05 / *04:01 / *01:01 / *10:01', 'DRB1*04:05 is the main shared-epitope allele in East Asians'), '强直性脊柱炎 Ankylosing spondylitis': ('Ankylosing spondylitis', 'B*27', 'about 90% of patients carry B*27; its absence makes the diagnosis much less likely'), '白塞病 Behçet disease': ('Behçet disease', 'B*51', 'the main risk allele in East Asian and Mediterranean populations'), '发作性睡病 1 型 Narcolepsy': ('Narcolepsy type 1', 'DQB1*06:02', 'carried by nearly 100% of patients; absence makes type 1 narcolepsy much less likely, a statistical association, not a diagnostic exclusion'), '多发性硬化 Multiple sclerosis': ('Multiple sclerosis', 'DRB1*15:01', 'the main European risk allele, uncommon in East Asia'), 'Graves 病': ('Graves disease', 'B*46:01 / DPB1*05:01 (East Asian)', 'DPB1*05:01 associates with Graves disease in East Asians, small effect'), 'Vogt–小柳–原田病 VKH': ('Vogt-Koyanagi-Harada disease', 'DRB1*04:05', 'enriched in East Asian VKH patients with a high odds ratio, but the disease itself is rare'), 'IgA 肾病': ('IgA nephropathy', 'DQB1*06:02 / DRB1*04:05 and others, evidence inconsistent', 'reports do not agree; exploratory only'), '氨苯砜超敏综合征 Dapsone': ('Dapsone hypersensitivity syndrome', 'B*13:01', 'odds ratio around 20 in Chinese cohorts; screen before prescribing'), '别嘌醇 SJS/TEN Allopurinol': ('Allopurinol SJS/TEN', 'B*58:01', 'strong association in East Asians; a negative does not mean zero risk'), '卡马西平 SJS/TEN Carbamazepine': ('Carbamazepine SJS/TEN', 'B*15:02', "strong in Han Chinese; both the FDA and China's regulator advise screening"), '卡马西平药疹 Carbamazepine DRESS': ('Carbamazepine DRESS', 'A*31:01', 'linked to delayed hypersensitivity across several populations'), '阿巴卡韦超敏 Abacavir': ('Abacavir hypersensitivity', 'B*57:01', 'a very strong association; screening is mandatory before use'), '奈韦拉平 Nevirapine': ('Nevirapine hypersensitivity', 'B*35:05 / DRB1*01:01', ''), '氟氯西林肝损伤 Flucloxacillin': ('Flucloxacillin liver injury', 'B*57:01', ''), '拉帕替尼肝损伤 Lapatinib': ('Lapatinib liver injury', 'DQA1*02:01 / DRB1*07:01', ''), '甲巯咪唑粒缺 Methimazole': ('Methimazole agranulocytosis', 'B*38:02 / DRB1*08:03', 'reported in East Asian cohorts'), 'HIV 病程控制 HIV control': ('HIV disease control', 'B*57:01 / B*27:05 (protective)', 'carriers progress more slowly'), '乙肝疫苗无应答 HBV vaccine': ('Hepatitis B vaccine non-response', 'DRB1*07:01 / DQB1*02:01', 'a moderate association'), '鼻咽癌 Nasopharyngeal carcinoma': ('Nasopharyngeal carcinoma', 'A*02:07 (East Asian risk) / A*11:01 (protective)', 'A*02:07 associates with raised risk in southern Han studies, small effect')}
rows=[]
for cat,cond,cond_en,req,req_en,pos,untyp,note,note_en in CHECKS:
    e=EN.get(cond,(cond_en,req_en,note_en))
    carried = "unknown" if untyp else ("yes" if pos else "no")
    # when a rule could not be fully evaluated, say which gene was missing instead of
    # letting the note's population statement read as a clinical exclusion of the sample.
    nt = note if carried != "unknown" else (f"{note}（{untyp} 未分型，无法判定）" if note else f"{untyp} 未分型，无法判定")
    nt_en = e[2] if carried != "unknown" else (f"{e[2]} ({untyp} not typed; cannot evaluate)" if e[2] else f"{untyp} not typed; cannot evaluate")
    rows.append(dict(category=cat,condition=cond,condition_en=e[0],requirement=req,requirement_en=e[1],
                     carried=carried,untyped=untyp,note=nt,note_en=nt_en))
df=pd.DataFrame(rows); df.to_csv(f"{W}/hla_disease.tsv",sep="\t",index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",52)
for cat,g in df.groupby("category",sort=False):
    print(f"\n=== {cat}")
    print(g[["condition","requirement","carried"]].to_string(index=False))
print("\n阳性条目：")
print(df[df.carried=="yes"][["condition","requirement","note"]].to_string(index=False))
json.dump({"genotype":G,"checks":rows},open(f"{W}/hla_disease.json","w"),ensure_ascii=False)
