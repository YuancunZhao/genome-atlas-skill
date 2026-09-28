#!/bin/bash
# Re-call chrX with correct ploidy (PAR diploid, non-PAR haploid for a male) and pile up MT for heteroplasmy.
set -euo pipefail
source "$(dirname "$0")/env.sh"
mkdir -p $WGS/03_haplo $WGS/00_input
cd $WGS/00_input
# X 的倍性按**证据**决定，而不是按配置声明（H4）：配置说的是样本被声明成什么，它说不出这份 CRAM
# 里到底有没有 Y、X 是否真的测到了两条。判定依据由 _sex_evidence.py 从比对索引算出并留下记录；
# 证据不可用时回退到声明值，但明确告警——回退可以选择，沉默不可以。
# 复用的前提不只是"文件在"：它必须是**本样本**的证据。上一个样本留下的旧 sex_evidence.json
# 会把 chrX 倍性整整定错一轮，所以先校验 sample_id，不匹配就重算；SM 块里再校验一次，重算
# 失败也读不到别人的证据，而是按声明回退并留话。
_ev="$WGS/sex_evidence.json"
if [ -f "$_ev" ]; then
  _ev_sid=$(python3 -c "import json;print(json.load(open('$_ev')).get('sample_id',''))" 2>/dev/null || true)
  if [ "$_ev_sid" != "$SAMPLE" ]; then
    echo "note: $_ev does not carry this sample's id (found '${_ev_sid:-unreadable}', want '$SAMPLE'); recomputing" >&2
    python3 "$S/_sex_evidence.py" >/dev/null 2>&1 || true
  fi
else
  python3 "$S/_sex_evidence.py" >/dev/null 2>&1 || true
fi
SM=$(python3 - "$_ev" "${SEX:-male}" "$SAMPLE" <<'PY'
import json, sys
fallback = sys.argv[2][:1].upper() or "M"
try:
    d = json.load(open(sys.argv[1]))
except Exception:
    print(fallback, end=""); sys.exit(0)
if d.get("sample_id") not in (None, sys.argv[3]):
    print(fallback, end="")
    print(f"sex_evidence.json carries sample_id {d.get('sample_id')!r}, not this sample; "
          f"chrX ploidy taken from the declared sex ({sys.argv[2]})", file=sys.stderr)
    sys.exit(0)
if d.get("x_conflict"):
    print(fallback, end="")
    print(f"X and Y evidence conflict ({d.get('why')}); chrX ploidy taken from the declared sex "
          f"({sys.argv[2]})", file=sys.stderr)
    sys.exit(0)
inf = d.get("inferred")
print({"has_y": "M", "no_y": "F"}.get(inf, fallback), end="")
if inf not in ("has_y", "no_y"):
    print(f"chrX ploidy taken from the declared sex ({sys.argv[2]}); read-level evidence was {inf or d.get('state')}",
          file=sys.stderr)
PY
)
echo -e "$SAMPLE\t$SM" > sample_sex.txt
echo "chrX ploidy: $SM (evidence: $(python3 -c "import json;d=json.load(open('$_ev'));print(d.get('inferred') or d.get('state'))" 2>/dev/null || echo none))"
# X: parallel over 4 chunks
for r in X:1-40000000 X:40000001-80000000 X:80000001-120000000 X:120000001-155270560; do
  ( bcftools mpileup -f $REF -r $r -q "$MIN_MQ" -Q "$MIN_BQ" -d "$MPILEUP_MAX_DP" -a FORMAT/AD,FORMAT/DP -Ou "$CRAM" 2>/dev/null \
    | bcftools call -m -v --ploidy GRCh37 -S sample_sex.txt -Oz -o X.$(echo $r|tr ':-' '__').vcf.gz 2>/dev/null ) &
done
wait
bcftools concat -a -Oz -o X.recall.raw.vcf.gz X.X_*.vcf.gz 2>/dev/null || { for f in X.X_*.vcf.gz; do tabix -f $f; done; bcftools concat -a -Oz -o X.recall.raw.vcf.gz X.X_*.vcf.gz; }
tabix -f X.recall.raw.vcf.gz
# QUAL 与深度门槛同样来自配置：DP>=8 就是 MIN_DP，写死会让「改配置生效」在 X 召回上落空。
bcftools norm -f $REF -m -both -Ou X.recall.raw.vcf.gz | bcftools filter -i "QUAL>=$MIN_VQ && FORMAT/DP>=$MIN_DP" -Oz -o X.recall.vcf.gz
tabix -f X.recall.vcf.gz
rm -f X.X_*.vcf.gz*
echo "PAR1 het/hom:"; bcftools query -r X:60001-2699520 -f '[%GT]\n' X.recall.vcf.gz | sort | uniq -c
echo "nonPAR GTs:"; bcftools query -r X:2699521-154931043 -f '[%GT]\n' X.recall.vcf.gz | sort | uniq -c
# MT pileup: all positions, allele depths
bcftools mpileup -f $REF -r MT -q "$MIN_MQ" -Q "$MIN_BQ" -d "$MT_MAX_DP" -a FORMAT/AD -Ou "$CRAM" 2>/dev/null \
  | bcftools query -f '%POS\t%REF\t%ALT\t[%AD]\n' > $WGS/03_haplo/mt_pileup_ad.tsv
wc -l $WGS/03_haplo/mt_pileup_ad.tsv
echo RECALL_DONE
