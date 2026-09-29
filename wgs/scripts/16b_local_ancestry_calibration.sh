#!/bin/bash
# Calibration: score a held-out sample of the configured source panels with the same panels (AN3).
#
# Population list, holdout count, chromosomes and seed all come from the configuration. The holdout set
# is verified by *name* against the subset VCF: the previous version left `cols` empty after disabling the
# expression that computed it, so the awk subset carried a header and no samples at all, and nothing
# checked the result (a non-empty file was treated as success).
set -euo pipefail
source "$(dirname "$0")/env.sh"
# Scripts live in the repository; $WGS is the work tree. Resolve this file's own directory before cd.
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
W=$WGS/12_localanc; mkdir -p "$W" "$WGS/logs"; cd "$W"
R=$REF_DIR
JAVA=${JAVA:-java}
MEM_MB=$(( ${MEM_GB%.*} * 1000 ))

[ -s ref.panel ] || { echo "ref.panel is missing; run 16_local_ancestry.sh first" >&2; exit 1; }
[ -s "$R/all_phase3.psam" ] || { echo "no reference panel at $R/all_phase3.psam" >&2; exit 1; }

# 复审 AN0/AN5/H6 生命周期：开跑即作废上一轮的证据——旧的 ok 记录与旧 calib.* 是上一次尝试的
# 产物；本次若失败，不能留下它们冒充"本次校准成功"（17b 只认 state=ok 且样本相符的记录）。
# 校准是增强不是前置：本步失败不写 disabled 记录进 12_localanc，原始 LA（17/17b）照常交付，
# 失败证据在 run_info（failed）与 17b 的 calibration_state 里。
rm -f manifest.16b-la-calibration.json calib.*

# --- holdout selection: fixed seed, without replacement, min(config n, floor(group/5)) per population
python3 "$S/ancestry_data.py" --pick-holdout "$R/all_phase3.psam" \
  --pops "$CALIB_POPS" --n "$CALIB_N" --seed "$CALIB_SEED" --out holdout.ids || {
    echo "calibration holdout selection failed: state != ok (see above)" >&2; exit 2; }
[ -s holdout.ids ] || { echo "holdout selection wrote an empty list" >&2; exit 2; }
grep -vwFf holdout.ids ref.panel > ref.calib.panel
echo "holdout: $(wc -l < holdout.ids) individuals; reference panel rows now $(wc -l < ref.calib.panel)"
cut -f2 ref.calib.panel | sort | uniq -c

for c in ${CALIB_CHROMS//,/ }; do
  rm -f "ho.$c.vcf.gz" "ho.$c.vcf.gz.tbi"
  # Sample-subset mode first: it is the supported path when the header can be parsed. These 1000G VCFs
  # have minimal headers, so when bcftools refuses we cut the sample columns with awk -- computing the
  # column list for real, and failing loudly if it comes out empty.
  if bcftools view -S holdout.ids -Oz -o "ho.$c.vcf.gz" "$R/imp/ALL.chr$c.vcf.gz" 2>/dev/null; then
    tabix -f "ho.$c.vcf.gz"
    echo "chr$c: subset by bcftools -S"
  else
    cols=$(zcat "$R/imp/ALL.chr$c.vcf.gz" | grep -m1 "^#CHROM" | tr '\t' '\n' \
             | grep -nxFf holdout.ids | cut -d: -f1 | paste -sd, || true)
    [ -n "$cols" ] || { echo "chr$c: none of the holdout samples appear in the VCF header" >&2; exit 3; }
    zcat "$R/imp/ALL.chr$c.vcf.gz" \
      | awk -v c="$cols" 'BEGIN{n=split(c,C,","); for(i=1;i<=9;i++) k[i]=1; for(i=1;i<=n;i++) k[C[i]]=1} { s=$1; for(i=2;i<=NF;i++) if(k[i]) s=s "\t" $i; print s }' \
      | bgzip -@ 4 > "ho.$c.vcf.gz"
    tabix -f "ho.$c.vcf.gz"
    echo "chr$c: subset by awk column cut (bcftools refused the header)"
  fi
  # The subset must carry exactly the holdout set -- checked by name, not by size.
  python3 "$S/ancestry_data.py" --check-samples "ho.$c.vcf.gz" --expect holdout.ids \
    || { echo "chr$c: the subset does not carry exactly the holdout samples" >&2; exit 3; }

  "$JAVA" "-Xmx${MEM_GB%.*}g" -jar "$FLARE" ref="$R/imp/ALL.chr$c.vcf.gz" ref-panel=ref.calib.panel \
    gt="ho.$c.vcf.gz" map="$R/imp/maps/plink.chr$c.GRCh37.map" out="calib.$c" gen=100 \
    nthreads="$THREADS" > "$WGS/logs/flare_calib.$c.log" 2>&1
  rm -f "ho.$c.vcf.gz" "ho.$c.vcf.gz.tbi"
  echo "chr$c calibration done"
done
# 成功凭证（复审 AN0/AN5/H6）：17b 凭这份 state=ok 且样本相符的记录读 calib.*——只看文件存在
# 会把上一轮失败前的旧校准当成本次结果。写不进记录就视为失败，不留半成功状态。
python3 "$S/ancestry_data.py" --step-ok 16b-la-calibration \
  --out "$W/manifest.16b-la-calibration.json" --sample "$SAMPLE" \
  || { echo "calibration ran but its success record could not be written" >&2; exit 7; }
echo CALIB_DONE
