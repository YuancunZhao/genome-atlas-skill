#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/env.sh"
mkdir -p "$WGS/00_input" "$WGS/01_qc" "$WGS/logs"; cd "$WGS/00_input"
IN="${CRAM}"; EXT="${IN##*.}"
[ -e "$SAMPLE.cram" ] || ln -s "$IN" "$SAMPLE.cram"
case "$EXT" in
  cram) [ -e "$SAMPLE.cram.crai" ] || ln -sf "$IN.crai" "$SAMPLE.cram.crai" 2>/dev/null || true ;;
  bam)  [ -e "$SAMPLE.cram.bai"  ] || ln -sf "${IN%.bam}.bam.bai" "$SAMPLE.cram.bai" 2>/dev/null \
          || ln -sf "$IN.bai" "$SAMPLE.cram.bai" 2>/dev/null || true ;;
esac
if [ -z "${VENDOR_VCF:-}" ] && [ ! -f "$SAMPLE.call.vcf.gz" ]; then
  echo "no vendor VCF; calling from reads (slow)…"
  bcftools mpileup -f "$REF" -a FORMAT/AD,FORMAT/DP -Ou "$SAMPLE.cram" | bcftools call -m -v -Oz -o "$SAMPLE.call.vcf.gz"
  tabix -f "$SAMPLE.call.vcf.gz"
fi
# 参数语义（`mosdepth --help` 原文）：`-Q --mapq <mapq>  mapping quality threshold. reads with a quality
# less than this value are ignored`。所以这里接的是**比对质量** MIN_MQ，不是碱基质量——前一版接了
# MIN_BQ：接错的后果是"改 callable_min_mapq 却不影响深度"，而深度喂给可调用区间与 CNV 归一化，
# 偏差会一路传下去。mosdepth 没有 base-quality 阈值选项；MIN_BQ 用于 bcftools 的 -Q（02_recall_x_mt）。
mosdepth --no-per-base -t "$THREADS" -f "$REF" -Q "${MIN_MQ:?MIN_MQ not exported by env.sh}" --by 1000 --quantize 0:1:4:8:20:60: "$WGS/01_qc/depth" "$SAMPLE.cram"
tabix -f -p bed "$WGS/01_qc/depth.regions.bed.gz" || true
echo "QC_DONE -> $WGS/01_qc/depth.quantized.bed.gz"
