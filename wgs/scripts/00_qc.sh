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
mosdepth --no-per-base -t "$THREADS" -f "$REF" -Q 20 --by 1000 --quantize 0:1:4:8:20:60: "$WGS/01_qc/depth" "$SAMPLE.cram"
tabix -f -p bed "$WGS/01_qc/depth.regions.bed.gz" || true
echo "QC_DONE -> $WGS/01_qc/depth.quantized.bed.gz"
