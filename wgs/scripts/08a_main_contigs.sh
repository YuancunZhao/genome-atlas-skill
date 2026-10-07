#!/bin/bash
# Main-contig alignment feeding Delly (08d): samtools view over the main chromosomes of $CRAM, as an
# embedded-reference CRAM. Parameters come from the delivered pipeline's own records --
# `samtools view -T <b37> -O cram,embed_ref=2 -@16` is the surviving @PG chain, and -F 4 / -L are the
# documented reading (unmapped reads carry no coordinates, so they cannot contribute to SV calling;
# keeping them quintuples the file for nothing).
#
# Verified against the delivery (2026-09-28): this extract + 08d's `delly sr -g b37 -h 16`
# reproduced the delivered target.sv.bcf exactly -- 42,909/42,909 records identical on
# (CHROM,POS,SVTYPE,END), md5 differing only in the BCF ##fileDate line. Evidence:
# work/wgs/08_sv/delly_retest_20260928/RETEST_NOTES.md. The retest extract (46.8 GB) differs in
# representation from the 1.24 GB file the delivery left behind, so sizes are not expected to match;
# what was verified is the SV-calling content.
#
# The delivered tree's own norm.main.cram (no provenance marker, 1.24 GB) is the last surviving
# artifact of the original run, so this script never silently overwrites it: rm it first or set
# REBUILD_MAIN_CRAM=1 to replace it deliberately.
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
W=$WGS/08_sv/main_contigs; mkdir -p "$W" "$WGS/logs"; cd "$W"

[ -f "$CRAM" ] || { echo "input alignment not found: $CRAM" >&2; exit 1; }

# 主染色体列表来自参考索引本身，而不是写死 "1..22,X,Y,MT"：换参考时列表要跟着变。
FAI=${REF}.fai          # env.sh 导出的 REF 就是 FASTA 路径（不是目录）
[ -f "$FAI" ] || { echo "reference index not found: $FAI" >&2; exit 1; }
awk 'BEGIN{OFS="\t"} $1 ~ /^([0-9]+|X|Y|MT|M)$/ {print $1, 0, $2}' "$FAI" > main_contigs.bed
n=$(wc -l < main_contigs.bed)
[ "$n" -ge 24 ] || { echo "expected at least 24 main contigs, got $n" >&2; exit 1; }

# 覆盖保护（见头注释）：不带 provenance 标记的现存 norm.main.cram 视为交付原件，不默默替换。
if [ -s norm.main.cram ] && [ ! -f norm.main.cram.provenance.json ]; then
  echo "refusing to overwrite the existing norm.main.cram: it carries no provenance marker, so it is" >&2
  echo "most likely the delivered original from the run this pipeline reconstructs." >&2
  echo "To replace it deliberately: rm norm.main.cram && re-run, or set REBUILD_MAIN_CRAM=1." >&2
  [ "${REBUILD_MAIN_CRAM:-0}" = "1" ] || exit 0
  echo "REBUILD_MAIN_CRAM=1: overwriting anyway" >&2
fi

samtools view -T "$REF" -@ "${THREADS:-8}" -O cram,embed_ref=2 -F 4 -L main_contigs.bed \
  -o norm.main.cram "$CRAM" 2> "$WGS/logs/main_contigs.log"

# 产出即记录本次运行做了什么（生产者/命令/输入/时间）。“是否复现某次历史交付”是复测证据
# （头注释），不是每次运行都能重新主张的运行时断言——新样本的 provenance 只如实描述自身。
cat > norm.main.cram.provenance.json <<JSON
{
  "producer": "08a_main_contigs.sh",
  "command": "samtools view -T $REF -@ ${THREADS:-8} -O cram,embed_ref=2 -F 4 -L main_contigs.bed -o norm.main.cram $CRAM",
  "input_cram": "$CRAM",
  "bed": "main_contigs.bed",
  "contigs": $n,
  "date": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON
samtools index norm.main.cram
[ -s norm.main.cram ] || { echo "no output written" >&2; exit 1; }
echo "MAIN_CONTIGS_DONE -> $W/norm.main.cram ($n contigs)"
