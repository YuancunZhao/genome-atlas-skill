#!/bin/bash
# Main-contig alignment (H5): a hand-made intermediate that Delly's recovered command line depends on.
#
# Delly was invoked on 08_sv/main_contigs/norm.main.cram, not on the pipeline's own CRAM. That file is the
# alignment restricted to the reference's main chromosomes, which is what the head of delly.log is about:
# dozens of "BAM file chromosome <decoy> is NOT present in your reference file" warnings, because the
# original alignment carries hg19 decoy/alt contigs that b37 does not. Feeding Delly the unrestricted CRAM
# would either drop those reads noisily or mis-handle them; restricting first makes the input match the
# reference exactly.
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

# -C（CRAM）而不是 -b（BAM）：文件名是 .cram，交付件也是 CRAM。用 -b 会写出一份 BAM 却叫 .cram，
# 体积差 6 倍以上（实测 8.16 GB vs 交付 1.24 GB），而且下游按 CRAM 读时会出错。
samtools view -@ "${THREADS:-8}" -C -L main_contigs.bed -o norm.main.cram "$CRAM" 2> "$WGS/logs/main_contigs.log"
samtools index norm.main.cram
[ -s norm.main.cram ] || { echo "no output written" >&2; exit 1; }
echo "MAIN_CONTIGS_DONE -> $W/norm.main.cram ($n contigs)"
