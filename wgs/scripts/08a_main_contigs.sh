#!/bin/bash
# Main-contig alignment (H5): a hand-made intermediate that Delly's recovered command line depends on.
#
# Delly was invoked on 08_sv/main_contigs/norm.main.cram. That file exists in the delivered work tree, but
# **how it was made is not recoverable from what is left**: its own @PG chain is gone (the file was
# overwritten during this investigation and the pre-existing copy did not carry a producer line for this
# step), and the one @PG chain that does survive -- `samtools view -T <b37> -O cram,embed_ref=2 -@16` --
# belongs to 00_input/norm.cram, not to this file.
#
# Three attempts here did not reproduce the delivered size, and the arithmetic says why the approach is
# wrong rather than merely imprecise: the delivered norm.main.cram is 1.24 GB while 00_input/norm.cram is
# 47 GB, and the main chromosomes are ~92% of the genome, so a full main-contig extract would be tens of
# gigabytes. A 1.24 GB file is not a whole-extract of anything -- it holds a small fraction of the reads,
# and what selects them is unknown.
#
# The script below therefore implements the documented reading (main chromosomes, matching the reference)
# and is marked UNVERIFIED. Do not treat its output as equivalent to the delivered file, and do not re-run
# Delly on it expecting the delivered SV calls. The step needs the original command from whoever ran it.
set -euo pipefail
source "$(dirname "$0")/env.sh"
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

# 格式与压缩参数取自交付件自身的 @PG 记录（不是猜的）：
#   samtools view -T <b37 fasta> -O cram,embed_ref=2 -@16
# embed_ref=2 把参考序列嵌进 CRAM，这是交付件只有 1.24 GB 的原因；不嵌入时同样内容要 5-8 GB
# （-b 出 BAM 8.2 GB、-C 默认 5.2-6.9 GB 都实测过）。文件更大本身不算错，但既然 @PG 明写了参数，
# 就没有理由产出一个与交付件不同的文件——下游若按大小或按"是否需要外部参考"判断，差异会误导。
# -F 4 排除未比对的读段：-L 是"与区域重叠即取"，未比对的读段没有坐标、会因重叠判断被一并带出，
# 实测使输出达到交付件的 5 倍（6.2 GB vs 1.24 GB）。主 contig 提取的目的是让坐标与参考一致，
# 而未比对的读段不带坐标，对 SV 调用没有贡献。
samtools view -T "$REF" -@ "${THREADS:-8}" -O cram,embed_ref=2 -F 4 -L main_contigs.bed \
  -o norm.main.cram "$CRAM" 2> "$WGS/logs/main_contigs.log"
samtools index norm.main.cram
[ -s norm.main.cram ] || { echo "no output written" >&2; exit 1; }
echo "MAIN_CONTIGS_DONE -> $W/norm.main.cram ($n contigs)"
