#!/bin/bash
# HLA and KIR typing with T1K (H5): this product had no producer in the repository.
#
# The whole chain is recovered from the tool's own logs, which record every sub-command it ran
# (06_pgx/t1k/t1k_hla.log and t1k_kir.log, run-t1k v1.0.10-r263). The per-stage arguments are NOT
# symmetric and are copied verbatim from those SYSTEM CALL lines:
#
#   HLA  fastq-extractor -t 8 -f hlaidx/hlaidx_dna_seq.fa -o dayu_candidate -s 0.97 -1 R1 -2 R2
#        genotyper        -s 0.97 -o dayu -t 8 -f hlaidx/... -1 dayu_candidate_1.fq -2 ...
#        analyzer         -s 0.97 -o dayu -t 8 -f hlaidx/... -a dayu_allele.tsv -1 dayu_aligned_1.fa ...
#   KIR  fastq-extractor  (no extra arguments) -t 8 -f kiridx/kiridx_dna_seq.fa -o kir_candidate -1 R1 -2 R2
#        genotyper        -s 0.9 --relaxIntronAlign -o kir -t 8 -f kiridx/...
#        analyzer         -s 0.9 --relaxIntronAlign -o kir -t 8 -f kiridx/...
#
# An earlier draft generalized "-s 0.9 --relaxIntronAlign" onto the KIR extractor and dropped the
# HLA extractor's -s 0.97: the HLA candidate set ballooned from the delivered 5.8 MB to 39.5 MB
# and the KIR extractor returned zero candidates, leaving the genotyper nothing to type. The
# delivered parameters are not "tidied" into symmetry: changing them changes the calls.
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
W=$WGS/06_pgx/t1k; mkdir -p "$W" "$WGS/logs"; cd "$W"
T1K=${T1K_DIR:-$TOOLS/T1K}
PREFIX=${T1K_PREFIX:-dayu}
TH=${THREADS:-8}

for x in fastq-extractor genotyper analyzer; do
  [ -x "$T1K/$x" ] || { echo "T1K component missing: $T1K/$x" >&2; exit 1; }
done

# HLA/KIR 读段：交付的 hla_kir_R{1,2}.fq 是 364 MB 一对，来自比对文件的区域提取。
# 区域来源未在产物中留下记录；这里取 b37 上 MHC 与 KIR 两个位点的规范坐标（索引 FASTA 本身是
# 等位基因命名，不含坐标）。坐标按交付读段抽样核验：交付 FASTQ 的读段全部落在 MHC(6:29.9-33.5M)
# 与 KIR(19:54.0-55.4M)，而早期草稿里的 19:59-64M 与 6:162-171M 零命中/纯垃圾跨度，已弃用。
R1=hla_kir_R1.fq; R2=hla_kir_R2.fq
if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
  [ -f "$CRAM" ] || { echo "input alignment not found: $CRAM" >&2; exit 1; }
  # 区域按参考升序给出，且输出必须过一道 sort：多个区域交叉传给 samtools view 时输出
  # 会出现"染色体块不连续"（实测 6:… 19:… 6:… 的顺序会让 6:162M 的读段排在 chr19 之后），
  # samtools index 随即失败、整条链中断。sort 保证坐标连续，索引稳定建立。
  samtools view -@ "$TH" -h "$CRAM" 6:29900000-33500000 19:54000000-55400000 \
    | samtools sort -@ "$TH" -o hla_kir.bam -
  samtools index hla_kir.bam
  # fastq 前必须 collate：samtools fastq 按相邻配对，坐标排序输入下不相邻的对会被当作单端
  # 丢弃（实测 1,943,434 读段只剩 7,243 对，交付流程同样先 collate——见 HANDOFF §5 的
  # "T1K 要求严格配对，曾用 collate→fastq"）。collate 后只剩真单端（配对另一半在区域外）。
  samtools collate -@ "$TH" -O hla_kir.bam \
    | samtools fastq -@ "$TH" -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null -
fi

run_one(){
  local tag=$1 idx=$2 geno_args=$3 extract_args=$4
  "$T1K/fastq-extractor" -t "$TH" -f "$T1K/$idx" -o "${tag}_candidate" $extract_args -1 "$R1" -2 "$R2"
  "$T1K/genotyper" $geno_args -o "$tag" -t "$TH" -f "$T1K/$idx" \
      -1 "${tag}_candidate_1.fq" -2 "${tag}_candidate_2.fq"
  "$T1K/analyzer" $geno_args -o "$tag" -t "$TH" -f "$T1K/$idx" \
      -a "${tag}_allele.tsv" -1 "${tag}_aligned_1.fa" -2 "${tag}_aligned_2.fa"
  echo "$tag: $(wc -l < "${tag}_genotype.tsv" 2>/dev/null || echo 0) genotype rows"
}

run_one "$PREFIX" hlaidx/hlaidx_dna_seq.fa "-s 0.97" "-s 0.97" > "$WGS/logs/t1k_hla.log" 2>&1 || { cat "$WGS/logs/t1k_hla.log" >&2; exit 1; }
echo "T1K_HLA_DONE"
run_one kir kiridx/kiridx_dna_seq.fa "-s 0.9 --relaxIntronAlign" "" > "$WGS/logs/t1k_kir.log" 2>&1 || { cat "$WGS/logs/t1k_kir.log" >&2; exit 1; }
echo "T1K_KIR_DONE"

for f in "${PREFIX}_genotype.tsv" kir_genotype.tsv; do
  [ -s "$f" ] || { echo "T1K produced no $f" >&2; exit 1; }
done
echo "T1K_DONE -> $W/{${PREFIX},kir}_genotype.tsv"
