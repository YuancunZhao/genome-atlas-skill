#!/bin/bash
# HLA and KIR typing with T1K (H5): this product had no producer in the repository.
#
# The whole chain is recovered from the tool's own logs, which record every sub-command it ran
# (06_pgx/t1k/t1k_hla.log and t1k_kir.log, run-t1k v1.0.10-r263):
#
#   HLA  fastq-extractor -t 8 -f hlaidx/hlaidx_dna_seq.fa -o dayu_candidate -s 0.97 -1 R1 -2 R2
#        genotyper     -s 0.97 -o dayu -t 8 -f hlaidx/hlaidx_dna_seq.fa -1 dayu_candidate_1.fq -2 ...
#        analyzer      -s 0.97 -o dayu -t 8 -f hlaidx/... -a dayu_allele.tsv -1 dayu_aligned_1.fa -2 ...
#   KIR  the same three with -s 0.9 --relaxIntronAlign and kiridx/kiridx_dna_seq.fa, output prefix kir
#
# Note the asymmetry, which is copied rather than smoothed over: HLA uses -s 0.97, KIR uses 0.9 with
# --relaxIntronAlign. Those are the delivered parameters; "tidying" them into one value would change the
# calls. The sample prefix in the delivered outputs is dayu.
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
  samtools fastq -@ "$TH" -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null hla_kir.bam
fi

run_one(){
  local tag=$1 idx=$2 score=$3 extra=$4
  "$T1K/fastq-extractor" -t "$TH" -f "$T1K/$idx" -o "${tag}_candidate" $extra -1 "$R1" -2 "$R2"
  "$T1K/genotyper" -s "$score" $extra -o "$tag" -t "$TH" -f "$T1K/$idx" \
      -1 "${tag}_candidate_1.fq" -2 "${tag}_candidate_2.fq"
  "$T1K/analyzer" -s "$score" $extra -o "$tag" -t "$TH" -f "$T1K/$idx" \
      -a "${tag}_allele.tsv" -1 "${tag}_aligned_1.fa" -2 "${tag}_aligned_2.fa"
  echo "$tag: $(wc -l < "${tag}_genotype.tsv" 2>/dev/null || echo 0) genotype rows"
}

run_one "$PREFIX" hlaidx/hlaidx_dna_seq.fa 0.97 "" > "$WGS/logs/t1k_hla.log" 2>&1 || { cat "$WGS/logs/t1k_hla.log" >&2; exit 1; }
echo "T1K_HLA_DONE"
run_one kir       kiridx/kiridx_dna_seq.fa 0.9 "--relaxIntronAlign" > "$WGS/logs/t1k_kir.log" 2>&1 || { cat "$WGS/logs/t1k_kir.log" >&2; exit 1; }
echo "T1K_KIR_DONE"

for f in "${PREFIX}_genotype.tsv" kir_genotype.tsv; do
  [ -s "$f" ] || { echo "T1K produced no $f" >&2; exit 1; }
done
echo "T1K_DONE -> $W/{${PREFIX},kir}_genotype.tsv"
