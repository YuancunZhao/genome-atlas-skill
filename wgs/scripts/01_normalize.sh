#!/bin/bash
# Normalize the delivered VCF: split multiallelics, left-align, keep PASS, drop decoy contigs.
set -euo pipefail
source "$(dirname "$0")/env.sh"
mkdir -p "$WGS/00_input"
cd $WGS/00_input
# Main contigs are derived from the VCF header (1-22/X/Y/MT with or without a chr prefix), not
# hardcoded: a bare list breaks on chr-prefixed inputs, and scaffolds/decoys must fall through.
# Region names keep the file's own spelling for bcftools -r; matching elsewhere canonicalizes.
MAIN_LIST=$WGS/00_input/main_contigs.txt
bcftools view -h "$VCF_RAW" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' \
  | awk '{bare=$1; sub(/^chr/,"",bare); u=toupper(bare); if (u=="M") u="MT";
          if ((bare ~ /^[0-9]+$/ && bare+0>=1 && bare+0<=22) || u=="X" || u=="Y" || u=="MT") print $1}' \
  | sort -u > "$MAIN_LIST"
[ -s "$MAIN_LIST" ] || { echo "no main contigs (1-22/X/Y/MT, chr-prefixed or bare) in the header of $VCF_RAW" >&2; exit 1; }
MAIN=$(paste -sd, "$MAIN_LIST")
# 1. callable BED: depth>=8 & MQ>=20 (from mosdepth quantized), main contigs only. Positive match
# on the canonical bare name (chrM==MT), so the old negative /^GL/ filter that silently kept
# chrUn_*/random scaffolds on chr-prefixed inputs is gone, and bed/VCF naming may differ.
zcat $WGS/01_qc/depth.quantized.bed.gz \
  | awk -F'\t' 'NR==FNR{b=$1; sub(/^chr/,"",b); u=toupper(b); if (u=="M") u="MT"; m[u]=1; next}
              {b=$1; sub(/^chr/,"",b); u=toupper(b); if (u=="M") u="MT";
               if (m[u] && ($4=="8:20"||$4=="20:60"||$4=="60:inf")) print}' "$MAIN_LIST" - \
  | bedtools merge -i - > callable.bed
awk '{s+=$3-$2} END{printf "callable bp: %d\n", s}' callable.bed
[ -s callable.bed ] || { echo "callable.bed is empty: no depth/MQ-passing main-contig regions -- check the mosdepth input and its contig naming" >&2; exit 1; }
# 2. normalized full VCF (all filters kept, FILTER column retained) on main contigs
bcftools view -r $MAIN "$VCF_RAW" -Ou \
 | bcftools norm -f $REF -m -both -c w --threads $THREADS -Oz -o target.norm.vcf.gz 2> norm.log
tabix -f -p vcf target.norm.vcf.gz
# 3. PASS-only
bcftools view -f PASS --threads $THREADS -Oz -o target.pass.vcf.gz target.norm.vcf.gz
tabix -f -p vcf target.pass.vcf.gz
# 4. stats
bcftools stats -F $REF -s - target.norm.vcf.gz > $WGS/01_qc/stats.norm.txt
bcftools stats -F $REF -s - target.pass.vcf.gz > $WGS/01_qc/stats.pass.txt
tail -3 norm.log
echo NORMALIZE_DONE
