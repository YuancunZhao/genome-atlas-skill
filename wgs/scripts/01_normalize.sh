#!/bin/bash
# Normalize the delivered VCF: split multiallelics, left-align, keep PASS, drop decoy contigs.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd $WGS/00_input
# 1. callable BED: depth>=8 & MQ>=20 (from mosdepth quantized), main contigs only
zcat $WGS/01_qc/depth.quantized.bed.gz | awk -F'\t' '$1!~/^GL/ && ($4=="8:20"||$4=="20:60"||$4=="60:inf")' \
  | bedtools merge -i - > callable.bed
awk '{s+=$3-$2} END{printf "callable bp: %d\n", s}' callable.bed
# 2. normalized full VCF (all filters kept, FILTER column retained) on main contigs
MAIN=$(seq 1 22 | tr '\n' ',')X,Y,MT
# delivered VCFs frequently carry a "chr" prefix (and chrM); -r 1,2,...,X,Y,MT then matches
# nothing and the view below writes an empty VCF without any error. Rename the contigs first.
IN="$VCF_RAW"
bcftools view -h "$IN" > header.txt
if awk -F'[<>,]' '/^##contig=<ID=chr/{f=1} END{exit !f}' header.txt; then
  echo "chr-prefixed contigs: renaming to $MAIN"
  awk -F'[<>,]' '/^##contig=<ID=chr/{n=substr($3,4); if(n=="M") n="MT"; print $3"\t"n}' header.txt > chr_rename.txt
  bcftools annotate --rename-chrs chr_rename.txt -Oz -o input.renamed.vcf.gz "$IN"
  tabix -f -p vcf input.renamed.vcf.gz
  IN=input.renamed.vcf.gz
fi
bcftools view -r $MAIN "$IN" -Ou \
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
