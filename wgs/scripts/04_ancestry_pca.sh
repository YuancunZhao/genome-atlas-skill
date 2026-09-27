#!/bin/bash
# 1000G PCA on LD-pruned common autosomal SNVs; project the target (complete genotype set) with --score.
set -euo pipefail
source "$(dirname "$0")/env.sh"
W=$WGS/04_ancestry; mkdir -p $W; cd $W
R=$KG_PFILE
# Candidate sites come from the reference panel alone: autosomal biallelic SNVs with global MAF>=0.05.
# The target used to select the training sites (--extract target.present.ids), which made the reference
# space depend on the sample and dropped panel sites the sample could not be called at; the target's own
# gaps are now reported as coverage instead. MEM is given in GB in the config.
MEM_MB=$(( ${MEM_GB%.*} * 1000 ))
$PLINK2 --pfile $R --autosome --snps-only just-acgt --max-alleles 2 --maf 0.05 \
  --rm-dup exclude-all --make-pgen --out kg.common --threads $THREADS --memory $MEM_MB >/dev/null
# drop strand-ambiguous not needed (same reference build, same alleles); thin then LD-prune
$PLINK2 --pfile kg.common --bp-space 2000 --indep-pairwise 200 50 0.2 --out prune --threads $THREADS --memory $MEM_MB >/dev/null
$PLINK2 --pfile kg.common --extract prune.prune.in --freq --pca 10 allele-wts --out kg.pca --threads $THREADS --memory $MEM_MB >/dev/null
for s in kg:kg.common target:$WGS/02_complete/$SAMPLE.1kg; do n=${s%%:*}; f=${s#*:}
  $PLINK2 --pfile $f --extract prune.prune.in --read-freq kg.pca.afreq \
     --score kg.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
     --score-col-nums 6-15 --out $n.proj --threads $THREADS --memory $MEM_MB >/dev/null
done
echo "pruned SNPs: $(wc -l < prune.prune.in)"; head -2 target.proj.sscore
# Regional PCA. The super-population comes from the configuration (SUPERPOP); the products are named
# after the scope, not "eas", so a non-East-Asian panel is not described by the wrong file name.
SCOPE=regional
{ printf '#FID\tIID\n'; awk -v sp="$SUPERPOP" 'NR>1 && $4==sp{print $1"\t"$2}' $R.psam; } > $SCOPE.ids
$PLINK2 --pfile kg.common --keep $SCOPE.ids --maf 0.05 --bp-space 2000 --indep-pairwise 200 50 0.2 --out prune.$SCOPE --threads $THREADS --memory $MEM_MB >/dev/null
$PLINK2 --pfile kg.common --keep $SCOPE.ids --extract prune.$SCOPE.prune.in --freq --pca 10 allele-wts --out $SCOPE.pca --threads $THREADS --memory $MEM_MB >/dev/null
$PLINK2 --pfile kg.common --keep $SCOPE.ids --extract prune.$SCOPE.prune.in --read-freq $SCOPE.pca.afreq --score $SCOPE.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out $SCOPE.proj --threads $THREADS --memory $MEM_MB >/dev/null
$PLINK2 --pfile $WGS/02_complete/$SAMPLE.1kg --extract prune.$SCOPE.prune.in --read-freq $SCOPE.pca.afreq --score $SCOPE.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out $SCOPE.target.proj --threads $THREADS --memory $MEM_MB >/dev/null
echo "$SUPERPOP pruned SNPs: $(wc -l < prune.$SCOPE.prune.in)"
echo PCA_DONE
