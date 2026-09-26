#!/bin/bash
# Annotate PASS VCF with ClinVar (pos+ref+alt), Ensembl consequences (bcftools csq), 1000G EAS AF.
set -euo pipefail
source "$(dirname "$0")/env.sh"
A=$REF_DIR/annot; W=$WGS/05_clinvar; mkdir -p $W; cd $W
CV=$REF_DIR/clinvar_grch37.vcf.gz
[ -f $CV.tbi ] || tabix -p vcf $CV
# 1. EAS AF annotation file straight from the 1000G panel, whose INFO column already carries
# EAS_AF and AF. The two .afreq files this used to read are never produced by anything in the
# repository, and the process substitutions hid that: paste then wrote an empty eas_af.tsv.gz
# and every variant silently lost its frequencies.
awk -F'\t' 'BEGIN{OFS="\t"} $1!~/^#/ { e="."; a="."; n=split($8,kv,";"); for(i=1;i<=n;i++){ if(kv[i]~/^EAS_AF=/) e=substr(kv[i],8); else if(kv[i]~/^AF=/) a=substr(kv[i],4) } print $1,$2,$4,$5,e,a }' $REF_DIR/all_phase3.pvar | bgzip -@4 > eas_af.tsv.gz
tabix -f -s1 -b2 -e2 eas_af.tsv.gz
printf '##INFO=<ID=EAS_AF,Number=1,Type=Float,Description="1000G phase3 EAS alt allele frequency">\n##INFO=<ID=ALL_AF,Number=1,Type=Float,Description="1000G phase3 global alt allele frequency">\n' > eas_af.hdr
# 2. ClinVar + EAS AF + csq
bcftools annotate -a $CV -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/GENEINFO,INFO/CLNSIGCONF,INFO/CLNVC,INFO/ALLELEID --threads 8 -Ou $WGS/00_input/target.pass.vcf.gz \
 | bcftools annotate -a eas_af.tsv.gz -h eas_af.hdr -c CHROM,POS,REF,ALT,INFO/EAS_AF,INFO/ALL_AF -Ou \
 | bcftools csq -f $REF -g $A/Homo_sapiens.GRCh37.87.gff3.gz -p a --ncsq 32 -l --threads 8 -Oz -o target.pass.annot.vcf.gz
tabix -f -p vcf target.pass.annot.vcf.gz
# 3. tables
bcftools query -i 'INFO/CLNSIG!=""' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/GENEINFO\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\t%INFO/CLNSIGCONF\t%INFO/EAS_AF\t%INFO/ALL_AF\t%INFO/BCSQ\t[%GT\t%DP\t%GQ\t%AD]\n' target.pass.annot.vcf.gz > clinvar_all_hits.tsv
bcftools query -i 'INFO/BCSQ~"stop_gained" || INFO/BCSQ~"frameshift" || INFO/BCSQ~"splice_acceptor" || INFO/BCSQ~"splice_donor" || INFO/BCSQ~"start_lost"' \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/EAS_AF\t%INFO/ALL_AF\t%INFO/CLNSIG\t%INFO/BCSQ\t[%GT\t%DP\t%GQ\t%AD]\n' target.pass.annot.vcf.gz > lof_all.tsv
wc -l clinvar_all_hits.tsv lof_all.tsv
echo ANNOT_DONE
