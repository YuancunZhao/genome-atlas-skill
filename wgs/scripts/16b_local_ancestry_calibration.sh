#!/bin/bash
# Calibration: score 20 CHB and 20 CHS reference individuals with the same north/south panels, holding them out.
set -euo pipefail
source "$(dirname "$0")/env.sh"
JAVA=$TOOLS/env/bin/java
W=$WGS/12_localanc; cd $W
R=$REF_DIR
awk 'NR>1 && $5=="CHB"{print $1}' $R/all_phase3.psam | head -20 > holdout.ids
awk 'NR>1 && $5=="CHS"{print $1}' $R/all_phase3.psam | head -20 >> holdout.ids
grep -vwFf holdout.ids ref.panel > ref.calib.panel
cut -f2 ref.calib.panel | sort | uniq -c
for c in 1 2 6 22; do
  # the 1000G phase3 VCFs have a minimal header (no contigs, GT only) and bcftools refuses
  # sample-subset mode on them ("Undefined tags in the header"), so cut the sample columns with awk
  cols= || true  # grep -m1 closes the pipe; zcat dies SIGPIPE and pipefail would abort the script$(zcat $R/imp/ALL.chr$c.vcf.gz | grep -m1 "^#CHROM" | tr '\t' '\n' | grep -nxFf holdout.ids | cut -d: -f1 | paste -sd,)
  zcat $R/imp/ALL.chr$c.vcf.gz | awk -v c="$cols" 'BEGIN{n=split(c,C,","); for(i=1;i<=9;i++) k[i]=1; for(i=1;i<=n;i++) k[C[i]]=1} { s=$1; for(i=2;i<=NF;i++) if(k[i]) s=s "\t" $i; print s }' | bgzip -@ 4 > ho.$c.vcf.gz
  tabix -f ho.$c.vcf.gz
  $JAVA -Xmx40g -jar $FLARE ref=$R/imp/ALL.chr$c.vcf.gz ref-panel=ref.calib.panel gt=ho.$c.vcf.gz \
    map=$R/imp/maps/plink.chr$c.GRCh37.map out=calib.$c gen=100 nthreads=12 > $WGS/logs/flare_calib.$c.log 2>&1
  rm -f ho.$c.vcf.gz*
  echo "chr$c calibration done"
done
echo CALIB_DONE
