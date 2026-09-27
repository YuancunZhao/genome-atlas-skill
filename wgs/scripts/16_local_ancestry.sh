#!/bin/bash
# Local ancestry along the genome with FLARE 0.6. Source and control panels come from the configuration
# (local_ancestry_a/b/control and their labels); nothing here assumes East Asia or a particular set of
# continental controls. The control output is called a control panel and is not asserted to be a noise
# floor by itself (7.3).
set -euo pipefail
source "$(dirname "$0")/env.sh"
JAVA=${JAVA:-java}
W=$WGS/12_localanc; mkdir -p $W; cd $W
R=$REF_DIR
# Panel map: sample <tab> panel. Source labels come from the configuration; each control population
# keeps its own name instead of being folded into a hard-coded "European"/"SouthAsian" (AN3). Panels
# are used as FLARE ancestries, so their order here is the AN index order.
awk -v a="$LA_A" -v b="$LA_B" -v ctrl="$LA_CTRL" -v la="$LA_LAB_A" -v lb="$LA_LAB_B" '
  BEGIN{ n=split(ctrl,C,","); for(i=1;i<=n;i++) if(C[i]!="") ctrl_lab[C[i]]=C[i];
         na=split(a,A,","); for(i=1;i<=na;i++) if(A[i]!="") map[A[i]]=la;
         nb=split(b,B,","); for(i=1;i<=nb;i++) if(B[i]!="") map[B[i]]=lb }
  NR>1 { p=map[$5]; if(p=="") p=ctrl_lab[$5]; if(p!="") print $1"\t"p }' $R/all_phase3.psam > ref.panel
cut -f2 ref.panel | sort | uniq -c
for c in $(seq 1 22); do
  [ -f $WGS/10_phase/ph.$c.vcf.gz ] || { echo "chr$c not phased yet, skipping"; continue; }
  $JAVA "-Xmx${MEM_GB%.*}g" -jar $FLARE ref=$R/imp/ALL.chr$c.vcf.gz ref-panel=ref.panel gt=$WGS/10_phase/ph.$c.vcf.gz \
    map=$R/imp/maps/plink.chr$c.GRCh37.map out=la.$c probs=true gen=100 nthreads=$THREADS > $WGS/logs/flare.$c.log 2>&1
  echo "chr$c done: $(bcftools view -H la.$c.anc.vcf.gz 2>/dev/null | wc -l) markers"
done
echo FLARE_DONE
