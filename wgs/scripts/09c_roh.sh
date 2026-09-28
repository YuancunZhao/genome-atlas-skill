#!/bin/bash
# Runs of homozygosity (H5): this product previously had no producer in the repository -- it was computed
# by hand, so the chain could not be re-run and a missing file surfaced late as an empty table.
#
# The parameters are not chosen here: all seven are read off the delivered run's own log
# (09_misc/plink_roh/roh_plink.log), so re-running reproduces the same definition of a run rather than a
# similar-looking one.
#
# It must be PLINK 1.9, not 2.x: run-recovery on the server showed a plink2 invocation in the sibling
# roh.log that failed with "Unrecognized flag ('--homozyg')" -- plink2 has no such flag. The delivered
# product is roh_plink.hom, which is 1.9 output. Using plink2 here would have looked right and produced
# nothing.
#
# Input: 02_complete/target.1kg -- the same 1000G-augmented set the ROH was originally computed on.
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
BASE=${WGS}/02_complete/target.1kg
W=$WGS/09_misc/plink_roh; mkdir -p "$W" "$WGS/logs"
PL=${PLINK19:-$TOOLS/plink19/plink}

[ -x "$PL" ] || { echo "PLINK 1.9 not executable: $PL (ROH needs --homozyg, which plink2 does not have)" >&2; exit 1; }
[ -f "${BASE}.pgen" ] || { echo "missing input ${BASE}.pgen (run 03_complete_set.py first)" >&2; exit 1; }

cd "$W"
# PLINK 1.9 读不了 pgen，而 02_complete 是 pfile。交付时这个目录里有一份转换好的 target.1kg.{bed,bim,fam}
# （就是当初手工转换留下的），这里在缺失时用同一份 pfile 重新生成，而不是要求用户自己准备二进制格式。
if [ ! -f target.1kg.bed ]; then
  P2=${PLINK2:-$TOOLS/plink2}
  [ -x "$P2" ] || { echo "plink2 needed to convert the pfile to bed: $P2" >&2; exit 1; }
  "$P2" --pfile "$BASE" --make-bed --out target.1kg >/dev/null 2>&1 \
    || { echo "could not convert $BASE to bed" >&2; exit 1; }
  echo "converted $(basename "$BASE") to target.1kg.{bed,bim,fam}"
fi

"$PL" --bfile target.1kg \
  --homozyg \
  --homozyg-density 50 \
  --homozyg-gap 100 \
  --homozyg-kb 1000 \
  --homozyg-snp 50 \
  --homozyg-window-het 1 \
  --homozyg-window-missing 5 \
  --homozyg-window-snp 50 \
  --out roh_plink 2>&1 | tee "$WGS/logs/roh_plink.log"

[ -s roh_plink.hom ] || { echo "PLINK produced no ROH output" >&2; exit 1; }
wc -l < roh_plink.hom | sed 's/^/roh segments: /'
echo "ROH_DONE -> $W/roh_plink.hom"
