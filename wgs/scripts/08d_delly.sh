#!/bin/bash
# Structural variants (H5): this product had no producer; the command below is the delivered run's own,
# recovered verbatim from the first line of 08_sv/delly/delly.log:
#
#   delly sr -g <ref>/b37/human_g1k_v37.fasta -o <work>/08_sv/delly/target.sv.bcf -h 16 <work>/08_sv/main_contigs/norm.main.cram
#
# "sr" is the split-read plus paired-end caller, -h is the thread count, and the input is the main-contig
# alignment produced by 08a -- see that script for why the unrestricted CRAM is not used.
#
# NOTE: the input is currently UNVERIFIED -- see the header of 08a_main_contigs.sh. Delly's own command
# line is recovered verbatim, but the file it was pointed at cannot be reproduced from what remains,
# so re-running this produces SV calls on a different input than the delivered ones.
#
# This step takes hours on a whole genome (the delivered run went 01:44 to 03:37), so it is wired into
# run_all.sh but not part of any quick check. It fails loudly rather than leaving a stale bcf behind.
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
W=$WGS/08_sv/delly; mkdir -p "$W" "$WGS/logs"; cd "$W"
IN=${MAIN_CRAM:-$WGS/08_sv/main_contigs/norm.main.cram}
DELLY=${DELLY:-$TOOLS/delly}

[ -x "$DELLY" ] || { echo "delly not executable: $DELLY" >&2; exit 1; }
[ -f "$IN" ] || { echo "main-contig alignment not found: $IN (run 08a_main_contigs.sh first)" >&2; exit 1; }

# 复审 P0：不要默认在未验证的输入上产出一份"看起来像交付"的 SV 调用集。08a 生成的文件带 provenance
# 标记；看到 verified_against_delivered=false 就停下来，除非调用者明确接受这个差别。
if [ -f "$IN.provenance.json" ]; then
  if grep -q '"verified_against_delivered": *false' "$IN.provenance.json" 2>/dev/null; then
    if [ "${ALLOW_UNVERIFIED_INPUT:-0}" != "1" ]; then
      echo "ERROR: $IN was produced by an unverified reconstruction (see $IN.provenance.json)." >&2
      echo "Running Delly on it yields structural variants that are NOT the delivered ones." >&2
      echo "Set ALLOW_UNVERIFIED_INPUT=1 if that is understood and intended." >&2
      exit 2
    fi
    echo "WARNING: running on an unverified input (ALLOW_UNVERIFIED_INPUT=1); output is not the delivered SV set" >&2
  fi
fi

"$DELLY" sr -g "$REF" -o target.sv.bcf -h "${THREADS:-8}" "$IN" > "$WGS/logs/delly.log" 2>&1
[ -s target.sv.bcf ] || { echo "delly produced no target.sv.bcf (see logs/delly.log)" >&2; exit 1; }
bcftools index -f target.sv.bcf 2>/dev/null || true
echo "DELLY_DONE -> $W/target.sv.bcf"
