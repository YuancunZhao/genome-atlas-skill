#!/bin/bash
# Structural variants (H5): `delly sr` on the main-contig alignment from 08a. The command is the
# delivered run's own, recovered verbatim from the first line of 08_sv/delly/delly.log:
#
#   delly sr -g <ref>/b37/human_g1k_v37.fasta -o <work>/08_sv/delly/target.sv.bcf -h 16 <work>/08_sv/main_contigs/norm.main.cram
#
# "sr" is the split-read plus paired-end caller, -h is the thread count, and the input is the
# main-contig alignment produced by 08a -- see that script for why the unrestricted CRAM is not used.
# Verified 2026-09-28 (see 08a_main_contigs.sh header): 08a's extract + this command reproduced the
# delivered target.sv.bcf exactly (42,909/42,909 records; md5 differs only in ##fileDate).
#
# This step takes hours on a whole genome (the delivered run went 01:44 to 03:37), so it is wired
# into run_all.sh but not part of any quick check. It fails loudly rather than leaving a stale bcf
# behind. The delivered tree's target.sv.bcf is a delivery artifact: like 08a, this script refuses
# to silently overwrite a bcf that carries no provenance marker of its own -- rm it first or set
# DELLY_OVERWRITE=1 to replace it deliberately.
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
W=$WGS/08_sv/delly; mkdir -p "$W" "$WGS/logs"; cd "$W"
IN=${MAIN_CRAM:-$WGS/08_sv/main_contigs/norm.main.cram}
DELLY=${DELLY:-$TOOLS/delly}

[ -x "$DELLY" ] || { echo "delly not executable: $DELLY" >&2; exit 1; }
[ -f "$IN" ] || { echo "main-contig alignment not found: $IN (run 08a_main_contigs.sh first)" >&2; exit 1; }

# 覆盖保护（见头注释）：不带 provenance 的现存 target.sv.bcf 视为交付原件——留下并说明，而不是默默重算。
if [ -s target.sv.bcf ] && [ ! -f target.sv.bcf.provenance.json ]; then
  echo "refusing to overwrite the existing target.sv.bcf: it carries no provenance marker, so it is" >&2
  echo "most likely the delivered original. To replace it deliberately: rm target.sv.bcf && re-run," >&2
  echo "or set DELLY_OVERWRITE=1." >&2
  [ "${DELLY_OVERWRITE:-0}" = "1" ] || exit 0
  echo "DELLY_OVERWRITE=1: overwriting anyway" >&2
fi

"$DELLY" sr -g "$REF" -o target.sv.bcf -h "${THREADS:-8}" "$IN" > "$WGS/logs/delly.log" 2>&1
[ -s target.sv.bcf ] || { echo "delly produced no target.sv.bcf (see logs/delly.log)" >&2; exit 1; }
bcftools index -f target.sv.bcf 2>/dev/null || true

# 产出即记录本次运行（同 08a：运行事实与“复现历史交付”分开）。
cat > target.sv.bcf.provenance.json <<JSON
{
  "producer": "08d_delly.sh",
  "command": "delly sr -g $REF -o target.sv.bcf -h ${THREADS:-8} $IN",
  "input_cram": "$IN",
  "input_provenance": $([ -f "$IN.provenance.json" ] && printf '"%s"' "$IN.provenance.json" || echo null),
  "date": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON
echo "DELLY_DONE -> $W/target.sv.bcf"
