#!/bin/bash
# SMN1/SMN2 copy number (H5): this product had no producer in the repository.
#
# Recovered from the delivered run rather than reconstructed: 08_sv/smn/manifest.txt lists exactly the
# input that was used (work/wgs/00_input/norm.cram) and smn.log records the tool processing sample "norm",
# which is the stem of that filename. The invocation follows the tool's own documented form with the same
# manifest/genome/prefix/outDir/threads/reference arguments Cyrius takes -- both are from the same authors
# and ship the same interface.
#
# The M1 report conclusion depends on this product, so a missing step here is not cosmetic.
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
W=$WGS/08_sv/smn; mkdir -p "$W" "$WGS/logs"; cd "$W"
SMN=${SMNCALLER:-$TOOLS/SMNCopyNumberCaller/smn_caller.py}
PY=${SMNPY:-$TOOLS/env/bin/python}

[ -f "$SMN" ] || { echo "SMNCopyNumberCaller not found: $SMN" >&2; exit 1; }
[ -f "$CRAM" ] || { echo "input alignment not found: $CRAM" >&2; exit 1; }

echo "$CRAM" > manifest.txt          # 与交付件同样的形态：一行一个绝对路径
"$PY" "$SMN" --manifest manifest.txt --genome 37 --prefix target --outDir . --threads "${THREADS:-8}" \
  --reference "$REF" > "$WGS/logs/smn.log" 2>&1

[ -s target.tsv ] || { echo "SMNCopyNumberCaller produced no target.tsv (see logs/smn.log)" >&2; exit 1; }
echo "SMN_DONE -> $W/target.tsv + target.json"
