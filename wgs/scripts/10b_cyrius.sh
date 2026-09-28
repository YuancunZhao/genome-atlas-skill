#!/bin/bash
# CYP2D6 star alleles (H5): this product had no producer in the repository, and 11_pgx_extra plus
# 10_pharmcat both read it -- PharmCAT does not resolve CYP2D6 on its own, which is why a separate caller
# is needed at all.
#
# Recovered the same way as SMN: 06_pgx/cyrius/manifest.txt names the input (00_input/norm.cram) and
# cyrius.log records the sample being processed as "norm". Arguments follow the tool's documented form.
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
W=$WGS/06_pgx/cyrius; mkdir -p "$W" "$WGS/logs"; cd "$W"
# 入口是 star_caller.py（不是 cyrius.py）：Cyrius 的仓库把它放在顶层，其余是库模块。
CY=${CYRIUS:-$TOOLS/Cyrius/star_caller.py}
PY=${CYRIUSPY:-$TOOLS/env/bin/python}

[ -f "$CY" ] || { echo "Cyrius entry point not found: $CY" >&2; exit 1; }
[ -f "$CRAM" ] || { echo "input alignment not found: $CRAM" >&2; exit 1; }

echo "$CRAM" > manifest.txt
"$PY" "$CY" --manifest manifest.txt --genome 37 --prefix target --outDir . --threads "${THREADS:-8}" \
  --reference "$REF" > "$WGS/logs/cyrius.log" 2>&1

[ -s target.tsv ] || { echo "Cyrius produced no target.tsv (see logs/cyrius.log)" >&2; exit 1; }
echo "CYRIUS_DONE -> $W/target.tsv + target.json"
