#!/bin/bash
# Repeat expansions (H5): this product had no producer in the repository.
#
# Recovered from the delivered run's own log (08_sv/eh/eh.log) rather than reconstructed:
#   [Starting ExpansionHunter v5.0.0]
#   [Loading variant catalog from disk .../tools/eh/variant_catalog_grch37.json]
#   [Running sample analysis in seeking mode]
# The mode matters: seeking mode is what lets ExpansionHunter find the repeat without being told a
# region per locus, and the recovered invocation ran with it.
#
# The region file is an input, not something derived here. It holds 37 windows of the catalog's loci
# +/-1000 bp; the delivered bed has more windows than the catalog has entries, so it was maintained
# separately and is now kept in panel/eh_loci_regions.bed alongside the other reference panels. Deriving
# it from today's catalog would silently drop whatever the extra windows were for.
set -euo pipefail
source "$(dirname "$0")/env.sh"
S=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
W=$WGS/08_sv/eh; mkdir -p "$W" "$WGS/logs"; cd "$W"
BED=${EH_REGIONS_BED:-$(cd "$S/.." && pwd)/panel/eh_loci_regions.bed}
CAT=${EH_CATALOG:-$TOOLS/eh/variant_catalog_grch37.json}
# 二进制在解压出的版本目录里（ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter）。
# 用 glob 找而不是写死版本号：换版本时不用改脚本，但找不到就明确失败。
if [ -n "${EXPANSIONHUNTER:-}" ]; then
  EH=$EXPANSIONHUNTER
else
  EH=$(ls -1 "$TOOLS"/eh/ExpansionHunter-*/bin/ExpansionHunter 2>/dev/null | head -1)
fi

[ -x "$EH" ] || { echo "ExpansionHunter not executable: $EH" >&2; exit 1; }
[ -f "$CAT" ] || { echo "variant catalog not found: $CAT" >&2; exit 1; }
[ -f "$BED" ] || { echo "region list not found: $BED" >&2; exit 1; }

# 只提取目标区域：全 CRAM 送给 ExpansionHunter 没有意义，也会慢到不可用。
if [ ! -s eh_regions.bam ] || [ eh_regions.bam -ot "$CRAM" ]; then
  samtools view -@ "${THREADS:-8}" -b -L "$BED" "$CRAM" > eh_regions.bam
  samtools index eh_regions.bam
fi

# --sex 决定性染色体位点（AR、FMR1、AFF2 在 X 上）的期望拷贝数。来自证据而非声明。
#
# 这一处与交付产物**不一致，且交付的那版与证据不符**：交付的 target.vcf 把 AR 写成 1/1:27/27
# （二倍体），而本样本的 Y/常染色体深度比 0.188 明确指向男性，男性 X 非 PAR 区应报单倍体。
# 本脚本按证据报单倍体，因此重跑会在 X 连锁位点上与交付产物不同——这是**修正**而不是漂移。
# 需要与交付完全一致的输出时，设 EXPANSIONHUNTER_SEX=female 可复现旧行为。
if [ -n "${EXPANSIONHUNTER_SEX:-}" ]; then
  SEXARG=$EXPANSIONHUNTER_SEX          # 复现旧行为的逃生口，见上面的说明
else
  SEXARG=$(python3 -c "
import sys; sys.path.insert(0, '$S')
import _sex_evidence as se
print(se.evidence_sex(log=lambda *a: None)[0])" 2>/dev/null || echo "${SEX:-female}")
fi
echo "ExpansionHunter --sex $SEXARG"

"$EH" --reference "$REF" --variant-catalog "$CAT" --reads eh_regions.bam \
  --output-prefix target --sex "$SEXARG" --analysis-mode seeking \
  > "$WGS/logs/eh.log" 2>&1

[ -s target.vcf ] || { echo "ExpansionHunter produced no target.vcf (see logs/eh.log)" >&2; exit 1; }

# eh_summary.tsv：报告与 30 用的紧凑表（每 locus 一行的 unit/genotype/max_allele）。
python3 - "$W" <<'PY'
import csv, sys, pathlib
w = pathlib.Path(sys.argv[1])
rows = []
with open(w / "target.vcf", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        # genotype 列是重复次数（REPCN，如 "27/27"；男性 X 连锁位点为单倍型 "27"），不是 GT 的
        # 等位基因索引（"1/1"）。30 用 max_allele 对逐位点阈值，写成索引会把每个位点都变成 1。
        # EH 原生 REPCN 两位点时顺序不定（实测 10/2），交付的 eh_summary 一律升序——按交付格式排序。
        _m = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
        _rep = _m.get("REPCN", "")
        al = [int(x) for x in _rep.replace("|", "/").split("/") if x.isdigit()]
        gt = "/".join(str(x) for x in sorted(al)) if al else _rep
        rows.append((info.get("REPID", f[2]), info.get("RU", ""), gt,
                     str(max(al)) if al else ""))
# 按 locus 名排序。交付的 eh_summary.tsv 就是字母序（AFF2, AR, ATN1, ATXN1 ...），而 ExpansionHunter
# 的 VCF 是基因组顺序；两者不同，按 VCF 顺序写会得到一个内容相同但逐行都不同的文件。
# 我先前以为"保持原始顺序"更好而删掉排序，实测与交付件逐行不符，遂恢复——以交付件为准。
rows.sort()
with open(w / "eh_summary.tsv", "w", newline="", encoding="utf-8") as out:
    wr = csv.writer(out, delimiter="\t")
    wr.writerow(["locus", "unit", "genotype", "max_allele"])
    wr.writerows(rows)
print(f"eh_summary.tsv: {len(rows)} loci")
PY
echo "EH_DONE -> $W/target.json + eh_summary.tsv"
