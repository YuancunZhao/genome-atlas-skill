#!/usr/bin/env bash
# Install every command-line tool the pipeline needs into a self-contained conda environment under work/tools/env,
# plus the Java and Python tools that are not on conda. Works on x86_64 and aarch64.
set -euo pipefail
source "$(dirname "$0")/../scripts/env.sh"
MAMBA=${MAMBA:-$(command -v mamba || command -v micromamba || command -v conda)}
[ -n "$MAMBA" ] || { echo "need mamba/micromamba/conda on PATH"; exit 1; }

# H6: same reasoning as setup/02 -- a truncated download leaves a file that later steps will happily use.
# This script cannot reuse that fetch() (it runs before the reference tree exists), so it carries its own.
# Jars and position files are checked for a plausible size and, for archives, for being readable.
_fetch_tool(){
  local url=$1 out=$2 min=${3:-10000}
  if [ -s "$out" ]; then
    local sz; sz=$(stat -c %s "$out" 2>/dev/null || echo 0)
    [ "$sz" -ge "$min" ] && return 0
    echo "WARNING: existing $out is only $sz bytes (< $min); re-downloading" >&2
  fi
  curl -fsSL --retry 3 -o "$out.part" "$url" || { echo "download failed: $url" >&2; rm -f "$out.part"; return 1; }
  local sz; sz=$(stat -c %s "$out.part" 2>/dev/null || echo 0)
  if [ "$sz" -lt "$min" ]; then
    echo "download too small ($sz < $min bytes): $url" >&2; rm -f "$out.part"; return 1
  fi
  case "$out" in *.zip) unzip -tq "$out.part" >/dev/null 2>&1 || { echo "not a valid zip: $url" >&2; rm -f "$out.part"; return 1; } ;; esac
  mv "$out.part" "$out"; echo "fetched $(basename "$out") ($sz bytes)"
}

# 1. core toolchain
"$MAMBA" create -y -p "$TOOLS/env" -c conda-forge -c bioconda \
  "samtools>=1.19" "bcftools>=1.19" htslib mosdepth delly bedtools expansionhunter t1k \
  "openjdk>=17" pysam cyvcf2 pandas numpy scipy pyarrow requests pyyaml pyliftover matplotlib pip
# WhatsHap pins an older Python, so it gets its own environment
"$MAMBA" create -y -p "$TOOLS/whatshap" -c conda-forge -c bioconda whatshap "python=3.10"

# 2. plink2
mkdir -p "$TOOLS"
if [ ! -x "$TOOLS/plink2" ]; then
  ARCH=$(uname -m)
  if [ "$ARCH" = "x86_64" ]; then
    _fetch_tool https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241222.zip "$TOOLS/plink2.zip" 5000000
    unzip -o -q "$TOOLS/plink2.zip" -d "$TOOLS" && chmod +x "$TOOLS/plink2" && rm -f "$TOOLS/plink2.zip"
  else
    echo "plink2 has no aarch64 build: run setup/00b_build_plink2_arm64.sh"
  fi
fi

# 3. Java and Python tools that are not packaged
mkdir -p "$TOOLS/pharmcat" "$TOOLS/flare"
PC=${PHARMCAT_VERSION:-3.4.0}
B=https://github.com/PharmGKB/PharmCAT/releases/download/v$PC
for f in pharmcat-$PC-all.jar pharmcat-preprocessor-$PC.tar.gz pharmcat_positions_$PC.vcf.bgz pharmcat_positions_$PC.vcf.bgz.csi; do
  [ -f "$TOOLS/pharmcat/$f" ] || _fetch_tool "$B/$f" "$TOOLS/pharmcat/$f" 10000
done
tar xzf "$TOOLS/pharmcat/pharmcat-preprocessor-$PC.tar.gz" -C "$TOOLS/pharmcat"
# picard 的 URL 要先从 GitHub API 取；给 _fetch_tool 补上输出路径与下限（少了参数它会当成 out 为空）
[ -f "$TOOLS/picard.jar" ] || _fetch_tool \
  "$(curl -fsSL https://api.github.com/repos/broadinstitute/picard/releases/latest | grep -o 'https[^"]*picard.jar' | head -1)" \
  "$TOOLS/picard.jar" 5000000
[ -f "$TOOLS/flare/flare.jar" ] || _fetch_tool https://faculty.washington.edu/browning/flare.jar "$TOOLS/flare/flare.jar" 100000
[ -d "$TOOLS/Cyrius" ] || git clone -q https://github.com/Illumina/Cyrius.git "$TOOLS/Cyrius"
[ -d "$TOOLS/SMNCopyNumberCaller" ] || git clone -q https://github.com/Illumina/SMNCopyNumberCaller.git "$TOOLS/SMNCopyNumberCaller"
"$TOOLS/env/bin/pip" install -q -r "$TOOLS/Cyrius/requirements.txt" -r "$TOOLS/pharmcat/preprocessor/requirements.txt"
# haplogrep3: the upstream tags are v-prefixed (verified: v3.2.2 ships haplogrep3-3.2.2-linux.zip;
# today's latest is v3.3.2 whose assets are .tar.gz), so the old latest-download + pinned-filename
# URL no longer resolves. Pin the tag explicitly and verify the archive.
HV=${HAPLOGREP3_VERSION:-3.2.2}
if [ ! -x "$TOOLS/haplogrep3" ]; then
  _fetch_tool "https://github.com/genepi/haplogrep3/releases/download/v$HV/haplogrep3-$HV-linux.zip" "$TOOLS/haplogrep3.zip" 10000000
  unzip -o -q "$TOOLS/haplogrep3.zip" -d "$TOOLS" && chmod +x "$TOOLS/haplogrep3" && rm -f "$TOOLS/haplogrep3.zip"
fi
# 06_mtdna.py runs `classify --tree phylotree-rcrs@17.2`; haplogrep3 resolves tree IDs only against
# the trees/ directory next to its own jar (a direct path in --tree is rejected, verified 2026-10-07),
# so a fresh install must fetch the tree through the tool itself, then prove it actually landed:
# a missing/partial tree would otherwise surface only as a classify failure much later.
H3TREE="$TOOLS/trees/phylotree-rcrs/17.2"
[ -f "$H3TREE/tree.yaml" ] || "$TOOLS/haplogrep3" install-tree phylotree-rcrs@17.2
[ -f "$H3TREE/tree.yaml" ] && [ -s "$H3TREE/rcrs.fasta" ] || {
  echo "haplogrep3 tree phylotree-rcrs@17.2 missing tree.yaml/rcrs.fasta under $H3TREE after install-tree" >&2
  exit 1
}
# --- Optional cross-check tools (never part of the pipeline; §4: two tools agreeing is not
# independent validation, it only rules out a copied-wrong table). Installed into the same conda env
# as everything else, so one environment reproduces the whole tool set.
if [ "${WITH_CROSSCHECK:-1}" = "1" ]; then
  # Yleaf: an independent Y-haplogroup caller. Ships its own hg19/hg38/t2t tables (449 MB) and can use
  # several published trees (yfull, yfull_v10, ftdna, isogg), which is what makes it useful as a check
  # against step 05. Verified here as Yleaf 4.1.4.
  if [ ! -x "$TOOLS/env/bin/Yleaf" ]; then
    "$TOOLS/env/bin/pip" install -q "git+https://github.com/genid/Yleaf.git" \
      || echo "  note: Yleaf not installed (no network?); rerun with WITH_CROSSCHECK=1 when reachable"
  fi
  if [ -x "$TOOLS/env/bin/Yleaf" ]; then
    echo "  Yleaf (cross-check only) -> $TOOLS/env/bin/Yleaf"
    # Yleaf ships its position tables but downloads a ~3 GB reference genome on first use. This project
    # already has hg19 -- the same FASTA every other step was verified against -- so point Yleaf at it:
    # no second copy, no download, no chance of the check running on a different reference.
    YCFG=$(ls "$TOOLS"/env/lib/python*/site-packages/yleaf/config.txt 2>/dev/null | head -1)
    if [ -n "$YCFG" ] && ! grep -q "^full hg19 genome fasta location = /" "$YCFG"; then
      sed -i "s|^full hg19 genome fasta location = .*|full hg19 genome fasta location = $TOOLS/../data/ref/b37/human_g1k_v37.fasta|" "$YCFG"
      echo "  Yleaf: hg19 reference pointed at the pipeline's own FASTA (avoids the 3 GB download)"
    fi
  fi
fi

echo "tools installed under $TOOLS"
