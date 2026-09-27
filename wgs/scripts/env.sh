# source this before any shell step: reads config.yaml through wgsconfig.py so shell and Python agree.
# usage:  source scripts/env.sh
# BASH_SOURCE is bash-only; when sourced from zsh it is unset and $0 is the sourced file
_here="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
eval "$(python3 - "$_here" <<'PY'
import sys, pathlib, shlex
sys.path.insert(0, sys.argv[1])
import wgsconfig as c
out = {
  "PROJ": str(c.P), "WGS": str(c.W), "REF_DIR": str(c.REF), "TOOLS": str(c.TOOLS), "PGSDIR": str(c.PGS),
  "SAMPLE": c.SAMPLE, "SEX": c.SEX, "THREADS": str(c.THREADS), "MEM_GB": str(c.MEM_GB),
  "CRAM": c.READS, "YBAM": c.Y_READS, "VENDOR_VCF": c.VENDOR_VCF,
  "VCF_RAW": c.VCF_RAW,
  "REF": c.FASTA, "REF_PATH": c.FASTA, "REF38": c.FASTA38,
  "CHAIN": c.CHAIN_19_38, "CHAIN_BACK": c.CHAIN_38_19,
  "KG_PFILE": c.KG_PFILE, "KG_VCF_DIR": c.KG_VCF_DIR, "CLINVAR": c.CLINVAR, "GFF3": c.GFF3,
  "PRS_REF": c.PRS_REF, "AADR": c.AADR, "YTREE": c.YTREE,
  "PLINK2": c.PLINK2, "PICARD": c.PICARD, "PHARMCAT_DIR": c.PHARMCAT_DIR, "FLARE": c.FLARE,
  "JAVA": c.JAVA, "WHATSHAP": c.WHATSHAP,
  "MIN_DP": str(c.MIN_DP), "MIN_MQ": str(c.MIN_MQ), "MIN_BQ": str(c.MIN_BQ),
  # Local-ancestry sources and labels come from the validated options, not from the module defaults:
  # an unconfigured sample exports empty values and LOCAL_ENABLED=0, instead of inheriting a panel it
  # never asked for. SUPERPOP/SUBPOPS keep their previous meaning for PRS and other existing consumers.
  "LA_A": ",".join(c.OPT["la_a"]), "LA_B": ",".join(c.OPT["la_b"]), "LA_CTRL": ",".join(c.OPT["la_control"]),
  "LA_LAB_A": c.OPT["la_labels"][0] if len(c.OPT["la_labels"]) > 0 else "",
  "LA_LAB_B": c.OPT["la_labels"][1] if len(c.OPT["la_labels"]) > 1 else "",
  "AXIS": ",".join(c.OPT["axis_pops"]),
  "SUPERPOP": c.SUPERPOP, "SUBPOPS": ",".join(c.SUBPOPS),
  # validated switches and thresholds: shell steps and Python read the same values (§7 AN0)
  "REGIONAL_ENABLED": "1" if c.REGIONAL_ENABLED else "0",
  "AADR_ENABLED": "1" if c.AADR_ENABLED else "0",
  "LOCAL_ENABLED": "1" if c.LOCAL_ENABLED else "0",
  "MIN_CR_MODERN": str(c.MIN_CR_MODERN), "MIN_CR_TARGET": str(c.MIN_CR_TARGET),
  "MIN_CR_ANCIENT": str(c.MIN_CR_ANCIENT), "MIN_PROJECTION_SNPS": str(c.MIN_PROJECTION_SNPS),
  "MIN_GROUP_N": str(c.MIN_GROUP_N), "CALIB_N": str(c.CALIB_N), "CALIB_SEED": str(c.CALIB_SEED),
  "CALIB_POPS": ",".join(c.CALIB_POPS), "CALIB_CHROMS": ",".join(c.CALIB_CHROMS),
  "AADR_MODERN": ",".join(c.AADR_MODERN), "AADR_ANCIENT_PREFIX": ",".join(c.AADR_ANCIENT_PREFIX),
}
for k, v in out.items():
    print(f"export {k}={shlex.quote(v)}")
PY
)"
unset _here
# tools installed by setup/00_install_tools.sh live in a conda env; put it first if it exists
if [ -d "$TOOLS/env/bin" ]; then export PATH="$TOOLS/env/bin:$PATH"; fi
