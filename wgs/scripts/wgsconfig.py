"""Shared configuration for the genome-atlas-wgs pipeline.

Reads config.yaml at the repository root (copy config.example.yaml first). Nothing personal lives in the
repository: raw reads, references, tools and every output live under the work directory, which is git-ignored.

Layout created under WORK:
    data/raw    the CRAM/BAM/FASTQ and any delivered VCF, plus symlinks without spaces in the name
    data/ref    reference genomes and panels (setup/download_references.sh)
    data/pgs    PGS Catalog scoring files
    tools       plink2, PharmCAT, Picard, Cyrius, FLARE, SMNCopyNumberCaller, T1K indexes
    wgs/NN_*    one directory per analysis step
    report      the finished HTML
"""
import os, pathlib

try:
    import yaml
except ImportError:
    yaml = None

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _config_path():
    """Which config file to read. $WGS_CONFIG wins when set (absolute, or relative to the pipeline
    root); otherwise config.yaml at the root. One config per sample/run, no shared mutable defaults."""
    env = os.environ.get("WGS_CONFIG", "").strip()
    if not env:
        return ROOT / "config.yaml"
    p = pathlib.Path(env).expanduser()
    return (p if p.is_absolute() else ROOT / p).resolve()


def _load():
    f = _config_path()
    if not f.exists():
        return {}
    if yaml:
        return yaml.safe_load(open(f)) or {}
    out = {}
    for line in open(f):
        line = line.split("#")[0].strip()
        if ":" in line:
            k, v = line.split(":", 1)
            out[k.strip()] = v.strip().strip('"').strip("'")
    return out


CFG = _load()
_g = lambda k, d: CFG.get(k, os.environ.get(k.upper(), d))

SAMPLE = str(_g("sample_id", "SAMPLE"))
NAME_ZH = str(_g("name_zh", SAMPLE))
NAME_EN = str(_g("name_en", SAMPLE))
SEX = str(_g("sex", "male")).lower()

# P is the work root: every path below hangs off it, and every script writes only inside it.
# Relative work/ref/tools dirs are anchored at the pipeline root (this wgs/ directory), never
# at the caller's CWD: resolving against the CWD made every output path depend on which
# directory a step happened to be launched from.
def _anchored(v):
    p = pathlib.Path(v)
    return (p if p.is_absolute() else ROOT / p).expanduser().resolve()

P = _anchored(_g("work_dir", ROOT / "work"))
W = P / "wgs"                      # analysis outputs, one directory per step
DATA = P / "data"
RAW = DATA / "raw"
REF = _anchored(_g("ref_dir", DATA / "ref"))
TOOLS = _anchored(_g("tools_dir", P / "tools"))
PGS = DATA / "pgs"
REPORT = P / "report"
for d in (DATA, RAW, REF, TOOLS, PGS, W, REPORT):
    d.mkdir(parents=True, exist_ok=True)

# ---- input reads. Give either a CRAM/BAM or a pair of FASTQ files; setup/01_align.sh turns FASTQ into a CRAM.
READS = str(_g("reads", RAW / f"{SAMPLE}.cram"))          # CRAM or BAM, coordinate-sorted and indexed
FASTQ1 = str(_g("fastq1", ""))                            # optional, only used by setup/01_align.sh
FASTQ2 = str(_g("fastq2", ""))
VENDOR_VCF = str(_g("vendor_vcf", ""))                    # optional: a VCF delivered with the reads
Y_READS = str(_g("y_reads", ""))                          # optional: a separate Y-only BAM

# ---- references
FASTA = str(_g("fasta", REF / "b37" / "human_g1k_v37.fasta"))
FASTA38 = str(_g("fasta38", REF / "hg38" / "Homo_sapiens_assembly38.fasta"))
CHAIN_19_38 = str(_g("chain_hg19_hg38", REF / "chain" / "hg19ToHg38.over.chain.gz"))
CHAIN_38_19 = str(_g("chain_hg38_hg19", REF / "chain" / "hg38ToHg19.over.chain.gz"))
KG_PFILE = str(_g("kg_pfile", REF / "all_phase3"))        # 1000G phase3 plink2 bundle, GRCh37
KG_VCF_DIR = str(_g("kg_vcf_dir", REF / "imp"))           # per-chromosome 1000G VCFs + genetic maps + beagle.jar
CLINVAR = str(_g("clinvar", REF / "clinvar_grch37.vcf.gz"))
GFF3 = str(_g("gff3", REF / "annot" / "Homo_sapiens.GRCh37.87.gff3.gz"))
GNOMAD_CONSTRAINT = str(_g("gnomad_constraint", REF / "annot" / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"))
AADR = str(_g("aadr_prefix", REF / "aadr" / "v66.p1_HO.aadr.patch.PUB"))
SPRIME_DIR = str(_g("sprime_dir", REF / "archaic" / "mendeley_data"))
YTREE = str(_g("ytree_dir", REF / "ytree"))
PRS_REF = str(_g("prs_ref", REF / "prs" / "kg_all"))      # 1000G subset with chrom:pos IDs, built by setup/03

# ---- tools
PLINK2 = str(_g("plink2", TOOLS / "plink2"))
HAPLOGREP3 = str(_g("haplogrep3", TOOLS / "haplogrep3"))
PHARMCAT_DIR = str(_g("pharmcat_dir", TOOLS / "pharmcat"))
PICARD = str(_g("picard", TOOLS / "picard.jar"))
CYRIUS = str(_g("cyrius_dir", TOOLS / "Cyrius"))
SMN_CALLER = str(_g("smn_dir", TOOLS / "SMNCopyNumberCaller"))
FLARE = str(_g("flare", TOOLS / "flare" / "flare.jar"))
JAVA = str(_g("java", "java"))                            # needs Java 17+ for PharmCAT and FLARE
WHATSHAP = str(_g("whatshap", "whatshap"))
THREADS = str(_g("threads", 8))
MEM_GB = str(_g("mem_gb", 30))

# ---- analysis parameters
import ancestry_data as _ad  # normalize_build: 唯一的 build 校验点
BUILD = _ad.normalize_build(_g("build", "GRCh37"))
# 校验只有一处（ancestry_data.normalize_build）：注释里说"只支持 GRCh37"而代码不校验，等于没限制 ——
# `build: GRCh38` 会被静默接受，而 FASTA/chain/1000G/AADR/注释全是 GRCh37，位点按错误坐标系解释。
MIN_DP = int(_g("callable_min_depth", 8))                 # callable mask: minimum depth
MIN_MQ = int(_g("callable_min_mapq", 20))                 # callable mask: minimum mapping quality
MIN_BQ = int(_g("min_base_quality", 20))                  # read depth/allele counting: minimum BASE quality
                                                          # (mosdepth -Q; distinct from mapping quality above)
MEAN_DEPTH = _g("mean_depth", None)                       # filled in by 01_qc; used to normalise CNV depth ratios

# ---- reference populations (1000 Genomes labels). Defaults suit an East Asian sample.
SUPERPOP = str(_g("ref_superpop", "EAS"))
SUBPOPS = list(CFG.get("ref_subpops", ["CHB", "CHS"]))            # second percentile column, and the PCA highlight
LA_NORTH = list(CFG.get("local_ancestry_a", ["CHB", "JPT"]))      # local-ancestry panel A
LA_SOUTH = list(CFG.get("local_ancestry_b", ["CDX", "KHV"]))      # local-ancestry panel B
LA_CONTROL = list(CFG.get("local_ancestry_control", ["CEU", "GBR", "GIH", "PJL"]))  # noise-floor panels
LA_LABELS = list(CFG.get("local_ancestry_labels", ["NorthEA", "SouthEA"]))
AXIS = list(CFG.get("axis_pops", ["CHS", "CHB"]))                 # per-chromosome axis: 0 = first, 1 = second
AADR_MODERN = list(CFG.get("aadr_modern", []))                    # extra present-day AADR groups to keep
AADR_ANCIENT_PREFIX = list(CFG.get("aadr_ancient_prefix", []))    # ancient AADR group-ID prefixes to project

# ---- ancestry switches and validated parameters (AN0). These come from read_options(), which decides
# "explicitly configured" from the keys written in config.yaml itself -- a shared default is never
# mistaken for a user choice. The constants above keep their previous meaning and defaults so that
# PRS and every existing consumer behave exactly as before.
import sys as _sys
_sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
OPT = _ad.read_options(CFG)
REGIONAL_ENABLED = OPT["regional_enabled"]
AADR_ENABLED = OPT["aadr_enabled"]
LOCAL_ENABLED = OPT["local_enabled"]
MIN_CR_MODERN = OPT["min_call_rate_modern"]
MIN_CR_TARGET = OPT["min_call_rate_target"]
MIN_CR_ANCIENT = OPT["min_call_rate_ancient"]
MIN_PROJECTION_SNPS = OPT["min_projection_snps"]
MIN_GROUP_N = OPT["min_group_n"]
CALIB_POPS = OPT["calibration_pops"]
CALIB_CHROMS = OPT["calibration_chroms"]
CALIB_N = OPT["calibration_n"]
CALIB_SEED = OPT["calibration_seed"]
AADR_ANNOTATION = OPT["aadr_annotation"]
LINEAGE_HISTORY_FILE = OPT["lineage_history_file"]
ANALYSIS_ID = str(_g("analysis_id", ""))                          # optional; empty = derived from inputs

def rel(p):
    """Path relative to the work root, for printing."""
    try:
        return str(pathlib.Path(p).resolve().relative_to(P))
    except Exception:
        return str(p)

# raw VCF to normalize: delivered VCF if any, else a VCF called from reads by 00_qc.sh
VCF_RAW = str(VENDOR_VCF) if VENDOR_VCF else str(W / "00_input" / f"{SAMPLE}.call.vcf.gz")
