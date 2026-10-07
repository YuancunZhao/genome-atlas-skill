"""旧缺陷 #8c 回归：对着**生产 01_normalize.sh** 子进程断言，不对着复制品。

旧实现把主 contig 写死成裸名 `1..22,X,Y,MT`（bcftools view -r 在 chr 前缀输入上直接选错
contig），callable bed 的 decoy 过滤是负向 `$1!~/^GL/`（chr 前缀输入下 chrUn_*/random
scaffold 全部漏进 callable.bed），且 bed 为空时只打印 0 继续跑。本测试用真实脚本 + 合成
VCF 头/mosdepth bed + bcftools/bedtools/tabix/zcat 桩验证三件事：主 contig 从 VCF 头推导
（chr 前缀与裸名都接受、区域串沿用文件自身拼法）、bed 正向匹配（chrM 与 MT 视为同一
canonical 名，混合命名也能对上）、无主 contig / 空 callable bed 时响亮失败。

依赖仅 stdlib（wgsconfig/env.sh 无 yaml 也能解析本夹具的标量 config），本地与服务器都实跑。
"""
import os, pathlib, subprocess, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "01_normalize.sh"

CHR_HEADER = "".join(
    f"##contig=<ID={c},length=100000>\n"
    for c in ("chr1", "chr2", "chrX", "chrY", "chrM", "GL000250.2", "chrUn_gl000220")
) + "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1000\t.\tA\tC\t.\tPASS\t.\n"

BARE_HEADER = "".join(
    f"##contig=<ID={c},length=100000>\n"
    for c in ("1", "2", "X", "Y", "MT", "GL000250.2", "hs37d5")
) + "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1000\t.\tA\tC\t.\tPASS\t.\n"

SCAFFOLD_ONLY_HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=scaffold_1,length=5000>\n"
    "##contig=<ID=chrUn_gl000220,length=5000>\n##contig=<ID=GL000250.2,length=5000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nscaffold_1\t100\t.\tA\tC\t.\tPASS\t.\n")

# chr 前缀床：主 contig（chr1/chrX 保留）、深度不够（2:8 丢）、非主 contig（chrUn/GL 丢）、
# 裸名 MT 行——VCF 头里是 chrM，canonical 化后必须仍能对上（混合命名兼容）。
CHR_BED = ("chr1\t100\t200\t8:20\n"
           "chr1\t500\t600\t2:8\n"
           "chrX\t1000\t1500\t60:inf\n"
           "chrUn_gl000220\t10\t20\t8:20\n"
           "GL000250.2\t10\t20\t8:20\n"
           "MT\t30\t40\t8:20\n")

GL_ONLY_BED = "GL000250.2\t10\t20\t8:20\nchrUn_gl000220\t30\t40\t60:inf\n"


def _fixture(td, header, bed):
    root = pathlib.Path(td)
    ref, raw, tools, bin_ = root / "ref", root / "raw", root / "tools", root / "bin"
    for d in (ref, raw, tools, bin_):
        d.mkdir(parents=True, exist_ok=True)
    vcf = raw / "in.vcf"
    vcf.write_text(header, encoding="utf-8")
    (ref / "fake.fasta").write_text(">1\nACGT\n", encoding="utf-8")
    cfg = root / "config.yaml"
    cfg.write_text(
        f"sample_id: NORMSAMPLE\n"
        f"work_dir: {root / 'work'}\n"
        f"fasta: {ref / 'fake.fasta'}\n"
        f"vendor_vcf: {vcf}\n"
        f"tools_dir: {tools}\n"
        f"threads: 1\n", encoding="utf-8")
    qc = root / "work/wgs/01_qc"          # 真实流程由 00_qc 创建；测试里直接预置
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "depth.quantized.bed.gz").write_text(bed, encoding="utf-8")  # zcat 桩只 cat，不真解压
    (bin_ / "zcat").write_text('#!/bin/bash\nexec cat "$@"\n', encoding="utf-8")
    (bin_ / "bedtools").write_text("#!/bin/bash\nexec cat\n", encoding="utf-8")  # merge 桩：直通
    (bin_ / "tabix").write_text("#!/bin/bash\nexit 0\n", encoding="utf-8")
    args_log = root / "bcftools_args.log"
    # bcftools 桩：view -h 吐头行、view -r 记录区域串并吐数据行、view -f PASS 拷贝、norm 收
    # stdin 写 -o、stats 出一行。值取型旗标（-r/-o/-f/-F/--threads）的下一个参数不当定位参。
    (bin_ / "bcftools").write_text(
        "#!/bin/bash\n"
        'log="${BCFTOOLS_ARGS_LOG:?}"\n'
        'printf "%s\\n" "$*" >> "$log"\n'
        'cmd=$1; shift\n'
        'prev=""; out=""; regions=""; hflag=0; pos=()\n'
        'for a in "$@"; do\n'
        '  if [ -n "$prev" ]; then\n'
        '    case "$prev" in -r) regions="$a";; -o) out="$a";; -f|-F|--threads) :;; esac\n'
        '    prev=""; continue\n'
        '  fi\n'
        '  case "$a" in\n'
        '    -h) hflag=1;;\n'
        '    -r|-o|-f|-F|--threads) prev="$a";;\n'
        '    -*) :;;\n'
        '    *) pos+=("$a");;\n'
        '  esac\n'
        'done\n'
        'in=""; for p in "${pos[@]}"; do in="$p"; done\n'
        'if [ "$cmd" = view ] && [ "$hflag" = 1 ]; then grep "^#" "$in"; exit 0; fi\n'
        'if [ "$cmd" = view ] && [ -n "$regions" ]; then\n'
        '  printf "REGIONS\\t%s\\n" "$regions" >> "$log"\n'
        '  grep -v "^#" "$in" || true; exit 0\n'
        "fi\n"
        'if [ "$cmd" = view ]; then cat "$in" > "$out"; exit 0; fi\n'
        'if [ "$cmd" = norm ]; then cat > "$out"; exit 0; fi\n'
        'if [ "$cmd" = stats ]; then echo "STATS $in"; exit 0; fi\n'
        "exit 0\n", encoding="utf-8")
    for f in bin_.iterdir():
        f.chmod(0o755)
    return cfg, bin_, args_log, root / "work/wgs/00_input"


def _run(cfg, bin_dir, args_log):
    return subprocess.run(
        ["bash", str(SCRIPT)], capture_output=True, text=True,
        env={**os.environ, "WGS_CONFIG": str(cfg), "BCFTOOLS_ARGS_LOG": str(args_log),
             "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"},
        timeout=120)


def _regions(args_log):
    for line in args_log.read_text(encoding="utf-8").splitlines():
        if line.startswith("REGIONS\t"):
            return set(line.split("\t", 1)[1].split(","))
    return None


class TestMainContigsDerivedFromHeader(unittest.TestCase):

    def test_chr_prefixed_input_uses_the_files_own_names(self):
        """核心回归：chr 前缀输入时区域串必须来自头（chr1,chr2,chrM,chrX,chrY），不是写死的
        裸名全表——旧代码在这里把 22 条裸名一股脑塞给 -r。"""
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir, args_log, out_dir = _fixture(td, CHR_HEADER, CHR_BED)
            r = _run(cfg, bin_dir, args_log)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("NORMALIZE_DONE", r.stdout)
            self.assertEqual(_regions(args_log), {"chr1", "chr2", "chrM", "chrX", "chrY"})
            self.assertTrue((out_dir / "target.norm.vcf.gz").exists())

    def test_bare_input_still_selects_only_contigs_the_header_declares(self):
        """裸名输入照常工作，且不再把头里不存在的 3..22 塞进区域串。"""
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir, args_log, _ = _fixture(td, BARE_HEADER, CHR_BED.replace("chr", ""))
            r = _run(cfg, bin_dir, args_log)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(_regions(args_log), {"1", "2", "X", "Y", "MT"})

    def test_no_main_contigs_fails_loudly_before_any_selection(self):
        """只有 scaffold/decoy 的头：必须在选区之前失败，不能静默产出空集后续崩。"""
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir, args_log, _ = _fixture(td, SCAFFOLD_ONLY_HEADER, CHR_BED)
            r = _run(cfg, bin_dir, args_log)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("no main contigs", r.stderr)
            self.assertIsNone(_regions(args_log), "失败必须发生在 bcftools view -r 之前")


class TestCallableBedPositiveFilter(unittest.TestCase):

    def test_keeps_main_contigs_only_and_matches_across_naming(self):
        """正向过滤：丢 chrUn/GL 与深度不够的行；VCF 头写 chrM 而床写 MT 也能对上。旧代码的
        负向 /^GL/ 过滤会把 chrUn_gl000220 行留在 callable.bed 里。"""
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir, args_log, out_dir = _fixture(td, CHR_HEADER, CHR_BED)
            r = _run(cfg, bin_dir, args_log)
            self.assertEqual(r.returncode, 0, r.stderr)
            bed = (out_dir / "callable.bed").read_text(encoding="utf-8")
            self.assertIn("chr1\t100\t200", bed)
            self.assertIn("chrX\t1000\t1500", bed)
            self.assertIn("MT\t30\t40", bed, "床里的裸名 MT 必须与头里的 chrM canonical 对上")
            self.assertNotIn("chrUn", bed)
            self.assertNotIn("GL000250", bed)
            self.assertNotIn("500\t600", bed, "2:8 深度不够的行必须丢")

    def test_empty_callable_bed_fails_loudly(self):
        """主 contig 上没有任何过深度/MQ 门槛的区域：旧代码打印 0 继续跑，必须改为失败。"""
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir, args_log, _ = _fixture(td, CHR_HEADER, GL_ONLY_BED)
            r = _run(cfg, bin_dir, args_log)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("callable.bed is empty", r.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
