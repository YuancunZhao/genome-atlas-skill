"""血型标志位点的 no-call 基因型不得伪装成参考纯合——对**生产 23** 断言。

缺陷（血型核查 2026-10-13）：`23_bloodgroups.py` 用 `count("1")` 数 ALT 等位，GT 为 "./."（或
"0/."）时数出 0 个，二倍体分支把它渲染成 `REF/REF`——缺失呼叫被报告成参考纯合。本文件用
真实脚本 + 合成 work 树（pvar/VCF/callable/FASTA 全部自造，位点来自真实面板）复现：
杂合位点照常 G/A；./. 位点必须输出"no call"，而不是 C/C；VCF 无记录且 callable 的位点
按约定输出 FASTA 参考纯合并标注 hom-ref (callable)。

23 需要 pandas/pysam 与 PATH 里的 bcftools、samtools；缺任一整文件跳过（本地跳过、服务器
执行，与 test_30_states_production 同一约定）。
"""
import importlib.util, json, os, pathlib, shutil, subprocess, sys, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "23_bloodgroups.py"
_SKIP = ""
for _mod in ("pandas", "pysam"):
    if importlib.util.find_spec(_mod) is None:
        _SKIP = f"{_mod} not available (the production script needs it)"
        break
if not _SKIP and (shutil.which("bcftools") is None or shutil.which("samtools") is None):
    _SKIP = "bcftools/samtools not on PATH (the production script shells out to them)"


def _write_config(td, kg, fasta):
    cfg = {"sample_id": "TESTSAMPLE", "work_dir": str(pathlib.Path(td) / "work"),
           "kg_pfile": str(kg), "fasta": str(fasta)}
    p = pathlib.Path(td) / "config.yaml"
    p.write_text("".join(f"{k}: {json.dumps(v)}\n" for k, v in cfg.items()), encoding="utf-8")
    return p


def _write_pvar(kg):
    # 位点取自真实 panel/blood_groups.tsv；坐标是合成的（pvar↔VCF↔FASTA 三方一致即可）
    kg.parent.mkdir(parents=True, exist_ok=True)
    kg.with_suffix(".pvar").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "18\t100\trs1058396\tG\tA\n"     # VCF 有 0|1 记录 → G/A
        "7\t200\trs8176058\tC\tT\n"       # VCF 有 ./. 记录 → 必须是 no call，不是 C/C
        "1\t300\trs12075\tG\tA\n",        # VCF 无记录、callable → FASTA 参考纯合 G/G
        encoding="utf-8")


def _write_vcf(out):
    out.parent.mkdir(parents=True, exist_ok=True)
    text = ("##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=7>\n##contig=<ID=18>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTSAMPLE\n"
            "18\t100\t.\tG\tA\t50\tPASS\t.\tGT:DP:AD\t0|1:30:15,15\n"
            "7\t200\t.\tC\tT\t50\tLowQual\t.\tGT:DP:AD\t./.:5:0,0\n")
    r = subprocess.run(["bcftools", "view", "-Oz", "-o", str(out)], input=text,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    subprocess.run(["bcftools", "index", str(out)], check=True, capture_output=True)


def _write_fasta(fasta):
    fasta.parent.mkdir(parents=True, exist_ok=True)
    def _seq(n, at, base):
        s = ["A"] * n
        s[at - 1] = base
        return "".join(s)
    seqs = {"1": _seq(400, 300, "G"), "7": _seq(300, 200, "C"), "18": _seq(200, 100, "G")}
    lines = []
    for name, s in seqs.items():
        lines.append(f">{name}")
        lines += [s[i:i + 60] for i in range(0, len(s), 60)]
    fasta.write_text("\n".join(lines) + "\n", encoding="utf-8")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True, capture_output=True)


@unittest.skipIf(_SKIP, _SKIP)
class TestNoCallGenotypeIsNotReferenceHomozygote(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        td = pathlib.Path(self._td.name)
        self.cfg = _write_config(td, td / "work/data/ref/testkg", td / "work/data/ref/fa/test.fa")
        _write_pvar(td / "work/data/ref/testkg")
        _write_vcf(td / "work/wgs/00_input/target.norm.vcf.gz")
        (td / "work/wgs/00_input/callable.bed").parent.mkdir(parents=True, exist_ok=True)
        (td / "work/wgs/00_input/callable.bed").write_text("1\t0\t1000\n7\t0\t1000\n18\t0\t1000\n",
                                                            encoding="utf-8")
        _write_fasta(td / "work/data/ref/fa/test.fa")
        self.out = td / "work/wgs/16_panels/bloodgroups.tsv"

    def tearDown(self):
        self._td.cleanup()

    def _run_and_rows(self):
        r = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                           env={**os.environ, "WGS_CONFIG": str(self.cfg)}, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(self.out.exists(), "23 did not write bloodgroups.tsv")
        return {l.split("\t")[1]: l.split("\t") for l in
                self.out.read_text(encoding="utf-8").strip().splitlines()[1:]}

    def test_no_call_gt_is_reported_as_no_call(self):
        row = self._run_and_rows()["rs8176058"]
        self.assertIn("no call", row[4],
                      f"a ./. genotype must not masquerade as reference homozygote (got {row[4]!r})")
        self.assertNotEqual("C/C", row[4])

    def test_het_and_callable_absence_paths_are_unchanged(self):
        rows = self._run_and_rows()
        self.assertEqual("G/A", rows["rs1058396"][4], "a 0|1 record still renders as REF/ALT")
        self.assertEqual("G/G", rows["rs12075"][4], "absence + callable still renders the FASTA base")
        self.assertEqual("hom-ref (callable)", rows["rs12075"][5])


if __name__ == "__main__":
    unittest.main()
