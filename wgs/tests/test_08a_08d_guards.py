"""旧缺陷 #17（H5 残余）回归：对着**生产 08a/08d 脚本**子进程断言，不对着复制品。

08a 曾把每次产出永久标记 verified_against_delivered:false，08d 见标记即 exit 2——把“本次
运行做了什么”与“是否复现某次历史交付”耦在一起：全新目录的正常流程也必断（复审 H5 点名）。
而 2026-09-28 复测已证明该链完全可复现（delly BCF 42,909/42,909 条一致，md5 仅差 ##fileDate，
见 work/wgs/08_sv/delly_retest_20260928/RETEST_NOTES.md）。修复后：provenance 只记录本次运行
事实；交付原件靠输出侧覆盖保护（无 provenance 标记的现存产物不默默重算）。

替身：samtools view/index、bcftools index 走 PATH 桩；delly 桩放在 tools_dir（08d 按配置找
$TOOLS/delly）。config/wgsconfig/env.sh 全是真实生产路径；依赖仅 stdlib，本地与服务器都实跑。
"""
import json, os, pathlib, subprocess, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
A08A = REPO / "scripts" / "08a_main_contigs.sh"
A08D = REPO / "scripts" / "08d_delly.sh"


def _fixture(td, delly_stub=True):
    """临时 work 根：config + 假参考(.fai 25 contigs) + 假 CRAM + PATH 桩。返回 (cfg, bin_dir)。"""
    root = pathlib.Path(td)
    ref = root / "ref"; raw = root / "raw"; tools = root / "tools"; bin_ = root / "bin"
    for d in (ref, raw, tools, bin_):
        d.mkdir(parents=True, exist_ok=True)
    contigs = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    (ref / "fake.fasta").write_text(">1\nACGT\n", encoding="utf-8")
    (ref / "fake.fasta.fai").write_text(
        "".join(f"{c}\t1000\t6\t60\t61\n" for c in contigs), encoding="utf-8")
    (raw / "in.cram").write_bytes(b"FAKE-INPUT-CRAM")
    cfg = root / "config.yaml"
    cfg.write_text(
        f"sample_id: GUARDSAMPLE\n"
        f"work_dir: {root / 'work'}\n"
        f"fasta: {ref / 'fake.fasta'}\n"
        f"reads: {raw / 'in.cram'}\n"
        f"tools_dir: {tools}\n"
        f"threads: 2\n", encoding="utf-8")
    (bin_ / "samtools").write_text(
        "#!/bin/bash\n"
        'cmd=$1; shift\n'
        'prev=""; out=""\n'
        'for a in "$@"; do [ "$prev" = "-o" ] && out="$a"; prev="$a"; done\n'
        'if [ "$cmd" = view ]; then printf "FAKE-CRAM" > "$out"; exit 0; fi\n'
        'if [ "$cmd" = index ]; then f=$1; [ -n "$out" ] && f="$out"; touch "$f.crai" 2>/dev/null; exit 0; fi\n'
        "exit 0\n", encoding="utf-8")
    (bin_ / "bcftools").write_text("#!/bin/bash\nexit 0\n", encoding="utf-8")
    if delly_stub:
        (tools / "delly").write_text(
            "#!/bin/bash\n"
            'prev=""; out="target.sv.bcf"\n'
            'for a in "$@"; do [ "$prev" = "-o" ] && out="$a"; prev="$a"; done\n'
            'printf "FAKE-BCF" > "$out"\n'
            "exit 0\n", encoding="utf-8")
        (tools / "delly").chmod(0o755)
    for f in bin_.iterdir():
        f.chmod(0o755)
    return cfg, bin_


def _run(script, cfg, bin_dir, extra_env=None):
    env = {**os.environ, "WGS_CONFIG": str(cfg),
           "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"}
    env.update(extra_env or {})
    return subprocess.run(["bash", str(script)], capture_output=True, text=True,
                          env=env, timeout=120)


def _w(td):  # work/wgs root
    return pathlib.Path(td) / "work" / "wgs"


class TestFreshDirectoryChain(unittest.TestCase):
    """核心回归：全新目录 08a→08d 一路跑通；provenance 是运行事实，不再携带交付比对标记。"""

    def test_chain_runs_and_provenance_records_run_facts(self):
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir = _fixture(td)
            r = _run(A08A, cfg, bin_dir)
            self.assertEqual(r.returncode, 0, r.stderr)
            mc = _w(td) / "08_sv/main_contigs/norm.main.cram"
            self.assertEqual(mc.read_bytes(), b"FAKE-CRAM")
            prov = json.loads((mc.parent / "norm.main.cram.provenance.json")
                              .read_text(encoding="utf-8"))
            self.assertEqual(prov["producer"], "08a_main_contigs.sh")
            self.assertIn("samtools view", prov["command"])
            self.assertGreaterEqual(prov["contigs"], 24)
            # 旧缺陷本体：此键把新样本流程与历史交付复现耦合，08d 见 false 即断链
            self.assertNotIn("verified_against_delivered", prov,
                             "provenance 记录本次运行事实，不携带交付比对断言")
            r = _run(A08D, cfg, bin_dir)
            self.assertEqual(r.returncode, 0, r.stderr)
            bcf = _w(td) / "08_sv/delly/target.sv.bcf"
            self.assertEqual(bcf.read_bytes(), b"FAKE-BCF")
            dprov = json.loads((bcf.parent / "target.sv.bcf.provenance.json")
                               .read_text(encoding="utf-8"))
            self.assertEqual(dprov["producer"], "08d_delly.sh")
            self.assertIn("delly sr", dprov["command"])


class TestDeliveredArtifactProtection(unittest.TestCase):
    """守卫仍在，但改为输出侧：无 provenance 标记的现存产物（交付原件的指纹）不默默重算。"""

    def test_08a_keeps_unmarked_existing_cram(self):
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir = _fixture(td)
            mc_dir = _w(td) / "08_sv/main_contigs"; mc_dir.mkdir(parents=True)
            (mc_dir / "norm.main.cram").write_bytes(b"DELIVERED-1.24G")
            r = _run(A08A, cfg, bin_dir)
            self.assertEqual(r.returncode, 0, "保护性跳过不是失败")
            self.assertIn("refusing to overwrite", r.stderr)
            self.assertEqual((mc_dir / "norm.main.cram").read_bytes(), b"DELIVERED-1.24G")
            self.assertFalse((mc_dir / "norm.main.cram.provenance.json").exists())

    def test_08d_keeps_unmarked_existing_bcf(self):
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir = _fixture(td)
            self.assertEqual(_run(A08A, cfg, bin_dir).returncode, 0)
            d_dir = _w(td) / "08_sv/delly"; d_dir.mkdir(parents=True, exist_ok=True)
            (d_dir / "target.sv.bcf").write_bytes(b"DELIVERED-BCF")
            r = _run(A08D, cfg, bin_dir)
            self.assertEqual(r.returncode, 0, "保护性跳过不是失败")
            self.assertIn("refusing to overwrite", r.stderr)
            self.assertEqual((d_dir / "target.sv.bcf").read_bytes(), b"DELIVERED-BCF")
            self.assertFalse((d_dir / "logs/delly.log").exists(), "delly 不应被调用")

    def test_08d_overwrites_when_provenance_marker_present(self):
        """带 provenance 的现存 bcf 是本管线自己产的：正常重算覆盖。"""
        with tempfile.TemporaryDirectory() as td:
            cfg, bin_dir = _fixture(td)
            self.assertEqual(_run(A08A, cfg, bin_dir).returncode, 0)
            d_dir = _w(td) / "08_sv/delly"; d_dir.mkdir(parents=True, exist_ok=True)
            (d_dir / "target.sv.bcf").write_bytes(b"OLD-RUN-BCF")
            (d_dir / "target.sv.bcf.provenance.json").write_text(
                json.dumps({"producer": "08d_delly.sh", "date": "2026-01-01T00:00:00Z"}),
                encoding="utf-8")
            r = _run(A08D, cfg, bin_dir)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual((d_dir / "target.sv.bcf").read_bytes(), b"FAKE-BCF")


if __name__ == "__main__":
    unittest.main()
