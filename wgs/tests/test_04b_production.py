"""AN2 验收（04b）：断言对着**生产 04b 脚本**跑，不对着测试里的复制品。

子进程执行 scripts/04b_ancestry_summary.py：$WGS_CONFIG 指向临时 config.yaml，04_ancestry 与
02_complete 的产物在临时目录合成。覆盖走两条真实路径：预写 qc.<tag>.{ref,target}.smiss（缓存命中，
不调 plink2），或用假 plink2 替身从 --out 生成 smiss（计算路径）。解释器没有 pandas 时跳过。
"""
import importlib.util, json, os, pathlib, subprocess, sys, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "04b_ancestry_summary.py"
if importlib.util.find_spec("pandas") is None:
    _SKIP = "pandas not available (the production script imports it)"
else:
    _SKIP = ""

S_MISS_HDR = "#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS"
KG_HDR = ("#FID\tIID\tSuperPop\tPopulation\tALLELE_CT\tNAMED_ALLELE_DOSAGE_SUM\t"
          "PC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG")
TGT_HDR = "#FID\tIID\tALLELE_CT\tNAMED_ALLELE_DOSAGE_SUM\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG"

# 四个参考个体：CHS.1/CHS.2 合格；CDX.1 缺失 60% 不合格；CDX.2 合格但组内只剩 1 人。
KG_ROWS = [
    ("CHS", "CHS.1", "EAS", "CHS", 0.10, 0.01, 0.02, 0.03),
    ("CHS", "CHS.2", "EAS", "CHS", 0.12, 0.02, 0.01, 0.04),
    ("CDX", "CDX.1", "EAS", "CDX", 0.30, 0.05, 0.06, 0.07),
    ("CDX", "CDX.2", "EAS", "CDX", 0.32, 0.06, 0.05, 0.08),
]
REF_SMISS_GOOD = [
    "CHS\tCHS.1\t0\t100000\t0", "CHS\tCHS.2\t0\t100000\t0",
    "CDX\tCDX.1\t60000\t100000\t0.6", "CDX\tCDX.2\t0\t100000\t0",
]
TGT_SMISS_GOOD = "T\tTESTSAMPLE\t0\t100000\t0"
TGT_SMISS_ALL_MISSING = "T\tTESTSAMPLE\t100000\t100000\t1"


def _write_config(td, **over):
    cfg = {
        "sample_id": "TESTSAMPLE",
        "work_dir": str(pathlib.Path(td) / "work"),
        "ancestry_min_call_rate_modern": 0.95,
        "ancestry_min_call_rate_target": 0.95,
        "ancestry_min_projection_snps": 10000,
        "ancestry_min_group_n": 2,
    }
    cfg.update(over)
    p = pathlib.Path(td) / "config.yaml"
    p.write_text("".join(f"{k}: {v}\n" for k, v in cfg.items()), encoding="utf-8")
    return p


def _write_space(td, base="global", eas_files=False):
    """合成 04_ancestry 的一个参考空间 + 02_complete 目标 pfile 占位。"""
    w = pathlib.Path(td) / "work" / "wgs" / "04_ancestry"
    w.mkdir(parents=True, exist_ok=True)
    kg = w / (f"{base}.proj.sscore" if base != "global" else "kg.proj.sscore")
    tgt = w / (f"{base}.target.proj.sscore" if base != "global" else "target.proj.sscore")
    prune = w / (f"prune.{base}.prune.in" if base != "global" else "prune.prune.in")
    kg.write_text("\n".join([KG_HDR] + [
        f"{fid}\t{iid}\t{sp}\t{pop}\t100000\t100000\t" + "\t".join(f"{v:.4f}" for v in pcs)
        for fid, iid, sp, pop, *pcs in KG_ROWS]) + "\n", encoding="utf-8")
    tgt.write_text(TGT_HDR + "\nTESTSAMPLE\tTESTSAMPLE\t100000\t100000\t0.0000\t0.0000\t0.0000\t0.0000\n",
                   encoding="utf-8")
    prune.write_text("rs1\nrs2\n", encoding="utf-8")
    for e in ("pgen", "psam", "pvar"):
        (w / f"kg.common.{e}").touch()
    c = pathlib.Path(td) / "work" / "wgs" / "02_complete"
    c.mkdir(parents=True, exist_ok=True)
    for e in ("pgen", "psam", "pvar"):
        (c / f"TESTSAMPLE.1kg.{e}").touch()
    return w


def _pre_smiss(w, tag, ref_rows, tgt_row):
    (w / f"qc.{tag}.ref.smiss").write_text("\n".join([S_MISS_HDR] + ref_rows) + "\n", encoding="utf-8")
    (w / f"qc.{tag}.target.smiss").write_text(S_MISS_HDR + "\n" + tgt_row + "\n", encoding="utf-8")


def _run_04b(cfg):
    return subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                          env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=120)


def _summary(w):
    return json.loads((w / "summary.json").read_text(encoding="utf-8"))


@unittest.skipIf(_SKIP, _SKIP)
class TestCoverageIsReal(unittest.TestCase):
    """1000G 没有例外：覆盖与门槛和 09b 同一套（最终位点上的真实 smiss）。"""

    def test_gate_excludes_low_coverage_references(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_space(td)
            _pre_smiss(w, "global", REF_SMISS_GOOD, TGT_SMISS_GOOD)
            r = _run_04b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            s = _summary(w)
            self.assertEqual(s["state"], "ok")
            a = s["analyses"][0]
            self.assertEqual(a["counts"]["excluded"], 1, "CDX.1 覆盖 0.4 必须被排除")
            self.assertEqual(a["target"]["call_rate"], 1.0)
            by = {g["group_id"]: g for g in a["groups"]}
            self.assertEqual(by["CHS"]["n"], 2)
            self.assertEqual(by["CHS"]["rank"], 1)
            self.assertEqual(by["CDX"]["n"], 1, "CDX 剩 1 人：不进默认排名")
            self.assertIsNone(by["CDX"]["rank"])

    def test_missingness_is_computed_over_the_final_sites(self):
        """缓存缺失时 04b 自己跑 plink2（--extract 该空间的 prune set）并把结果接到同一门槛。"""
        with tempfile.TemporaryDirectory() as td:
            stub = pathlib.Path(td) / "fake_plink2"
            stub.write_text("#!/bin/bash\n"
                            "out=\"\"; prev=\"\"\n"
                            "for a in \"$@\"; do [ \"$prev\" = \"--out\" ] && out=\"$a\"; prev=\"$a\"; done\n"
                            "case \"$out\" in\n"
                            "  *.ref) printf '" + S_MISS_HDR + "\\n" + "\\n".join(REF_SMISS_GOOD) + "' > \"$out.smiss\" ;;\n"
                            "  *.target) printf '" + S_MISS_HDR + "\\n" + TGT_SMISS_GOOD + "' > \"$out.smiss\" ;;\n"
                            "esac\n"
                            "exit 0\n", encoding="utf-8")
            stub.chmod(0o755)
            cfg = _write_config(td, plink2=str(stub))
            w = _write_space(td)
            r = _run_04b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue((w / "qc.global.ref.smiss").exists(), "plink2 替身应被调用并落盘 smiss")
            a = _summary(w)["analyses"][0]
            self.assertEqual(a["counts"]["excluded"], 1)


@unittest.skipIf(_SKIP, _SKIP)
class TestTargetGate(unittest.TestCase):
    """目标在最终位点集上覆盖不足：该空间不出排名，summary/manifest 记 unavailable 与原因。"""

    def test_all_missing_target_blocks_the_space(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_space(td)
            _pre_smiss(w, "global", REF_SMISS_GOOD, TGT_SMISS_ALL_MISSING)
            r = _run_04b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)      # 04b 仍退出 0：状态进 manifest
            s = _summary(w)
            self.assertEqual(s["state"], "unavailable")
            self.assertEqual(s["reason_code"], "target_below_coverage_gate")
            self.assertEqual(s["analyses"], [])
            man = json.loads((w / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(man["state"], "unavailable")


@unittest.skipIf(_SKIP, _SKIP)
class TestLegacyEasIdentity(unittest.TestCase):
    """旧 eas.* 文件只在配置的区域就是东亚时才可当区域结果；别的区域配置宁可缺不可错。"""

    def test_non_eas_scope_rejects_legacy_files(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td, ref_superpop="SAS")
            w = _write_space(td)                       # global 完整
            _write_space(td, base="eas")               # 只有旧 eas.*，没有 regional.*
            _pre_smiss(w, "global", REF_SMISS_GOOD, TGT_SMISS_GOOD)
            _pre_smiss(w, "eas", REF_SMISS_GOOD, TGT_SMISS_GOOD)
            r = _run_04b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            s = _summary(w)
            self.assertEqual([a["analysis_id"] for a in s["analyses"]], ["kg-global"])
            self.assertIn("legacy eas", (w / "summary.txt").read_text(encoding="utf-8"))

    def test_eas_scope_accepts_legacy_files(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td, ref_superpop="EAS")
            w = _write_space(td)
            _write_space(td, base="eas")
            _pre_smiss(w, "global", REF_SMISS_GOOD, TGT_SMISS_GOOD)
            _pre_smiss(w, "eas", REF_SMISS_GOOD, TGT_SMISS_GOOD)
            r = _run_04b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            s = _summary(w)
            self.assertEqual([a["analysis_id"] for a in s["analyses"]], ["kg-global", "kg-eas"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
