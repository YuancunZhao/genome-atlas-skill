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
        "plink2": str(pathlib.Path(td) / "fake_plink2"),   # _write_space 落盘的替身
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
    """合成 04_ancestry 的一个参考空间 + 02_complete 目标 pfile 占位 + 假 plink2 替身。

    04b 现在**每次重算** qc.<tag>.{ref,target}.smiss（复审 AN2/H6：存在性缓存挡不住换目标/换
    prune 集），所以缺失率经替身从 --out 落盘：替身复制 <td>/stub.{ref,target}.smiss（_pre_smiss
    写入本次想要的行）并 touch plink_called 留证。"""
    td = pathlib.Path(td)
    w = td / "work" / "wgs" / "04_ancestry"
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
    c = td / "work" / "wgs" / "02_complete"
    c.mkdir(parents=True, exist_ok=True)
    for e in ("pgen", "psam", "pvar"):
        (c / f"TESTSAMPLE.1kg.{e}").touch()
    (td / "stub.ref.smiss").write_text("\n".join([S_MISS_HDR] + REF_SMISS_GOOD) + "\n", encoding="utf-8")
    (td / "stub.target.smiss").write_text(S_MISS_HDR + "\n" + TGT_SMISS_GOOD + "\n", encoding="utf-8")
    stub = td / "fake_plink2"
    stub.write_text("#!/bin/bash\n"
                    "out=\"\"; prev=\"\"\n"
                    "for a in \"$@\"; do [ \"$prev\" = \"--out\" ] && out=\"$a\"; prev=\"$a\"; done\n"
                    "case \"$out\" in\n"
                    f"  *.ref) cp \"{td}/stub.ref.smiss\" \"$out.smiss\" ;;\n"
                    f"  *.target) cp \"{td}/stub.target.smiss\" \"$out.smiss\" ;;\n"
                    "esac\n"
                    f"touch \"{td}/plink_called\"\n"
                    "exit 0\n", encoding="utf-8")
    stub.chmod(0o755)
    return w


def _pre_smiss(w, tag, ref_rows, tgt_row):
    """每个测试想要的缺失率行：写给替身的 sidecar（04b 每次重算），同时落一份 qc.<tag>.*——
    那份旧文件正是"重算必须覆盖存在性缓存"回归里的 stale 候选。"""
    td = w.parents[2]
    (td / "stub.ref.smiss").write_text("\n".join([S_MISS_HDR] + ref_rows) + "\n", encoding="utf-8")
    (td / "stub.target.smiss").write_text(S_MISS_HDR + "\n" + tgt_row + "\n", encoding="utf-8")
    (w / f"qc.{tag}.ref.smiss").write_text("\n".join([S_MISS_HDR] + ref_rows) + "\n", encoding="utf-8")
    (w / f"qc.{tag}.target.smiss").write_text(S_MISS_HDR + "\n" + tgt_row + "\n", encoding="utf-8")


def _run_04b(cfg):
    return subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                          env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=120)


def _summary(w):
    return json.loads((w / "summary.json").read_text(encoding="utf-8"))


@unittest.skipIf(_SKIP, _SKIP)
class TestCoverageIsReal(unittest.TestCase):
    def test_stale_qc_smiss_is_overwritten_not_reused(self):
        """复审 AN2/H6：qc.<tag>.*.smiss 不再存在性缓存。预置一份旧文件（上一轮 prune 集算的，
        CDX.1 数字全错）后运行：plink2 必须真的被调用，旧文件被本次结果覆盖，门槛用新数字。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_space(td)
            stale = [S_MISS_HDR,
                     "CHS\tCHS.1\t100000\t100000\t1", "CHS\tCHS.2\t100000\t100000\t1",
                     "CDX\tCDX.1\t0\t100000\t0", "CDX\tCDX.2\t0\t100000\t0"]
            (w / "qc.global.ref.smiss").write_text("\n".join(stale) + "\n", encoding="utf-8")
            r = _run_04b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue((pathlib.Path(td) / "plink_called").exists(),
                            "qc smiss 已存在也必须重算（plink2 被真正调用）")
            fresh = (w / "qc.global.ref.smiss").read_text(encoding="utf-8")
            self.assertNotIn("CHS.1\t100000", fresh, "旧 stale smiss 必须被覆盖")
            a = _summary(w)["analyses"][0]
            by = {g["group_id"]: g for g in a["groups"]}
            self.assertEqual(by["CHS"]["n"], 2, "门槛读的是本次重算的数字（CHS 合格），不是旧缓存（全缺失）")
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


@unittest.skipIf(_SKIP, _SKIP)
class TestFailureInvalidatesOldOk(unittest.TestCase):
    """复审 §3.2 P0（H6/AN0/AN5 失败生命周期）：成功→同 ID 改输入→中途崩溃→单独准入。

    04b 此前只在成功结尾写 manifest；计算中途崩溃时旧 manifest=ok 与旧 summary.json 原样
    保留，30 的 analysis_state() 照旧准入。现在 04b 开工先写失效记录。"""

    def test_midrun_crash_invalidates_previous_ok(self):
        sys.path.insert(0, str(REPO / "scripts"))
        import ancestry_data as ad
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_space(td)
            r = _run_04b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            state, _, doc = ad.analysis_state(
                w, {"sample_id": "TESTSAMPLE"},
                owning_steps=("04-ancestry-pca", "04b-ancestry-summary"))
            self.assertEqual(state, "ok")
            self.assertIsNotNone(doc)
            # 同 ID 换输入：kg.proj.sscore 换成读不了的东西（目录），复现计算中途的未处理崩溃
            (w / "kg.proj.sscore").unlink()
            (w / "kg.proj.sscore").mkdir()
            r2 = _run_04b(cfg)
            self.assertNotEqual(r2.returncode, 0, "输入损坏必须以非零退出终止，不得静默产出")
            # 单独准入（30 的调用形状，global-only：04c 禁用记录同场也不得干扰这条结论）
            ad.write_manifest(w / "manifest.04c-per-chromosome-axis.json",
                              ad.disabled_manifest("TESTSAMPLE", "04c-per-chromosome-axis",
                                                   "regional_axis_not_configured"))
            state, reason, doc = ad.analysis_state(
                w, {"sample_id": "TESTSAMPLE"},
                owning_steps=("04-ancestry-pca", "04b-ancestry-summary"))
            self.assertEqual((state, reason), ("unavailable", "run_begun_not_published"))
            self.assertIsNone(doc, "开工失效后，上一轮的 summary.json 不得继续交付")


if __name__ == "__main__":
    unittest.main(verbosity=2)
