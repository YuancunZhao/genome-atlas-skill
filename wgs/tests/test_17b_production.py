"""AN3 验收（17b）：断言对着**生产 17b 脚本**跑。

子进程执行 scripts/17b_local_ancestry_calibrated.py：$WGS_CONFIG 指向临时 config.yaml，
12_localanc 的 FLARE 产物在临时目录按真实格式合成（la.<c>.global.anc.gz 为 SAMPLE+各面板列的
TSV，.model 带 "# list of ancestries"，psam 提供 holdout 的群体归属）。无 fai 时长度回退内置
hg19 表（GRCh37 任务合法回退）。解释器没有 pandas 时跳过。
"""
import gzip, importlib.util, json, os, pathlib, subprocess, sys, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "17b_local_ancestry_calibrated.py"
if importlib.util.find_spec("pandas") is None:
    _SKIP = "pandas not available (the production script imports it)"
else:
    _SKIP = ""

# 面板列顺序故意不按"北在前"：European/SouthEA/SouthAsian/NorthEA 与真实产物一致，
# 标签换序测试据此校验 17b 是否跟着**配置**的 A/B 走而不是跟着列顺序走。
PANELS = ["European", "SouthEA", "SouthAsian", "NorthEA"]
A, B = "NorthEA", "SouthEA"          # local_ancestry_labels 配置的两个来源面板


def _anc_row(sample, north, south):
    return {"SAMPLE": sample, "European": 0.005, "SouthEA": south, "SouthAsian": 0.02, "NorthEA": north}


def _write_config(td, **over):
    cfg = {
        "sample_id": "TESTSAMPLE",
        "work_dir": str(pathlib.Path(td) / "work"),
        "kg_pfile": str(pathlib.Path(td) / "work/data/ref/all_phase3"),
        "local_ancestry_a": "[CHB]", "local_ancestry_b": "[CDX]",
        "local_ancestry_labels": "[NorthEA, SouthEA]",
    }
    cfg.update(over)
    p = pathlib.Path(td) / "config.yaml"
    p.write_text("".join(f"{k}: {v}\n" for k, v in cfg.items()), encoding="utf-8")
    return p


def _write_anc(w, prefix, chrom, rows):
    cols = "\t".join(["SAMPLE"] + PANELS)
    body = "\n".join("\t".join(str(r[c]) for c in ["SAMPLE"] + PANELS) for r in rows)
    with gzip.open(w / f"{prefix}.{chrom}.global.anc.gz", "wt", encoding="utf-8") as fh:
        fh.write(cols + "\n" + body + "\n")
    (w / f"{prefix}.{chrom}.model").write_text(
        "# list of ancestries\n" + "\t".join(PANELS) + "\n\n# list of reference panels\n"
        + "\t".join(PANELS) + "\n", encoding="utf-8")


def _write_la(td, target_rows_by_chrom):
    """目标 la.* 产物 + 17 预置的 local_ancestry.json + 空 psam（无校准时 17b 也会读它）。"""
    w = pathlib.Path(td) / "work" / "wgs" / "12_localanc"
    w.mkdir(parents=True, exist_ok=True)
    for c, rows in target_rows_by_chrom.items():
        _write_anc(w, "la", c, rows)
    (w / "local_ancestry.json").write_text(json.dumps({
        "state": "ok", "reason_code": "", "panels": [], "global": [], "per_chrom": [],
        "calibration": [], "missing_chroms": [], "aggregation_method": "marker_counts_for_global",
    }), encoding="utf-8")
    ref = pathlib.Path(td) / "work" / "data" / "ref"
    ref.mkdir(parents=True, exist_ok=True)
    (ref / "all_phase3.psam").write_text("#IID\tPopulation\nNONE\tNONE\n", encoding="utf-8")
    return w


def _write_calib(td, rows_by_chrom, psam_lines):
    """校准产物 + 真实 psam（holdout 的群体归属来自这里）。"""
    w = pathlib.Path(td) / "work" / "wgs" / "12_localanc"
    for c, rows in rows_by_chrom.items():
        _write_anc(w, "calib", c, rows)
    ref = pathlib.Path(td) / "work" / "data" / "ref"
    (ref / "all_phase3.psam").write_text("#IID\tPopulation\n" + psam_lines, encoding="utf-8")
    return w


def _run_17b(cfg):
    return subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                          env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=120)


def _json(w):
    return json.loads((w / "local_ancestry.json").read_text(encoding="utf-8"))


@unittest.skipIf(_SKIP, _SKIP)
class TestRawDeliverableWithoutCalibration(unittest.TestCase):
    """无校准产物：原始 LA 独立交付，校准状态明确 not_run（而不是留一个分不清含义的空列表）。"""

    def test_no_calib_still_delivers_and_marks_not_run(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_la(td, {"1": [_anc_row("TESTSAMPLE", 0.60, 0.35)],
                                "2": [_anc_row("TESTSAMPLE", 0.58, 0.37)]})
            r = _run_17b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue((w / "dayu_global.tsv").exists(), "原始 LA 的全局 posterior 照常写盘")
            d = _json(w)
            self.assertEqual(d["calibration"], [])
            self.assertEqual(d["calibration_state"], "not_run")
            self.assertIn("16b", d["calibration_reason"])


@unittest.skipIf(_SKIP, _SKIP)
class TestCalibrationStatistics(unittest.TestCase):
    """列序不是标签、SD=0 不算 z、口径写明并回写同集合目标值（HANDOFF AN3 的四条验收）。"""

    def test_labels_follow_the_config_not_the_column_order(self):
        """面板列序 European 在前而配置的 A=NorthEA 在末位：统计必须取 NorthEA。若按 anc[0] 取，
        north_mean 会变成 European 的 0.005——数值看着也合理，含义完全错。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_la(td, {"1": [_anc_row("TESTSAMPLE", 0.60, 0.35)],
                            "2": [_anc_row("TESTSAMPLE", 0.58, 0.37)]})
            w = _write_calib(td,
                {"1": [_anc_row("N1", 0.95, 0.03), _anc_row("N2", 0.93, 0.05)],
                 "2": [_anc_row("N1", 0.95, 0.03), _anc_row("N2", 0.93, 0.05)]},
                "N1\tP_NORTH\nN2\tP_NORTH\n")
            r = _run_17b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            cal = {c["population"]: c for c in _json(w)["calibration"]}
            self.assertAlmostEqual(cal["P_NORTH"]["north_mean"], 0.94, places=6,
                                   msg="north_mean 必须是 NorthEA 列的均值，不是列序第一的 European")
            self.assertEqual(_json(w)["calibration_panels"], ["NorthEA", "SouthEA"])

    def test_zero_and_undefined_sd_give_no_zscore(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_la(td, {"1": [_anc_row("TESTSAMPLE", 0.60, 0.35)],
                            "2": [_anc_row("TESTSAMPLE", 0.58, 0.37)]})
            w = _write_calib(td,
                {"1": [_anc_row("G1", 0.90, 0.06), _anc_row("G2", 0.80, 0.16),   # 正常群体
                       _anc_row("F1", 0.50, 0.45), _anc_row("F2", 0.50, 0.45),    # 完全同值 → sd=0
                       _anc_row("S1", 0.70, 0.25)],                                # 单个体 → sd=NaN
                 "2": [_anc_row("G1", 0.92, 0.04), _anc_row("G2", 0.78, 0.18),
                       _anc_row("F1", 0.50, 0.45), _anc_row("F2", 0.50, 0.45),
                       _anc_row("S1", 0.70, 0.25)]},
                "G1\tP_GOOD\nG2\tP_GOOD\nF1\tP_FLAT\nF2\tP_FLAT\nS1\tP_SINGLE\n")
            r = _run_17b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("sd on NorthEA", r.stdout, "正常群体照常给 z 分数")
            self.assertEqual(r.stdout.count("no z-score"), 2, "sd=0 与单个体都不给 z 分数")
            cal = {c["population"]: c for c in _json(w)["calibration"]}
            self.assertEqual(cal["P_FLAT"]["north_sd"], 0.0)
            self.assertIsNone(cal["P_SINGLE"]["north_sd"], "单个体的 sd 未定义 → null，不是 NaN 字面量")

    def test_calibration_uses_chromosomes_common_to_everyone(self):
        """缺染色体时统一到交集：C2 缺 chr1 → 比较只在 chr2 上做，目标值也取 chr2（同分母）。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_la(td, {"1": [_anc_row("TESTSAMPLE", 0.90, 0.08)],   # 目标 chr1 偏北
                            "2": [_anc_row("TESTSAMPLE", 0.60, 0.35)]})  # chr2 居中
            w = _write_calib(td,
                {"1": [_anc_row("C1", 0.95, 0.03)],
                 "2": [_anc_row("C1", 0.95, 0.03), _anc_row("C2", 0.93, 0.05)]},
                "C1\tP_A\nC2\tP_A\n")
            r = _run_17b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _json(w)
            self.assertEqual(d["calibration_chroms"], ["2"], "C2 缺 chr1，比较只落在共同的 chr2")
            self.assertAlmostEqual(d["calibration"][0]["target_north"], 0.60, places=6,
                                   msg="目标值也限制在同一染色体集合上")

    def test_empty_intersection_records_the_state(self):
        """两个 holdout 各缺对方那条染色体 → 交集为空：明确记 no_common_chroms，不做比较。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_la(td, {"1": [_anc_row("TESTSAMPLE", 0.60, 0.35)],
                            "2": [_anc_row("TESTSAMPLE", 0.58, 0.37)]})
            w = _write_calib(td,
                {"1": [_anc_row("X1", 0.90, 0.06)],
                 "2": [_anc_row("X2", 0.80, 0.16)]},
                "X1\tP_A\nX2\tP_B\n")
            r = _run_17b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _json(w)
            self.assertEqual(d["calibration_state"], "no_common_chroms")
            self.assertEqual(d["calibration"], [])
            self.assertIn("no chromosome is shared", r.stderr)

    def test_calibrated_global_declares_its_caliber(self):
        """17 的 global 是 marker 计数，17b 的全局是 posterior 长度加权：口径必须写进 JSON。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_la(td, {"1": [_anc_row("TESTSAMPLE", 0.60, 0.35)],
                            "2": [_anc_row("TESTSAMPLE", 0.58, 0.37)]})
            w = _write_calib(td,
                {"1": [_anc_row("N1", 0.95, 0.03)],
                 "2": [_anc_row("N1", 0.95, 0.03)]},
                "N1\tP_NORTH\n")
            r = _run_17b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _json(w)
            self.assertEqual(d["calibrated_global_caliber"], "flare_posterior_length_weighted")
            g = {row["panel_id"]: row["value"] for row in d["calibrated_global"]}
            self.assertAlmostEqual(g["NorthEA"], 0.590113, places=4,
                                   msg="长度加权（chr1 249.25 Mb / chr2 243.20 Mb）的全局 posterior")


if __name__ == "__main__":
    unittest.main(verbosity=2)
