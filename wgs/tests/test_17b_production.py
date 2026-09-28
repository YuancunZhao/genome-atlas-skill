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


if __name__ == "__main__":
    unittest.main(verbosity=2)
