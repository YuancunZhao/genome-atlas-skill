"""AN2 验收：断言对着**生产 09b 脚本**跑，不对着测试里的复制品。

子进程执行 scripts/09b_aadr_summary.py：$WGS_CONFIG 指向临时 config.yaml，work 根与 11_aadr
产物全部在临时目录合成，缺失率文件预写为 aadr.pruned.smiss（plink2 因此不会被调用）。HANDOFF
要求"验收必须对生产09b断言"——辅助函数绿灯不能替代这里。生产脚本 import pandas，当前解释器没有
pandas 时整文件跳过（本地套件跳过、服务器套件执行，与 test_f3_stats 同一约定）。
"""
import importlib.util, json, os, pathlib, subprocess, sys, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "09b_aadr_summary.py"
if importlib.util.find_spec("pandas") is None:
    _SKIP = "pandas not available (the production script imports it)"
else:
    _SKIP = ""

# 合成面板：目标 + 两个现代个体（同组凑 n>=2）+ 一个古代个体。PC 值任意但互不相同。
PANEL = [
    # (iid, label, kind, source_population_id, PC1, PC2, PC3, PC4)
    ("TESTSAMPLE", "Target",  "target",  "target",     0.00, 0.00, 0.00, 0.00),
    ("HAN1",       "Han",     "modern",  "Han",        0.10, 0.01, 0.02, 0.03),
    ("HAN2",       "Han",     "modern",  "Han",        0.12, 0.02, 0.01, 0.04),
    ("AM1",        "China_Am","ancient", "China_Am",   0.30, 0.05, 0.06, 0.07),
]


def _write_config(td, **over):
    cfg = {
        "sample_id": "TESTSAMPLE",
        "work_dir": str(pathlib.Path(td) / "work"),
        "aadr_prefix": str(pathlib.Path(td) / "work/data/ref/aadr/panel"),
        "ancestry_min_call_rate_modern": 0.95,
        "ancestry_min_call_rate_target": 0.95,
        "aadr_min_call_rate_ancient": 0.50,
        "ancestry_min_projection_snps": 10000,
        "ancestry_min_group_n": 2,
    }
    cfg.update(over)
    p = pathlib.Path(td) / "config.yaml"
    p.write_text("".join(f"{k}: {v}\n" for k, v in cfg.items()), encoding="utf-8")
    return p


def _write_panel(td, smiss_rows):
    """合成 11_aadr：sscore + samples + 预计算 smiss。smiss_rows: {(iid): (MISSING_CT, OBS_CT)}。"""
    w = pathlib.Path(td) / "work" / "wgs" / "11_aadr"
    w.mkdir(parents=True, exist_ok=True)
    hdr = "#FID\tIID\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG\tPHENO\tC1\tC2\tC3\tC4"
    lines = [f"{lab}\t{iid}\t" + "\t".join(f"{v:.4f}" for v in pcs) + "\tNA\t0\t0\t0\t0"
             for iid, lab, kind, pop, *pcs in PANEL]
    (w / "proj.sscore").write_text("\n".join([hdr] + lines) + "\n", encoding="utf-8")
    (w / "samples.tsv").write_text(
        "iid\tlabel\tkind\tsource_population_id\tdate_mean_bp\n" +
        "".join(f"{iid}\t{lab}\t{kind}\t{pop}\t{1000 if kind == 'ancient' else ''}\n"
                for iid, lab, kind, pop, *_ in PANEL), encoding="utf-8")
    smiss = ["#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS"]
    for iid, (mc, oc) in smiss_rows.items():
        smiss.append(f"X\t{iid}\t{mc}\t{oc}\t{mc / oc:.6f}")
    (w / "aadr.pruned.smiss").write_text("\n".join(smiss) + "\n", encoding="utf-8")
    (w / "prune.prune.in").write_text("rs1\nrs2\n", encoding="utf-8")   # 存在即可：smiss 已预写
    return w


def _run_09b(cfg):
    return subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                          env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=120)


def _summary(w):
    return json.loads((w / "summary.json").read_text(encoding="utf-8"))


@unittest.skipIf(_SKIP, _SKIP)
class TestCoverageCaliber(unittest.TestCase):
    """n_called 必须是 OBS_CT−MISSING_CT：plink2 的 OBS_CT 是分母，不是已调用数。"""

    def test_n_called_subtracts_missing(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_panel(td, {                      # AM1 缺 30000/130000，已调用应为 100000
                "TESTSAMPLE": (130, 130000),
                "HAN1": (0, 130000), "HAN2": (650, 130000),
                "AM1": (30000, 130000),
            })
            r = _run_09b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            recs = {x["iid"]: x for x in _summary(w)["records"]}
            self.assertEqual(recs["AM1"]["n_called_snps"], 100000)
            self.assertAlmostEqual(recs["AM1"]["call_rate"], 1 - 30000 / 130000, places=5)
            self.assertEqual(recs["HAN1"]["n_called_snps"], 130000,
                             "零缺失时 n_called 才等于 OBS_CT")


@unittest.skipIf(_SKIP, _SKIP)
class TestTargetGate(unittest.TestCase):
    """目标不过覆盖门槛就没有排名：合成"旧有限投影 + 目标在最终位点集上全缺失"不能再 state=ok。"""

    def test_fully_missing_target_yields_no_ranking(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_panel(td, {
                "TESTSAMPLE": (130000, 130000),        # 目标在剪枝集上 100% 缺失
                "HAN1": (0, 130000), "HAN2": (650, 130000),
                "AM1": (30000, 130000),
            })
            r = _run_09b(cfg)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("target_below_coverage_gate", r.stderr + r.stdout)
            self.assertFalse((w / "summary.json").exists(), "不过门槛不得写出带排名的 summary")
            self.assertFalse((w / "near_modern.tsv").exists())

    def test_target_absent_from_smiss_fails_the_gate(self):
        """目标不在缺失率文件里 → 覆盖未知，按不过门槛处理，而不是默认通过。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_panel(td, {                   # smiss 里没有 TESTSAMPLE 这一行
                "HAN1": (0, 130000), "HAN2": (650, 130000), "AM1": (30000, 130000),
            })
            r = _run_09b(cfg)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("target_below_coverage_gate", r.stderr + r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
