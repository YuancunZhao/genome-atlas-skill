"""AN2 验收：断言对着**生产 09b 脚本**跑，不对着测试里的复制品。

子进程执行 scripts/09b_aadr_summary.py：$WGS_CONFIG 指向临时 config.yaml，work 根与 11_aadr
产物全部在临时目录合成。缺失率经**假 plink2 替身**从 --out 落盘（复审 AN2/H6 后 09b 每次重算
aadr.pruned.smiss，不再吃存在性缓存），替身把每个测试想要的行复制到 $out.smiss 并留调用凭证。
HANDOFF 要求"验收必须对生产09b断言"——辅助函数绿灯不能替代这里。生产脚本 import pandas，当前
解释器没有 pandas 时整文件跳过（本地套件跳过、服务器套件执行，与 test_f3_stats 同一约定）。
"""
import hashlib, importlib.util, json, os, pathlib, subprocess, sys, tempfile, unittest

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
        "plink2": str(pathlib.Path(td) / "fake_plink2"),   # _write_panel 落盘的替身
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


def _write_panel(td, smiss_rows, preexisting_stale_smiss=None):
    """合成 11_aadr：sscore + samples + bed/bim/fam + 假 plink2 替身。

    09b 现在**每次重算** aadr.pruned.smiss（复审 AN2/H6：存在性缓存挡不住换目标/换 prune 集），
    所以缺失率不再预写为成品，而是经替身 plink2 从 --out 落盘——替身把本测试想要的行写到
    <td>/stub.smiss 再复制过去，并 touch plink_called 留证。preexisting_stale_smiss 用来预置
    一份**旧** smiss：重算必须覆盖它，否则那份旧文件就是"换输入后命中缓存"的复现。
    smiss_rows: {(iid): (MISSING_CT, OBS_CT)}。"""
    td = pathlib.Path(td)
    w = td / "work" / "wgs" / "11_aadr"
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
    (td / "stub.smiss").write_text("\n".join(smiss) + "\n", encoding="utf-8")
    stub = td / "fake_plink2"
    stub.write_text("#!/bin/bash\n"
                    "out=\"\"; prev=\"\"\n"
                    "for a in \"$@\"; do [ \"$prev\" = \"--out\" ] && out=\"$a\"; prev=\"$a\"; done\n"
                    f"cp \"{td}/stub.smiss\" \"$out.smiss\" && touch \"{td}/plink_called\"\n"
                    "exit 0\n", encoding="utf-8")
    stub.chmod(0o755)
    if preexisting_stale_smiss is not None:
        (w / "aadr.pruned.smiss").write_text(preexisting_stale_smiss, encoding="utf-8")
    (w / "aadr.bed").write_text("", encoding="utf-8")   # 存在即可：替身 plink2 不读内容
    (w / "aadr.bim").write_text("", encoding="utf-8")
    (w / "aadr.fam").write_text("", encoding="utf-8")
    (w / "prune.prune.in").write_text("rs1\nrs2\n", encoding="utf-8")
    return w


def _run_09b(cfg):
    return subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                          env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=120)


def _summary(w):
    return json.loads((w / "summary.json").read_text(encoding="utf-8"))


@unittest.skipIf(_SKIP, _SKIP)
class TestSmissIsRecomputed(unittest.TestCase):
    """复审 AN2/H6：aadr.pruned.smiss 是**每次重算**的，存在性缓存已删。预置一份旧 smiss
    （上一轮 prune 集算出的、数字全错）后运行：plink2 必须真的被调用，旧文件必须被覆盖，
    门槛用的是本次替换身落盘的数字——"换 prune 集/换目标后命中旧缓存"从此不复现。"""

    def test_stale_smiss_is_overwritten_not_reused(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            # 替身给出：目标 0 缺失/130000（过门槛），其余个体也全过
            rows = {"TESTSAMPLE": (0, 130000), "HAN1": (100, 130000),
                    "HAN2": (200, 130000), "AM1": (60000, 130000)}
            stale = ("#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS\n"
                     "X\tTESTSAMPLE\t130000\t130000\t1.000000\n")   # 旧 prune 集：目标全缺失
            w = _write_panel(td, rows, preexisting_stale_smiss=stale)
            r = _run_09b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue((pathlib.Path(td) / "plink_called").exists(),
                            "smiss 已存在也必须重算（plink2 被真正调用）")
            self.assertNotIn("1.000000", (w / "aadr.pruned.smiss").read_text(encoding="utf-8"),
                             "旧的 stale smiss 必须被本次结果覆盖")
            t = _summary(w)["target"]
            self.assertEqual(t["n_called_snps"], 130000, "门槛读的是本次重算的数字，不是旧缓存")
            # manifest 记 prune 指纹（复审 AN0/AN5/H6）：30 据此拒绝"换 prune 集后的旧 summary"
            man = json.loads((w / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(man["parameters"]["prune_sites"], 2)
            self.assertEqual(man["parameters"]["prune_sha"],
                             hashlib.sha256((w / "prune.prune.in").read_bytes()).hexdigest()[:12])


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


@unittest.skipIf(_SKIP, _SKIP)
class TestKindAnnotations(unittest.TestCase):
    """modern/ancient 连续两筛互不污染：每条记录带的是自己 kind 的判定。"""

    def test_both_kinds_keep_their_own_verdicts(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            # HAN2 覆盖 0.846 < 0.95 → modern 判不合格；AM1 覆盖 0.769 >= 0.50 → ancient 判合格。
            # 09b 先筛 modern 再筛 ancient：修复前第二遍把合格的 HAN1 改成 other_kind。
            w = _write_panel(td, {
                "TESTSAMPLE": (130, 130000),
                "HAN1": (0, 130000), "HAN2": (20000, 130000), "AM1": (30000, 130000),
            })
            r = _run_09b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            recs = {x["iid"]: x for x in _summary(w)["records"]}
            self.assertTrue(recs["HAN1"]["eligible"])
            self.assertIsNone(recs["HAN1"]["exclusion_reason"])
            self.assertFalse(recs["HAN2"]["eligible"])
            self.assertEqual(recs["HAN2"]["exclusion_reason"], "low_call_rate")
            self.assertTrue(recs["AM1"]["eligible"], "ancient 记录带 ancient 轮的判定")
            # 目标不是任何一轮判定的对象：不带合格标注，也不计入被排除
            self.assertNotIn("eligible", recs["TESTSAMPLE"])
            counts = _summary(w)["counts"]
            self.assertEqual(counts["excluded"], 1)


@unittest.skipIf(_SKIP, _SKIP)
class TestGroupingBySite(unittest.TestCase):
    """同组不同遗址分开成两组；同一 iid 出现两行只算一个人。"""

    def _panel_with_sites(self, td):
        w = _write_panel(td, {
            "TESTSAMPLE": (130, 130000),
            "HAN1": (0, 130000), "HAN2": (650, 130000), "AM1": (30000, 130000),
        })
        # HAN1/HAN2 同属 Han 群体，但位于两个遗址：修复前被合成一个 n=2 的"群体"，
        # location_id 取第一个地点，还挤进默认排名。
        (w / "reference_metadata.tsv").write_text(
            "record_id\tlocation_id\tlocality\tlatitude\tlongitude\tdate_min_bp\tdate_max_bp\n"
            "HAN1\tlocA\tSiteA\t30\t120\t\t\n"
            "HAN2\tlocB\tSiteB\t31\t121\t\t\n"
            "AM1\tlocC\tSiteC\t32\t122\t500\t2000\n", encoding="utf-8")
        return w

    def test_same_group_at_two_sites_is_two_groups(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = self._panel_with_sites(td)
            r = _run_09b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            groups = _summary(w)["groups"]["modern"]
            self.assertEqual({g["group_id"] for g in groups}, {"Han|locA", "Han|locB"},
                             "同一群体在两个遗址是两组，不能平均到一起")
            self.assertTrue(all(g["n"] == 1 for g in groups), "两遗址各 1 人，没有 n=2 的合并组")
            self.assertTrue(all(g["rank"] is None for g in groups), "n=1 的组不进默认排名")

    def test_a_duplicated_iid_counts_once(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = self._panel_with_sites(td)
            with open(w / "proj.sscore", "a", encoding="utf-8") as fh:   # HAN1 再来一行
                fh.write("Han\tHAN1\t0.1100\t0.0200\t0.0300\t0.0400\tNA\t0\t0\t0\t0\n")
            r = _run_09b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            counts = _summary(w)["counts"]
            self.assertEqual(counts["selected"], 4, "同一 iid 两行只算一个人")
            self.assertIn("duplicated iid", r.stderr)


@unittest.skipIf(_SKIP, _SKIP)
class TestTargetBlockCarriesRealCoverage(unittest.TestCase):
    """summary 的 target 块要带合并 smiss 后的覆盖数字，不是 samples.tsv 的提取阶段陈旧值。"""

    def test_target_block_reads_the_merged_frame(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            w = _write_panel(td, {
                "TESTSAMPLE": (650, 130000),           # 目标真实覆盖 129350 / 0.995
                "HAN1": (0, 130000), "HAN2": (650, 130000), "AM1": (30000, 130000),
            })
            # 旧格式 samples.tsv 自带提取阶段的 call_rate（目标恰为 1.0）：smiss 合并重建 s 之后，
            # target 块若仍读合并前的旧帧，就会原样回显这个 1.0，而 n_called_snps 整列缺失得 None。
            (w / "samples.tsv").write_text(
                "iid\tlabel\tkind\tdate\tcall_rate\n"
                "TESTSAMPLE\tTarget\ttarget\t\t1.0\n"
                "HAN1\tHan\tmodern\t\t1.0\n"
                "HAN2\tHan\tmodern\t\t1.0\n"
                "AM1\tChina_Am\tancient\t1000\t0.8\n", encoding="utf-8")
            r = _run_09b(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            t = _summary(w)["target"]
            self.assertEqual(t["n_called_snps"], 130000 - 650)
            self.assertAlmostEqual(t["call_rate"], 1 - 650 / 130000, places=5)
            self.assertNotEqual(t["call_rate"], 1.0, "不得回显 samples.tsv 的提取阶段值")


if __name__ == "__main__":
    unittest.main(verbosity=2)
