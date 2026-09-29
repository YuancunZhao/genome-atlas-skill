"""AN5 的行为测试：报告侧契约（HANDOFF §7 的最小验收）。

这些用例只依赖 ancestry_data 的纯函数——构建器 30/31 在 import 时就会读数据、建目录，不能拿来
当测试夹具；被它们调用的逻辑因此都放在这里。
"""
import json, pathlib, subprocess, sys, tempfile, unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import ancestry_data as ad  # noqa: E402


class TestPayloadHygiene(unittest.TestCase):
    """载荷不能把页面弄坏，也不能不是合法 JSON。"""

    EVIL = "</script><img src=x onerror=alert(1)>"

    def test_an_injected_name_stays_text(self):
        # §7：姓名注入时只显示文本——载荷里不允许出现裸的闭合序列，否则浏览器提前结束脚本块。
        lit = ad.js_string_literal(self.EVIL)
        self.assertNotIn("</", lit)
        self.assertEqual(json.loads(lit), self.EVIL, "转义后数据必须一模一样")

    def test_nan_and_infinity_become_null(self):
        txt = ad.json_text({"a": float("nan"), "b": [float("inf"), float("-inf"), 1.0], "c": "x"})
        self.assertNotIn("NaN", txt)
        self.assertNotIn("Infinity", txt)
        self.assertEqual(json.loads(txt), {"a": None, "b": [None, None, 1.0], "c": "x"})

    def test_chunked_payload_round_trips(self):
        # 模拟 31 的分块：切块 → 各自写成 JS 字面量 → 拼回 → 解析，数据必须完全一致
        payload = ad.json_text({"x": "y" * 200, "z": self.EVIL})
        chunks = [payload[i:i + 16] for i in range(0, len(payload), 16)]
        joined = "".join(json.loads(ad.js_string_literal(c)) for c in chunks)
        self.assertEqual(json.loads(joined), {"x": "y" * 200, "z": self.EVIL})

    def test_numpy_scalars_become_native(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("numpy not available")
        txt = ad.json_text({"n": np.int64(7), "f": np.float64(1.5), "nan": np.float64("nan")})
        self.assertEqual(json.loads(txt), {"n": 7, "f": 1.5, "nan": None})


class TestStaleManifestRejected(unittest.TestCase):
    """§7：伪造旧 sample_id 的 manifest 必须被拒绝，而不是拿来出报告。"""

    def _dir(self, td, sample, state="ok", summary='{"counts": {"selected": 1}}'):
        d = pathlib.Path(td)
        ad.write_manifest(d / "manifest.json", ad.build_manifest(sample, "aadr-human-origins", state=state))
        if summary is not None:
            (d / "summary.json").write_text(summary, encoding="utf-8")
        return d

    def test_forged_sample_id_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A")
            state, reason, doc = ad.analysis_state(d, {"sample_id": "SAMPLE_B"})
            self.assertEqual((state, reason), ("unavailable", "stale_result"))
            self.assertIsNone(doc, "指纹不符时不能把结果交给调用方")
            self.assertEqual(ad.analysis_state(d, {"sample_id": "SAMPLE_A"})[0], "ok")

    def test_missing_manifest_is_not_a_result(self):
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A")
            (d / "manifest.json").unlink()
            self.assertEqual(ad.analysis_state(d, {"sample_id": "SAMPLE_A"})[1], "missing_manifest")

    def test_unreadable_and_absent_results(self):
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A", summary="{not json")
            self.assertEqual(ad.analysis_state(d, {"sample_id": "SAMPLE_A"})[1], "unreadable_result")
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A", summary=None)
            self.assertEqual(ad.analysis_state(d, {"sample_id": "SAMPLE_A"})[1], "missing_result")

    def test_disabled_state_is_reported_not_hidden(self):
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A", state="disabled")
            state, reason, doc = ad.analysis_state(d, {"sample_id": "SAMPLE_A"})
            self.assertEqual(state, "disabled")
            self.assertIsNotNone(doc, "禁用也要能读出结果文件（里面写着原因）")

    def test_step_level_disabled_manifest_overrides_a_healthy_directory_one(self):
        """run_all 禁用某步时写 manifest.<step>.json；它必须压过目录里上一次运行留下的
        manifest.json——否则 A→B→A 换回配置后，旧结果照旧以 ok 进报告（复审 AN0+AN5+H6）。"""
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A")            # 健康的 manifest.json + summary.json（上次运行的）
            ad.write_manifest(d / "manifest.09b-aadr-summary.json",
                              ad.disabled_manifest("SAMPLE_A", "09b-aadr-summary", "aadr_not_configured"))
            state, reason, doc = ad.analysis_state(d, {"sample_id": "SAMPLE_A"})
            self.assertEqual((state, reason), ("disabled", "aadr_not_configured"))
            self.assertIsNone(doc, "禁用时不能把旧 summary.json 当结果交给调用方")
        # 步骤级记录全部是 ok 时，回落到目录级 manifest：不因存在 manifest.*.json 而误判
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A")
            ad.write_manifest(d / "manifest.09b-aadr-summary.json",
                              ad.build_manifest("SAMPLE_A", "09b-aadr-summary"))
            self.assertEqual(ad.analysis_state(d, {"sample_id": "SAMPLE_A"})[0], "ok")

    def test_disable_then_rerun_success_ends_the_disabled_state(self):
        """复审 AN0/AN5/H6 生命周期：禁用→重新启用并成功产出后，状态必须回到 ok。
        生产者成功时通过 clear_step_manifests 撤下**自己的**步骤级记录——不撤别人的。"""
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A")
            ad.write_manifest(d / "manifest.09b-aadr-summary.json",
                              ad.disabled_manifest("SAMPLE_A", "09b-aadr-summary", "aadr_not_configured"))
            ad.write_manifest(d / "manifest.08-aadr-extract.json",
                              ad.disabled_manifest("SAMPLE_A", "08-aadr-extract", "aadr_not_configured"))
            self.assertEqual(ad.analysis_state(d, {"sample_id": "SAMPLE_A"})[0], "disabled")
            # 09b 本次成功：只撤自己的记录；08 的禁用证据必须留下（不粗暴清场）
            removed = ad.clear_step_manifests(d, "09b-aadr-summary")
            self.assertEqual(removed, ["manifest.09b-aadr-summary.json"])
            self.assertTrue((d / "manifest.08-aadr-extract.json").exists(),
                            "clearing one step's record must not delete another step's evidence")
            # 08 仍有一条 disabled 记录压着目录——照实报告，直到 08 也成功撤下
            self.assertEqual(ad.analysis_state(d, {"sample_id": "SAMPLE_A"})[0], "disabled")
            ad.clear_step_manifests(d, "08-aadr-extract")
            state, reason, doc = ad.analysis_state(d, {"sample_id": "SAMPLE_A"})
            self.assertEqual((state, doc), ("ok", {"counts": {"selected": 1}}))

    def test_clear_step_is_a_noop_without_a_record(self):
        with tempfile.TemporaryDirectory() as td:
            self.assertEqual(ad.clear_step_manifests(td, "17b-la-calibrated"), [])

    def test_expected_parameters_gate_config_changes(self):
        """复审 AN0/AN5/H6：30 的准入要比 sample_id 更多——同一样本换了门槛/prune 集/参考子集，
        旧 manifest 必须判 stale_result。parameters 里**缺键**同样 stale：没记录不等于一致。"""
        with tempfile.TemporaryDirectory() as td:
            d = self._dir(td, "SAMPLE_A")
            m = ad.read_manifest(d / "manifest.json")
            m["parameters"] = {"min_group_n": 2, "prune_sha": "abc123def456"}
            ad.write_manifest(d / "manifest.json", m)
            # 完全一致 → ok
            self.assertEqual(ad.analysis_state(
                d, {"sample_id": "SAMPLE_A"},
                expected_parameters={"min_group_n": 2, "prune_sha": "abc123def456"})[0], "ok")
            # 换门槛 → stale_result（旧结果是另一套参数算的）
            self.assertEqual(ad.analysis_state(
                d, {"sample_id": "SAMPLE_A"}, expected_parameters={"min_group_n": 5})[1], "stale_result")
            # 换 prune 集 → stale_result
            self.assertEqual(ad.analysis_state(
                d, {"sample_id": "SAMPLE_A"}, expected_parameters={"prune_sha": "ffffffffffff"})[1],
                "stale_result")
            # manifest 根本没记录这个键 → stale_result，不能当"一致"
            self.assertEqual(ad.analysis_state(
                d, {"sample_id": "SAMPLE_A"}, expected_parameters={"reference_release": "v66"})[1],
                "stale_result")
            # 数据依赖键不传就不比：missing_chroms 由本次运行决定，不是准入条件
            self.assertEqual(ad.analysis_state(d, {"sample_id": "SAMPLE_A"})[0], "ok")

    def test_manifest_matches_rejects_missing_expected_keys(self):
        """复审 AN0/AN5/H6：此前非 REQUIRED 键缺失会跳过——没有 reference_release 的 manifest
        也能通过按 reference_release 的比对。现在任何 expected 键缺失都拒绝。"""
        bare = {"sample_id": "A", "state": "ok"}
        self.assertFalse(ad.manifest_matches(bare, {"sample_id": "A", "reference_release": "v66"}),
                        "a manifest that never recorded the reference must not pass a reference check")
        full = {**bare, "reference_release": "v66"}
        self.assertTrue(ad.manifest_matches(full, {"sample_id": "A", "reference_release": "v66"}))

    def test_run_all_clear_cli_matches_the_documented_invocation(self):
        """run_all 用 `ancestry_data.py --clear-step <id> --dir <dir>` 撤记录；这个入口本身
        要按 run_all 的写法被驱动一次（生产路径），而不是只测函数。"""
        with tempfile.TemporaryDirectory() as td:
            d = pathlib.Path(td)
            ad.write_manifest(d / "manifest.17b-la-calibrated.json",
                              ad.disabled_manifest("S1", "17b-la-calibrated", "x"))
            r = subprocess.run(
                [sys.executable, str(pathlib.Path(__file__).resolve().parents[1]
                                     / "scripts" / "ancestry_data.py"),
                 "--clear-step", "17b-la-calibrated", "--dir", str(d)],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertFalse((d / "manifest.17b-la-calibrated.json").exists())
            # 幂等：再次清除同样退出 0
            r2 = subprocess.run(
                [sys.executable, str(pathlib.Path(__file__).resolve().parents[1]
                                     / "scripts" / "ancestry_data.py"),
                 "--clear-step", "17b-la-calibrated", "--dir", str(d)],
                capture_output=True, text=True)
            self.assertEqual(r2.returncode, 0, r2.stderr)


class TestEligibleSetIsShared(unittest.TestCase):
    """§7：同一 Analysis 进入所有视图的合格集合完全相同。"""

    RECS = [
        {"record_id": "a", "kind": "ancient", "call_rate": 0.9, "n_called_snps": 60000, "distance_to_target": 0.01},
        {"record_id": "b", "kind": "ancient", "call_rate": 0.4, "n_called_snps": 60000, "distance_to_target": 0.001},
        {"record_id": "c", "kind": "modern", "call_rate": 0.99, "n_called_snps": 60000, "distance_to_target": 0.02},
    ]

    def test_repeated_calls_give_the_same_ids_in_the_same_order(self):
        ids = lambda: [r["record_id"] for r in ad.eligible_records(
            [dict(x) for x in self.RECS], kind="ancient", min_rate=0.5, min_snps=10000)]
        self.assertEqual(ids(), ids())
        self.assertEqual(ids(), ["a"], "低覆盖的 b 被排除，而它恰是最接近的那个")

    def test_kind_filter_keeps_views_apart(self):
        anc = ad.eligible_records([dict(x) for x in self.RECS], kind="ancient")
        mod = ad.eligible_records([dict(x) for x in self.RECS], kind="modern")
        self.assertEqual({r["record_id"] for r in anc} & {r["record_id"] for r in mod}, set())


class TestDegradedModes(unittest.TestCase):
    """§7：有效现代 + 禁用古代、有 mt 缺 Y，都要能产出成对的双语文案。"""

    def test_modern_only_without_ancient(self):
        copy = ad.ancestry_copy({"ho_affinity": [{"label": "Han", "kind": "modern", "n": 5, "d": 0.02}]})
        self.assertIn("现代人群", copy["c_anc"][0])
        self.assertNotIn("古代群体", copy["c_anc"][0])
        for k, v in copy.items():
            self.assertEqual(len(v), 2, k)
            self.assertTrue(all(isinstance(x, str) and x.strip() for x in v), k)

    def test_no_projection_says_so(self):
        copy = ad.ancestry_copy({})
        self.assertIn("没有", copy["c_anc"][0])
        self.assertIn("no usable", copy["c_anc"][1].lower())

    def test_missing_y_is_explained_not_invented(self):
        copy = ad.ancestry_copy({"y_state": "unavailable",
                                 "mt": {"hg": "A13", "quality": 0.92}})
        self.assertIn("没有", copy["c_y"][0])
        self.assertNotIn("属于", copy["c_y"][0], "没有判定时不能写'属于某支系'")

    def test_generated_copy_never_claims_identity(self):
        copy = ad.ancestry_copy({"ho_affinity": [{"label": "Han_Chongqing", "kind": "modern", "n": 3, "d": 0.01}],
                                 "ho_near_individual": [{"iid": "X.SG", "d": 0.003}]})
        self.assertIn("不由此做亲缘或族群推断", copy["n_anc"][0])
        self.assertIn("no kinship or ethnic inference", copy["n_anc"][1])


if __name__ == "__main__":
    unittest.main(verbosity=2)
