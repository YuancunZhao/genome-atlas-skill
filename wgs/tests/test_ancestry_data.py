"""AN0 的行为测试：配置校验（read_options）与产物 manifest（读写、指纹匹配）。

约束（§7 AN0）：测试**不读真实 config.yaml、不创建目录、不导入 wgsconfig**——wgsconfig 在 import
时就会 mkdir 工作目录，导入它会让测试产生副作用。ancestry_data 因此是纯模块：只接收 dict，返回 dict。
"""
import json, os, pathlib, subprocess, sys, tempfile, unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import ancestry_data as ad  # noqa: E402


class TestReadOptions(unittest.TestCase):
    """read_options(cfg) -> dict：把配置翻译成开关与已校验参数，非法输入一律 ValueError。"""

    def test_everything_off_by_default(self):
        o = ad.read_options({})
        self.assertFalse(o["regional_enabled"])
        self.assertFalse(o["aadr_enabled"])
        self.assertFalse(o["local_enabled"])
        # 默认值仍然是可用的参数，只是开关关着
        self.assertEqual(o["min_group_n"], 2)
        self.assertEqual(o["min_projection_snps"], 10000)
        self.assertAlmostEqual(o["min_call_rate_ancient"], 0.50)
        self.assertAlmostEqual(o["min_call_rate_modern"], 0.95)
        self.assertAlmostEqual(o["min_call_rate_target"], 0.95)

    def test_explicit_superpop_turns_regional_on(self):
        self.assertFalse(ad.read_options({"ref_superpop": ""})["regional_enabled"])
        self.assertTrue(ad.read_options({"ref_superpop": "EUR"})["regional_enabled"])
        # 键存在但值为空 = 明确关掉，不回退默认
        self.assertEqual(ad.read_options({"ref_superpop": ""})["ref_superpop"], "")

    def test_aadr_needs_a_selector(self):
        self.assertFalse(ad.read_options({"aadr_modern": []})["aadr_enabled"])
        self.assertTrue(ad.read_options({"aadr_modern": ["Han"]})["aadr_enabled"])
        self.assertTrue(ad.read_options({"aadr_ancient_prefix": ["China_"]})["aadr_enabled"])

    def test_local_ancestry_needs_two_disjoint_sources(self):
        self.assertFalse(ad.read_options({})["local_enabled"])
        self.assertFalse(ad.read_options({"local_ancestry_a": ["CHB"]})["local_enabled"])
        ok = ad.read_options({"local_ancestry_a": ["CHB", "JPT"], "local_ancestry_b": ["CDX"]})
        self.assertTrue(ok["local_enabled"])
        # 控制面板可以为空，标签默认按来源名生成
        self.assertEqual(ok["la_control"], [])
        self.assertEqual(len(ok["la_labels"]), 2)

    def test_string_instead_of_list_is_rejected(self):
        for key in ("ref_subpops", "axis_pops", "local_ancestry_a", "local_ancestry_b",
                    "local_ancestry_control", "local_ancestry_labels", "aadr_modern",
                    "aadr_ancient_prefix", "local_ancestry_calibration_pops",
                    "local_ancestry_calibration_chroms"):
            with self.assertRaises(ValueError, msg=key):
                ad.read_options({key: "CHB"})  # 字符串会被 list() 拆成字符，必须挡住

    def test_overlapping_or_duplicate_local_ancestry_is_rejected(self):
        with self.assertRaises(ValueError):
            ad.read_options({"local_ancestry_a": ["CHB"], "local_ancestry_b": ["CHB"]})
        with self.assertRaises(ValueError):
            ad.read_options({"local_ancestry_a": ["CHB"], "local_ancestry_b": ["CDX"],
                             "local_ancestry_labels": ["X", "X"]})
        with self.assertRaises(ValueError):
            ad.read_options({"local_ancestry_a": ["CHB", "CHB"], "local_ancestry_b": ["CDX"]})

    def test_unsupported_build_fails_before_any_work(self):
        self.assertEqual(ad.read_options({})["build"], "GRCh37")
        with self.assertRaises(ValueError):
            ad.read_options({"build": "GRCh38"})

    def test_thresholds_are_range_checked(self):
        for bad in (0, -0.1, 1.5, "half"):
            with self.assertRaises(ValueError, msg=repr(bad)):
                ad.read_options({"ancestry_min_call_rate_modern": bad})
        for bad in (0, -3):
            with self.assertRaises(ValueError, msg=repr(bad)):
                ad.read_options({"ancestry_min_group_n": bad})
        for bad in (0, "many"):
            with self.assertRaises(ValueError, msg=repr(bad)):
                ad.read_options({"ancestry_min_projection_snps": bad})

    def test_calibration_defaults_and_validation(self):
        o = ad.read_options({"local_ancestry_a": ["CHB"], "local_ancestry_b": ["CDX"]})
        self.assertEqual(o["calibration_chroms"], ["1", "2", "6", "22"])
        self.assertEqual(o["calibration_n"], 20)
        self.assertEqual(o["calibration_seed"], 1)
        # 未配置轴群体时校准群体跟随来源面板
        self.assertEqual(o["calibration_pops"], ["CHB", "CDX"])
        with self.assertRaises(ValueError):
            ad.read_options({"local_ancestry_calibration_n": 0})

    def test_empty_ancestor_prefix_string_is_rejected(self):
        with self.assertRaises(ValueError):
            ad.read_options({"aadr_ancient_prefix": ["China_", ""]})
        with self.assertRaises(ValueError):
            ad.read_options({"aadr_ancient_prefix": [""]})

    def test_display_names_do_not_touch_the_identity_keys(self):
        o = ad.read_options({"sample_id": "S1", "name_zh": "名字 带空格", "name_en": "Name B"})
        self.assertEqual(o["sample_id"], "S1")
        self.assertEqual(o["name_zh"], "名字 带空格")
        self.assertEqual(o["name_en"], "Name B")
        # 显示名不得参与内部键：命名空间化后的目标键仍然只由 sample_id 决定
        self.assertEqual(ad.target_key("S1"), ad.target_key("S1"))
        self.assertNotEqual(ad.target_key("S1"), ad.target_key("S2"))

    def test_partial_config_never_raises_on_missing_keys(self):
        o = ad.read_options({"threads": "4", "mem_gb": "16"})
        self.assertEqual(o["threads"], "4")
        self.assertEqual(o["mem_gb"], "16")


class TestManifest(unittest.TestCase):
    """manifest 是最小事实记录：写、读、比对。比对失败必须让调用方拒绝旧产物。"""

    def _m(self, **kw):
        base = {"schema_version": 1, "sample_id": "S1", "analysis_id": "aadr-modern-v1",
                "state": "ok", "reason_code": "", "reference_release": "v66.p1",
                "build": "GRCh37", "parameters": {"min_call_rate": 0.95},
                "tool_versions": {"plink2": "v2.0"}, "input_fingerprints": {"aadr": "sha256:x"},
                "outputs": ["11_aadr/summary.json"]}
        base.update(kw)
        return base

    def test_sample_id_mismatch_rejects(self):
        self.assertFalse(ad.manifest_matches({"sample_id": "A"}, {"sample_id": "B"}))
        self.assertTrue(ad.manifest_matches({"sample_id": "A"}, {"sample_id": "A"}))

    def test_fingerprint_change_rejects(self):
        old, expected = self._m(), self._m()
        expected["input_fingerprints"] = {"aadr": "sha256:y"}
        self.assertFalse(ad.manifest_matches(old, expected))

    def test_parameter_change_rejects(self):
        old, expected = self._m(), self._m()
        expected["parameters"] = {"min_call_rate": 0.5}
        self.assertFalse(ad.manifest_matches(old, expected))

    def test_reference_release_and_build_change_reject(self):
        a, b = self._m(), self._m(reference_release="v54")
        self.assertFalse(ad.manifest_matches(a, b))
        self.assertFalse(ad.manifest_matches(self._m(), self._m(build="GRCh38")))

    def test_incomplete_or_unfinished_manifest_rejects(self):
        self.assertFalse(ad.manifest_matches(None, self._m()))
        self.assertFalse(ad.manifest_matches({}, self._m()))
        self.assertFalse(ad.manifest_matches(self._m(state="failed"), self._m()))
        self.assertFalse(ad.manifest_matches(self._m(state="unavailable"), self._m()))
        self.assertFalse(ad.manifest_matches(self._m(schema_version=2), self._m()))
        self.assertFalse(ad.manifest_matches(self._m(analysis_id="other"), self._m()))

    def test_state_vocabulary_is_enforced(self):
        with self.assertRaises(ValueError):
            ad.write_manifest(pathlib.Path(tempfile.mkdtemp()) / "m.json", self._m(state="finished"))

    def test_write_is_atomic_and_read_roundtrips(self):
        with tempfile.TemporaryDirectory() as td:
            p = pathlib.Path(td) / "out" / "manifest.json"
            ad.write_manifest(p, self._m())
            got = ad.read_manifest(p)
            self.assertEqual(got["sample_id"], "S1")
            self.assertEqual(got["state"], "ok")
            self.assertEqual(ad.read_manifest(pathlib.Path(td) / "missing.json"), None)
            # 临时文件不得留在目录里
            self.assertEqual([f.name for f in p.parent.iterdir()], ["manifest.json"])

    def test_disabled_analyses_still_produce_a_manifest(self):
        """§7：显式禁用也要写 manifest，让 30 能区分"没跑"和"跑了但禁用"。"""
        with tempfile.TemporaryDirectory() as td:
            p = pathlib.Path(td) / "manifest.json"
            ad.write_manifest(p, ad.disabled_manifest("S1", "regional-pca", "disabled_by_config"))
            got = ad.read_manifest(p)
            self.assertEqual(got["state"], "disabled")
            self.assertEqual(got["reason_code"], "disabled_by_config")
            self.assertFalse(ad.manifest_matches(got, self._m()))


class TestConfigIsolation(unittest.TestCase):
    """§7 AN0：两份配置各用各的 work_dir，A→B→A 不读对方的缓存；指纹变了旧 manifest 不被接受。"""

    SCRIPT = ("import sys; sys.path.insert(0, %r); import wgsconfig as c; "
              "print(c.SAMPLE, c.P, c.REGIONAL_ENABLED)")

    def _run(self, cfg_path, scripts_dir):
        env = {**os.environ, "WGS_CONFIG": str(cfg_path)}
        r = subprocess.run([sys.executable, "-c", self.SCRIPT % str(scripts_dir)],
                           capture_output=True, text=True, env=env, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr[-500:])
        return r.stdout.strip().split()

    def test_two_configs_are_isolated_and_reusable(self):
        scripts = pathlib.Path(__file__).resolve().parents[1] / "scripts"
        with tempfile.TemporaryDirectory() as td:
            td = pathlib.Path(td)
            a, b = td / "a.yaml", td / "b.yaml"
            a.write_text(f"sample_id: SAMPLE_A\nwork_dir: {td/'workA'}\n", encoding="utf-8")
            b.write_text(f"sample_id: SAMPLE_B\nwork_dir: {td/'workB'}\nref_superpop: EUR\n", encoding="utf-8")
            sa, pa, ra = self._run(a, scripts)
            sb, pb, rb = self._run(b, scripts)
            sa2, pa2, _ = self._run(a, scripts)
            self.assertNotEqual(pa, pb, "两份配置必须落在不同 work_dir")
            self.assertEqual((sa, pa), (sa2, pa2), "A→B→A 必须回到同一目标与同一缓存目录")
            self.assertEqual(sa, "SAMPLE_A")
            self.assertEqual(ra, "False", "A 没写 ref_superpop，区域分析不得自己打开")
            self.assertEqual(rb, "True", "B 显式写了 ref_superpop")

    def test_bad_config_fails_before_work_happens(self):
        scripts = pathlib.Path(__file__).resolve().parents[1] / "scripts"
        with tempfile.TemporaryDirectory() as td:
            bad = pathlib.Path(td) / "bad.yaml"
            bad.write_text(f"sample_id: S\nwork_dir: {pathlib.Path(td) / 'w'}\nbuild: GRCh38\n", encoding="utf-8")
            env = {**os.environ, "WGS_CONFIG": str(bad)}
            r = subprocess.run([sys.executable, "-c", self.SCRIPT % str(scripts)],
                               capture_output=True, text=True, env=env, timeout=120)
            self.assertNotEqual(r.returncode, 0, "不支持的 build 必须在计算前失败")
            self.assertIn("GRCh38", r.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
