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


# ─────────────────────────────────────────────────── AN1: metadata, dosage, locations

ANNO_COLS = {
    "gid": 'Genetic ID (suffices: ".DG" is a high coverage shotgun genome; ".AG" Agilent 1240K)',
    "master": "Persistent Genetic ID",
    "iid": "Individual ID",
    "group": "Group ID",
    "locality": "Locality",
    "political": "Political Entity",
    "lat": "Latitude",
    "lon": "Longitude",
    "dmean": "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date]",
    "dsd": "Date standard deviation in BP [OxCal sigma]",
    "draw": "Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age",
    "dbasis": "Method for Determining Date; unless otherwise specified, calibrations use 95.4% intervals",
    "suffices": "Suffices (indicating data types used for sources which can be a subset of that in bam)",
    "yhg": "Y haplogroup in terminal mutation notation automatically called based on Y-full 12.03",
    "yhg_isogg": "Y haplogroup  in ISOGG notation automatically called based on Yfull 12.03",
    "mthg": "mtDNA haplogroup if >2x or published",
    "pub": "Publication abbreviation",
}


def anno_row(**kw):
    """一行 .anno 数据，列名与真实文件同前缀（真实列名很长，必须靠前缀别名适配）。"""
    row = {v: ".." for v in ANNO_COLS.values()}
    for k, v in kw.items():
        row[ANNO_COLS[k]] = v
    return row


class TestNormalizeMetadata(unittest.TestCase):
    """normalize_metadata 把 .anno 行变成 7.3 的统一记录：列名适配只在这里发生。"""

    def test_aliases_and_missing_optional_columns_become_null(self):
        rec, = ad.normalize_metadata([anno_row(gid="Loschbour.AG", master="Loschbour", iid="Loschbour",
                                              group="Luxembourg_Loschbour_Mesolithic",
                                              locality="Loschbour (Echternach)", political="Luxembourg",
                                              lat="49.81", lon="6.40", dmean="8025", dsd="64",
                                              suffices="AG", yhg="I-V6473", mthg="U5b1a",
                                              pub="MathiesonReichNature2018")],
                                     dataset="AADR", release="v66.p1")
        self.assertEqual(rec["record_id"], "Loschbour.AG")
        self.assertEqual(rec["individual_id"], "Loschbour")
        self.assertEqual(rec["source_population_id"], "Luxembourg_Loschbour_Mesolithic")
        self.assertEqual(rec["dataset"], "AADR")
        self.assertEqual(rec["reference_release"], "v66.p1")
        self.assertEqual(rec["genotype_representation"], "AG")
        self.assertEqual(rec["y_hg_raw"], "I-V6473")
        self.assertEqual(rec["mt_hg_raw"], "U5b1a")
        self.assertAlmostEqual(rec["latitude"], 49.81)
        self.assertAlmostEqual(rec["longitude"], 6.40)
        self.assertEqual(rec["location_precision"], "site")
        self.assertEqual(rec["date_mean_bp"], 8025)
        self.assertEqual(rec["date_sd_bp"], 64)
        # 没有提供的可选字段是 None，不是 ""、不是 0
        self.assertIsNone(rec["date_min_bp"])
        self.assertIsNone(rec["date_max_bp"])
        self.assertIsNone(rec["distance_to_target"])
        self.assertIsNone(rec["group_id"])

    def test_unknown_date_stays_null_and_never_becomes_zero(self):
        rec, = ad.normalize_metadata([anno_row(gid="X.SG", master="X", iid="X", group="G",
                                               lat="..", lon="..", dmean="..", dsd="..", draw="..")],
                                     dataset="AADR", release="v66")
        self.assertIsNone(rec["date_mean_bp"], "未知年代不得填 0 伪装现代")
        self.assertIsNone(rec["date_min_bp"])
        self.assertIsNone(rec["date_max_bp"])
        self.assertIsNone(rec["date_raw"])

    def test_contextual_range_gives_bounds_and_marks_the_basis(self):
        # 有区间就不是"未知"：允许取中点，但必须标明是区间中点（BP 基准 1950）
        rec, = ad.normalize_metadata([anno_row(gid="Y.AG", master="Y", iid="Y", group="G",
                                               dmean="..", dsd="..", draw="2500-1700 BCE")],
                                     dataset="AADR", release="v66")
        self.assertEqual(rec["date_raw"], "2500-1700 BCE")
        self.assertEqual((rec["date_min_bp"], rec["date_max_bp"]), (3650, 4450))
        self.assertEqual(rec["date_mean_bp"], 4050)
        self.assertEqual(rec["date_basis"], "contextual_range_midpoint")
        self.assertLessEqual(rec["date_min_bp"], rec["date_max_bp"])
        # 95.4% CI 里的 calBCE 区间同样解析，但基准仍是 1950
        ci, = ad.normalize_metadata([anno_row(gid="Z.AG", master="Z", iid="Z", group="G",
                                              dmean="7205", dsd="50", dbasis="Direct: IntCal20",
                                              draw="6221-5986 calBCE (7205+-50 BP, OxA-7738)")],
                                    dataset="AADR", release="v66")
        self.assertEqual((ci["date_min_bp"], ci["date_max_bp"]), (7936, 8171))
        self.assertEqual(ci["date_mean_bp"], 7205, "有实测均值时不得用区间中点覆盖它")
        self.assertEqual(ci["date_basis"], "Direct: IntCal20")

    def test_coordinates_out_of_range_or_empty_never_become_zero(self):
        rec, = ad.normalize_metadata([anno_row(gid="A.AG", master="A", iid="A", group="G",
                                               lat="..", lon="..")], dataset="AADR", release="v66")
        self.assertIsNone(rec["latitude"])
        self.assertIsNone(rec["longitude"])
        self.assertEqual(rec["location_precision"], "unknown")
        bad, = ad.normalize_metadata([anno_row(gid="B.AG", master="B", iid="B", group="G",
                                               lat="95", lon="200")], dataset="AADR", release="v66")
        self.assertIsNone(bad["latitude"], "越界纬度必须丢弃，不能保留假位置")
        self.assertIsNone(bad["longitude"])
        self.assertEqual(bad["location_precision"], "unknown")

    def test_missing_required_columns_fail_loudly(self):
        with self.assertRaises(ValueError):
            ad.normalize_metadata([{"Group ID": "G"}], dataset="AADR", release="v66")

    def test_master_id_deduplication_is_decided_by_the_caller(self):
        """同一 Master ID 的多种表示只算一个人：normalize 如实保留，去重规则单独可测。"""
        rows = [anno_row(gid="P.SG", master="P", iid="P", group="G", lat="1", lon="2", suffices="SG"),
                anno_row(gid="P.AG", master="P", iid="P", group="G", lat="1", lon="2", suffices="AG")]
        recs = ad.normalize_metadata(rows, dataset="AADR", release="v66")
        self.assertEqual(len(recs), 2, "两行都要保留")
        pick = ad.dedupe_by_master_id(recs)
        self.assertEqual(len(pick["kept"]), 1)
        self.assertEqual(len(pick["dropped"]), 1)
        self.assertEqual(pick["dropped"][0]["reason_code"], "duplicate_representation")
        # 推荐表示优先（SG/DG 优于 AG/TW）
        self.assertEqual(pick["kept"][0]["genotype_representation"], "SG")

    def test_kind_comes_from_the_panel_selection_not_from_a_date_cutoff(self):
        """§7：不用 500 BP 硬切。900 BP 的记录若被点名算古代，它仍是古代。"""
        rows = [anno_row(gid="R.SG", master="R", iid="R", group="China_Recent_IA", dmean="900"),
                anno_row(gid="M.HO", master="M", iid="M", group="Han", dmean="0")]
        recs = ad.normalize_metadata(rows, dataset="AADR", release="v66")
        out = ad.assign_kind(recs, ancient_prefixes=["China_"], modern_groups=["Han"])
        by_id = {r["record_id"]: r for r in out}
        self.assertEqual(by_id["R.SG"]["kind"], "ancient")
        self.assertEqual(by_id["M.HO"]["kind"], "modern")
        # 没被任何选择命中的记录是 unknown，不进入 PCA 集合
        other = ad.assign_kind(ad.normalize_metadata([anno_row(gid="Z.AG", master="Z", iid="Z", group="Zzz")],
                                                     dataset="AADR", release="v66"),
                               ancient_prefixes=["China_"], modern_groups=["Han"])
        self.assertEqual(other[0]["kind"], "unknown")

    def test_overlapping_modern_and_ancient_selection_is_rejected(self):
        with self.assertRaises(ValueError):
            ad.assign_kind([], ancient_prefixes=["Han"], modern_groups=["Han_sub"])

    def test_locations_table_is_validated(self):
        with tempfile.TemporaryDirectory() as td:
            f = pathlib.Path(td) / "ancestry_locations.tsv"
            cols = ["dataset", "source_id", "location_id", "label_zh", "label_en", "locality",
                    "latitude", "longitude", "precision", "source_url", "note"]
            data = [
                ("AADR", "China_Baligang_LN_Longshan", "CN-HA-Baligang", "河南 八里岗", "Baligang, Henan",
                 "Baligang (Dengzhou)", "32.7", "112.1", "site", "https://example.org/x", "corrected to Dengzhou"),
                ("AADR", "China_Tibet_Kangyu", "CN-XZ-Kangyu", "", "Kangyu", "Kangyu", "", "", "unknown", "", "no coords"),
            ]
            f.write_text("# header comment\n" + "\t".join(cols) + "\n"
                         + "\n".join("\t".join(r) for r in data) + "\n", encoding="utf-8")
            loc = ad.load_locations(f)
            self.assertEqual(len(loc), 2)
            a = loc["CN-HA-Baligang"]
            self.assertEqual(a["label_zh"], "河南 八里岗")
            self.assertAlmostEqual(a["latitude"], 32.7)
            # 无中文名时保留原始名称：不再有"必须汉字开头"的校验
            b = loc["CN-XZ-Kangyu"]
            self.assertEqual(b["label_zh"], "")
            self.assertIsNone(b["latitude"])
            self.assertEqual(b["precision"], "unknown")
        # 缺列必须报错，而不是静默
        with tempfile.TemporaryDirectory() as td:
            f = pathlib.Path(td) / "bad.tsv"
            f.write_text("dataset\tsource_id\nAADR\tx\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                ad.load_locations(f)

    def test_locations_do_not_merge_different_sites_of_one_group(self):
        """同一 Group ID 的三个遗址是三个地点，不得合成一个坐标。"""
        recs = ad.normalize_metadata([
            anno_row(gid="S.AG", master="S1", iid="S1", group="China_MLBA", locality="Baiyangcun", lat="25.9", lon="100.2"),
            anno_row(gid="S2.AG", master="S2", iid="S2", group="China_MLBA", locality="Guchengcun", lat="27.1", lon="101.0"),
        ], dataset="AADR", release="v66")
        keys = {ad.location_key_of(r) for r in recs}
        self.assertEqual(len(keys), 2, "同组不同遗址必须是不同地点")


class TestA1Dosage(unittest.TestCase):
    """a1_dosage：只数能确定的等位基因，其余一律 None（绝不猜链翻转或把缺失当参考纯合）。"""

    def test_textbook_cases(self):
        self.assertEqual(ad.a1_dosage("0|1", "A", "G", "G"), 1)
        self.assertEqual(ad.a1_dosage("0/0", "A", "G", "A"), 2)
        self.assertIsNone(ad.a1_dosage("0/.", "A", "G", "G"))
        self.assertIsNone(ad.a1_dosage("2/2", "A", "G", "G"))

    def test_only_identical_or_swapped_alleles_are_accepted(self):
        self.assertEqual(ad.a1_dosage("1|1", "A", "G", "G"), 2)
        self.assertEqual(ad.a1_dosage("1|1", "A", "G", "A"), 0)
        # a1 既不是 REF 也不是 ALT（例如链翻转过的第三等位）→ 缺失，不猜
        self.assertIsNone(ad.a1_dosage("0|1", "A", "G", "C"))
        # REF/ALT 互换：a1 等于 ALT 时计数必须翻转
        self.assertEqual(ad.a1_dosage("0|1", "G", "A", "G"), 1)
        self.assertEqual(ad.a1_dosage("0|1", "G", "A", "A"), 1)

    def test_missing_ploidy_and_missingness_are_reported_with_reasons(self):
        self.assertIsNone(ad.a1_dosage("./.", "A", "G", "G"))
        self.assertIsNone(ad.a1_dosage(".", "A", "G", "G"))
        d, why = ad.a1_dosage_with_reason("0|1|1", "A", "G", "G")
        self.assertIsNone(d); self.assertEqual(why, "not_diploid")
        d, why = ad.a1_dosage_with_reason("0/.", "A", "G", "G")
        self.assertEqual(why, "partial_missing")
        d, why = ad.a1_dosage_with_reason("./.", "A", "G", "G")
        self.assertEqual(why, "all_missing")
        d, why = ad.a1_dosage_with_reason("2/2", "A", "G", "G")
        self.assertEqual(why, "allele_not_observable")
        d, why = ad.a1_dosage_with_reason("0|1", "A", "G,T", "T")
        self.assertEqual(why, "not_biallelic")


if __name__ == "__main__":
    unittest.main(verbosity=2)
