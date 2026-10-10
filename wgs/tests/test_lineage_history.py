"""AN4 的行为测试：父母系的树关系匹配与历史证据组装。

纯函数测试：不读真实 YFull/PhyloTree 树、不跑 haplogrep3、不改报告数据。
"""
import json, pathlib, sys, tempfile, unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import lineage_history as lh  # noqa: E402

try:                       # 生产 05 子进程需要 pandas（wgsconfig 还要 yaml）
    import pandas  # noqa: F401
    _SKIP05 = ""
except ImportError:
    _SKIP05 = "pandas not available (the production 05 script imports it)"


class TestMatchLineage(unittest.TestCase):
    """match_lineage(query, record, parents) 回答的是"记录相对查询支系是什么关系"。"""

    PARENTS = {"mt:A": None, "mt:A1": "mt:A", "mt:A13": "mt:A", "mt:A13a": "mt:A13"}

    def test_plan_examples(self):
        self.assertEqual(lh.match_lineage("mt:A13", "mt:A13a", self.PARENTS), "descendant")
        self.assertEqual(lh.match_lineage("mt:A1", "mt:A13", self.PARENTS), "unrelated")
        self.assertEqual(lh.match_lineage("mt:A13", "mt:A", self.PARENTS), "ancestor")

    def test_exact_and_self(self):
        self.assertEqual(lh.match_lineage("mt:A13", "mt:A13", self.PARENTS), "exact")
        self.assertEqual(lh.match_lineage("mt:A13a", "mt:A13", self.PARENTS), "ancestor")

    def test_namespaces_never_mix(self):
        # Y 与 mt 的单倍群名可以重名（N 两支都有）；不区分命名空间就会把 Y 的记录挂到 mt 树上。
        parents = {"y:A13": "y:A", "mt:A13": "mt:A"}
        self.assertEqual(lh.match_lineage("mt:A13", "y:A13", parents), "unresolved")
        self.assertEqual(lh.match_lineage("y:A13", "y:A", parents), "ancestor")

    def test_missing_edges_are_unresolved_not_guessed(self):
        # 树版本对不上 / 没有树边：不能靠字符串前缀推断。
        parents = {"mt:A": None, "mt:A13": "mt:A"}
        self.assertEqual(lh.match_lineage("mt:A13", "mt:A13a", parents), "unresolved")
        self.assertEqual(lh.match_lineage("mt:A13", "mt:Z9", parents), "unresolved")
        self.assertEqual(lh.match_lineage("mt:Z9", "mt:A13", parents), "unresolved")

    def test_same_version_canonical_ids_stay_exact_without_edges(self):
        # §7：缺树边时，只有可信的同版本规范 ID 可以 exact；其余 unresolved。
        self.assertEqual(lh.match_lineage("mt:A13", "mt:A13", {}, same_tree=True), "exact")
        self.assertEqual(lh.match_lineage("mt:A13", "mt:A13", {}, same_tree=False), "unresolved")

    def test_prefix_matching_is_not_used(self):
        # startswith("A1") 会把 A13 误认成 A1 的下游。
        parents = {"mt:A": None, "mt:A1": "mt:A", "mt:A13": "mt:A"}
        self.assertEqual(lh.match_lineage("mt:A1", "mt:A13", parents), "unrelated")
        self.assertNotEqual(lh.match_lineage("mt:A1", "mt:A13", parents), "exact")


class TestCanonicalize(unittest.TestCase):
    """AADR 的标签来自 YFull 12.03，我们的树是 14.06.0：比较之前必须规范化，并说清做了什么。"""

    def test_placeholder_labels_are_missing_not_branches(self):
        # AADR 用这些字符串表示"没有标签"（实测 (female) 9409 行、(sex unknown) 1122 行）
        for v in ("", "..", "nan", "n/a", "n/a (female)", "n/a (sex unknown)", "N/A"):
            self.assertEqual(lh.canonicalize(v, "y", {}), ("", "missing"), v)

    def test_alias_maps_an_older_name_to_the_current_node(self):
        h = {"aliases": {"y:I-V6473": "y:I-V6473-新名"}}
        self.assertEqual(lh.canonicalize("I-V6473", "y", h), ("y:I-V6473-新名", "alias"))

    def test_node_present_in_the_tree(self):
        self.assertEqual(lh.canonicalize("N-CTS4714", "y", {}, known_nodes={"y:N-CTS4714"}),
                         ("y:N-CTS4714", "none"))

    def test_absent_node_is_flagged_as_version_mismatch(self):
        # 当前树里没有该节点：明确标注版本不匹配，而不是当成"无关支系"悄悄丢掉
        self.assertEqual(lh.canonicalize("I-V6473", "y", {}, known_nodes={"y:N-CTS4714"}),
                         ("y:I-V6473", "version_mismatch"))

    def test_without_a_node_set_it_says_unverified(self):
        self.assertEqual(lh.canonicalize("I-V6473", "y", {}), ("y:I-V6473", "unverified"))


class TestLineageResult(unittest.TestCase):
    """05/06 的结构化结果：判定字段、不可用状态、以及 caller 与人工复核的分工。"""

    def _y(self, **kw):
        base = {"kind": "y", "state": "ok", "reason_code": "", "reported_hg": "N-CTS4714",
                "conservative_hg": "N-M1845", "tree_source": "YFull", "tree_version": "v14.0",
                "call_quality": "manual_review", "supported_path": [
                    {"node": "N", "der": 420, "anc": 0, "na": 0}, {"node": "N-CTS4714", "der": 3, "anc": 0, "na": 2}],
                "uncertain_nodes": [], "route_review": "reviewed_by_hand", "conflicts": []}
        base.update(kw)
        return base

    def test_terminal_is_never_taken_as_proven_just_because_it_is_deepest(self):
        r = self._y(conservative_hg="N-M1845",
                    uncertain_nodes=[{"node": "N-CTS4714", "reason": "only 3 supporting sites"}])
        out = lh.lineage_summary(r)
        self.assertEqual(out["reported_hg"], "N-CTS4714")
        self.assertEqual(out["conservative_hg"], "N-M1845", "保守回退必须保留")
        self.assertEqual(out["uncertain_nodes"][0]["node"], "N-CTS4714")
        self.assertEqual(out["call_quality"], "manual_review", "人工复核与 caller 不能混为一谈")

    def test_unavailable_states_carry_a_reason(self):
        r = lh.unavailable_lineage("y", "no_haplogroup_call", "the call set carries no Y interval")
        self.assertEqual(r["state"], "unavailable")
        self.assertEqual(r["reason_code"], "no_haplogroup_call")
        self.assertEqual(r["supported_path"], [])
        self.assertIsNone(r["reported_hg"])
        summary = lh.lineage_summary(r)
        self.assertEqual(summary["state"], "unavailable")
        self.assertTrue(summary["reason_code"])

    def test_empty_path_does_not_crash_or_invent_a_terminal(self):
        r = self._y(supported_path=[], reported_hg=None)
        out = lh.lineage_summary(r)
        self.assertIsNone(out["terminal"], "没有路径就不能有末端支系")
        self.assertEqual(out["state"], "unavailable" if not r.get("reported_hg") else out["state"])

    def test_conflicting_manual_and_automatic_calls_are_kept_and_flagged(self):
        r = self._y(reported_hg="N-CTS4714", conservative_hg="N-M1845",
                    conflicts=[{"field": "y_hg", "automatic": "N-CTS4714", "manual": "N-M1845"}])
        out = lh.lineage_summary(r)
        self.assertEqual(len(out["conflicts"]), 1)
        self.assertFalse(out["strict_match_allowed"],
                         "有未解决的冲突时，严格匹配视图必须排除该记录")


class TestObservationsAndRoutes(unittest.TestCase):
    """历史记录与迁移路线：不套用常染色体门槛，不凭证据缺失造故事。"""
    ROWS = [
        {"record_id": "X.SG", "individual_id": "X", "master_id": "X", "dataset": "AADR",
         "record_release": "v66", "y_hg_raw": "N-CTS4714", "mt_hg_raw": "A13", "call_rate": 0.49,
         "locality": "Somewhere", "latitude": 30.0, "longitude": 110.0, "location_precision": "site",
         "date_mean_bp": 3300, "date_min_bp": 3000, "date_max_bp": 3600, "hg_source_tree": "YFull12.03"},
        {"record_id": "X.AG", "individual_id": "X", "master_id": "X", "dataset": "AADR",
         "record_release": "v66", "y_hg_raw": "N-CTS4714", "mt_hg_raw": "A13", "call_rate": 0.90,
         "locality": "Somewhere", "latitude": 30.0, "longitude": 110.0, "location_precision": "site",
         "date_mean_bp": 3300, "date_min_bp": 3000, "date_max_bp": 3600, "hg_source_tree": "YFull12.03"},
        {"record_id": "Y.SG", "individual_id": "Y", "master_id": "Y", "dataset": "AADR",
         "record_release": "v66", "y_hg_raw": "O-F438", "mt_hg_raw": "D4", "call_rate": 0.99,
         "locality": "Elsewhere", "latitude": 35.0, "longitude": 115.0, "location_precision": "region",
         "date_mean_bp": None, "date_min_bp": None, "date_max_bp": None, "hg_source_tree": "YFull12.03"},
    ]

    def test_ancient_records_below_the_autosomal_gate_still_count(self):
        # 0.49 覆盖在常染色体视图里会被排除；历史视图用的是来源对 Y/mt 的可用性，不是常染色体覆盖率。
        low = [dict(self.ROWS[0])]                       # 只有 X.SG（call_rate 0.49）
        obs = lh.lineage_observations(low, "mt:A13", parents={}, tree_kind="mt")
        self.assertEqual([o["record_id"] for o in obs], ["X.SG"],
                         "低覆盖的古人只要来源给了 mt 标签，就仍然进历史视图")

    def test_one_person_counts_once_across_representations(self):
        obs = lh.lineage_observations(self.ROWS, "mt:A13", parents={}, tree_kind="mt")
        self.assertEqual([o["record_id"] for o in obs], ["X.SG"],
                         "同一 Master ID 的 .SG/.AG 只算一次（.SG 优先）")

    def test_records_without_a_call_are_skipped_and_unknown_dates_stay_null(self):
        obs = lh.lineage_observations(self.ROWS, "mt:D4", parents={}, tree_kind="mt")
        self.assertEqual([o["record_id"] for o in obs], ["Y.SG"])
        self.assertIsNone(obs[0]["date_range"]["mean"], "未知年代保持 null，不填 0")
        self.assertEqual(obs[0]["precision"], "region")

    def test_placeholder_haplogroups_never_match_a_branch(self):
        rows = [{"record_id": "F.AG", "master_id": "F", "y_hg_raw": "n/a (female)", "call_rate": 0.9},
                {"record_id": "U.SG", "master_id": "U", "y_hg_raw": "n/a (sex unknown)", "call_rate": 0.9},
                {"record_id": "R.SG", "master_id": "R", "y_hg_raw": "R", "call_rate": 0.9}]
        obs = lh.lineage_observations(rows, "y:R", parents={}, tree_kind="y")
        self.assertEqual([o["record_id"] for o in obs], ["R.SG"])

    def test_observations_carry_the_required_fields(self):
        o = lh.lineage_observations(self.ROWS, "mt:A13", parents={}, tree_kind="mt")[0]
        for k in ("record_id", "node_id", "relation", "locality", "coordinates", "precision",
                  "date_range", "date_basis", "call_source", "publication"):
            self.assertIn(k, o)

    def test_version_mismatched_label_never_claims_exact(self):
        """复审 AN4：标签不在当前树（version_mismatch）时，与查询字符串相等不等于版本等价。

        旧行为：canonicalize 记 version_mismatch 只进计数，match_lineage 仍按 same_tree=True
        判 exact——不兼容树标签可以冒充同版本精确匹配。现在一律 unresolved，不进默认视图。
        """
        rows = [{"record_id": "V.SG", "master_id": "V", "dataset": "AADR",
                 "y_hg_raw": "N-CTS4714", "hg_source_tree": "YFull12.03", "call_rate": 0.9,
                 "locality": "Somewhere"}]
        stats = {}
        obs = lh.lineage_observations(rows, "y:N-CTS4714", parents={}, tree_kind="y",
                                      known_nodes={"y:N-M1845"}, stats=stats)
        self.assertEqual(stats.get("version_mismatch"), 1)
        self.assertEqual(obs, [], "未验证版本等价的记录不得冒充 exact 进默认视图")
        # 对照：同一标签在当前树节点集里（note=none）时，字符串相等才可以说 exact
        ok_stats = {}
        ok = lh.lineage_observations(rows, "y:N-CTS4714", parents={}, tree_kind="y",
                                     known_nodes={"y:N-CTS4714"}, stats=ok_stats)
        self.assertEqual(ok_stats.get("none"), 1)
        self.assertEqual([o["relation"] for o in ok], ["exact"])

    def test_weak_chain_above_the_terminal_forces_a_step_back(self):
        """真实的 12 / 1 / 0 / 1 / 5 形态：末端那个 5 只是弱链的末尾，不是独立证据。"""
        path = [("N-Z4762", 132, 0, 1), ("N-F2905", 57, 0, 2), ("N-CTS12473", 119, 0, 1),
                ("N-M1845", 12, 0, 0), ("N-M1928", 1, 0, 0), ("N-Y125475", 0, 0, 0),
                ("N-M1793", 1, 0, 0), ("N-CTS4714", 5, 0, 0)]
        self.assertEqual(lh.conservative_from_path(path, solid=5, tail=4), "N-M1845")

    def test_solid_terminal_is_kept(self):
        path = [("N", 420, 0, 0), ("N-M1845", 12, 0, 0), ("N-A", 9, 0, 0), ("N-B", 7, 0, 0)]
        self.assertEqual(lh.conservative_from_path(path, solid=5, tail=4), "N-B")

    def test_empty_and_single_level_paths(self):
        self.assertIsNone(lh.conservative_from_path([], solid=5, tail=4))
        self.assertEqual(lh.conservative_from_path([("N", 420, 0, 0)], solid=5, tail=4), "N")
        self.assertIsNone(lh.conservative_from_path([("N-A", 1, 0, 0)], solid=5, tail=4),
                          "唯一一级就弱，且上面没有别的：保守落点为空而不是硬报它")

    def test_reviewed_call_beats_any_automatic_rule(self):
        """复核记录优先于自动规则——但只在**它属于这个样本、且绑定对得上**时（复审 AN4 + AN4-a）。

        旧版本不带样本就返回结论，那会让一个样本的人工判断顺着支系名传播给所有同支系样本；
        AN4-a 再加两层绑定：记录写明它核对过的末端 SNP 表（evidence_sha256）与项目树版本
        （tree_version），调用方传本次运行的对应值，全部一致才生效。
        """
        hist = lh.load_history({"schema_version": 1, "reviewed_calls": {
            "y:N-CTS4714": {"conservative_hg": "N-M1845", "reason": "1-5 sites per level below",
                            "sample_id": "S1", "tree_version": "14.06.0",
                            "evidence_sha256": "a" * 64}}})
        got = lh.reviewed_call(hist, "y", "N-CTS4714", sample_id="S1",
                               tree_version="14.06.0", evidence_sha="a" * 64)
        self.assertEqual(got, ("N-M1845", "1-5 sites per level below"))
        self.assertEqual(lh.reviewed_call(hist, "y", "y:N-CTS4714", sample_id="S1",
                                          tree_version="14.06.0", evidence_sha="a" * 64)[0],
                         lh.reviewed_call(hist, "y", "N-CTS4714", sample_id="S1",
                                          tree_version="14.06.0", evidence_sha="a" * 64)[0],
                         "带/不带 kind 前缀的键等价")
        self.assertEqual(lh.reviewed_call(hist, "mt", "A13", sample_id="S1",
                                          tree_version="14.06.0", evidence_sha="a" * 64), (None, None),
                         "别的支系不受影响")
        self.assertEqual(lh.reviewed_call({}, "y", "N-CTS4714", sample_id="S1",
                                          tree_version="14.06.0", evidence_sha="a" * 64), (None, None))

    def test_a_review_belongs_to_one_sample_not_to_the_branch(self):
        """同支系的**另一个**样本不得继承本样本的复核结论——这是 AN4 的核心。"""
        hist = lh.load_history({"schema_version": 1, "reviewed_calls": {
            "y:N-CTS4714": {"conservative_hg": "N-M1845", "reason": "reviewed for S1", "sample_id": "S1",
                            "tree_version": "14.06.0", "evidence_sha256": "a" * 64}}})
        said = []
        got = lh.reviewed_call(hist, "y", "N-CTS4714", sample_id="S2",
                               tree_version="14.06.0", evidence_sha="a" * 64, log=said.append)
        self.assertEqual(got, (None, None), "S2 调出同一支系也不得套用 S1 的结论")
        self.assertTrue(any("belongs to S1" in m for m in said), said)
        self.assertEqual(lh.reviewed_call(hist, "y", "N-CTS4714", sample_id=None,
                                          tree_version="14.06.0", evidence_sha="a" * 64), (None, None),
                         "未指明样本时不套用")

    def test_an_unmarked_review_is_not_a_generic_rule(self):
        """没有 sample_id 的旧记录不得被当成通用规则——它无法与"通用规则"区分，而正是这种混同要修。"""
        hist = lh.load_history({"schema_version": 1, "reviewed_calls": {
            "y:N-CTS4714": {"conservative_hg": "N-M1845", "reason": "no owner recorded"}}})
        said = []
        self.assertEqual(lh.reviewed_call(hist, "y", "N-CTS4714", sample_id="S1", log=said.append),
                         (None, None))
        self.assertTrue(any("no sample_id" in m for m in said), said)

    def test_no_routes_without_evidence(self):
        hist = lh.load_history({"schema_version": 1, "tree_source": "YFull", "tree_version": "v14.0",
                                "parents": {}, "aliases": {}, "routes": []})
        out = lh.history_view("mt:A13", hist, [], parents={})
        self.assertEqual(out["routes"], [])
        self.assertEqual(out["history_state"], "distribution_only")
        self.assertTrue(out["history_reason_code"])

    def test_routes_require_a_tree_version_and_waypoints(self):
        hist = lh.load_history({"schema_version": 1, "routes": [
            {"route_id": "r1", "applicable_node": "mt:A13", "tree_version": "v17.2",
             "waypoints": [{"locality": "A", "latitude": 1, "longitude": 2}], "time_range": [1000, 2000],
             "sources": ["doi:10.0/x"]}]})
        self.assertEqual(len(hist["routes"]), 1)
        ok = lh.history_view("mt:A13", hist, [], parents={})
        self.assertEqual(len(ok["routes"]), 1, "有来源的路线应匹配到查询支系")
        self.assertEqual(lh.history_view("mt:Z9", hist, [], parents={})["routes"], [])


class TestObservationSourceTiers(unittest.TestCase):
    """AN4 查询扩展：观测来源分三层（完整元数据 > 已配置 .anno > 09b 筛选 summary），
    且"没查过"绝不写成"数据集没有记录"。

    旧行为：CLI 只读 11_aadr/summary.json（08/09 的 PCA 筛选子集），AADR 关闭时缺查询或沿用
    旧 summary；来源缺失时 history_view 照样报 distribution_only/no_records_in_this_dataset——
    把"从来没查过"冒充"查过、没有"。这些用例在旧代码上全部失败。
    """

    METADATA_COLS = ("record_id\tindividual_id\tmaster_id\tdataset\treference_release\t"
                     "genotype_representation\ty_hg_raw\tmt_hg_raw\tlocality\tlatitude\tlongitude\t"
                     "location_precision\tdate_mean_bp\tdate_min_bp\tdate_max_bp\thg_source_tree\t"
                     "publication\n")

    def _yard(self, td):
        yard = pathlib.Path(td) / "03_haplo"
        yard.mkdir(parents=True, exist_ok=True)
        (yard / "y_result.json").write_text(json.dumps({
            "kind": "y", "state": "ok", "reported_hg": "N-CTS4714",
            "supported_path": [{"node": "N-CTS4714", "der": 5, "anc": 0, "na": 0}]}), encoding="utf-8")
        return yard

    def _write_metadata(self, td, rows):
        p = pathlib.Path(td) / "reference_metadata.tsv"
        p.write_text(self.METADATA_COLS + "".join("\t".join(str(c) for c in r) + "\n" for r in rows),
                     encoding="utf-8")
        return p

    def test_full_metadata_is_preferred_over_the_screened_summary(self):
        """完整元数据里的记录（09b 筛选会丢的那种）必须能被查到——这正是"未读完整元数据"要修的。"""
        with tempfile.TemporaryDirectory() as td:
            yard = self._yard(td)
            # 两行：一行在 summary 里也有；一行低覆盖、只存在于完整元数据（09b 的 PCA 筛选会丢）
            meta = self._write_metadata(td, [
                ("Keep.SG", "Keep", "Keep", "AADR", "v66", "SG", "N-CTS4714", "", "Site", "30.0", "110.0",
                 "site", "3300", "3000", "3600", "YFull12.03", "pub"),
                ("ScreenedOut.DG", "ScreenedOut", "ScreenedOut", "AADR", "v66", "DG", "N-CTS4714", "",
                 "Other", "31.0", "111.0", "region", "", "", "", "YFull12.03", "pub"),
            ])
            out = pathlib.Path(td) / "lineage_history.json"
            rc = lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1",
                          "--history", "", "--metadata", str(meta)])
            self.assertEqual(rc, 0)
            doc = json.loads(out.read_text(encoding="utf-8"))
            q = doc["y"]["history"]["query"]
            self.assertEqual(q["source_kind"], "reference_metadata")
            self.assertEqual(q["n_rows"], 2)
            self.assertEqual(q["state"], "ok")
            self.assertEqual(sorted(o["record_id"] for o in doc["y"]["history"]["observations"]),
                             ["Keep.SG", "ScreenedOut.DG"],
                             "完整元数据里被 09b 筛掉的记录也要能查到")
            # 数值列从 TSV 字符串转回来了（模板要拿它们做投影）
            obs = [o for o in doc["y"]["history"]["observations"] if o["record_id"] == "Keep.SG"][0]
            self.assertEqual(obs["coordinates"], {"latitude": 30.0, "longitude": 110.0})
            self.assertEqual(obs["date_range"]["mean"], 3300)

    def test_anno_tier_queries_when_aadr_pca_is_off(self):
        """AADR 关闭（08/09b 没跑、没有 reference_metadata.tsv）时，已配置 .anno 仍可查询。"""
        with tempfile.TemporaryDirectory() as td:
            yard = self._yard(td)
            anno = pathlib.Path(td) / "panel.anno"
            anno.write_text(
                "Genetic ID\tPersistent Genetic ID\tIndividual ID\tLocality\tLatitude\tLongitude\t"
                "Date mean in BP\tDate standard deviation in BP\tFull Date\tMethod for Determining Date\t"
                "Suffices\tY haplogroup in terminal\tmtDNA haplogroup\tPublication abbreviation\n"
                "I1.SG\tI1\tI1\tSite\t30\t110\t3300\t150\t3000-3600 calBCE\tradiocarbon\tSG\t"
                "N-CTS4714\t\t pub\n",
                encoding="utf-8")
            out = pathlib.Path(td) / "lineage_history.json"
            rc = lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1",
                          "--history", "", "--anno", str(anno)])
            self.assertEqual(rc, 0)
            doc = json.loads(out.read_text(encoding="utf-8"))
            q = doc["y"]["history"]["query"]
            self.assertEqual((q["source_kind"], q["state"]), ("anno", "ok"))
            self.assertEqual(q["n_rows"], 1)
            self.assertEqual([o["record_id"] for o in doc["y"]["history"]["observations"]], ["I1.SG"])

    def test_anno_missing_required_columns_is_unreadable_not_an_empty_dataset(self):
        """.anno 缺必需列：如实报"来源不可读"，不是静默空数据集。"""
        with tempfile.TemporaryDirectory() as td:
            yard = self._yard(td)
            anno = pathlib.Path(td) / "broken.anno"
            anno.write_text("Some\tOther\tColumns\nx\ty\tz\n", encoding="utf-8")
            out = pathlib.Path(td) / "lineage_history.json"
            rc = lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1",
                          "--history", "", "--anno", str(anno)])
            self.assertEqual(rc, 0)
            h = json.loads(out.read_text(encoding="utf-8"))["y"]["history"]
            self.assertEqual(h["query"]["state"], "unreadable")
            self.assertEqual(h["history_state"], "unavailable")
            self.assertEqual(h["history_reason_code"], "history_source_unreadable")
            self.assertEqual(h["observations"], [])

    def test_no_source_at_all_is_not_available_not_no_records(self):
        """没有任何可查来源：历史资料不可用。旧行为报 no_records_in_this_dataset——把没查过写成没有记录。"""
        with tempfile.TemporaryDirectory() as td:
            yard = self._yard(td)                    # 没有 11_aadr/，没有 .anno，没有 summary
            out = pathlib.Path(td) / "lineage_history.json"
            rc = lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1", "--history", ""])
            self.assertEqual(rc, 0)
            h = json.loads(out.read_text(encoding="utf-8"))["y"]["history"]
            self.assertEqual(h["query"]["state"], "not_available")
            self.assertEqual(h["history_state"], "unavailable")
            self.assertEqual(h["history_reason_code"], "history_source_not_available")
            self.assertNotEqual(h["history_reason_code"], "no_records_in_this_dataset",
                                "没查过不等于数据集没有记录")
            man = json.loads((out.parent / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(man["parameters"]["source_kind"], "none")
            self.assertEqual(man["parameters"]["query_state"], "not_available")

    def test_summary_tier_is_labelled_as_the_screened_subset(self):
        """显式 --rows（09b summary）仍可用，但 query 必须说明它是筛选子集，不是完整数据集。"""
        with tempfile.TemporaryDirectory() as td:
            yard = self._yard(td)
            summary = pathlib.Path(td) / "summary.json"
            summary.write_text(json.dumps({"records": [
                {"record_id": "K.SG", "master_id": "K", "y_hg_raw": "N-CTS4714",
                 "locality": "Site", "latitude": 30.0, "longitude": 110.0}]}), encoding="utf-8")
            out = pathlib.Path(td) / "lineage_history.json"
            rc = lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1",
                          "--history", "", "--rows", str(summary)])
            self.assertEqual(rc, 0)
            h = json.loads(out.read_text(encoding="utf-8"))["y"]["history"]
            self.assertEqual(h["query"]["source_kind"], "summary")
            self.assertIn("screened", h["query"]["note"])
            self.assertEqual([o["record_id"] for o in h["observations"]], ["K.SG"])


class TestReviewBindingAndThresholds(unittest.TestCase):
    """AN4-a：人工复核绑定输入证据/当前树版本；solid/tail 不再从当前样本泛化；每样本复核与
    共享 panel 分离。旧代码上这些用例全部失败。"""

    REC = {"conservative_hg": "N-M1845", "reason": "resolution limit below N-M1845",
           "sample_id": "S1", "tree_version": "14.06.0", "evidence_sha256": "a" * 64}

    def _hist(self, **kw):
        rec = dict(self.REC); rec.update(kw)
        return lh.load_history({"schema_version": 1, "reviewed_calls": {"y:N-CTS4714": rec}})

    def test_review_recorded_under_a_different_tree_is_not_applied(self):
        """复核是对着某个树版本做的：树换了，结论不能原样照搬。"""
        said = []
        got = lh.reviewed_call(self._hist(), "y", "N-CTS4714", sample_id="S1",
                               tree_version="15.01.0", evidence_sha="a" * 64, log=said.append)
        self.assertEqual(got, (None, None))
        self.assertTrue(any("recorded under tree 14.06.0" in m for m in said), said)

    def test_review_without_a_recorded_tree_version_is_not_applied(self):
        said = []
        got = lh.reviewed_call(self._hist(tree_version=""), "y", "N-CTS4714", sample_id="S1",
                               tree_version="14.06.0", evidence_sha="a" * 64, log=said.append)
        self.assertEqual(got, (None, None), "没写树版本的记录没有绑定，不套用")
        self.assertTrue(any("unrecorded" in m for m in said), said)

    def test_caller_that_cannot_state_the_tree_cannot_use_the_review(self):
        """调用方不知道自己跑的是哪个树版本（tree_version=None）：无法核对绑定，不套用。"""
        got = lh.reviewed_call(self._hist(), "y", "N-CTS4714", sample_id="S1",
                               tree_version=None, evidence_sha="a" * 64)
        self.assertEqual(got, (None, None))

    def test_review_bound_to_different_evidence_is_not_applied(self):
        """复核核对的是某一份末端 SNP 表：输入证据换了（重新测序/换 BAM），结论不能照搬。"""
        said = []
        got = lh.reviewed_call(self._hist(), "y", "N-CTS4714", sample_id="S1",
                               tree_version="14.06.0", evidence_sha="b" * 64, log=said.append)
        self.assertEqual(got, (None, None))
        self.assertTrue(any("bound to different caller evidence" in m for m in said), said)
        said2 = []
        self.assertEqual(lh.reviewed_call(self._hist(evidence_sha256=""), "y", "N-CTS4714",
                                          sample_id="S1", tree_version="14.06.0",
                                          evidence_sha="a" * 64, log=said2.append), (None, None),
                         "记录没写证据指纹 = 没有绑定，不套用")
        self.assertEqual(lh.reviewed_call(self._hist(), "y", "N-CTS4714", sample_id="S1",
                                          tree_version="14.06.0", evidence_sha=None), (None, None),
                         "调用方没有证据指纹 = 无法核对，不套用")

    def test_load_review_reads_a_per_sample_file_and_tolerates_absence(self):
        """复核文件是样本私有的 work 产物：load_review 读它，文件不存在/坏 JSON 都是空复核，不阻断。"""
        with tempfile.TemporaryDirectory() as td:
            p = pathlib.Path(td) / "lineage_review.json"
            p.write_text(json.dumps({"reviewed_calls": {"y:N": dict(self.REC)}}), encoding="utf-8")
            self.assertIn("y:N", lh.load_review(p)["reviewed_calls"])
            self.assertEqual(lh.load_review(pathlib.Path(td) / "absent.json")["reviewed_calls"], {})
            bad = pathlib.Path(td) / "bad.json"
            bad.write_text("{not json", encoding="utf-8")
            self.assertEqual(lh.load_review(bad)["reviewed_calls"], {})

    def test_conservative_rule_has_no_sample_derived_defaults(self):
        """solid/tail 必须显式传入（来自本样本 config）：函数不藏 5/4 默认——那是从首个样本的
        12/1/0/1/5 路径归纳的，把它当默认等于把一个样本的分辨率极限写进所有样本。"""
        import inspect
        params = list(inspect.signature(lh.conservative_from_path).parameters.values())
        for p in params[1:]:                    # path 之后的 solid/tail 都不得有默认值
            self.assertIs(p.default, inspect.Parameter.empty,
                         f"{p.name} 必须显式配置，不设默认")
        path = [("N-M1845", 12, 0, 0), ("N-CTS4714", 5, 0, 0)]
        self.assertIsNone(lh.conservative_from_path(path, None, 4),
                          "未配置 solid：规则不生效，不硬造保守落点")
        self.assertIsNone(lh.conservative_from_path(path, 5, None),
                          "未配置 tail：规则不生效")

    def test_step05_reads_the_per_sample_review_not_the_shared_panel(self):
        """接线：05 读 work 目录的复核文件（配置可覆盖），不再读共享 panel；阈值来自 wgsconfig。"""
        src = (pathlib.Path(__file__).resolve().parents[1] / "scripts" / "05_y_haplogroup.py").read_text(
            encoding="utf-8")
        self.assertNotIn('"panel" /', src, "05 不得再读共享 panel")
        self.assertIn("lineage_review.json", src)
        self.assertIn("LINEAGE_REVIEW_FILE", src)
        self.assertIn("SOLID = LINEAGE_SOLID", src)
        self.assertIn("TAIL = LINEAGE_TAIL", src)
        self.assertNotIn("SOLID = 5", src, "样本归纳出的阈值不得硬编码")
        self.assertIn("tree_version=TREE_VERSION", src, "复核必须按当前树版本核对")
        self.assertIn("evidence_sha=_EVID_SHA", src, "复核必须按末端 SNP 表指纹核对")

    def test_shared_panel_carries_no_private_review(self):
        """共享 panel 的 reviewed_calls 必须为空：私人样本的复核住在样本自己的 work 目录。"""
        panel = json.loads((pathlib.Path(__file__).resolve().parents[1] / "panel" /
                            "lineage_history.json").read_text(encoding="utf-8"))
        self.assertEqual(panel.get("reviewed_calls"), {})
        self.assertTrue(str(panel.get("reviewed_calls_note", "")).strip(),
                       "panel 要说明复核去哪了，避免后人再把私人记录写回来")


class TestRunAllWiring(unittest.TestCase):
    """接线断言：run_all 用配置的 LINEAGE_HISTORY_FILE，观测来源给完整元数据/.anno，
    不再写死 panel 路径、不再只喂 09b 的 summary。"""

    REPO = pathlib.Path(__file__).resolve().parents[1]

    def _09d_block(self):
        text = (self.REPO / "run_all.sh").read_text(encoding="utf-8")
        i = text.find("step 09d")
        self.assertGreater(i, 0, "run_all.sh must still run step 09d")
        return text[i:text.find("\n", text.find("--sample", i))]

    def test_run_all_uses_config_and_full_metadata_not_hardcoded_paths(self):
        block = self._09d_block()
        self.assertIn('--history "$LINEAGE_HISTORY_FILE"', block)
        self.assertNotIn("panel/lineage_history.json", block, "history 面板路径来自配置，不写死")
        self.assertNotIn("--rows", block, "run_all 不再把 09b 筛选 summary 当观测来源")
        run_all = (self.REPO / "run_all.sh").read_text(encoding="utf-8")
        self.assertIn("--metadata $WGS/11_aadr/reference_metadata.tsv", run_all)
        self.assertIn("--anno $AADR_ANNO", run_all)

    def test_env_sh_exports_the_lineage_paths(self):
        env = (self.REPO / "scripts" / "env.sh").read_text(encoding="utf-8")
        self.assertIn('"LINEAGE_HISTORY_FILE"', env)
        self.assertIn('"AADR_ANNO"', env)

    def test_wgsconfig_falls_back_to_the_shipped_panel(self):
        """config 未覆盖时回到仓库自带 panel（通用面板，不是私人数据）；覆盖则用覆盖值。"""
        src = (self.REPO / "scripts" / "wgsconfig.py").read_text(encoding="utf-8")
        self.assertIn('OPT["lineage_history_file"] or str(ROOT / "panel" / "lineage_history.json")', src)


class TestCliWritesAdmissibleManifest(unittest.TestCase):
    """复审 H6/AN5 残留：09d 的 lineage_history.json 此前没有 manifest，30 只能按"文件存在"
    读取。现在 CLI 写 analysis_id=09d-lineage-history 的 manifest，30 走 analysis_state 准入。
    """

    def _yard(self, td):
        yard = pathlib.Path(td) / "03_haplo"
        yard.mkdir(parents=True)
        (yard / "y_result.json").write_text(json.dumps({
            "kind": "y", "state": "ok", "reported_hg": "N-CTS4714", "conservative_hg": "N-M1845",
            "supported_path": [{"node": "N-CTS4714", "der": 5, "anc": 0, "na": 0}]}), encoding="utf-8")
        (yard / "mt_result.json").write_text(json.dumps({
            "kind": "mt", "state": "ok", "reported_hg": "A13",
            "supported_path": [{"node": "A13", "der": 3, "anc": 0, "na": 0}]}), encoding="utf-8")
        return yard

    def test_cli_writes_manifest_that_analysis_state_admits(self):
        import ancestry_data as ad
        with tempfile.TemporaryDirectory() as td:
            yard = self._yard(td)
            out = pathlib.Path(td) / "lineage_history.json"
            rc = lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1",
                          "--history", ""])
            self.assertEqual(rc, 0)
            self.assertTrue(out.exists())
            man = json.loads((out.parent / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(man["analysis_id"], "09d-lineage-history")
            self.assertEqual(man["state"], "ok")
            self.assertEqual(man["sample_id"], "S1")
            self.assertEqual(man["outputs"], ["lineage_history.json"])
            # 与 30 相同的准入调用：state+sample_id+analysis_id 全部一致才交出 doc
            state, reason, doc = ad.analysis_state(
                out.parent, {"sample_id": "S1", "analysis_id": "09d-lineage-history"},
                names=("lineage_history.json",))
            self.assertEqual(state, "ok", reason)
            self.assertEqual(doc["y"]["reported_hg"], "N-CTS4714")
            # 换样本的旧 manifest 必须被判 stale——存在性读取要修掉的正是这个
            other = ad.analysis_state(out.parent, {"sample_id": "OTHER", "analysis_id": "09d-lineage-history"},
                                       names=("lineage_history.json",))
            self.assertEqual(other[0], "unavailable")
            self.assertEqual(other[1], "stale_result")

    def test_manifest_binds_history_content_not_path(self):
        """复审 §3.2 P0 指纹绑定：09d 的 manifest 记 history 文件的内容 sha（不是路径）；
        换面板证据内容后 30 的准入形状判 stale，还原后重新 ok（A→B→A）。"""
        import ancestry_data as ad
        with tempfile.TemporaryDirectory() as td:
            yard = self._yard(td)
            hist = pathlib.Path(td) / "panel_history.json"
            hist.write_text(json.dumps({"parents": {}}), encoding="utf-8")
            out = pathlib.Path(td) / "lineage_history.json"
            rc = lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1",
                          "--history", str(hist)])
            self.assertEqual(rc, 0)
            man = json.loads((out.parent / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(man["parameters"]["history_sha"], ad.file_sha(hist))
            admit = lambda: ad.analysis_state(
                out.parent, {"sample_id": "S1", "analysis_id": "09d-lineage-history"},
                names=("lineage_history.json",),
                expected_parameters={"history_sha": ad.file_sha(hist)})
            self.assertEqual(admit()[0], "ok")
            orig = hist.read_bytes()
            hist.write_text(json.dumps({"parents": {}, "v": 2}), encoding="utf-8")
            self.assertEqual(admit()[:2], ("unavailable", "stale_result"))
            hist.write_bytes(orig)
            self.assertEqual(admit()[0], "ok")

    def test_midrun_crash_invalidates_previous_ok(self):
        """复审 §3.2 P0 失败生命周期：09d 开工先失效旧 ok。

        成功一轮后再跑一轮、--history 指向损坏 JSON（未处理崩溃路径）：开工失效记录必须已经
        写下，旧 ok manifest 与旧 lineage_history.json 不得继续充当本次结果。"""
        import ancestry_data as ad
        with tempfile.TemporaryDirectory() as td:
            yard = self._yard(td)
            out = pathlib.Path(td) / "lineage_history.json"
            rc = lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1",
                          "--history", ""])
            self.assertEqual(rc, 0)
            state, _, doc = ad.analysis_state(
                out.parent, {"sample_id": "S1", "analysis_id": "09d-lineage-history"},
                names=("lineage_history.json",))
            self.assertEqual(state, "ok")
            self.assertIsNotNone(doc)
            # 第二轮：history 文件损坏 → json.loads 在失效记录之后崩溃
            bad = pathlib.Path(td) / "bad_history.json"
            bad.write_text("{not json", encoding="utf-8")
            with self.assertRaises(json.JSONDecodeError):
                lh._cli(["--yard", str(yard), "--out", str(out), "--sample", "S1",
                         "--history", str(bad)])
            state, reason, doc = ad.analysis_state(
                out.parent, {"sample_id": "S1", "analysis_id": "09d-lineage-history"},
                names=("lineage_history.json",))
            self.assertEqual((state, reason), ("unavailable", "run_begun_not_published"))
            self.assertIsNone(doc, "开工失效后，上一轮的 lineage_history.json 不得继续交付")


@unittest.skipIf(_SKIP05, _SKIP05)
class Test05EmptyPathProduction(unittest.TestCase):
    """复审 §3.2 P0（AN4/H4）：05 在空 Y 路径下曾先读 path[-1]——IndexError、exit 1、不写
    y_result.json，run_all 停在这一步。三种生产路径：无 Y 证据 / 全祖先 / 有效 Y。
    子进程跑真实 05；pileup 预先落盘，绕开 samtools/CRAM。"""

    SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "05_y_haplogroup.py"

    def _run05(self, td, pileup_text, preexisting_ok=False):
        import os
        import subprocess
        W = pathlib.Path(td) / "work" / "wgs" / "03_haplo"
        W.mkdir(parents=True)
        yt = pathlib.Path(td) / "ytree"
        yt.mkdir()
        (yt / "snps_hg19.csv").write_text(
            "Name,start,allele_anc,allele_der\nSNP1,1000,A,G\n", encoding="utf-8")
        (yt / "current_tree.json").write_text(json.dumps(
            {"id": "root", "snps": "", "children": [{"id": "N-TEST1", "snps": "SNP1"}]}),
            encoding="utf-8")
        cv = pathlib.Path(td) / "work" / "data" / "ref" / "ytree"
        cv.mkdir(parents=True)
        (cv / "current_version.txt").write_text("14.06.0\n", encoding="utf-8")
        (W / "y_pileup.tsv").write_text(pileup_text, encoding="utf-8")
        if preexisting_ok:   # 上一轮的成功必须被本次的 unavailable 覆盖
            (W / "y_result.json").write_text(json.dumps(
                {"kind": "y", "state": "ok", "reported_hg": "N-OLD"}), encoding="utf-8")
            (W / "y_terminal_snps.tsv").write_text("branch\tsnp\nN-OLD\tX\n", encoding="utf-8")
        cfg = pathlib.Path(td) / "config.yaml"
        cfg.write_text("".join([
            "sample_id: TESTSAMPLE\n",
            f"work_dir: {json.dumps(str(pathlib.Path(td) / 'work'))}\n",
            f"ytree_dir: {json.dumps(str(yt))}\n"]), encoding="utf-8")
        r = subprocess.run([sys.executable, str(self.SCRIPT)], capture_output=True, text=True,
                           env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=120)
        return r, W

    def test_no_y_evidence_and_all_ancestral_deliver_unavailable(self):
        """无 pileup 行 / 全祖先：exit 0、y_result.json=unavailable（no_supported_path），
        上一轮的 ok 结果与旧末端 SNP 表被本次覆盖，不残留 N-OLD。"""
        for label, pileup in (("no pileup rows", ""), ("all ancestral", "Y\t1000\tA\t5\t.....\tIIIII\n")):
            with tempfile.TemporaryDirectory() as td:
                r, W = self._run05(td, pileup, preexisting_ok=True)
                self.assertEqual(r.returncode, 0, f"{label}: {r.stderr[-400:]}")
                res = json.loads((W / "y_result.json").read_text(encoding="utf-8"))
                self.assertEqual(res["state"], "unavailable", label)
                self.assertEqual(res["reason_code"], "no_supported_path", label)
                self.assertNotEqual(res.get("reported_hg"), "N-OLD", f"{label}: 旧成功不得残留")
                self.assertEqual((W / "y_terminal_snps.tsv").read_text().splitlines(),
                                 ["branch\tsnp\tpos_hg19\tanc\tder\tstate\tdepth\tn_anc\tn_der"],
                                 f"{label}: 旧末端 SNP 表必须清空")
                self.assertIn("no supported branch", (W / "y_haplogroup_yfull.txt").read_text())

    def test_valid_y_still_works(self):
        """有效 Y：路径非空，正常 ok——空路径处理不得破坏既有生产路径。"""
        with tempfile.TemporaryDirectory() as td:
            r, W = self._run05(td, "Y\t1000\tA\t5\tGGGGG\tIIIII\n")
            self.assertEqual(r.returncode, 0, r.stderr)
            res = json.loads((W / "y_result.json").read_text(encoding="utf-8"))
            self.assertEqual(res["state"], "ok")
            self.assertEqual(res["reported_hg"], "N-TEST1")
