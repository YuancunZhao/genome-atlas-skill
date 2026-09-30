"""AN4 的行为测试：父母系的树关系匹配与历史证据组装。

纯函数测试：不读真实 YFull/PhyloTree 树、不跑 haplogrep3、不改报告数据。
"""
import json, pathlib, sys, tempfile, unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import lineage_history as lh  # noqa: E402


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
        """复核记录优先于自动规则——但只在**它属于这个样本**时（复审 AN4）。

        旧版本不带样本就返回结论，那会让一个样本的人工判断顺着支系名传播给所有同支系样本。
        """
        hist = lh.load_history({"schema_version": 1, "reviewed_calls": {
            "y:N-CTS4714": {"conservative_hg": "N-M1845", "reason": "1-5 sites per level below",
                            "sample_id": "S1"}}})
        self.assertEqual(lh.reviewed_call(hist, "y", "N-CTS4714", sample_id="S1"),
                         ("N-M1845", "1-5 sites per level below"))
        self.assertEqual(lh.reviewed_call(hist, "y", "N-CTS4714", sample_id="S1")[0],
                         lh.reviewed_call(hist, "y", "y:N-CTS4714", sample_id="S1")[0],
                         "带/不带 kind 前缀的键等价")
        self.assertEqual(lh.reviewed_call(hist, "mt", "A13", sample_id="S1"), (None, None), "别的支系不受影响")
        self.assertEqual(lh.reviewed_call({}, "y", "N-CTS4714", sample_id="S1"), (None, None))

    def test_a_review_belongs_to_one_sample_not_to_the_branch(self):
        """同支系的**另一个**样本不得继承本样本的复核结论——这是 AN4 的核心。"""
        hist = lh.load_history({"schema_version": 1, "reviewed_calls": {
            "y:N-CTS4714": {"conservative_hg": "N-M1845", "reason": "reviewed for S1", "sample_id": "S1"}}})
        said = []
        got = lh.reviewed_call(hist, "y", "N-CTS4714", sample_id="S2", log=said.append)
        self.assertEqual(got, (None, None), "S2 调出同一支系也不得套用 S1 的结论")
        self.assertTrue(any("belongs to S1" in m for m in said), said)
        self.assertEqual(lh.reviewed_call(hist, "y", "N-CTS4714", sample_id=None), (None, None),
                         "未指明样本时不套用")

    def test_an_unmarked_review_is_not_a_generic_rule(self):
        """没有 sample_id 的旧记录不得被当成通用规则——它无法与"通用规则"区分，而正是这种混同要修。"""
        hist = lh.load_history({"schema_version": 1, "reviewed_calls": {
            "y:N-CTS4714": {"conservative_hg": "N-M1845", "reason": "no owner recorded"}}})
        said = []
        self.assertEqual(lh.reviewed_call(hist, "y", "N-CTS4714", sample_id="S1", log=said.append), (None, None))
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
