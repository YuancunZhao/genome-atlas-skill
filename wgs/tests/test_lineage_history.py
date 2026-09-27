"""AN4 的行为测试：父母系的树关系匹配与历史证据组装。

纯函数测试：不读真实 YFull/PhyloTree 树、不跑 haplogrep3、不改报告数据。
"""
import pathlib, sys, unittest

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

    def test_observations_carry_the_required_fields(self):
        o = lh.lineage_observations(self.ROWS, "mt:A13", parents={}, tree_kind="mt")[0]
        for k in ("record_id", "node_id", "relation", "locality", "coordinates", "precision",
                  "date_range", "date_basis", "call_source", "publication"):
            self.assertIn(k, o)

    def test_reviewed_call_beats_any_automatic_rule(self):
        hist = lh.load_history({"schema_version": 1, "reviewed_calls": {
            "y:N-CTS4714": {"conservative_hg": "N-M1845", "reason": "1-5 sites per level below"}}})
        self.assertEqual(lh.reviewed_call(hist, "y", "N-CTS4714"), ("N-M1845", "1-5 sites per level below"))
        self.assertEqual(lh.reviewed_call(hist, "y", "N-CTS4714")[0], lh.reviewed_call(hist, "y", "y:N-CTS4714")[0])
        self.assertEqual(lh.reviewed_call(hist, "mt", "A13"), (None, None), "别的支系不受影响")
        self.assertEqual(lh.reviewed_call({}, "y", "N-CTS4714"), (None, None))

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
