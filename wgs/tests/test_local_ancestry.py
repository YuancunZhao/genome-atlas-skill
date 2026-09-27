"""AN3 的行为测试：校准 holdout 选择与局部祖源的通用化。

约束：纯函数测试，不读真实参考、不跑 FLARE。真实的 VCF 子集与 FLARE 集成由 16b 在服务器上验证。
"""
import pathlib, sys, unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import ancestry_data as ad  # noqa: E402


class TestPickHoldout(unittest.TestCase):
    """按固定 seed 无放回抽样：人数 = min(配置 n, floor(群体人数/5))，样本不足则校准不可用。"""

    def _pops(self, n_chb=100, n_chs=40, n_ceu=3):
        rows = []
        for i in range(n_chb):
            rows.append({"iid": f"CHB{i:04d}", "population": "CHB"})
        for i in range(n_chs):
            rows.append({"iid": f"CHS{i:04d}", "population": "CHS"})
        for i in range(n_ceu):
            rows.append({"iid": f"CEU{i:04d}", "population": "CEU"})
        return rows

    def test_counts_are_min_of_config_and_a_fifth_of_the_group(self):
        out = ad.pick_holdout(self._pops(), pops=["CHB", "CHS"], n=20, seed=1)
        self.assertEqual(out["state"], "ok")
        self.assertEqual(out["per_pop"]["CHB"], 20, "100 人 → min(20, 20)")
        self.assertEqual(out["per_pop"]["CHS"], 8, "40 人 → min(20, floor(40/5))")
        self.assertEqual(len(out["ids"]), 28)

    def test_same_seed_same_choice_and_a_different_seed_differs(self):
        a = ad.pick_holdout(self._pops(), pops=["CHB"], n=20, seed=1)
        b = ad.pick_holdout(self._pops(), pops=["CHB"], n=20, seed=1)
        c = ad.pick_holdout(self._pops(), pops=["CHB"], n=20, seed=2)
        self.assertEqual(a["ids"], b["ids"], "同 seed 必须可复现")
        self.assertNotEqual(a["ids"], c["ids"], "换 seed 应换样本")

    def test_small_groups_make_calibration_unavailable_not_silently_small(self):
        out = ad.pick_holdout(self._pops(n_chb=100, n_chs=4, n_ceu=3), pops=["CHB", "CHS"], n=20, seed=1)
        # CHS 4 人 → floor(4/5)=0，不足 2 人 holdout：校准不可用，而不是拿 0 个人硬算
        self.assertEqual(out["state"], "unavailable")
        self.assertEqual(out["reason_code"], "insufficient_holdout")
        self.assertIn("CHS", out["detail"])
        # 剩余参考少于 2 人也一样
        out2 = ad.pick_holdout(self._pops(n_chb=100, n_chs=40), pops=["CHB", "CHS"], n=20, seed=1)
        self.assertEqual(out2["state"], "ok")
        self.assertNotIn("CEU", out2["per_pop"], "未配置的群体不得出现在校准里")

    def test_unknown_population_is_reported(self):
        out = ad.pick_holdout(self._pops(), pops=["CHB", "PEL"], n=20, seed=1)
        self.assertEqual(out["state"], "unavailable")
        self.assertEqual(out["reason_code"], "unknown_population")
        self.assertIn("PEL", out["detail"])

    def test_holdout_ids_are_disjoint_from_the_remaining_reference(self):
        rows = self._pops()
        out = ad.pick_holdout(rows, pops=["CHB", "CHS"], n=20, seed=1)
        remaining = [r["iid"] for r in rows if r["iid"] not in set(out["ids"])]
        self.assertTrue(set(out["ids"]).isdisjoint(set(remaining)))
        self.assertEqual(len(remaining) + len(out["ids"]), len(rows))

    def test_no_duplicates_and_original_order_is_irrelevant(self):
        rows = self._pops()
        a = ad.pick_holdout(rows, pops=["CHB"], n=20, seed=7)
        b = ad.pick_holdout(list(reversed(rows)), pops=["CHB"], n=20, seed=7)
        self.assertEqual(len(set(a["ids"])), len(a["ids"]))
        self.assertEqual(sorted(a["ids"]), sorted(b["ids"]), "抽样只依赖 iid 集合与 seed")
