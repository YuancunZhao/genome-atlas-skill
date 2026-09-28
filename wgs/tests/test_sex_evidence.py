"""性染色体构成必须由**数据**判定，而不是由配置声明（H4）。

配置说的是样本被声明成什么，它说不出这份数据里到底测到了什么。信 `sex: male` 就对女性样本算 Y
统计，会得到一堆看着正常、实际没有意义的数字。判不了的时候要说判不了，不能默认成"没有 Y"。
"""
import pathlib, sys, tempfile, unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import _sex_evidence as se  # noqa: E402


def _write_summary(rows, regions=True):
    d = pathlib.Path(tempfile.mkdtemp())
    p = d / "depth.mosdepth.summary.txt"
    body = ["chrom\tlength\tbases\tmean\tmin\tmax"]
    for name, mean in rows:
        body.append(f"{name}\t1000\t100000\t{mean}\t0\t200")
        if regions:
            body.append(f"{name}_region\t1000\t100000\t{mean}\t0\t200")
    p.write_text("\n".join(body) + "\n", encoding="utf-8")
    return p


AUTOS = [(str(i), 50.0) for i in range(1, 23)]


class TestSexEvidence(unittest.TestCase):
    def test_male_shape_is_detected(self):
        """X 半量 + Y 约半数深度 = 有 Y。真实样本是 0.552 / 0.188，这里用同形状的数据。"""
        p = _write_summary(AUTOS + [("X", 26.65), ("Y", 9.1)])
        means, reason = se._depth_means(p)
        self.assertEqual(reason, "")
        rec = se.classify_depths(means, "male")
        self.assertEqual(rec["inferred"], "has_y")
        self.assertTrue(rec["agrees_with_declared"])
        self.assertAlmostEqual(rec["evidence"]["x_ratio"], 0.533, places=2)

    def test_female_shape_is_detected_and_conflict_is_reported(self):
        p = _write_summary(AUTOS + [("X", 52.0), ("Y", 0.05)])
        means, _ = se._depth_means(p)
        rec = se.classify_depths(means, "female")
        self.assertEqual(rec["inferred"], "no_y")
        # 同一份数据配上一个矛盾的声明：必须报冲突，而不是迁就声明
        conflict = se.classify_depths(means, "male")
        self.assertFalse(conflict["agrees_with_declared"])

    def test_region_rows_are_ignored(self):
        p = _write_summary(AUTOS + [("X", 50.0), ("Y", 25.0)], regions=True)
        means, _ = se._depth_means(p)
        self.assertEqual(len(means), 24)          # 22 + X + Y，_region 行已丢弃

    def test_chr_prefix_is_normalised(self):
        p = _write_summary([(f"chr{i}", 50.0) for i in range(1, 23)] + [("chrX", 25.0), ("chrY", 24.0)])
        means, _ = se._depth_means(p)
        self.assertIn("Y", means)
        self.assertEqual(se.classify_depths(means, "male")["inferred"], "has_y")

    def test_missing_summary_says_so(self):
        rec = se.classify_depths(None, "male", "no_depth_summary")
        self.assertEqual(rec["state"], "unavailable")
        self.assertEqual(rec["reason_code"], "no_depth_summary")
        self.assertIsNone(rec["evidence"])

    def test_summary_without_sex_chromosomes_is_unavailable(self):
        p = _write_summary(AUTOS)
        means, _ = se._depth_means(p)
        rec = se.classify_depths(means, "male")
        self.assertEqual(rec["state"], "unavailable")
        self.assertEqual(rec["reason_code"], "no_sex_chromosomes_in_summary")

    def test_unreadable_and_empty_summaries_are_reported(self):
        d = pathlib.Path(tempfile.mkdtemp())
        empty = d / "e.txt"
        empty.write_text("chrom\tlength\tbases\tmean\tmin\tmax\n", encoding="utf-8")
        self.assertEqual(se._depth_means(empty)[1], "depth_summary_empty")
        self.assertEqual(se._depth_means(d / "nope.txt")[1], "no_depth_summary")

    def test_unexpected_columns_are_rejected_not_guessed(self):
        """列名变了要明确拒绝：按位置猜"第 4 列是 mean"会在格式变化后给出错数字。"""
        d = pathlib.Path(tempfile.mkdtemp())
        bad = d / "b.txt"
        bad.write_text("chrom\tlength\tbases\tdepth\tmin\tmax\n1\t1000\t50000\t50\t0\t100\n", encoding="utf-8")
        self.assertEqual(se._depth_means(bad)[1], "depth_summary_unexpected_columns")


class TestEvidenceBoundaries(unittest.TestCase):
    """H4 复审补的三条边界：低深度 Y 不断言、X 与 Y 矛盾不替人下核型结论、判不了要留话。"""

    def test_shallow_y_depth_is_not_asserted_as_present(self):
        """比值达标但 Y 绝对深度不足：2× 文库的噪声就能把 0.15× 推成 0.075 的"达标"比值，
        这种"有 Y"降级为 ambiguous，由 evidence_sex 走声明回退。"""
        autos = [(str(i), 2.0) for i in range(1, 23)]      # 常染色体 2× 的浅文库
        p = _write_summary(autos + [("X", 1.1), ("Y", 0.15)])   # yr=0.075 >= 0.05，但 Y 只有 0.15×
        means, _ = se._depth_means(p)
        rec = se.classify_depths(means, "male")
        self.assertEqual(rec["inferred"], "ambiguous")
        self.assertIn("absolute floor", rec["why"])
        logs = []
        sex, rec2 = se.evidence_sex(declared="male", summary=p, log=logs.append)
        self.assertEqual(sex, "male")                       # 回退声明，不是把 ambiguous 当 male
        self.assertTrue(any("falling back" in m for m in logs), logs)

    def test_diploid_x_with_y_evidence_is_a_conflict_not_a_call(self):
        """有 Y 而 X 深度像两条：证据自相矛盾。保留两侧证据并标记 x_conflict，
        evidence_sex 回退声明——Y 一条证据不能当核型结论发布。"""
        p = _write_summary(AUTOS + [("X", 48.0), ("Y", 9.1)])    # xr=0.96 >= 0.90，yr 达标
        means, _ = se._depth_means(p)
        rec = se.classify_depths(means, "male")
        self.assertEqual(rec["inferred"], "has_y")          # Y 侧读数保留
        self.assertTrue(rec["x_conflict"])
        self.assertIn("conflicting X evidence", rec["why"])
        logs = []
        sex, rec2 = se.evidence_sex(declared="female", summary=p, log=logs.append)
        self.assertEqual(sex, "female")                     # 回退声明，而不是由 Y 定 male
        self.assertTrue(any("X and Y evidence conflict" in m for m in logs), logs)

    def test_haploid_x_without_y_evidence_is_a_conflict_too(self):
        """无 Y 而 X 深度像一条：同样的矛盾，同样回退。X_SINGLE_MAX 原先定义了却从未参与判定。"""
        p = _write_summary(AUTOS + [("X", 26.0), ("Y", 0.05)])   # xr=0.52 <= 0.75，yr 达标无 Y
        means, _ = se._depth_means(p)
        rec = se.classify_depths(means, "female")
        self.assertEqual(rec["inferred"], "no_y")
        self.assertTrue(rec["x_conflict"])
        logs = []
        sex, _ = se.evidence_sex(declared="male", summary=p, log=logs.append)
        self.assertEqual(sex, "male")
        self.assertTrue(any("X and Y evidence conflict" in m for m in logs), logs)

    def test_normal_shapes_carry_no_conflict_flag(self):
        """常规男性/女性形状不能被新门槛误伤：夹具即真实样本的形状（0.533 / 1.04）。"""
        male = _write_summary(AUTOS + [("X", 26.65), ("Y", 9.1)])
        female = _write_summary(AUTOS + [("X", 52.0), ("Y", 0.05)])
        for p in (male, female):
            means, _ = se._depth_means(p)
            rec = se.classify_depths(means, None)
            self.assertFalse(rec["x_conflict"])
            self.assertNotIn("absolute floor", rec["why"])


class TestEffectiveSex(unittest.TestCase):
    """下游怎么得到"该按哪种性别处理"：证据优先，回退可以但必须留话。"""

    def test_evidence_wins_over_the_declaration(self):
        p = _write_summary(AUTOS + [("X", 52.0), ("Y", 0.05)])   # 数据说无 Y
        logs = []
        sex, rec = se.evidence_sex(declared="male", summary=p, log=logs.append)
        self.assertEqual(sex, "female", "证据必须压过声明")
        self.assertFalse(rec["agrees_with_declared"])
        self.assertTrue(any("disagrees" in m for m in logs), logs)

    def test_male_shape_yields_male(self):
        p = _write_summary(AUTOS + [("X", 26.65), ("Y", 9.1)])
        sex, rec = se.evidence_sex(declared="male", summary=p, log=lambda *_: None)
        self.assertEqual(sex, "male")
        self.assertTrue(rec["agrees_with_declared"])

    def test_unavailable_evidence_falls_back_loudly(self):
        """回退到声明值是可以的（否则缺深度汇总的环境会连带停掉 Y 分析），但必须说明。"""
        logs = []
        sex, rec = se.evidence_sex(declared="female", summary="/nonexistent.txt", log=logs.append)
        self.assertEqual(sex, "female")
        self.assertEqual(rec["state"], "unavailable")
        self.assertTrue(any("unavailable" in m and "falling back" in m for m in logs), logs)

    def test_used_sex_is_recorded_on_the_result(self):
        p = _write_summary(AUTOS + [("X", 26.65), ("Y", 9.1)])
        _, rec = se.evidence_sex(declared="male", summary=p, log=lambda *_: None)
        self.assertEqual(rec["used_sex"], "male")


if __name__ == "__main__":
    unittest.main(verbosity=2)
