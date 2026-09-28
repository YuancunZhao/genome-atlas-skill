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


if __name__ == "__main__":
    unittest.main(verbosity=2)
