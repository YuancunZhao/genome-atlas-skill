"""链外产物必须能被点名、被检查，缺了要报错而不是变成空结果（H5）。

这些产物在服务器上是手工跑出来的：仓库里没有产生它们的脚本。风险不在于"手工"本身，而在于
**缺失时下游会静默降级**——一张空表、一个悄悄少掉的小节，看起来和"这个样本没有该发现"一样。
"""
import pathlib, sys, tempfile, unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import _artifacts as art  # noqa: E402


class TestArtifactContract(unittest.TestCase):
    def _root(self, present=(), empty=()):
        d = pathlib.Path(tempfile.mkdtemp())
        for key in list(present) + list(empty):
            p = d / art.ARTIFACTS[key]["path"]
            if key in empty:
                p.mkdir(parents=True, exist_ok=True)      # 存在但为空：也算没有
            else:
                p.mkdir(parents=True, exist_ok=True)
                (p / "x.txt").write_text("data", encoding="utf-8")
        return d

    def test_all_products_are_named_with_a_consequence(self):
        for key, spec in art.ARTIFACTS.items():
            self.assertTrue(spec["path"], key)
            self.assertIn("consumers", spec, key)
            self.assertTrue(spec["consumers"], key)
            self.assertTrue(spec["note"], key)

    def test_manual_products_say_they_have_no_producer(self):
        """不编造生产者命令：猜一个 Delly 调用去替代实际用的那个，会让数字在无人被告知的情况下变化。"""
        manual = [k for k, v in art.ARTIFACTS.items() if v["produced_by_hand"]]
        self.assertTrue(manual)
        for k in manual:
            self.assertIsNone(art.ARTIFACTS[k]["producer"], k)

    def test_check_reports_missing_rather_than_passing(self):
        root = self._root(present=["delly_sv"])
        ok, rep = art.check(["delly_sv", "smn"], root=root)
        self.assertFalse(ok)
        self.assertEqual(rep["missing"], ["smn"])

    def test_empty_directory_counts_as_missing(self):
        """一个空目录不该被当成"有产物"：那正是静默降级的样子。"""
        root = self._root(empty=["smn"])
        ok, rep = art.check(["smn"], root=root)
        self.assertFalse(ok)
        self.assertIn("smn", rep["missing"])

    def test_require_raises_with_the_consequence_and_the_gap(self):
        root = self._root(present=[])
        with self.assertRaises(FileNotFoundError) as cm:
            art.require(["smn"], root=root)
        msg = str(cm.exception)
        self.assertIn("smn", msg)
        self.assertIn("SMN1/SMN2", msg)                    # 说清缺了会怎样
        self.assertIn("run by hand on the server", msg)     # 说清为什么链里没有

    def test_require_passes_when_present(self):
        root = self._root(present=["cyrius"])
        rep = art.require(["cyrius"], root=root)
        self.assertEqual(rep["missing"], [])

    def test_unknown_product_key_is_reported_not_ignored(self):
        ok, rep = art.check(["not_a_product"], root=self._root())
        self.assertFalse(ok)
        self.assertIn("not_a_product", rep["missing"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
