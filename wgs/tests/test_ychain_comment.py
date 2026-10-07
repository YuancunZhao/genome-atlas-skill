"""旧缺陷 #18 收尾回归：Y 链卡片注释不得再点名示例样本的支系。

上游模板的 ychain 注释写死 "YFull, from O-F438 down"——那是示例样本（O 系）的路径。代码自
3e212d1/d6069b1 起已按 D.ypath 通用行走，32_check 也在成品层拦 O-M122/O-F438 字面量，但注释
仍在向维护者宣称这段是"O-F438 往下"的专属逻辑。注释即契约：换成按树根通用行走的描述，
并用源级断言锁住（成品模板是唯一生产路径，node 不需要）。
"""
import pathlib, unittest

JS = pathlib.Path(__file__).resolve().parents[1] / "templates" / "report_script.js"


class TestYChainComment(unittest.TestCase):
    def test_section_comment_is_haplogroup_agnostic(self):
        src = JS.read_text(encoding="utf-8")
        m = [l for l in src.splitlines() if "01 Y chain" in l]
        self.assertTrue(m, "the ychain section header comment must exist")
        # 任何"从某支系往下"的样本特定措辞都说明这段被当成专属逻辑——路径由 05 从树根走出
        self.assertNotIn("from O-", src, "comment must not name the example sample's clade")
        self.assertIn("tree root", m[0])


if __name__ == "__main__":
    unittest.main()
