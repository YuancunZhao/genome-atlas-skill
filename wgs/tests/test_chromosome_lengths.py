"""染色体长度必须来自参考基因组的索引，而不是代码里的内置表。

内置表只对 hg19/GRCh37 正确。换参考版本后，用旧长度去算百分比会给出偏大的数字——图看着合理，
含义全错。这些断言同时钉住"回退时会告警"，避免它悄悄退回去。
"""
import pathlib, sys, tempfile, unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import ancestry_data as ad  # noqa: E402


class TestChromosomeLengths(unittest.TestCase):
    def _fai(self, rows):
        d = tempfile.mkdtemp()
        p = pathlib.Path(d) / "ref.fasta.fai"
        p.write_text("".join(f"{n}\t{l}\t0\t60\t61\n" for n, l in rows), encoding="utf-8")
        return str(p)

    def test_reads_the_index_and_strips_the_chr_prefix(self):
        p = self._fai([("chr1", 248956422), ("chr22", 50818468)])
        got = ad.chromosome_lengths([p], fallback={"1": 249250621}, log=lambda *_: None)
        self.assertEqual(got["1"], 248956422)
        self.assertEqual(got["22"], 50818468)
        self.assertNotIn("chr1", got)

    def test_prefers_the_first_readable_candidate(self):
        a = self._fai([("1", 111)])
        b = self._fai([("1", 222)])
        got = ad.chromosome_lengths(["/nonexistent.fai", a, b], fallback={"1": 999}, log=lambda *_: None)
        self.assertEqual(got["1"], 111)

    def test_scaffolds_are_excluded(self):
        p = self._fai([("1", 100), ("GL000191.1", 106433), ("KI270728.1", 2750178), ("X", 156040895)])
        got = ad.chromosome_lengths([p], fallback=None, log=lambda *_: None)
        self.assertEqual(sorted(got), ["1", "X"])

    def test_fallback_is_announced_not_silent(self):
        """没有索引时必须开口：静默回退到 hg19 表正是换版本后出错的方式。"""
        said = []
        got = ad.chromosome_lengths([], fallback={"1": 249250621}, log=said.append)
        self.assertEqual(got, {"1": 249250621})
        self.assertTrue(any("built-in hg19" in m and "WRONG" in m for m in said), said)

    def test_blank_candidates_are_skipped_quietly(self):
        """配置里没有这一项时传进来的是空串；Path('') 等于当前目录，会冒出一句假的读不了告警。"""
        said = []
        p = self._fai([("1", 100)])
        got = ad.chromosome_lengths(["", p], fallback=None, log=said.append)
        self.assertEqual(got, {"1": 100})
        self.assertEqual([m for m in said if "could not read" in m], [], said)

    def test_malformed_lines_do_not_raise(self):
        p = self._fai([("1", 100)])
        with open(p, "a", encoding="utf-8") as fh:
            fh.write("garbage\n\n2\tnotanumber\n")
        got = ad.chromosome_lengths([p], fallback=None, log=lambda *_: None)
        self.assertEqual(got, {"1": 100})


if __name__ == "__main__":
    unittest.main(verbosity=2)
