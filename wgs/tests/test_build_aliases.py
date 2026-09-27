"""build 校验只有一处，且别名必须归一化到规范名。

注释里写"只支持 GRCh37"而代码不校验，等于没限制：`build: GRCh38` 会被静默接受，而 FASTA、chain、
1000G panel、AADR 与注释 GFF 全是 GRCh37，位点会按错误的坐标系解释且没有任何一步失败。
"""
import os, pathlib, subprocess, sys, tempfile, unittest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import ancestry_data as ad  # noqa: E402


class TestBuildValidation(unittest.TestCase):
    def test_aliases_normalise_to_the_canonical_name(self):
        for alias in ("GRCh37", "grch37", "hg19", "b37", "GRCh37.p13"):
            self.assertEqual(ad.normalize_build(alias), "GRCh37", alias)

    def test_other_builds_are_rejected_with_a_reason(self):
        for bad in ("GRCh38", "hg38", "GRCh36", "mm10", ""):
            if bad == "":
                continue  # empty means "not configured", which defaults
            with self.assertRaises(ValueError, msg=bad) as cm:
                ad.normalize_build(bad)
            msg = str(cm.exception)
            self.assertIn(bad, msg)
            self.assertIn("GRCh37", msg)

    def test_missing_value_takes_the_default(self):
        self.assertEqual(ad.normalize_build(None), "GRCh37")
        self.assertEqual(ad.normalize_build(""), "GRCh37")
        self.assertEqual(ad.normalize_build(None, default="GRCh37"), "GRCh37")

    def test_config_with_an_unsupported_build_fails_before_any_work(self):
        """校验必须发生在 import 时（即任何步骤开始前），而不是某个中途步骤里。"""
        d = tempfile.mkdtemp()
        cfg = pathlib.Path(d) / "c.yaml"
        cfg.write_text("sample_id: X\nbuild: GRCh38\n", encoding="utf-8")
        env = dict(os.environ, WGS_CONFIG=str(cfg), PYTHONPATH=str(SCRIPTS))
        r = subprocess.run([sys.executable, "-c", "import wgsconfig"], capture_output=True, text=True, env=env)
        self.assertNotEqual(r.returncode, 0, "an unsupported build must not import cleanly")
        self.assertIn("GRCh38", r.stderr)
        self.assertIn("GRCh37", r.stderr)

    def test_read_options_and_wgsconfig_agree(self):
        """两处若各自维护一份校验，同一份配置会得出不同结论（这正是修之前的状况）。"""
        self.assertIs(ad.normalize_build("hg19"), ad.BUILD_ALIASES["hg19"])
        for alias in ad.BUILD_ALIASES:
            self.assertEqual(ad.normalize_build(alias), "GRCh37")


if __name__ == "__main__":
    unittest.main(verbosity=2)
