"""Unit tests for the pure f3 math in 28_f3_stats.py (jackknife, block coding, allele-checked
frequency merge, matrix assembly). Importing the module must not touch wgsconfig: that creates
the work directories, a side effect tests cannot have.

numpy/pandas live on the analysis server (the local suite is stdlib-only, like 30's), so the
tests skip themselves when those imports are unavailable."""
import importlib.util
import pathlib
import unittest

try:
    import numpy as np
    import pandas as pd
    _spec = importlib.util.spec_from_file_location(
        "f3_stats", pathlib.Path(__file__).resolve().parents[1] / "scripts" / "28_f3_stats.py")
    f3 = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(f3)
except ImportError:            # local runner: everything below needs the analysis environment
    np = pd = f3 = None

_NEEDS_NP = unittest.skipIf(f3 is None, "numpy/pandas not installed here (runs on the server)")


def _afreq(rows):
    """rows: (ID, CHROM, POS, REF, ALT, freq) -> plink2 .afreq-shaped frame."""
    return pd.DataFrame(rows, columns=["ID", "CHROM", "POS", "REF", "ALT", "ALT_FREQS"])


@_NEEDS_NP
class TestJackknife(unittest.TestCase):
    def test_two_blocks_by_hand(self):
        # blocks {0:[1,2], 1:[3,4]}: theta 2.5, leave-one-block-out means 1.5 and 3.5,
        # se = sqrt((nblk-1)/nblk * var_around_loo_mean) = sqrt(1/2 * (1+1)) = 1
        theta, se, z = f3.jackknife(np.array([1.0, 2, 3, 4]), np.array([0, 0, 1, 1]), 2)
        self.assertAlmostEqual(theta, 2.5)
        self.assertAlmostEqual(se, 1.0)
        self.assertAlmostEqual(z, 2.5)

    def test_zero_statistic_gives_nan_z_not_crash(self):
        theta, se, z = f3.jackknife(np.zeros(6), np.array([0, 0, 1, 1, 2, 2]), 3)
        self.assertEqual((theta, se), (0.0, 0.0))
        self.assertTrue(np.isnan(z))  # Z is meaningless for an exactly-zero statistic


@_NEEDS_NP
class TestBlockCodes(unittest.TestCase):
    def test_blocks_follow_chrom_and_5mb(self):
        code = f3.block_codes(["1", "1", "1", "2", "chrX"], [1, 6_000_001, 5, 5, 5])
        self.assertEqual(len(set(code)), 4)          # 1_0, 1_1, 2_0, X_0
        self.assertEqual(code[0], code[2])          # chr1 pos 1 and 5 share the first block
        self.assertNotEqual(code[0], code[1])       # crossing 5 Mb starts a new block
        self.assertNotEqual(code[3], code[4])       # different chromosomes never share a block


@_NEEDS_NP
class TestMergeFreqs(unittest.TestCase):
    def _pos(self):
        return pd.DataFrame({"CHROM_N": ["1", "1", "2"], "POS": [100, 200, 50]},
                            index=pd.Index(["s1", "s2", "s3"], name="ID"))

    def test_allele_flip_is_dropped_and_counted(self):
        a = _afreq([("s1", "1", 100, "A", "G", 0.1), ("s2", "1", 200, "C", "T", 0.2),
                    ("s3", "2", 50, "G", "A", 0.3)])
        b = _afreq([("s1", "1", 100, "A", "G", 0.1), ("s2", "1", 200, "T", "C", 0.8),  # flipped
                    ("s3", "2", 50, "G", "A", 0.3)])
        df, nbad = f3.merge_freqs({"a": a, "b": b}, self._pos())
        self.assertEqual(nbad, 1)
        self.assertEqual(list(df.index), ["s1", "s3"])   # s2 dropped, not mixed

    def test_missing_site_in_one_source_is_dropped(self):
        a = _afreq([("s1", "1", 100, "A", "G", 0.1), ("s2", "1", 200, "C", "T", 0.2),
                    ("s3", "2", 50, "G", "A", 0.3)])
        b = _afreq([("s1", "1", 100, "A", "G", 0.1), ("s3", "2", 50, "G", "A", 0.3)])
        df, nbad = f3.merge_freqs({"a": a, "b": b}, self._pos())
        self.assertEqual((nbad, len(df)), (0, 2))        # s2 has no b frequency -> dropna removes it

    def test_sorted_by_chrom_then_position(self):
        a = _afreq([("s3", "2", 50, "G", "A", 0.3), ("s1", "1", 100, "A", "G", 0.1),
                    ("s2", "1", 200, "C", "T", 0.2)])
        df, _ = f3.merge_freqs({"a": a}, self._pos())
        self.assertEqual(list(df.index), ["s1", "s2", "s3"])


@_NEEDS_NP
class TestF3Matrix(unittest.TestCase):
    def _df(self):
        # Six sites in three 2-site blocks, shaped like real data: the outgroup sits near the
        # ancestral allele and both target and references are shifted derived, so outgroup-f3
        # rises as a profile group moves closer to the target (P1 is the near one, P2 far).
        rows = {
            "og":     [0.05, 0.95, 0.10, 0.90, 0.15, 0.85],
            "target": [0.60, 0.40, 0.55, 0.45, 0.52, 0.48],
            "P1":     [0.58, 0.42, 0.53, 0.47, 0.50, 0.50],
            "P2":     [0.30, 0.70, 0.35, 0.65, 0.40, 0.60],
            "poolA":  [0.44, 0.56, 0.44, 0.56, 0.46, 0.54],
            "poolB":  [0.36, 0.64, 0.41, 0.59, 0.43, 0.57],
        }
        df = pd.DataFrame(rows)
        df["CHROM_N"] = ["1"] * 6
        df["POS"] = [1, 2, 6_000_001, 6_000_002, 12_000_001, 12_000_002]
        return df

    def test_ranking_contrasts_admixture(self):
        res = f3.f3_matrix(self._df(), "og", "target", ["P1", "P2"], ["poolA", "poolB"])
        self.assertEqual(res["sites"], 6)
        self.assertEqual(res["blocks"], 3)
        self.assertEqual([r["set"] for r in res["ranked"]], ["P1", "P2"])   # desc by f3
        self.assertEqual([(c["a"], c["b"]) for c in res["contrasts"]], [("P1", "P2")])
        self.assertGreater(res["ranked"][0]["f3"], res["ranked"][1]["f3"])  # P1 shares more drift
        self.assertGreater(res["ranked"][0]["se"], 0)
        # admixture f3(target; poolA, poolB) = mean (t-A)(t-B): per-site products are
        # .0384 .0384 .0154 .0154 .0054 .0054 -> mean 0.019733
        adm = res["admixture"]
        self.assertEqual([(r["a"], r["b"]) for r in adm], [("poolA", "poolB")])
        self.assertAlmostEqual(adm[0]["f3"], 0.019733, places=5)
        self.assertGreater(adm[0]["f3"], 0)

    def test_ancient_mode_has_no_ranking(self):
        res = f3.f3_matrix(self._df(), None, "target", [], ["poolA", "poolB"])
        self.assertEqual(res["ranked"], [])
        self.assertEqual(res["contrasts"], [])
        self.assertEqual(len(res["admixture"]), 1)

    def test_contrast_uses_same_site_difference(self):
        # og and target constant, P2 = P1 + 0.1: each ranked row's per-site statistic still
        # varies across blocks (P1's block means differ), so both rows have SE > 0, but the
        # per-site difference is the constant (og-target)*(P2-P1) -> contrast SE 0. A naive
        # combination of the two rows' SEs would report noise the difference does not have.
        df = pd.DataFrame({"og": [0.05] * 4, "target": [0.55] * 4,
                           "P1": [0.50, 0.52, 0.60, 0.52]})
        df["P2"] = df["P1"] + 0.10
        df["CHROM_N"] = ["1"] * 4
        df["POS"] = [1, 2, 6_000_001, 6_000_002]
        res = f3.f3_matrix(df, "og", "target", ["P1", "P2"], [])
        self.assertGreater(res["ranked"][0]["se"], 0)
        self.assertGreater(res["ranked"][1]["se"], 0)
        self.assertAlmostEqual(res["contrasts"][0]["f3"], 0.05, places=12)
        self.assertLess(res["contrasts"][0]["se"], 1e-9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
