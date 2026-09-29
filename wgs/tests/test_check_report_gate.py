"""Regression tests for the 32_check_report gate itself.

32 is the acceptance checker; when it demands things unconditionally, every configuration
that legitimately lacks them fails as a false positive -- which is exactly what happened to
ho_affinity on non-AADR builds (AN7). The module name starts with a digit, so it is loaded
via importlib; main() stays untouched, the tests drive the pure check functions and inspect
the FAIL accumulator."""
import importlib.util
import pathlib
import unittest

_spec = importlib.util.spec_from_file_location(
    "check_report", pathlib.Path(__file__).resolve().parents[1] / "scripts" / "32_check_report.py")
chk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(chk)


class TestHoAffinityGate(unittest.TestCase):
    def setUp(self):
        chk.FAIL.clear()

    def _sections(self, status):
        return {"sections": [{"id": "aadr", "status": status, "name_zh": "a", "name_en": "a"}],
                "ho_affinity": []}

    def test_empty_affinity_is_fine_when_aadr_unavailable(self):
        chk.check_shapes(self._sections("unavailable"))
        self.assertFalse([f for f in chk.FAIL if "ho_affinity" in f],
                         "a build without the ancient projection must not fail on ho_affinity")

    def test_empty_affinity_fails_when_aadr_claims_ok(self):
        chk.check_shapes(self._sections("ok"))
        self.assertTrue(any("ho_affinity" in f for f in chk.FAIL),
                        "aadr ok with no rows means the figure has nothing to draw")


class TestSourceChecksFollowTheReport(unittest.TestCase):
    """AN7: the source checks must read the outputs of the work tree the report came from,
    not the repository's own work/ directory. check_sv_source is exercised on a temp tree
    because it is pure file IO; check_naming needs bcftools and is covered by real runs."""

    def setUp(self):
        chk.FAIL.clear()
        import tempfile
        self._td = tempfile.TemporaryDirectory()
        self.run_dir = pathlib.Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def _tsv(self, genes):
        d = self.run_dir / "08_sv"
        d.mkdir(parents=True, exist_ok=True)
        (d / "sv_filtered.tsv").write_text(
            "chrom\tpos\twhole_gene_del\n1\t1000\t" + (",".join(genes) if genes else "") + "\n",
            encoding="utf-8")

    def test_stray_gene_fails_against_the_reports_own_tree(self):
        self._tsv(["REAL1"])                                   # the report's own source has REAL1 only
        D = {"sv_gene_dels": [{"chrom": "1", "pos": 1, "gene": "GHOST", "frac": 0.5}]}
        chk.check_sv_source(self.run_dir, D)
        self.assertTrue(any("GHOST" in f for f in chk.FAIL),
                        "a gene absent from this tree's whole_gene_del column must be flagged")

    def test_matching_tree_passes_and_missing_tree_is_skipped(self):
        self._tsv(["REAL1"])
        chk.check_sv_source(self.run_dir, {"sv_gene_dels": [{"chrom": "1", "pos": 1, "gene": "REAL1", "frac": 0.5}]})
        self.assertFalse(chk.FAIL)
        chk.FAIL.clear()
        empty = self.run_dir / "elsewhere"                    # no sv_filtered.tsv there: not an error
        empty.mkdir()
        chk.check_sv_source(empty, {"sv_gene_dels": [{"chrom": "1", "pos": 1, "gene": "REAL1", "frac": 0.5}]})
        self.assertFalse(chk.FAIL)


class TestF3Gate(unittest.TestCase):
    """W-T2: the f3 payload may carry explicit nulls (unavailable: zero-SE / non-finite),
    and must never carry NaN/Infinity -- which json.loads hands back as floats, so an
    isinstance-only gate waves them through as numbers."""

    def setUp(self):
        chk.FAIL.clear()

    def _row(self, **over):
        row = {"set": "P1", "n": 85, "label_zh": "P1", "label_en": "P1",
               "f3": 0.01, "se": 0.001, "z": 10.0}
        row.update(over)
        return row

    def _D(self, rows):
        return {"f3": {"modern": {"outgroup_f3": rows, "contrasts": []},
                       "ancient": {"admixture": []}}}

    def test_null_triples_are_unavailable_not_failures(self):
        chk.check_shapes(self._D([self._row(f3=None, se=None, z=None)]))
        self.assertFalse([f for f in chk.FAIL if "f3" in f],
                         "an unavailable row is a rendered state, not a build error")

    def test_point_estimate_with_null_z_is_valid(self):
        chk.check_shapes(self._D([self._row(z=None)]))
        self.assertFalse([f for f in chk.FAIL if "f3" in f],
                         "se==0 yields (f3, 0, null): the point estimate is still quotable")

    def test_nan_f3_fails_the_old_isinstance_hole(self):
        chk.check_shapes(self._D([self._row(f3=float("nan"))]))
        self.assertTrue(any("is neither a finite number nor null" in f for f in chk.FAIL),
                        "NaN parses as float -- the gate must reject it explicitly")

    def test_infinity_and_z_without_estimate_fail(self):
        chk.check_shapes(self._D([self._row(se=float("inf"))]))
        self.assertTrue(any("is neither a finite number nor null" in f for f in chk.FAIL))
        chk.FAIL.clear()
        chk.check_shapes(self._D([self._row(f3=None, se=None, z=3.0)]))
        self.assertTrue(any("a Z without its f3/SE" in f for f in chk.FAIL),
                        "a verdict number with no estimate behind it is not a statistic")


class TestNodeGate(unittest.TestCase):
    def test_missing_node_fails_instead_of_passing(self):
        chk.FAIL.clear()
        real_which = chk.shutil.which
        chk.shutil.which = lambda name: None          # simulate a machine without node
        try:
            chk.check_runtime(None)
        finally:
            chk.shutil.which = real_which
        self.assertTrue(any("node is not installed" in f for f in chk.FAIL),
                        "a skipped runtime check must not count as passed")


if __name__ == "__main__":
    unittest.main(verbosity=2)
