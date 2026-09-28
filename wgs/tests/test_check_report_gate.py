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


if __name__ == "__main__":
    unittest.main(verbosity=2)
