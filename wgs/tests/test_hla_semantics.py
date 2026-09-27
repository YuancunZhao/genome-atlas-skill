"""Regression for the HLA/KIR semantics shared by steps 24 and 25.

Run: python3 wgs/tests/test_hla_semantics.py   (standard library only, no sample data needed).
Covers the two defects the inline code shipped: rsplit truncating two-field alleles to one
field (so no Bw4 whitelist entry ever matched), and Bw6 being the default for anything the
whitelist missed -- including the split B*15 family and genes the run never typed.
"""
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
from hla_semantics import a311, bw_epitope, c_ligand, two_fields  # noqa: E402


class TestTwoFields(unittest.TestCase):
    def test_three_field_reduced(self):
        self.assertEqual(two_fields("46:01:01"), "46:01")

    def test_two_field_preserved(self):
        # rsplit turned 44:02 into 44, silently unmatching every whitelist entry
        self.assertEqual(two_fields("44:02"), "44:02")
        self.assertEqual(two_fields("51:01"), "51:01")

    def test_no_fields(self):
        self.assertEqual(two_fields("46"), "46")


class TestBwEpitope(unittest.TestCase):
    def test_bwd_families(self):
        self.assertEqual(bw_epitope("44:02"), "Bw4")
        self.assertEqual(bw_epitope("27:05"), "Bw4")
        self.assertEqual(bw_epitope("57:01"), "Bw4")

    def test_bw6_families(self):
        self.assertEqual(bw_epitope("46:01"), "Bw6")
        self.assertEqual(bw_epitope("35:01"), "Bw6")

    def test_b15_split_is_unknown_not_guessed(self):
        self.assertEqual(bw_epitope("15:10"), "Bw4")
        self.assertEqual(bw_epitope("15:01"), "unknown (B*15 split family)")


class TestCLigand(unittest.TestCase):
    def test_groups(self):
        self.assertEqual(c_ligand("01:02"), "C1")
        self.assertEqual(c_ligand("07:02"), "C1")
        self.assertEqual(c_ligand("02:07"), "C2")
        self.assertEqual(c_ligand("99:99"), "?")


class TestA311(unittest.TestCase):
    def test_selection(self):
        self.assertEqual(a311(["24:20", "03:01", "11:01"]), ["03:01", "11:01"])
        self.assertEqual(a311(["02:07"]), [])


if __name__ == "__main__":
    unittest.main()
