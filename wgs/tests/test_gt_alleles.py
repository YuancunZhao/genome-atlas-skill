"""Regression for the GT parser used by 11_pgx_extra / 27_candidate_genes.

Run: python3 wgs/tests/test_gt_alleles.py   (standard library only, no sample data needed).
The cases mirror what '.count("1")' got wrong: a full or partial no-call, and any genotype
involving the second or later ALT, all collapsed into a reference homozygote.
"""
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
from gt_alleles import gt_alleles  # noqa: E402


class TestGtAlleles(unittest.TestCase):
    def test_clean_diploid(self):
        self.assertEqual(gt_alleles("0/0", "A", "T"), ["A", "A"])
        self.assertEqual(gt_alleles("0|1", "A", "T"), ["A", "T"])
        self.assertEqual(gt_alleles("1/1", "A", "T"), ["T", "T"])

    def test_full_no_call(self):
        self.assertIsNone(gt_alleles("./.", "A", "T"))
        self.assertIsNone(gt_alleles(".|.", "A", "T"))
        self.assertIsNone(gt_alleles(".", "A", "T"))  # haploid no-call

    def test_partial_no_call(self):
        # count("1") read both of these as 0 alt alleles -> printed ref/ref
        self.assertEqual(gt_alleles("0/.", "A", "T"), ["A", None])
        self.assertEqual(gt_alleles("./1", "A", "T"), [None, "T"])

    def test_multiallelic(self):
        # count("1") read 2/2 as 0 -> printed a hom-ref; the true genotype is G/G
        self.assertEqual(gt_alleles("2/2", "A", "T,G"), ["G", "G"])
        self.assertEqual(gt_alleles("1/2", "A", "T,G"), ["T", "G"])
        self.assertEqual(gt_alleles("0/2", "A", "T,G"), ["A", "G"])

    def test_allele_index_beyond_record(self):
        self.assertEqual(gt_alleles("3/3", "A", "T,G"), ["<alt3>", "<alt3>"])

    def test_haploid(self):
        self.assertEqual(gt_alleles("0", "A", "T"), ["A"])
        self.assertEqual(gt_alleles("1", "A", "T"), ["T"])

    def test_malformed_shows_raw(self):
        self.assertEqual(gt_alleles("T/T", "A", "T"), ["GT:T/T"])


if __name__ == "__main__":
    unittest.main()
