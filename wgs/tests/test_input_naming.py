"""H3 regression: consumers of step-01 outputs read target.* first.

Step 01 writes target.norm/pass.vcf.gz; several consumers historically read
{SAMPLE}.* names that only existed via hand-made symlinks. These tests pin the
unified naming so a revert fails here instead of silently reading a stale link.
"""
import pathlib
import unittest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"


class InputNaming(unittest.TestCase):
    def test_norm_vcf_consumers_read_target_first(self):
        for name in ("03_complete_set.py", "11_pgx_extra.py", "21_somatic.py",
                     "23_bloodgroups.py", "27_candidate_genes.py"):
            src = (SCRIPTS / name).read_text()
            self.assertIn("00_input/target.norm.vcf.gz", src,
                          f"{name} lost the target.norm path")
            self.assertLess(src.index("00_input/target.norm.vcf.gz"),
                            src.index("00_input/{SAMPLE}.norm.vcf.gz"),
                            f"{name}: the {{SAMPLE}}-literal path is used before target.*")

    def test_pass_vcf_consumer_reads_target_first(self):
        src = (SCRIPTS / "06_mtdna.py").read_text()
        self.assertIn("00_input/target.pass.vcf.gz", src)
        self.assertLess(src.index("00_input/target.pass.vcf.gz"),
                        src.index("00_input/{SAMPLE}.pass.vcf.gz"),
                        "06_mtdna: the SAMPLE-literal path is used before target.*")

    def test_cram_readers_use_configured_reads(self):
        for name in ("05_y_haplogroup.py", "21_somatic.py"):
            src = (SCRIPTS / name).read_text()
            self.assertIn("CRAM = READS", src,
                          f"{name} should read the configured CRAM (READS), not a hand-linked {{SAMPLE}}.cram")


if __name__ == "__main__":
    unittest.main()
