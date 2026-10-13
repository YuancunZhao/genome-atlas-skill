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


class TestKirGate(unittest.TestCase):
    """M3/H2: gene rows and kir_summary must agree; unknown (not run) is not absent."""

    def setUp(self):
        chk.FAIL.clear()

    def _D(self, rows, summary):
        return {"kir": rows, "kir_summary": summary}

    def test_rows_with_real_summary_pass(self):
        D = self._D([{"gene": "KIR2DL1", "present": True, "call": "KIR2DL1*003"}],
                    {"haplotype": "Bx (at least one B haplotype)",
                     "present": "KIR2DL1", "C_ligands": "C1", "B_epitopes": "Bw6/Bw6"})
        chk.check_shapes(D)
        self.assertFalse([f for f in chk.FAIL if "kir" in f.lower()])

    def test_rows_without_summary_fail(self):
        chk.check_shapes(self._D([{"gene": "KIR2DL1", "present": True}], None))
        self.assertTrue(any("kir_summary is missing" in f for f in chk.FAIL),
                        "the haplotype line must not fall back to an invented example")

    def test_summary_claim_without_rows_is_stale(self):
        D = self._D([], {"haplotype": "AA (two A haplotypes)"})
        chk.check_shapes(D)
        self.assertTrue(any("but the gene table is empty" in f for f in chk.FAIL),
                        "a haplotype call with no typing rows is another run's output")

    def test_unavailable_summary_with_empty_rows_is_consistent(self):
        D = self._D([], {"haplotype": "unavailable (T1K was not run)"})
        chk.check_shapes(D)
        self.assertFalse([f for f in chk.FAIL if "kir" in f.lower()])


class TestNodeGate(unittest.TestCase):
    def test_missing_node_fails_instead_of_passing(self):
        chk.FAIL.clear()
        real_which = chk.shutil.which
        chk.shutil.which = lambda name: None          # simulate a machine without node
        try:
            chk.check_runtime(None, {})    # (html, D): the blood render check rides along since 2026-10-13
        finally:
            chk.shutil.which = real_which
        self.assertTrue(any("node is not installed" in f for f in chk.FAIL),
                        "a skipped runtime check must not count as passed")


if __name__ == "__main__":
    unittest.main(verbosity=2)


class TestBloodMarkerGate(unittest.TestCase):
    """血型证据表绑定（血型核查 2026-10-13）：叙述表来自人工文案，30 另注入 23 的逐位点
    测序输出 BLOODMK，31/模板把它渲染成第二张表。本门要守住三件事：数据与 sections 状态
    互洽、node 运行时该表逐行渲染（报告因此永远携带本样本的血型证据）、三方接线齐备。"""

    def setUp(self):
        chk.FAIL.clear()

    @staticmethod
    def _rows():
        return [{"system": "Kidd (JK)", "rsid": "rs1058396", "variant": "SLC14A1 838G>A",
                 "locus": "18:43319519", "genotype": "G/A", "evidence": "PASS DP=77 AD=36,37",
                 "meaning": "ref=JK*A, alt=JK*B"},
                {"system": "Vel", "rsid": "rs566629828", "variant": "SMIM1 64_80del",
                 "locus": "-", "genotype": "not in 1000G panel", "evidence": "",
                 "meaning": ""}]

    def test_data_without_section_row_fails_and_vice_versa(self):
        chk.check_blood_markers({"blood_markers": self._rows(),
                                 "sections": [{"id": "kir", "status": "ok", "name_zh": "k", "name_en": "k"}]})
        self.assertTrue(any("blood" in f for f in chk.FAIL),
                        "blood_markers without a sections.blood row must be flagged")
        chk.FAIL.clear()
        chk.check_blood_markers({"blood_markers": [],
                                 "sections": [{"id": "blood", "status": "ok", "name_zh": "b", "name_en": "b"}]})
        self.assertTrue(any("blood" in f for f in chk.FAIL),
                        "a blood section claiming ok with no marker rows is stale output")

    def test_consistent_data_and_section_pass(self):
        chk.check_blood_markers({"blood_markers": self._rows(),
                                 "sections": [{"id": "blood", "status": "ok", "name_zh": "b", "name_en": "b"}]})
        self.assertFalse([f for f in chk.FAIL if "blood" in f])
        chk.FAIL.clear()
        chk.check_blood_markers({"blood_markers": [],
                                 "sections": [{"id": "blood", "status": "unavailable", "code": "missing_output",
                                               "reason_zh": "r", "reason_en": "r", "name_zh": "b", "name_en": "b"}]})
        self.assertFalse([f for f in chk.FAIL if "blood" in f])
        chk.FAIL.clear()
        chk.check_blood_markers({})   # an older report_data.json without the key at all
        self.assertFalse([f for f in chk.FAIL if "blood" in f])

    def test_rendered_table_must_carry_every_row(self):
        rows = self._rows()
        good = ("<tr><td>Kidd (JK)</td><td>rs1058396</td><td>SLC14A1 838G>A</td><td>G/A</td>"
                "<td>18:43319519 · PASS DP=77 AD=36,37</td><td>ref=JK*A, alt=JK*B</td></tr>"
                "<tr><td>Vel</td><td>rs566629828</td><td>SMIM1 64_80del</td>"
                "<td>not in 1000G panel</td><td></td><td></td></tr>")
        chk.check_bloodmk_render(good, rows)
        self.assertFalse([f for f in chk.FAIL if "blood" in f])
        chk.FAIL.clear()
        stripped = good.replace("rs566629828", "")     # one marker vanished from the render
        chk.check_bloodmk_render(stripped, rows)
        self.assertTrue(any("rs566629828" in f for f in chk.FAIL))
        chk.FAIL.clear()
        chk.check_bloodmk_render("", rows)             # rows exist but the table never rendered
        self.assertTrue(any("blood" in f for f in chk.FAIL))
        chk.FAIL.clear()
        chk.check_bloodmk_render("", [])               # no rows: an empty table is the honest state
        chk.check_bloodmk_render(good, [])             # rows rendered with none in the data: stale payload
        self.assertTrue(any("blood" in f for f in chk.FAIL))

    def test_report_wiring_binds_the_marker_table_end_to_end(self):
        repo = pathlib.Path(__file__).resolve().parents[1]
        js = (repo / "templates/report_script.js").read_text(encoding="utf-8")
        self.assertIn("__BLOODMK__", js, "the payload declaration line is missing from the template JS")
        self.assertIn("tb_blood_mk", js, "deepRender no longer fills the marker table")
        body = (repo / "templates/report_body.html").read_text(encoding="utf-8")
        self.assertIn('id="tb_blood_mk"', body, "the marker table has no element to render into")
        s31 = (repo / "scripts/31_html_report.py").read_text(encoding="utf-8")
        self.assertIn('__BLOODMK__', s31, "31 does not inject the marker payload")
        self.assertIn('D.get("blood_markers")', s31, "31 must tolerate an older report_data.json")


class TestPgxMachineTableGate(unittest.TestCase):
    """修-3（叙述审计 2026-10-13）：PGx 叙述表来自人工文案，PharmCAT/Cyrius 的机器输出
    （30 注入的 D.pgx_pharmcat）从未被渲染。与血型同一方案：31 以 PGXMK 载荷注入、模板
    在叙述表下渲染机器表、32 在 node 运行时逐行核对——换样本重跑后叙述与机器表的矛盾
    当场露出，而不是无声沿用旧样本的双倍型。"""

    def setUp(self):
        chk.FAIL.clear()

    @staticmethod
    def _rows():
        return [{"gene": "CYP2C19", "diplotype": "*1/*2", "phenotype": "Intermediate Metabolizer",
                 "call_source": "MATCHER", "n_warnings": 0},
                {"gene": "CYP2D6", "diplotype": "*1/*1", "phenotype": "Normal Metabolizer",
                 "call_source": "OUTSIDE", "n_warnings": 0}]

    @staticmethod
    def _good():
        return ("<tr><td>CYP2C19</td><td>*1/*2</td><td>Intermediate Metabolizer</td>"
                "<td>MATCHER</td></tr><tr><td>CYP2D6</td><td>*1/*1</td>"
                "<td>Normal Metabolizer</td><td>OUTSIDE</td></tr>")

    def test_rendered_table_must_carry_every_row(self):
        chk.check_pgxmk_render(self._good(), self._rows())
        self.assertFalse([f for f in chk.FAIL if "pgx" in f])
        chk.FAIL.clear()
        chk.check_pgxmk_render(self._good().replace("CYP2D6", "CYP2D9"), self._rows())
        self.assertTrue(any("CYP2D6" in f for f in chk.FAIL))
        chk.FAIL.clear()
        chk.check_pgxmk_render("", self._rows())
        self.assertTrue(any("pgx" in f for f in chk.FAIL))
        chk.FAIL.clear()
        chk.check_pgxmk_render("", [])
        self.assertFalse([f for f in chk.FAIL if "pgx" in f])
        chk.FAIL.clear()
        chk.check_pgxmk_render(self._good(), [])
        self.assertTrue(any("pgx" in f for f in chk.FAIL))

    def test_report_wiring_binds_the_machine_pgx_table(self):
        repo = pathlib.Path(__file__).resolve().parents[1]
        js = (repo / "templates/report_script.js").read_text(encoding="utf-8")
        self.assertIn("__PGXMK__", js, "the payload declaration line is missing from the template JS")
        self.assertIn("tb_pgx_mk", js, "renderPgx no longer fills the machine table")
        body = (repo / "templates/report_body.html").read_text(encoding="utf-8")
        self.assertIn('id="tb_pgx_mk"', body, "the machine table has no element to render into")
        s31 = (repo / "scripts/31_html_report.py").read_text(encoding="utf-8")
        self.assertIn("__PGXMK__", s31, "31 does not inject the machine pgx payload")
        self.assertIn('D.get("pgx_pharmcat")', s31, "31 must tolerate an older report_data.json")


class TestLofMachineTableGate(unittest.TestCase):
    """修-4（叙述审计 2026-10-13）：findings 列表是人工解读（基因+基因型+读段级叙述），
    其机器事实集是 05c 的罕见 LoF 终表。30 逐行注入 D.lof_rare、模板在 findings 卡下渲染
    证据表、32 逐行核对——叙述声称的携带者若不在机器表里，当场露出。"""

    def setUp(self):
        chk.FAIL.clear()

    @staticmethod
    def _rows():
        return [{"gene": "SLC26A2", "variant": "5:149360024 TG>T", "csq": "frameshift",
                 "zyg": "het", "clnsig": "Pathogenic/Likely_pathogenic"},
                {"gene": "GENE2", "variant": "1:2000 A>G", "csq": "frameshift",
                 "zyg": "hom/hemi", "clnsig": "Pathogenic/Likely_pathogenic"}]

    @staticmethod
    def _good():
        return ("<tr><td>SLC26A2</td><td>5:149360024 TG>T</td><td>frameshift</td>"
                "<td>het</td><td>Pathogenic/Likely_pathogenic</td></tr>"
                "<tr><td>GENE2</td><td>1:2000 A>G</td><td>frameshift</td>"
                "<td>hom/hemi</td><td>Pathogenic/Likely_pathogenic</td></tr>")

    def test_rendered_table_must_carry_every_row(self):
        chk.check_lof_mk_render(self._good(), self._rows())
        self.assertFalse([f for f in chk.FAIL if "rare LoF" in f or "lof" in f])
        chk.FAIL.clear()
        chk.check_lof_mk_render(self._good().replace("SLC26A2", "SLC26A3"), self._rows())
        self.assertTrue(any("SLC26A2" in f for f in chk.FAIL))
        chk.FAIL.clear()
        chk.check_lof_mk_render("", self._rows())
        self.assertTrue(any("rare LoF" in f for f in chk.FAIL))
        chk.FAIL.clear()
        chk.check_lof_mk_render("", [])
        self.assertFalse([f for f in chk.FAIL if "rare LoF" in f or "lof" in f])
        chk.FAIL.clear()
        chk.check_lof_mk_render(self._good(), [])
        self.assertTrue(any("rare LoF" in f for f in chk.FAIL))

    def test_report_wiring_binds_the_lof_evidence_table(self):
        repo = pathlib.Path(__file__).resolve().parents[1]
        js = (repo / "templates/report_script.js").read_text(encoding="utf-8")
        self.assertIn("__LOFMK__", js)
        self.assertIn("tb_lof_mk", js)
        body = (repo / "templates/report_body.html").read_text(encoding="utf-8")
        self.assertIn('id="tb_lof_mk"', body)
        s31 = (repo / "scripts/31_html_report.py").read_text(encoding="utf-8")
        self.assertIn("__LOFMK__", s31)
        self.assertIn('D.get("lof_rare")', s31)
