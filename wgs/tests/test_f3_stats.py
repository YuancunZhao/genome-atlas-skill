"""Unit tests for the pure f3 math in 28_f3_stats.py (jackknife, block coding, allele-checked
frequency merge with per-site allele counts, finite-sample-corrected matrix assembly, the
None-safe rounding used by the payload). Importing the module must not touch wgsconfig: that
creates the work directories, a side effect tests cannot have.

numpy/pandas live on the analysis server (the local suite is stdlib-only, like 30's), so the
tests skip themselves when those imports are unavailable."""
import importlib.util
import json
import math
import pathlib
import tempfile
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

BIG_N = 10 ** 6        # panel-sized allele counts: correction term ~2.5e-7, invisible at 1e-5


def _afreq(rows):
    """rows: (ID, CHROM, POS, REF, ALT, freq, obs_ct) -> plink2 .afreq-shaped frame.

    OBS_CT is the called-allele count plink2 really emits; the correction reads it, so every
    fixture must carry it exactly like production input does."""
    return pd.DataFrame(rows, columns=["ID", "CHROM", "POS", "REF", "ALT", "ALT_FREQS", "OBS_CT"])


def _with_n(df, n=BIG_N):
    """Attach N:<set> columns (the production frame shape merge_freqs returns)."""
    for c in [c for c in df.columns if c not in ("CHROM_N", "POS", "REF", "ALT")
              and not c.startswith("N:")]:
        df[f"N:{c}"] = n
    return df


@_NEEDS_NP
class TestJackknife(unittest.TestCase):
    def test_two_blocks_by_hand(self):
        # blocks {0:[1,2], 1:[3,4]}: theta 2.5, leave-one-block-out means 1.5 and 3.5,
        # se = sqrt((nblk-1)/nblk * var_around_loo_mean) = sqrt(1/2 * (1+1)) = 1
        theta, se, z = f3.jackknife(np.array([1.0, 2, 3, 4]), np.array([0, 0, 1, 1]), 2)
        self.assertAlmostEqual(theta, 2.5)
        self.assertAlmostEqual(se, 1.0)
        self.assertAlmostEqual(z, 2.5)

    def test_zero_statistic_gives_null_z_not_crash(self):
        theta, se, z = f3.jackknife(np.zeros(6), np.array([0, 0, 1, 1, 2, 2]), 3)
        self.assertEqual((theta, se), (0.0, 0.0))
        self.assertIsNone(z)  # Z is meaningless for an exactly-zero statistic: null, not inf

    def test_empty_or_nonfinite_inputs_are_unavailable(self):
        self.assertEqual(f3.jackknife(np.array([]), np.array([], dtype=int), 1), (None, None, None))
        self.assertEqual(f3.jackknife(np.array([1.0, np.inf]), np.array([0, 1]), 2), (None, None, None))
        self.assertEqual(f3.jackknife(np.array([1.0]), np.array([0]), 0), (None, None, None))


@_NEEDS_NP
class TestNum(unittest.TestCase):
    def test_none_and_nonfinite_stay_none(self):
        self.assertIsNone(f3._num(None))
        self.assertIsNone(f3._num(float("nan")))
        self.assertIsNone(f3._num(float("inf")))
        self.assertEqual(f3._num(0.1234567), 0.123457)
        self.assertEqual(f3._num(1.0, 2), 1.0)

    def test_num_output_is_json_safe(self):
        vals = [f3._num(v) for v in (None, float("nan"), float("inf"), 0.5, -0.25)]
        json.dumps(vals, allow_nan=False)  # must not raise: payload paths all go through _num


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
        a = _afreq([("s1", "1", 100, "A", "G", 0.1, 100), ("s2", "1", 200, "C", "T", 0.2, 100),
                    ("s3", "2", 50, "G", "A", 0.3, 100)])
        b = _afreq([("s1", "1", 100, "A", "G", 0.1, 200), ("s2", "1", 200, "T", "C", 0.8, 200),  # flipped
                    ("s3", "2", 50, "G", "A", 0.3, 200)])
        df, nbad, nlow = f3.merge_freqs({"a": a, "b": b}, self._pos())
        self.assertEqual((nbad, nlow), (1, 0))
        self.assertEqual(list(df.index), ["s1", "s3"])   # s2 dropped, not mixed

    def test_missing_site_in_one_source_is_dropped(self):
        a = _afreq([("s1", "1", 100, "A", "G", 0.1, 100), ("s2", "1", 200, "C", "T", 0.2, 100),
                    ("s3", "2", 50, "G", "A", 0.3, 100)])
        b = _afreq([("s1", "1", 100, "A", "G", 0.1, 200), ("s3", "2", 50, "G", "A", 0.3, 200)])
        df, nbad, nlow = f3.merge_freqs({"a": a, "b": b}, self._pos())
        self.assertEqual((nbad, nlow, len(df)), (0, 0, 2))  # s2 has no b frequency -> dropna

    def test_sorted_by_chrom_then_position(self):
        a = _afreq([("s3", "2", 50, "G", "A", 0.3, 100), ("s1", "1", 100, "A", "G", 0.1, 100),
                    ("s2", "1", 200, "C", "T", 0.2, 100)])
        df, _, _ = f3.merge_freqs({"a": a}, self._pos())
        self.assertEqual(list(df.index), ["s1", "s2", "s3"])

    def test_obs_ct_is_carried_as_n_columns(self):
        # Dropping OBS_CT was the uncorrected-estimator bug: the merged frame must hand each
        # set's per-site allele count to f3_matrix as N:<set>.
        a = _afreq([("s1", "1", 100, "A", "G", 0.1, 100), ("s2", "1", 200, "C", "T", 0.2, 300)])
        b = _afreq([("s1", "1", 100, "A", "G", 0.1, 50), ("s2", "1", 200, "C", "T", 0.2, 60)])
        df, _, _ = f3.merge_freqs({"a": a, "b": b}, self._pos().iloc[:2])
        self.assertEqual(list(df["N:a"]), [100, 300])
        self.assertEqual(list(df["N:b"]), [50, 60])

    def test_site_below_min_alleles_is_dropped_and_counted(self):
        # nA < 2 cannot estimate its sampling variance; such a site must leave the panel
        # loudly (counted), not silently bias every statistic computed over it.
        a = _afreq([("s1", "1", 100, "A", "G", 0.1, 1),    # single called allele
                    ("s2", "1", 200, "C", "T", 0.2, 100)])
        b = _afreq([("s1", "1", 100, "A", "G", 0.1, 50), ("s2", "1", 200, "C", "T", 0.2, 60)])
        df, nbad, nlow = f3.merge_freqs({"a": a, "b": b}, self._pos().iloc[:2])
        self.assertEqual((nbad, nlow), (0, 1))
        self.assertEqual(list(df.index), ["s2"])


@_NEEDS_NP
class TestF3Matrix(unittest.TestCase):
    def _df(self):
        # Six sites in three 2-site blocks, shaped like real data: the outgroup sits near the
        # ancestral allele and both target and references are shifted derived, so outgroup-f3
        # rises as a profile group moves closer to the target (P1 is the near one, P2 far).
        # Allele counts are panel-sized, so the correction term (~2.5e-7/site) is invisible
        # at the 1e-5 the expectations below are stated to.
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
        return _with_n(df)

    def test_ranking_contrasts_admixture(self):
        res = f3.f3_matrix(self._df(), "og", "target", ["P1", "P2"], ["poolA", "poolB"])
        self.assertEqual(res["sites"], 6)
        self.assertEqual(res["blocks"], 3)
        self.assertEqual([r["set"] for r in res["ranked"]], ["P1", "P2"])   # desc by f3
        self.assertEqual([(c["a"], c["b"]) for c in res["contrasts"]], [("P1", "P2")])
        self.assertGreater(res["ranked"][0]["f3"], res["ranked"][1]["f3"])  # P1 shares more drift
        self.assertGreater(res["ranked"][0]["se"], 0)
        # corrected admixture f3(target; poolA, poolB): per-site products .0384 .0384 .0154
        # .0154 .0054 .0054, correction ~2.5e-7 -> 0.019733 to 5 places
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
        res = f3.f3_matrix(_with_n(df), "og", "target", ["P1", "P2"], [])
        self.assertGreater(res["ranked"][0]["se"], 0)
        self.assertGreater(res["ranked"][1]["se"], 0)
        self.assertAlmostEqual(res["contrasts"][0]["f3"], 0.05, places=12)
        self.assertLess(res["contrasts"][0]["se"], 1e-9)

    def test_finite_sample_correction_removes_diploid_target_bias(self):
        # The review's analytic zero: both sources exact 0.5, target dosages 0/1/1/2 over four
        # sites (frequencies 0, .5, .5, 1; a single diploid, nA=2). True f3 is 0 at every
        # site, but the raw product mean is +0.125. Per-site correction pa(1-pa)/(nA-1)
        # gives .25-.25-.25+.25 -> corrected f3 exactly 0 (all values exact in binary).
        df = pd.DataFrame({"target": [0.0, 0.5, 0.5, 1.0],
                           "poolA": [0.5] * 4, "poolB": [0.5] * 4})
        df["CHROM_N"] = ["1"] * 4
        df["POS"] = [1, 2, 6_000_001, 6_000_002]
        _with_n(df)
        df["N:target"] = 2                      # one diploid: the worst case for the bias
        res = f3.f3_matrix(df, None, "target", [], ["poolA", "poolB"])
        self.assertEqual(res["admixture"][0]["f3"], 0.0)
        # uncorrected mean would be exactly 0.125 -- regression guard on the estimator itself
        raw = ((df.target - .5) ** 2).mean()
        self.assertEqual(raw, 0.125)
        # each block's corrected values net to zero -> SE exactly 0 -> Z must be null
        # (the old code emitted z = 0/0 = Infinity here, not "no verdict")
        self.assertEqual(res["admixture"][0]["se"], 0.0)
        self.assertIsNone(res["admixture"][0]["z"])

    def test_constant_statistic_zero_se_yields_null_z(self):
        # Everything constant -> per-site statistic constant -> SE exactly 0: the point
        # estimate stays, Z is reported as null (unavailable), never Infinity.
        df = pd.DataFrame({"target": [0.5] * 4, "poolA": [0.5] * 4, "poolB": [0.5] * 4})
        df["CHROM_N"] = ["1"] * 4
        df["POS"] = [1, 2, 6_000_001, 6_000_002]
        res = f3.f3_matrix(_with_n(df), None, "target", [], ["poolA", "poolB"])
        self.assertAlmostEqual(res["admixture"][0]["f3"], -0.25 / (BIG_N - 1), places=12)
        self.assertEqual(res["admixture"][0]["se"], 0.0)
        self.assertIsNone(res["admixture"][0]["z"])
        # and the payload path can serialize it without NaN/Infinity
        json.dumps([f3._num(v) for v in
                    (res["admixture"][0]["f3"], res["admixture"][0]["se"],
                     res["admixture"][0]["z"])], allow_nan=False)


@_NEEDS_NP
class TestWriteTsv(unittest.TestCase):
    def test_unavailable_rows_write_na(self):
        with tempfile.TemporaryDirectory() as td:
            p = pathlib.Path(td) / "out.tsv"
            f3._write_tsv(p, [{"test": "ok", "f3": 0.5, "SE": 0.1, "Z": 5.0},
                              {"test": "degenerate", "f3": None, "SE": 0.0, "Z": None}])
            lines = p.read_text().strip().split("\n")
        self.assertEqual(lines[0].split("\t"), ["test", "f3", "SE", "Z"])
        self.assertEqual(lines[2].split("\t"), ["degenerate", "NA", "0.000000", "NA"])


@_NEEDS_NP
class TestRegionsFromLocations(unittest.TestCase):
    """地点表迁移：28 的南北区域池改从统一人工表 ancestry_locations.tsv 的 note token 取。

    旧行为读已删除的独立区域表的 region 列；迁移后 region=north|south 藏在 note 里，
    其余（无 token / region=unclassified / 1000G 采样地行）一律 ''——池分组本来就跳过 ''，
    与旧表空 region 行为一致。取图测试用真实统一表的行形状。"""

    def test_note_token_extraction(self):
        loc = pd.DataFrame([
            {"source_id": "China_Baligang_BA_EasternZhou",
             "note": "省份=河南；region=north；八里岗遗址在河南淅川，淮河以北的河南境"},
            {"source_id": "China_Taiwan_Han", "note": "省份=台湾；非 China_ 前缀不参与南北池"},   # 无 token
            {"source_id": "China_IA", "note": "region=unclassified；跨区域聚合标签，区域不适用"},
            {"source_id": "KHV", "note": "superpop=EAS；描述含采样城市"},                          # 1000G 行
            {"source_id": "China_Baoj", "note": "省份=陕西；region=south；宝鸡在秦岭—淮河以南"},
        ])
        self.assertEqual(f3.regions_from_locations(loc),
                         {"China_Baligang_BA_EasternZhou": "north",
                          "China_Taiwan_Han": "", "China_IA": "", "KHV": "",
                          "China_Baoj": "south"})

    def test_real_unified_table_maps_the_94_aadr_rows(self):
        """对真实统一表整表跑一遍：有 region= 的恰为北 42 + 南 30，与迁移前逐行一致。"""
        panel = pathlib.Path(__file__).resolve().parents[1] / "panel" / "ancestry_locations.tsv"
        loc = pd.read_csv(panel, sep="\t", comment="#", dtype=str, keep_default_na=False)
        reg = f3.regions_from_locations(loc)
        aadr = [k for k in reg if loc.set_index("source_id").loc[k, "dataset"] == "AADR"]
        north = sum(1 for k in aadr if reg[k] == "north")
        south = sum(1 for k in aadr if reg[k] == "south")
        self.assertEqual((north, south), (42, 30), "迁移前 94 行的区域判定必须原样保留")
        # 统一表里的 1000G 行永远不带 region token：采样地不是南北池成员
        kg = [k for k in reg if loc.set_index("source_id").loc[k, "dataset"] == "1000G"]
        self.assertTrue(kg and all(reg[k] == "" for k in kg))


class TestAncientRegionSource(unittest.TestCase):
    """源级接线（不依赖 numpy）：28 必须读统一表，且不得再引用任何已删除的旧表。"""

    def test_28_reads_the_unified_location_table(self):
        src = (pathlib.Path(__file__).resolve().parents[1] / "scripts" / "28_f3_stats.py"
               ).read_text(encoding="utf-8")
        self.assertIn('PANEL / "ancestry_locations.tsv"', src)
        self.assertNotIn("aadr_site_regions", src)
        self.assertNotIn("kg_population_locations", src)


@_NEEDS_NP
class Test28FailureInvalidatesOldOk(unittest.TestCase):
    """复审 §3.2 P0（H6/AN0/AN5 失败生命周期）——真实 28 子进程反例的回归。

    复现：28 因 f3_groups 的参考群不在 kg.common.psam 里 sys.exit(1)，旧 manifest=ok 与旧
    f3_stats.json 原样保留，30 的 analysis_state() 照旧准入。现在 28 开工先写失效记录；
    崩溃后旧 ok 不得再冒充本次结果。触发点在 plink2 调用之前，替身 plink2 不需要存在。"""

    def test_reference_group_missing_invalidates_previous_ok(self):
        import os
        import subprocess
        import sys
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
        import ancestry_data as ad
        with tempfile.TemporaryDirectory() as td:
            td = pathlib.Path(td)
            W = td / "work" / "wgs"
            d04 = W / "04_ancestry"
            f3dir = d04 / "f3"
            f3dir.mkdir(parents=True)
            # 上一轮的成功产物
            ad.write_manifest(f3dir / "manifest.json", ad.build_manifest(
                "TESTSAMPLE", "28-f3-stats", state="ok",
                parameters={"estimator": "site-mean-corrected", "block_mb": 5, "min_group_n": 2}))
            (f3dir / "f3_stats.json").write_text(
                json.dumps({"sample_id": "TESTSAMPLE", "modern": {"value": 0.01}}), encoding="utf-8")
            # 本次输入：modern 齐备但 psam 的群体不含 f3_groups 成员 → 真实复现的 sys.exit(1)
            (d04 / "prune.prune.in").write_text("rs1\n", encoding="utf-8")
            (d04 / "kg.common.psam").write_text(
                "#FID\tIID\tSuperPop\tPopulation\nFAM\tX\tZZZ\tZZZ\n", encoding="utf-8")
            (d04 / "kg.common.pvar").write_text(
                "#CHROM\tPOS\tID\tREF\tALT\n1\t1000\trs1\tA\tG\n", encoding="utf-8")
            (W / "02_complete").mkdir(parents=True)
            (W / "02_complete" / "TESTSAMPLE.1kg.pgen").write_bytes(b"")
            cfg = td / "config.yaml"
            cfg.write_text("".join([
                "sample_id: TESTSAMPLE\n",
                f"work_dir: {json.dumps(str(td / 'work'))}\n",
                f"plink2: {json.dumps(str(td / 'fake_plink2'))}\n",
            ]), encoding="utf-8")
            script = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "28_f3_stats.py"
            r = subprocess.run([sys.executable, str(script)], capture_output=True, text=True,
                               env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=120)
            self.assertNotEqual(r.returncode, 0, "参考群缺失必须以非零退出终止")
            self.assertIn("absent from kg.common.psam", r.stderr + r.stdout)
            # 单独准入（30 的调用形状）：旧 f3_stats.json 不得再作为本次结果交出
            state, reason, doc = ad.analysis_state(f3dir, {"sample_id": "TESTSAMPLE"},
                                                   names=("f3_stats.json",))
            self.assertEqual((state, reason), ("unavailable", "run_begun_not_published"))
            self.assertIsNone(doc)


if __name__ == "__main__":
    unittest.main(verbosity=2)
