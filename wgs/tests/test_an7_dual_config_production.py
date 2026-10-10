"""AN7 验收（复审"最终至少两份独立配置（含非EAS）的生产者→30→31→32"）。

两份干净独立配置各在各的空 work 目录里跑**真实脚本链**：
- 配置 A：ref_superpop=EUR + ref_subpops=[GBR]（非东亚区域）；
- 配置 B：无 ref_superpop（global-only）。
plink2 用一个替身（同时伺候 04 的剪枝/投影与 04b 的 smiss 重算），其余全是生产脚本：
04 → 04b → 30 → 31（文案走内置示例回退并带横幅）→ 32（含 node 运行时图检查）。
A 交付区域视图（pca_eas/near_eas/pop=EUR），B 按 global-only 交付（pca_eas=None +
reason=not_configured，旧 regional.* 不存在也不造）。32 必须 "all checks passed"。
"""
import importlib.util, json, os, pathlib, shutil, subprocess, sys, tempfile, unittest

import test_30_states_production as t30   # 复用同一套合成 work 树 fixtures

REPO = pathlib.Path(__file__).resolve().parents[1]
S04 = REPO / "scripts" / "04_ancestry_pca.sh"
S04B = REPO / "scripts" / "04b_ancestry_summary.py"
S30 = REPO / "scripts" / "30_build_report_data.py"
S31 = REPO / "scripts" / "31_html_report.py"
S32 = REPO / "scripts" / "32_check_report.py"

_SKIP = t30._SKIP
if not _SKIP and shutil.which("node") is None:
    _SKIP = "node not on PATH (32's runtime figure check shells out to it)"

# 参考面板：CHB/CHS/GBR 各 2 人——min_group_n 下也能形成可排名的组。
_REF_ROWS = [
    ("F1", "CHB.1", "EAS", "CHB"), ("F2", "CHB.2", "EAS", "CHB"),
    ("F3", "CHS.1", "EAS", "CHS"), ("F4", "CHS.2", "EAS", "CHS"),
    ("F5", "GBR.1", "EUR", "GBR"), ("F6", "GBR.2", "EUR", "GBR"),
]
_KG_SSCORE = "\n".join(["#FID\tIID\tSuperPop\tPopulation\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG"] + [
    f"{f}\t{i}\t{sp}\t{pop}\t{0.10 + 0.01 * k:.4f}\t{0.01 * k:.4f}\t0.0000\t0.0000"
    for k, (f, i, sp, pop) in enumerate(_REF_ROWS)]) + "\n"
_REG_SSCORE = ("#FID\tIID\tPopulation\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG\n"
               "F5\tGBR.1\tGBR\t-0.1000\t-0.0200\t0.0000\t0.0000\n"
               "F6\tGBR.2\tGBR\t-0.1100\t-0.0300\t0.0000\t0.0100\n")
_TGT_SSCORE = ("#FID\tIID\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG\n"
               "TEST\tTESTSAMPLE\t0.0900\t0.0100\t0.0000\t0.0000\n")
_SMISS_HDR = "#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS"
_SMISS_REF = "\n".join([_SMISS_HDR] + [f"{f}\t{i}\t0\t100000\t0" for f, i, _, _ in _REF_ROWS]) + "\n"
_SMISS_TGT = _SMISS_HDR + "\nTEST\tTESTSAMPLE\t0\t100000\t0\n"


def _plink_stub(td):
    """一个替身伺候 04（make-pgen/剪枝/PCA/投影）与 04b（smiss 重算）两类调用。"""
    stub = pathlib.Path(td) / "bin" / "plink2"
    stub.parent.mkdir(parents=True, exist_ok=True)
    stub.write_text(
        "#!/bin/bash\n"
        "out=\"\"; prev=\"\"\n"
        "for a in \"$@\"; do [ \"$prev\" = \"--out\" ] && out=\"$a\"; prev=\"$a\"; done\n"
        "base=\"${out##*/}\"\n"
        "case \"$*\" in\n"
        "  *--indep-pairwise*) printf 'rs1\\nrs2\\n' > \"$out.prune.in\"; : > \"$out.prune.out\"; exit 0;;\n"
        "esac\n"
        "case \"$base\" in\n"
        "  kg.common) for e in pgen psam pvar log; do : > \"$out.$e\"; done; exit 0;;\n"
        "  kg.pca|regional.pca) : > \"$out.afreq\"; : > \"$out.eigenvec.allele\"; : > \"$out.log\"; exit 0;;\n"
        "  kg.proj) printf '" + _KG_SSCORE.replace("'", "'\\''") + "' > \"$out.sscore\"; exit 0;;\n"
        "  regional.proj) printf '" + _REG_SSCORE.replace("'", "'\\''") + "' > \"$out.sscore\"; exit 0;;\n"
        "  target.proj|regional.target.proj) printf '" + _TGT_SSCORE.replace("'", "'\\''") + "' > \"$out.sscore\"; exit 0;;\n"
        "esac\n"
        "case \"$out\" in\n"
        "  *.ref) printf '" + _SMISS_REF.replace("'", "'\\''") + "' > \"$out.smiss\"; exit 0;;\n"
        "  *.target) printf '" + _SMISS_TGT.replace("'", "'\\''") + "' > \"$out.smiss\"; exit 0;;\n"
        "esac\n"
        ": > \"$out.log\"\n"
        "exit 0\n", encoding="utf-8")
    stub.chmod(0o755)
    return stub


_AADR_PARAMS = {"ancestry_min_call_rate_modern": 0.95, "aadr_min_call_rate_ancient": 0.50,
               "ancestry_min_projection_snps": 10000, "ancestry_min_group_n": 2}


def _aadR_products(W):
    """11_aadr 的最小合法产物：09b 形状的 summary + 准入匹配的 manifest + 投影表。"""
    a = W / "11_aadr"
    a.mkdir(parents=True, exist_ok=True)
    (a / "prune.prune.in").write_text("rs1\nrs2\n", encoding="utf-8")
    import hashlib
    prune_sha = hashlib.sha256((a / "prune.prune.in").read_bytes()).hexdigest()[:12]
    (a / "proj_annotated.tsv").write_text("\n".join([
        "iid\tlabel\tkind\tPC1_AVG\tPC2_AVG\tdate\td",
        "ANC1\tChina_X\tancient\t0.05\t0.01\t3000\t0.012",
        "MOD1\tHan\tmodern\t0.10\t0.02\t0\t0.031",
        "TESTSAMPLE\ttarget\ttarget\t0.09\t0.01\t0\t0.0"]) + "\n", encoding="utf-8")
    (a / "near_ancient.tsv").write_text(
        "iid\tlabel\td\tdate\nANC1\tChina_X\t0.012\t3000\n", encoding="utf-8")
    _rec = lambda rid, sid, kind, date: {"record_id": rid, "iid": rid, "master_id": rid,
                                         "dataset": "AADR", "reference_release": "v66",
                                         "kind": kind, "source_population_id": sid, "label": sid,
                                         "locality": "Somewhere", "latitude": 30.0, "longitude": 110.0,
                                         "location_precision": "region", "location_source": "anno",
                                         "date_mean_bp": date, "date_sd_bp": None,
                                         "date_min_bp": None, "date_max_bp": None,
                                         "date_basis": None, "date_raw": None,
                                         "y_hg_raw": "", "mt_hg_raw": "", "publication": None}
    # n>=2 on purpose: real 09b drops groups below min_group_n, and 32 (rightly) rejects an
    # affinity table that still carries such a row -- a 1-member group is not admissible product.
    _grp = lambda gid, kind, d, mid: {"group_id": gid, "label": gid, "kind": kind, "n": 2,
                                       "rank": 1, "small_group": False, "distance_mean": d,
                                       "member_ids": [mid]}
    (a / "summary.json").write_text(json.dumps({
        "analysis_id": "aadr-human-origins", "dataset": "AADR", "reference_release": "v66",
        "scope": "global", "state": "ok", "reason_code": "", "components": [1, 2],
        "metric": "mean_individual_pc_distance", "thresholds": {}, "counts": {"kept": 2},
        "target": {"sample_id": "TESTSAMPLE", "pcs": [0.09, 0.01], "call_rate": 1.0,
                   "n_called_snps": 100000},
        "records": [_rec("ANC1", "China_X", "ancient", 3000), _rec("MOD1", "Han", "modern", 0)],
        "groups": {"modern": [_grp("Han", "modern", 0.031, "MOD1")],
                   "ancient": [_grp("China_X", "ancient", 0.012, "ANC1")]},
        "sources": []}, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    (a / "manifest.json").write_text(json.dumps({
        "schema_version": "1", "sample_id": "TESTSAMPLE", "analysis_id": "09b-aadr-summary",
        "state": "ok", "reason_code": "", "build": "GRCh37", "reference_release": "v66",
        "parameters": {"min_call_rate_modern": _AADR_PARAMS["ancestry_min_call_rate_modern"],
                       "min_call_rate_ancient": _AADR_PARAMS["aadr_min_call_rate_ancient"],
                       "min_projection_snps": _AADR_PARAMS["ancestry_min_projection_snps"],
                       "min_group_n": _AADR_PARAMS["ancestry_min_group_n"],
                       "ancient_prefix": ["China_"], "modern_groups": ["Han"],
                       "prune_sha": prune_sha, "prune_sites": 2,
                       # 复审 §3.2 P0 指纹绑定：30 还比对 AADR 前缀/注释与 08/09 阶段产物内容。
                       # 本合成目录没有 proj.sscore/samples.tsv → sha 为 ""（30 的比对侧对缺失
                       # 文件同样算 ""）；前缀取该配置下 wgsconfig 的默认解析路径。
                       "aadr_prefix": str((a.parent.parent / "data" / "ref" / "aadr"
                                           / "v66.p1_HO.aadr.patch.PUB").resolve()),
                       "annotation": "", "proj_sha": "", "samples_sha": ""},
        "tool_versions": {}, "input_fingerprints": {}, "outputs": ["summary.json"]}) + "\n",
        encoding="utf-8")


def _figure_inputs(W):
    """32 的运行时图检查需要的最小输入：prs/behaviour/afam/hla_disease/candidate。"""
    pr = W / "07_prs"
    pr.mkdir(parents=True, exist_ok=True)
    (pr / "prs_wgs.tsv").write_text("\n".join([
        "trait\tscore\tpct_sub\tpct_EAS\tcoverage_pct\tz_vs_EAS",
        "Asthma\t0.5\t88.8\t91.2\t99.0\t1.5",
        "Height\t-0.2\t40.0\t35.5\t98.0\t-0.4"]) + "\n", encoding="utf-8")
    b = W / "20_behaviour"
    b.mkdir(parents=True, exist_ok=True)
    (b / "behaviour_prs.tsv").write_text("\n".join([
        "zh\ten\tpanel\tpct_EAS\tz\tcoverage",
        "近视\tMyopia\tEAS\t72.0\t0.8\t97",
        "失眠\tInsomnia\tEUR\t55.0\t0.1\t96"]) + "\n", encoding="utf-8")
    (b / "candidate_genes.tsv").write_text("\n".join([
        "gene\trsid\tvariant\tvariant_en\tgenotype\tpopular_claim\tpopular_claim_en\twhat_evidence_supports\tevidence_en",
        "DEC2\t-\t短睡眠\tshort sleeper\thet\t少睡\tless sleep\t对照研究\tassociation study"]) + "\n",
        encoding="utf-8")
    hd = W / "19_hla_disease"
    hd.mkdir(parents=True, exist_ok=True)
    (hd / "hla_disease.tsv").write_text("\n".join([
        "condition\tcondition_en\trequirement\trequirement_en\tcarried\tnote\tnote_en",
        "发作性睡病\tNarcolepsy\tDQB1*06:02\tDQB1*06:02\tunknown\t未分型\tnot typed"]) + "\n",
        encoding="utf-8")
    ar = W / "15_archaic"
    ar.mkdir(parents=True, exist_ok=True)
    (ar / "gene_families.tsv").write_text("\n".join([
        "family\tfamily_en\tpct\tcarried\ttotal\texamples",
        "嗅觉受体\tOlfactory receptors\t30.0\t3\t10\tOR4C1,OR4C2,OR4C3"]) + "\n", encoding="utf-8")
    (ar / "gene_summary.txt").write_text("n_genes\t123\nnull_mean\t45\n", encoding="utf-8")


def _build(td, regional):
    """合成 work 树（30 的其余输入沿用 t30 的 fixture）+ 替身 plink2 + 两份配置共用的骨架。"""
    over = {"plink2": str(_plink_stub(td)),
            "aadr_modern": ["Han"], "aadr_ancient_prefix": ["China_"], **_AADR_PARAMS}
    if regional:
        over.update({"ref_superpop": "EUR", "ref_subpops": ["GBR"]})
    cfg = t30._write_config(td, **over)
    W = t30._write_base(td, y=True, mt=True, regional=regional)
    ref = pathlib.Path(td) / "ref" / "all_phase3"
    ref.parent.mkdir(parents=True, exist_ok=True)
    (ref.parent / "all_phase3.psam").write_text(
        "#FID\tIID\tSEX\tSuperPop\tPopulation\n"
        + "\n".join(f"{f}\t{i}\t1\t{sp}\t{pop}" for f, i, sp, pop in _REF_ROWS) + "\n",
        encoding="utf-8")
    _cfg_patch(cfg, {"kg_pfile": str(ref)})
    c = W / "02_complete"
    c.mkdir(parents=True, exist_ok=True)
    for e in ("pgen", "psam", "pvar"):
        (c / f"TESTSAMPLE.1kg.{e}").touch()
    ytree = pathlib.Path(td) / "work" / "data" / "ref" / "ytree"
    ytree.mkdir(parents=True, exist_ok=True)
    (ytree / "current_version.txt").write_text("14.06.0\n", encoding="utf-8")
    _aadR_products(W)
    _figure_inputs(W)
    return cfg, W


def _cfg_patch(path, over):
    """改写既有 config.yaml 的键（保持 t30._write_config 的简单行格式）。"""
    lines = pathlib.Path(path).read_text(encoding="utf-8").splitlines()
    keep = [l for l in lines if not any(l.startswith(f"{k}:") for k in over)]
    keep += [f"{k}: {json.dumps(v)}" for k, v in over.items()]
    pathlib.Path(path).write_text("\n".join(keep) + "\n", encoding="utf-8")


def _run(cfg, script, *args):
    return subprocess.run([sys.executable, str(script), *args], capture_output=True, text=True,
                         env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=600)


@unittest.skipIf(_SKIP, _SKIP)
class TestDualConfigEndToEnd(unittest.TestCase):

    def _chain(self, td, regional):
        cfg, W = _build(td, regional)
        r04 = subprocess.run(["bash", str(S04)], capture_output=True, text=True,
                             env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=300)
        self.assertEqual(r04.returncode, 0, r04.stderr)
        r4b = _run(cfg, S04B)
        self.assertEqual(r4b.returncode, 0, r4b.stderr)
        r30 = _run(cfg, S30)
        self.assertEqual(r30.returncode, 0, r30.stderr)
        r31 = _run(cfg, S31)
        self.assertEqual(r31.returncode, 0, r31.stderr)
        html = pathlib.Path(td) / "work" / "report" / "report.html"
        self.assertTrue(html.exists(), "31 必须写出 report.html")
        r32 = subprocess.run([sys.executable, str(S32),
                              str(pathlib.Path(td) / "work/wgs/report_data.json"), str(html)],
                             capture_output=True, text=True, timeout=300)
        self.assertEqual(r32.returncode, 0, r32.stderr)
        self.assertIn("all checks passed", r32.stdout)
        return json.loads((pathlib.Path(td) / "work/wgs/report_data.json").read_text(encoding="utf-8"))

    def test_config_a_eur_regional_end_to_end(self):
        with tempfile.TemporaryDirectory() as td:
            d = self._chain(td, regional=True)
            self.assertIsNotNone(d["pca_eas"], "EUR 配置必须交付区域视图")
            self.assertEqual(d["pca_eas_reason"], "")
            self.assertEqual(d["pop"]["superpop"], "EUR")
            self.assertEqual(d["pop"]["n_super"], 2, "数的是 psam 里的 EUR 行")
            self.assertIsNotNone(d["near_eas"], "真实 04b 的 kg-regional 分析要进 near_eas")
            self.assertIn("GBR", d["near_eas"])
            self.assertIsNotNone(d["near_global"], "global 空间照常交付")

    def test_config_b_global_only_end_to_end(self):
        with tempfile.TemporaryDirectory() as td:
            d = self._chain(td, regional=False)
            self.assertIsNone(d["pca_eas"], "global-only 是合法形态，不造区域视图")
            self.assertEqual(d["pca_eas_reason"], "not_configured")
            self.assertIsNone(d["near_eas"])
            self.assertIsNotNone(d["near_global"], "global-only 的近邻来自真实 04b 的 kg-global")
            sec = {s["id"]: s for s in d["sections"]}
            self.assertEqual(sec["ancestry"]["status"], "ok")
            self.assertFalse(list((pathlib.Path(td) / "work/wgs/04_ancestry").glob("regional.*")),
                            "真实 04 在未配置时不得产出 regional.*")


if __name__ == "__main__":
    unittest.main(verbosity=2)
