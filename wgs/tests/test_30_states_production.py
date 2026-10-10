"""P1（复审 AN0/AN2/AN5）验收：global-only 与缺谱系源文件按状态交付，对**生产 30** 断言。

复审原话："30 在 manifest 准入前无条件读 regional.* 和 prune.regional，缺文件就失败，有旧文
件就读取……还在结构化空态之外必读 Y 文本/SNP 表、haplogrep3.txt 与 mt_heteroplasmy.tsv。"
本文件用真实 30_build_report_data.py（$WGS_CONFIG 指向临时 config.yaml，work 树全合成）跑
四种形态：仅 global 的新目录（无区域/无 Y/无 mt）、EUR 区域、缺 Y 有 mt、有 Y 无 mt——各自
按状态交付并 exit 0，绝不在读取处崩溃，也不拿旧文件冒充；另断言 near_* 来自 04b 的同一均
值口径（30 不再有第二套质心排名）与 global-only 下旧 regional.* 的 stale 免疫。

30 需要 pandas/numpy/yaml 与 PATH 里的 bcftools；缺任一整文件跳过（本地跳过、服务器执行，
与 test_09b_production 同一约定）。
"""
import hashlib, importlib.util, json, os, pathlib, shutil, subprocess, sys, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "30_build_report_data.py"
_SKIP = ""
for _mod in ("numpy", "pandas", "yaml"):
    if importlib.util.find_spec(_mod) is None:
        _SKIP = f"{_mod} not available (the production script needs it)"
        break
if not _SKIP and shutil.which("bcftools") is None:
    _SKIP = "bcftools not on PATH (the production script shells out to it)"


def _write_config(td, **over):
    cfg = {"sample_id": "TESTSAMPLE", "work_dir": str(pathlib.Path(td) / "work")}
    cfg.update(over)
    p = pathlib.Path(td) / "config.yaml"
    p.write_text("".join(f"{k}: {json.dumps(v)}\n" for k, v in cfg.items()), encoding="utf-8")
    return p


def _vcf_text():
    return ("##fileformat=VCFv4.2\n##contig=<ID=1>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTSAMPLE\n"
            "1\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\n")


def _gzip_vcf(out):
    r = subprocess.run(["bcftools", "view", "-Oz", "-o", str(out)], input=_vcf_text(),
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    subprocess.run(["bcftools", "index", str(out)], check=True, capture_output=True)


_SSCORE_HDR = "#FID\tIID\tSuperPop\tPopulation\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG"
_KG_ROWS = [
    "FAM1\tHAN1\tEAS\tCHB\t0.10\t0.01\t0.00\t0.00",
    "FAM2\tHAN2\tEAS\tCHS\t0.12\t0.02\t0.01\t0.00",
    "FAM3\tEUR1\tEUR\tGBR\t-0.10\t-0.02\t0.00\t0.01",
]
_TGT_SSCORE = ("#FID\tIID\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG\n"
               "TEST\tTESTSAMPLE\t0.09\t0.01\t0.00\t0.00\n")
_REG_SSCORE = ("#FID\tIID\tPopulation\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG\n"
               "FAM3\tEUR1\tGBR\t-0.10\t-0.02\t0.00\t0.01\n")

_GROUPS_GLOBAL = [
    {"group_id": "CHB", "label": "CHB", "kind": "modern", "n": 1, "rank": 1, "small_group": False,
     "distance_mean": 0.012, "member_ids": ["HAN1"]},
    {"group_id": "CHS", "label": "CHS", "kind": "modern", "n": 1, "rank": 2, "small_group": False,
     "distance_mean": 0.031, "member_ids": ["HAN2"]},
]
_GROUPS_REGIONAL = [
    {"group_id": "GBR", "label": "GBR", "kind": "modern", "n": 1, "rank": 1, "small_group": False,
     "distance_mean": 0.020, "member_ids": ["EUR1"]},
]


def _analysis(aid, scope, groups):
    return {"analysis_id": aid, "dataset": "1000G", "reference_release": "phase3", "scope": scope,
            "state": "ok", "reason_code": "", "components": [1, 2, 3, 4],
            "metric": "mean_individual_pc_distance", "thresholds": {}, "counts": {}, "target": {},
            "groups": groups, "nearest_individuals": {groups[0]["label"]: 2}}


_AADR_GROUP = {"group_id": "Han", "label": "Han", "kind": "modern", "n": 2, "rank": 1,
               "small_group": False, "distance_mean": 0.010, "member_ids": ["HAN1", "HAN2"]}
_AADR_ANCIENT_GROUP = {"group_id": "China_MLBA", "label": "China_MLBA", "kind": "ancient", "n": 2,
                       "rank": 1, "small_group": False, "distance_mean": 0.020,
                       "member_ids": ["ANC1", "ANC2"]}


def _add_aadr_module(td, aadr_cfg, groups):
    """11_aadr：09b 形状的 summary.json + 参数与配置一致的 ok manifest。

    没有真实的 08/09 产物文件——30 的期望指纹从**当前文件**重算，文件不存在即 ""/0，
    manifest 按同样值记录（准入比较的是两侧一致，不是文件存在）。groups 控制 modern/ancient
    两侧的群体汇总；modern 为空即复现"AADR 只选古代/现代侧全被筛掉"。"""
    w = pathlib.Path(td) / "work/wgs/11_aadr"
    w.mkdir(parents=True, exist_ok=True)
    (w / "summary.json").write_text(json.dumps({
        "schema_version": 1, "analysis_id": "aadr-human-origins", "dataset": "AADR",
        "reference_release": "v66.p1", "scope": "global", "state": "ok", "reason_code": "",
        "components": [1, 2, 3, 4], "metric": "mean_individual_pc_distance",
        "thresholds": {}, "counts": {}, "target": None, "records": [],
        "groups": groups, "sources": []}, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    (w / "manifest.json").write_text(json.dumps({
        "schema_version": "1", "sample_id": "TESTSAMPLE", "analysis_id": "09b-aadr-summary",
        "state": "ok", "reason_code": "", "build": "GRCh37", "reference_release": "v66.p1",
        "parameters": {"min_call_rate_modern": 0.95, "min_call_rate_ancient": 0.5,
                       "min_projection_snps": 10000, "min_group_n": 2,
                       "ancient_prefix": list(aadr_cfg.get("aadr_ancient_prefix", [])),
                       "modern_groups": list(aadr_cfg.get("aadr_modern", [])),
                       "prune_sha": "", "prune_sites": 0,
                       "aadr_prefix": str(aadr_cfg["aadr_prefix"]), "annotation": "",
                       "proj_sha": "", "samples_sha": ""},
        "tool_versions": {}, "input_fingerprints": {}, "outputs": []}) + "\n", encoding="utf-8")
    return w


def _write_base(td, y=True, mt=True, regional=False):
    """合成 30 的必读输入。y/mt/regional 决定对应形态的文件是否存在。"""
    W = pathlib.Path(td) / "work" / "wgs"
    (W / "01_qc").mkdir(parents=True)
    rows = [f"{c}\t1000000\t900000\t30.0" for c in [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]]
    (W / "01_qc/depth.mosdepth.summary.txt").write_text(
        "chrom\tlength\tbases\tmean\n" + "\n".join(rows) + "\ntotal\t2500000000\t2250000000\t30.0\n",
        encoding="utf-8")
    stats = ("SN\t0\tnumber of records:\t100\nSN\t0\tnumber of SNPs:\t90\nSN\t0\tnumber of indels:\t10\n"
             "TSTV\t0\t80\t20\t4.00\t4.00\t0.00\nPSC\t0\t0\t100\t30\t40\t30\t0\n")
    (W / "01_qc/stats.pass.txt").write_text(stats, encoding="utf-8")
    (W / "01_qc/stats.norm.txt").write_text(stats, encoding="utf-8")
    (W / "00_input").mkdir(parents=True)
    (W / "00_input/callable.bed").write_text("1\t1000\t2000\n", encoding="utf-8")
    _gzip_vcf(W / "00_input/target.pass.vcf.gz")
    # --- 谱系形态
    if y:
        (W / "03_haplo").mkdir(parents=True, exist_ok=True)
        (W / "03_haplo/y_haplogroup_yfull.txt").write_text(
            "  Z12345 der= 5 anc= 0 n/a= 1 formed= 5000 tmrca= 3000\n", encoding="utf-8")
        (W / "03_haplo/y_terminal_snps.tsv").write_text(
            "branch\tsnp\tdepth\tn_der\tn_anc\tstate\nZ12345\tZ12345\t20\t5\t0\tder\n", encoding="utf-8")
    if mt:
        (W / "03_haplo").mkdir(parents=True, exist_ok=True)
        (W / "03_haplo/haplogrep3.txt").write_text(
            "Haplogroup\tQuality\tFound_Polys\tRemaining_Polys\tNot_Found_Polys\n"
            "A1\t0.92\t195 8271\t143\t152 16390\n", encoding="utf-8")
        (W / "03_haplo/mt_heteroplasmy.tsv").write_text(
            "pos\talt\taf\tdepth\n310\tT\t0.02\t500\n", encoding="utf-8")
    # --- 祖源：global 空间 + （可选）regional 空间，04b 的 summary/manifest
    a = W / "04_ancestry"
    a.mkdir(parents=True)
    (a / "kg.proj.sscore").write_text("\n".join([_SSCORE_HDR] + _KG_ROWS) + "\n", encoding="utf-8")
    (a / "target.proj.sscore").write_text(_TGT_SSCORE, encoding="utf-8")
    (a / "prune.prune.in").write_text("rs1\nrs2\n", encoding="utf-8")
    (a / "summary.txt").write_text("audit log\n", encoding="utf-8")
    analyses = [_analysis("kg-global", "global", _GROUPS_GLOBAL)]
    regional_enabled = False
    if regional:
        regional_enabled = True
        (a / "regional.proj.sscore").write_text(_REG_SSCORE, encoding="utf-8")
        (a / "regional.target.proj.sscore").write_text(_TGT_SSCORE, encoding="utf-8")
        (a / "prune.regional.prune.in").write_text("rs1\nrs2\n", encoding="utf-8")
        analyses.append(_analysis("kg-regional", "regional", _GROUPS_REGIONAL))
    (a / "summary.json").write_text(json.dumps({
        "schema_version": 1, "default_analysis_id": "kg-global", "analyses": analyses,
        "local": {}, "state": "ok", "reason_code": ""}, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8")
    # 04b manifest 的参数必须与本次配置一致，30 的 expected_parameters 准入才放行。
    # 复审 §3.2 P0 指纹绑定：30 还比对阶段产物内容（sscore/prune 的 sha，与真实 04b 写侧
    # 同一算法）与 QC 门槛（本配置未覆盖 → wgsconfig 默认 0.95/0.95/10000/2）。
    superpop = "EUR" if regional else "EAS"
    _sha = lambda p: hashlib.sha256((a / p).read_bytes()).hexdigest()[:12]
    (a / "manifest.json").write_text(json.dumps({
        "schema_version": "1", "sample_id": "TESTSAMPLE", "analysis_id": "04b-ancestry-summary",
        "state": "ok", "reason_code": "", "build": "GRCh37", "reference_release": "",
        "parameters": {"regional_enabled": regional_enabled, "superpop": superpop,
                       "subpops": ["CHB", "CHS"],
                       "kg_sscore_sha": _sha("kg.proj.sscore"),
                       "target_sscore_sha": _sha("target.proj.sscore"),
                       "prune_sha": _sha("prune.prune.in"), "prune_sites": 2,
                       "min_call_rate_modern": 0.95, "min_call_rate_target": 0.95,
                       "min_projection_snps": 10000, "min_group_n": 2},
        "tool_versions": {}, "input_fingerprints": {}, "outputs": []}) + "\n", encoding="utf-8")
    # --- 其余必读
    ref = pathlib.Path(td) / "work" / "data" / "ref"
    ref.mkdir(parents=True)
    (ref / "all_phase3.psam").write_text(
        "#FID\tIID\tSEX\tSuperPop\tPopulation\nFAM1\tHAN1\t1\tEAS\tCHB\n"
        "FAM2\tHAN2\t2\tEAS\tCHS\nFAM3\tEUR1\t1\tEUR\tGBR\n", encoding="utf-8")
    c = W / "05_clinvar"
    c.mkdir(parents=True)
    (c / "clinvar_all_hits.tsv").write_text(
        "1\t1000\t12345\tA\tG\t50\tGENE:gene\tBenign\tno_assertion\tDisease\t"
        "some\t.\t.\t.\t0|1\t30\t45\t20\n", encoding="utf-8")
    (c / "lof_table.tsv").write_text("gene\tclass\tzyg\nGENE\tpLoF\thet\n", encoding="utf-8")
    (c / "lof_rare_final.tsv").write_text("gene\tzyg\toe_lof_upper\nGENE2\thom/hemi\t0.3\n", encoding="utf-8")
    pr = W / "07_prs"
    pr.mkdir(parents=True)
    (pr / "prs_wgs.tsv").write_text("score\tpct_sub\n0.5\t88.8\n", encoding="utf-8")
    return W


def _run_30(cfg):
    return subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                         env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=300)


def _data(td):
    return json.loads((pathlib.Path(td) / "work/wgs/report_data.json").read_text(encoding="utf-8"))


@unittest.skipIf(_SKIP, _SKIP)
class TestGlobalOnlyFresh(unittest.TestCase):
    """仅 global 的新目录：无区域/无 Y/无 mt——30 按状态交付并 exit 0。"""

    def test_global_only_directory_delivers_structured_states(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)                     # 无 ref_superpop → REGIONAL_ENABLED=False
            _write_base(td, y=False, mt=False)
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertIsNone(d["pca_eas"], "未配置区域时不得读 regional.*，也不能造空图")
            self.assertEqual(d["pca_eas_reason"], "not_configured")
            self.assertIsNone(d["anc"]["n_eas"])
            self.assertEqual(d["y_state"], "unavailable")
            self.assertEqual(d["y_snps"], [])
            self.assertEqual(d["mt"]["state"], "unavailable")
            self.assertEqual(d["mt_het"], [])
            self.assertEqual(d["near_global"], {"CHB": 0.012, "CHS": 0.031},
                             "near_global 来自 04b 的均值距离，不是 30 的质心排名")
            self.assertIsNone(d["near_eas"], "无区域分析时 near_eas 是结构化空态")
            self.assertEqual(d["knn_global"], {"CHB": 2}, "knn 转录 04b 记录的近个体人群计数")
            sec = {s["id"]: s for s in d["sections"]}
            self.assertEqual(sec["ancestry"]["status"], "ok",
                             "global-only 是合法交付形态，祖源节不该变 unavailable")
            self.assertIn("global-only", r.stderr)

    def test_stale_regional_products_are_ignored_when_disabled(self):
        """复审点名"有旧文件就读取"：换配置后留下的 regional.* 不得冒充本次结果。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_base(td, y=False, mt=False)
            a = pathlib.Path(td) / "work/wgs/04_ancestry"
            (a / "regional.proj.sscore").write_text(_REG_SSCORE, encoding="utf-8")
            (a / "regional.target.proj.sscore").write_text(_TGT_SSCORE, encoding="utf-8")
            (a / "prune.regional.prune.in").write_text("rs1\nrs2\n", encoding="utf-8")
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertIsNone(d["pca_eas"])
            self.assertEqual(d["pca_eas_reason"], "not_configured")
            self.assertIsNone(d["anc"]["n_eas"])

    def test_configured_regional_with_missing_products_is_not_a_crash(self):
        """配置了 ref_superpop 但 04 没跑：区域视图不可用（缺产物），30 照常交付。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td, **{"ref_superpop": "EUR"})
            _write_base(td, y=True, mt=True)            # 未写 regional.*，manifest 也不匹配？
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertIsNone(d["pca_eas"])
            self.assertEqual(d["pca_eas_reason"], "missing_products")
            # manifest 的 regional_enabled=False 与配置不符 → kg 分析整体 stale，near_* 也为空
            self.assertIsNone(d["near_global"])
            self.assertIn("unavailable", [s["status"] for s in d["sections"] if s["id"] == "ancestry"])


@unittest.skipIf(_SKIP, _SKIP)
class TestEurRegional(unittest.TestCase):
    """EUR 区域：pca_eas/近邻/标签按配置交付，百分位计数数的是配置的超群。"""

    def test_eur_regional_space_is_delivered(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td, **{"ref_superpop": "EUR"})
            _write_base(td, y=True, mt=True, regional=True)
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertIsNotNone(d["pca_eas"])
            self.assertEqual(d["pca_eas"]["pts"], [["GBR", -0.1, -0.02]])
            self.assertEqual(d["anc"]["n_eas"], 2)
            self.assertEqual(d["pop"]["superpop"], "EUR")
            self.assertEqual(d["pop"]["n_super"], 1, "数的是 EUR 行，不再写死 EAS")
            self.assertEqual(d["pop"]["n_sub"], 2)
            self.assertEqual(d["near_eas"], {"GBR": 0.02})
            self.assertEqual(d["near_global"], {"CHB": 0.012, "CHS": 0.031})


@unittest.skipIf(_SKIP, _SKIP)
class TestDefaultAnalysisSelection(unittest.TestCase):
    """复审 §3.2 P1 AN2/AN5：默认主分析按 §7.3 契约显式选择——有效现代 AADR → 已配置区域
    1000G → 全球 1000G。

    旧实现按数据顺序取第一个 state=ok：04b 的 spaces 把 global 排在 regional 前，于是
    AADR 不可用但区域已配置且 ok 时默认仍是 global；AADR 只选古代（现代侧没有群体汇总，
    无论因只配了古代前缀还是现代组全被门槛筛掉）时也冒称有效现代主空间。"""

    def _aadr_cfg(self, td, **over):
        cfg = {"aadr_prefix": str(pathlib.Path(td) / "work/data/ref/aadr/panel"),
               "aadr_modern": ["Han"], "aadr_ancient_prefix": ["China_"]}
        cfg.update(over)
        return cfg

    def test_valid_modern_aadr_wins_over_regional_and_global(self):
        """守卫：AADR 现代侧有群体汇总时仍是默认——修选择顺序不得把有效 AADR 挤掉。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td, **{"ref_superpop": "EUR"}, **self._aadr_cfg(td))
            _write_base(td, y=False, mt=False, regional=True)
            _add_aadr_module(td, self._aadr_cfg(td),
                             {"modern": [dict(_AADR_GROUP)], "ancient": []})
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(_data(td)["ancestry"]["default_analysis_id"], "aadr-human-origins")

    def test_regional_beats_global_when_aadr_is_not_delivered(self):
        """AADR 没有产物（未配置）时：已配置且 ok 的区域优先于全球——不是 04b 的书写顺序。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td, **{"ref_superpop": "EUR"})   # 无 aadr_* → 11_aadr 不存在
            _write_base(td, y=False, mt=False, regional=True)
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            anc = _data(td)["ancestry"]
            ids = {a["analysis_id"]: a["state"] for a in anc["analyses"]}
            self.assertEqual(ids, {"kg-global": "ok", "kg-regional": "ok"})
            self.assertEqual(anc["default_analysis_id"], "kg-regional",
                             "AADR 不可用时默认应是已配置区域，不是排在前的 global")

    def test_ancient_only_aadr_is_not_a_valid_modern_space(self):
        """AADR ok 但 groups.modern 为空（只配了古代前缀）：不得越过区域冒称默认现代空间。"""
        with tempfile.TemporaryDirectory() as td:
            ac = self._aadr_cfg(td, aadr_modern=[])
            cfg = _write_config(td, **{"ref_superpop": "EUR"}, **ac)
            _write_base(td, y=False, mt=False, regional=True)
            _add_aadr_module(td, ac, {"modern": [], "ancient": [dict(_AADR_ANCIENT_GROUP)]})
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            anc = _data(td)["ancestry"]
            self.assertEqual(anc["default_analysis_id"], "kg-regional",
                             "仅古代的 AADR 不是有效现代空间；区域与全球都在时选区域")

    def test_ancient_only_aadr_with_global_only_falls_to_global(self):
        """仅古代 AADR + 未配置区域：退到全球，仍不得选 AADR。"""
        with tempfile.TemporaryDirectory() as td:
            ac = self._aadr_cfg(td, aadr_modern=[])
            cfg = _write_config(td, **ac)                          # 无 ref_superpop → global-only
            _write_base(td, y=False, mt=False, regional=False)
            _add_aadr_module(td, ac, {"modern": [], "ancient": [dict(_AADR_ANCIENT_GROUP)]})
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(_data(td)["ancestry"]["default_analysis_id"], "kg-global")

    def test_global_only_without_aadr_keeps_global_default(self):
        """守卫：无 AADR、无区域的老形态默认仍是全球。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_base(td, y=False, mt=False, regional=False)
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(_data(td)["ancestry"]["default_analysis_id"], "kg-global")


@unittest.skipIf(_SKIP, _SKIP)
class TestLineageAbsenceShapes(unittest.TestCase):
    """缺 Y 有 mt / 有 Y 无 mt：各自按状态交付，缺失模块不带走整份报告。"""

    def test_no_y_with_mt(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_base(td, y=False, mt=True)
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertEqual(d["y_state"], "unavailable")
            self.assertEqual(d["y_terminal"], "")
            self.assertEqual(d["mt"]["hg"], "A1")
            self.assertEqual(len(d["mt_het"]), 1)

    def test_y_without_mt(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_base(td, y=True, mt=False)
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertEqual(d["y_state"], "ok")
            self.assertEqual(d["y_terminal"], "Z12345")
            self.assertEqual(d["mt"]["state"], "unavailable")
            self.assertEqual(d["mt_het"], [])


if __name__ == "__main__":
    unittest.main(verbosity=2)


@unittest.skipIf(_SKIP, _SKIP)
class TestOptionalProductAdmission(unittest.TestCase):
    """复审 H6/AN5 残留：f3_stats.json / lineage_history.json 此前按"文件存在"进场，绕过
    生产者的 manifest——禁用/失败/换样本的旧文件照样投递。现在 30 走 analysis_state 准入。"""

    @staticmethod
    def _write_unadmitted_products(W):
        f3 = W / "04_ancestry/f3"
        f3.mkdir(parents=True, exist_ok=True)
        (f3 / "f3_stats.json").write_text(json.dumps(
            {"sample_id": "TESTSAMPLE", "estimator": "site-mean-corrected", "block_mb": 5,
             "modern": {"value": 0.01, "se": 0.001}, "ancient": None}) + "\n", encoding="utf-8")
        (W / "03_haplo/lineage_history.json").write_text(json.dumps(
            {"y": {"state": "ok", "reported_hg": "Z9"}, "mt": None}) + "\n", encoding="utf-8")

    def test_files_without_manifests_are_not_delivered(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            self._write_unadmitted_products(_write_base(td, y=True, mt=True))
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertIsNone(d["f3"], "无 28 manifest 的 f3_stats.json 不得进场")
            self.assertEqual(d["lineages"], {"y": None, "mt": None})
            self.assertIn("f3 result not admitted", r.stderr)
            self.assertIn("lineage history not admitted", r.stderr)

    def test_files_with_matching_manifests_are_delivered(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            W = _write_base(td, y=True, mt=True)
            self._write_unadmitted_products(W)
            # 30 的 expected_parameters 与真实 28/09d 写侧一致：估计量串/块/门槛 + 分组面板与
            # prune 的内容 sha（复审 §3.2 P0 指纹绑定）；history_sha 是 LINEAGE_HISTORY_FILE
            # （默认 wgs/panel/lineage_history.json）的内容 sha。
            (W / "04_ancestry/f3/manifest.json").write_text(json.dumps({
                "schema_version": "1", "sample_id": "TESTSAMPLE", "analysis_id": "28-f3-stats",
                "state": "ok", "reason_code": "", "build": "GRCh37", "reference_release": "",
                "parameters": {
                    "estimator": ("site-mean[(pA-pB)(pA-pC) - pA(1-pA)/(nA-1)] "
                                  "+ weighted delete-one-block(5Mb) jackknife SE (unequal blocks)"),
                    "block_mb": 5, "min_group_n": 20,
                    "panel_sha": hashlib.sha256(
                        (REPO / "panel" / "f3_groups.tsv").read_bytes()).hexdigest()[:12],
                    "prune_sha": hashlib.sha256(
                        (W / "04_ancestry/prune.prune.in").read_bytes()).hexdigest()[:12]},
                "tool_versions": {}, "input_fingerprints": {}, "outputs": ["f3_stats.json"]}) + "\n",
                encoding="utf-8")
            (W / "03_haplo/manifest.json").write_text(json.dumps({
                "schema_version": "1", "sample_id": "TESTSAMPLE", "analysis_id": "09d-lineage-history",
                "state": "ok", "reason_code": "", "build": "GRCh37", "reference_release": "",
                # ytree_sha：本夹具的 work 树没有 data/ref/ytree/current_tree.json，30 比对侧
                # 用同一算法算出 ""（复审 §3.2 P0-2b 绑定树文件内容）。
                "parameters": {"history_sha": hashlib.sha256(
                    (REPO / "panel" / "lineage_history.json").read_bytes()).hexdigest()[:12],
                    "ytree_sha": ""},
                "tool_versions": {}, "input_fingerprints": {},
                "outputs": ["lineage_history.json"]}) + "\n", encoding="utf-8")
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertEqual(d["f3"]["estimator"], "site-mean-corrected")
            self.assertEqual(d["lineages"]["y"]["reported_hg"], "Z9")

    def test_stale_sample_manifest_is_rejected(self):
        """换样本留下的旧产物：文件与 manifest 都在，但 sample_id 不符 → 不得投递。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            W = _write_base(td, y=True, mt=True)
            self._write_unadmitted_products(W)
            (W / "04_ancestry/f3/manifest.json").write_text(json.dumps({
                "schema_version": "1", "sample_id": "OTHERSAMPLE", "analysis_id": "28-f3-stats",
                "state": "ok", "reason_code": "", "build": "GRCh37", "reference_release": "",
                "parameters": {}, "tool_versions": {}, "input_fingerprints": {},
                "outputs": ["f3_stats.json"]}) + "\n", encoding="utf-8")
            r = _run_30(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            d = _data(td)
            self.assertIsNone(d["f3"], "他样本的 f3 结果必须判 stale，不得按存在投递")
