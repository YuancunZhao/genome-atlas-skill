"""AN6 增强①回归：1000G 参考空间的人群采样点进入报告数据（panel + 30 组装）。

§3.2 P1 点名"1000G 没有 records"——kg 分析只有人群组，组上没有坐标，地图层无从下点。
修复分两层，都不碰 04b（其 manifest/参数不掺显示层）：
- panel/kg_population_locations.tsv：IGSR phase3 人群描述 + 采样地（region 级），描述含城市
  的用城市坐标、只到国家/区域的用区域质心并在 note 声明，离散人群坐标在采样地不在祖源地；
- 30 组装时把采样点附到 kg 分析的**人群组**（latitude/longitude/location_id/locality/
  location_precision/name_zh——与 AADR 记录同名字段，模板可用同一套几何与精度词表），并写
  group_location_counts。个体记录永远不落坐标：IGSR 不发布个体地理信息。

测试两层都对着生产件：面板文件直接解析校验；`_attach_kg_group_locations` 本体从 30 的源
文本抽出真跑（stdlib-only，本地/服务器都执行）；另加接线断言（30 必须真调用该函数并读该
panel），防"辅助函数留着、组装另写一套"。旧代码上这些用例全部失败（函数/文件/接线不存在）。
"""
import csv, io, pathlib, py_compile, re, subprocess, sys, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
PANEL = REPO / "panel" / "kg_population_locations.tsv"
B30 = REPO / "scripts" / "30_build_report_data.py"

PHASE3_POPS = {
    "CHB", "JPT", "CHS", "CDX", "KHV",          # EAS
    "CEU", "TSI", "FIN", "GBR", "IBS",           # EUR
    "YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB",   # AFR
    "MXL", "PUR", "CLM", "PEL",                  # AMR
    "GIH", "PJL", "BEB", "STU", "ITU",           # SAS
}
_COLS = ["dataset", "source_id", "location_id", "label_zh", "label_en",
         "locality", "latitude", "longitude", "precision", "source_url", "note"]


def _rows():
    lines = [ln for ln in PANEL.read_text(encoding="utf-8").splitlines()
             if ln.strip() and not ln.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))


def _fn():
    src = B30.read_text(encoding="utf-8")
    m = re.search(r"^def _attach_kg_group_locations\(analyses, loc_rows\):\n(?:.*\n)+?    return n_located\n",
                  src, re.M)
    if not m:
        raise AssertionError("_attach_kg_group_locations not found in the production script")
    ns = {}
    exec(m.group(0), ns)  # noqa: S102 -- 断言对象就是生产函数本体（stdlib-only，自包含）
    return ns["_attach_kg_group_locations"]


class TestPanelFile(unittest.TestCase):
    """面板文件自身可核查：26 人群齐全、坐标合法、精度一律 region、双语与采样地非空。"""

    def test_all_26_phase3_populations_present_once_with_valid_coords(self):
        rows = _rows()
        self.assertEqual(list(rows[0].keys()), _COLS, "列序应与表头声明一致")
        pops = [r["source_id"] for r in rows]
        self.assertEqual(set(pops), PHASE3_POPS, "phase3 的 26 个人群必须恰好各一行")
        self.assertEqual(len(pops), len(set(pops)), "人群不得重复")
        for r in rows:
            lat, lon = float(r["latitude"]), float(r["longitude"])
            self.assertTrue(-90 <= lat <= 90 and -180 <= lon <= 180, r["source_id"])
            self.assertEqual(r["precision"], "region",
                             f"{r['source_id']}: 人群采样点一律 region，不得冒充 site")
            self.assertEqual(r["location_id"], f"KG:{r['source_id']}")
            self.assertEqual(r["dataset"], "1000G")
            for col in ("label_zh", "label_en", "locality", "source_url"):
                self.assertTrue(r[col].strip(), f"{r['source_id']}.{col} 非空")
            self.assertIn("superpop=", r["note"], f"{r['source_id']}: note 应带 superpop")
            # 坐标的性质必须在 note 里声明：要么"描述含采样地"（城市级），要么"质心"（区域近似）；
            # 都不写就是来源不明的一个点。
            self.assertTrue("质心" in r["note"] or "描述含采样" in r["note"],
                            f"{r['source_id']}: note 必须写明坐标是采样地还是区域质心")

    def test_diaspora_populations_geocoded_at_sampling_site(self):
        by = {r["source_id"]: r for r in _rows()}
        for pop in ("STU", "ITU", "GIH", "MXL", "ASW"):
            self.assertIn("离散", by[pop]["note"],
                          f"{pop}: 离散人群的坐标在采样地，note 必须言明，避免读成祖源地")


class TestAttachFunction(unittest.TestCase):
    """抽取生产函数本体：只补不覆盖、不认识的人群如实计 unlocated、非 1000G 不碰。"""

    def test_attaches_to_matching_groups_and_counts_unlocated(self):
        f = _fn()
        an = [{"analysis_id": "kg-global", "dataset": "1000G", "groups": [
            {"group_id": "KHV", "label": "KHV"},
            {"group_id": "ZZZ", "label": "ZZZ"}]}]
        n = f(an, [{"source_id": "KHV", "location_id": "KG:KHV", "latitude": "10.8231",
                    "longitude": "106.6297", "locality": "Ho Chi Minh City, Vietnam",
                    "precision": "region", "label_zh": "京族（胡志明市）"}])
        self.assertEqual(n, 1)
        khv, zzz = an[0]["groups"]
        self.assertEqual(khv["latitude"], 10.8231)
        self.assertEqual(khv["longitude"], 106.6297)
        self.assertEqual(khv["location_id"], "KG:KHV")
        self.assertEqual(khv["location_precision"], "region")
        self.assertEqual(khv["locality"], "Ho Chi Minh City, Vietnam")
        self.assertEqual(khv["name_zh"], "京族（胡志明市）")
        self.assertNotIn("latitude", zzz, "面板里没有的人群保持无坐标，不造假")
        self.assertEqual(an[0]["group_location_counts"], {"located": 1, "unlocated": 1})

    def test_non_kg_analyses_untouched(self):
        f = _fn()
        an = [{"analysis_id": "aadr-human-origins", "dataset": "AADR",
               "groups": [{"group_id": "Han", "label": "Han"}]}]
        self.assertEqual(f(an, [{"source_id": "Han", "latitude": "1", "longitude": "2"}]), 0)
        self.assertNotIn("latitude", an[0]["groups"][0], "AADR 记录的坐标来自 .anno，不走本层")
        self.assertNotIn("group_location_counts", an[0])

    def test_existing_coordinates_not_overwritten(self):
        f = _fn()
        an = [{"dataset": "1000G", "groups": [{"group_id": "CHB", "latitude": 0.0, "longitude": 0.0}]}]
        f(an, [{"source_id": "CHB", "latitude": "39.9", "longitude": "116.4"}])
        self.assertEqual((an[0]["groups"][0]["latitude"], an[0]["groups"][0]["longitude"]), (0.0, 0.0),
                         "将来 04b 若自带坐标，本层只补不覆盖")

    def test_malformed_panel_rows_skipped(self):
        f = _fn()
        an = [{"dataset": "1000G", "groups": [{"group_id": "CHB"}]}]
        bad = [{"source_id": "CHB", "latitude": "", "longitude": ""},          # 空坐标
               {"source_id": "JPT", "latitude": "999", "longitude": "0"},      # 越界
               {"source_id": " ", "latitude": "1", "longitude": "1"},          # 空人群
               {"source_id": "YRI", "latitude": "7.3775", "longitude": "3.9470",
                "location_id": "KG:YRI", "locality": "Ibadan, Nigeria", "precision": "region"}]
        self.assertEqual(f(an, bad), 0, "YRI 在面板行里但分析组里没有它——定位的是组，不是面板行")
        self.assertNotIn("latitude", an[0]["groups"][0])
        self.assertEqual(an[0]["group_location_counts"], {"located": 0, "unlocated": 1})


class TestWiring(unittest.TestCase):
    """接线断言：30 必须真调用该函数并读该 panel（防死代码/防旁路），且能通过编译。"""

    def test_30_calls_the_function_with_the_panel(self):
        src = B30.read_text(encoding="utf-8")
        self.assertIn("_attach_kg_group_locations(_anc_analyses", src)
        self.assertIn("panel\" / \"kg_population_locations.tsv", src)
        py_compile.compile(str(B30), doraise=True, cfile="/tmp/_30_kgl.pyc")


if __name__ == "__main__":
    unittest.main(verbosity=2)
