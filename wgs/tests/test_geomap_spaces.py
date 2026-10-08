"""AN6 增强②回归：参考空间切换 + 1000G 人群采样点渲染（跑 report_script.js 里的真实纯函数）。

§7.4 AN6"1000G 位置/参考空间切换未做"：地理分布卡原先只画 default_analysis_id 那一个空间，
1000G 参考空间（kg-global/kg-regional）在报告里没有地图表达。本步把 1000G 分析的**人群组**
经 kgGroupPoints 变成采样点（一组一点，key=location_id，region 级，带 n 与组均值距离），并加
空间切换行（每个 state=ok 的分析一个按钮；切换只换已算好的结果，不重算统计）。

纯函数住在模板里不能 import（需要 DOM），沿用 node 抽取法；没有 node 时**跳过并说明**——按
§7，未执行不等于通过。另附源级断言：draw() 必须真消费 kgGroupPoints 与 SPACES[sa]（防纯
函数/切换器变死代码），1000G 空间不得出现无意义的古代时间轴，单空间报告不得冒出切换行。
"""
import pathlib, shutil, subprocess, unittest

JS = pathlib.Path(__file__).resolve().parents[1] / "templates" / "report_script.js"

HARNESS = r"""
const fs = require('fs');
const src = fs.readFileSync(process.argv[1], 'utf8');
const grab = (re) => { const m = src.match(re); if (!m) throw new Error('pattern not found: ' + re); return m[0]; };
const code = [
  grab(/const geoXY[\s\S]*?\n}/),
  grab(/const geoValid[\s\S]*?\n};/),
  grab(/function analysisGroups[\s\S]*?\n}/),
  grab(/function kgGroupPoints[\s\S]*?\n}/),
].join('\n').replace(/\bconst /g, 'var ');
eval(code);
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };
// —— kg 分析：两组已定位（带采样点字段，30 从 panel 附上）+ 一组无坐标。
const kg = {dataset: '1000G', analysis_id: 'kg-global', groups: [
  {group_id: 'KHV', label: 'KHV', kind: 'modern', n: 99, distance_mean: 0.0073,
   latitude: 10.8231, longitude: 106.6297, location_id: 'KG:KHV',
   locality: 'Ho Chi Minh City, Vietnam', location_precision: 'region', name_zh: '京族（胡志明市）'},
  {group_id: 'CHB', label: 'CHB', kind: 'modern', n: 103, distance_mean: 0.0121,
   latitude: 39.9042, longitude: 116.4074, location_id: 'KG:CHB', locality: 'Beijing, China'},
  {group_id: 'ZZZ', label: 'ZZZ', kind: 'modern', n: 5, distance_mean: 0.9, member_ids: ['x1']},
]};
const gp = kgGroupPoints(kg);
ok(gp.points.length === 2, 'located populations become points (got ' + gp.points.length + ')');
ok(gp.unlocated.length === 1 && gp.unlocated[0].group_id === 'ZZZ',
   'population without coordinates stays unlocated and is reported, never faked');
const khv = gp.points[0];
ok(khv.key === 'KG:KHV', 'selection key is the stable location_id');
ok(khv.label === 'KHV' && khv.n === 99 && Math.abs(khv.d - 0.0073) < 1e-12,
   'point carries the group identity, size and mean distance');
ok(khv.precision === 'region' && khv.locality === 'Ho Chi Minh City, Vietnam',
   'point carries region precision and sampling locality');
ok(khv.name_zh === '京族（胡志明市）', 'Chinese name survives for label display');
ok(khv.kind === 'modern', 'kind is carried so kindView keeps working');
ok(gp.points.every(p => geoValid(p.latitude, p.longitude)), 'all emitted points are placeable');
// —— groups 的 {modern,ancient} 字典形状（AADR/09b）同样可点化：字典归一化照走 analysisGroups。
const dictShape = kgGroupPoints({dataset: 'AADR', groups: {modern: [kg.groups[0]], ancient: []}});
ok(dictShape.points.length === 1 && dictShape.points[0].key === 'KG:KHV',
   'dict-shaped groups normalize through the same path');
// —— 生产路径消费（drawGeoMap 的教训：纯函数可以是死代码）。
ok(/kgGroupPoints\(a\)/.test(src), 'geomap draw() consumes kgGroupPoints in the 1000G space');
ok(/SPACES\[sa\]/.test(src) && /SPACES=A\.filter\(x=>String\(x\.state\|\|''\)==='ok'\)/.test(src),
   'space switching selects among the state=ok analyses via SPACES[sa]');
ok(/SPACES\.length>1/.test(src), 'the switch row only appears when more than one space exists');
ok(src.includes("kv!=='modern'&&!kgSpace"),
   'the ancient timeline is suppressed in the 1000G space (no dates there; an empty axis is not an honest view)');
ok(/sa=i;sel=null;/.test(src), 'switching space clears the selection (point keys differ across spaces)');
console.log(FAIL.length ? 'FAIL\n' + FAIL.join('\n') : 'PASS');
process.exit(FAIL.length ? 1 : 0);
"""


@unittest.skipIf(shutil.which("node") is None, "node not installed here (runs on the server)")
class TestGeomapSpaces(unittest.TestCase):
    def test_kg_points_and_space_switch_consumption(self):
        r = subprocess.run(["node", "-e", HARNESS, "--", str(JS)],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
