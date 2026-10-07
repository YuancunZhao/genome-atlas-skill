"""人工复核回归：地理分布卡取景随内容放大（跑 report_script.js 里的真实纯函数）。

用户复核成品报告：地图有内容的区域（本样本全在东亚）只占整幅世界底图的一角，点挤在一起
分不清。geomapViewport 按已落点包围盒取景：留边距、跨度下限 70×45（保区域上下文，不让
一两个点放大到街区级）、视口夹在世界边界内；无有效坐标或横跨全球回退全世界。纯函数住在
模板里不能 import（需要 DOM），沿用 node 抽取法；没有 node 时**跳过并说明**——按 §7，
未执行不等于通过。另附源级断言：底图定位必须走 use transform 而非 use 的 x/y——取景后
x/y 为负，个别引擎不应用负定位，底图落回原点、点悬大洋（人工复核第 3 条）。
"""
import pathlib, shutil, subprocess, unittest

JS = pathlib.Path(__file__).resolve().parents[1] / "templates" / "report_script.js"

HARNESS = r"""
const fs = require('fs');
const src = fs.readFileSync(process.argv[1], 'utf8');
const grab = (re) => { const m = src.match(re); if (!m) throw new Error('pattern not found: ' + re); return m[0]; };
const code = [
  grab(/const GEO_VB[\s\S]*?\n}/),
  grab(/const geoXY[\s\S]*?\n}/),
  grab(/const geoValid[\s\S]*?\n};/),
  grab(/function geomapViewport[\s\S]*?\n}/),
].join('\n').replace(/\bconst /g, 'var ');
eval(code);
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };
const inside = (vp, lat, lon) => { const q = geoXY(lat, lon);
  return q.x >= vp.x0 - 1e-9 && q.x <= vp.x0 + vp.w + 1e-9 && q.y >= vp.y0 - 1e-9 && q.y <= vp.y0 + vp.h + 1e-9; };
// —— 本样本形态：点集中在东亚（约经 75–135、纬 15–53）。视口应真正放大（小于全世界）、
//    装下全部点、并保住最小跨度（区域上下文）。
const eas = [].concat(...[75,100,135].map(lon => [15,34,53].map(lat => ({latitude: lat, longitude: lon}))));
const vp = geomapViewport(eas);
ok(vp.w < 360 && vp.h < 180, `east-asia cluster must zoom in (got ${vp.w}x${vp.h})`);
ok(vp.w >= 70 - 1e-9 && vp.h >= 45 - 1e-9, `span keeps regional context (>=70x45, got ${vp.w}x${vp.h})`);
eas.forEach(r => ok(inside(vp, r.latitude, r.longitude), `viewport must contain (${r.latitude},${r.longitude})`));
ok(vp.x0 >= 0 && vp.x0 + vp.w <= 360 + 1e-9 && vp.y0 >= 0 && vp.y0 + vp.h <= 180 + 1e-9,
   'viewport stays inside world bounds');
// —— 无点 / 全部缺坐标 → 全世界（没有内容可取景，不能拿几内亚湾 (0,0) 兜底）
const full = { x0: 0, y0: 0, w: 360, h: 180 };
ok(JSON.stringify(geomapViewport([])) === JSON.stringify(full), 'no rows -> whole world');
ok(JSON.stringify(geomapViewport([{latitude: null, longitude: 116}, {latitude: '', longitude: ''}]))
   === JSON.stringify(full), 'rows without usable coordinates -> whole world');
// —— 单点：最小跨度、不退化、包含该点（近世界边缘时夹回边界内）
const one = geomapViewport([{latitude: 40, longitude: 116}]);
ok(one.w >= 70 - 1e-9 && one.h >= 45 - 1e-9, 'single point keeps the minimum span');
ok(inside(one, 40, 116), 'single point is inside its viewport');
const edge = geomapViewport([{latitude: 70, longitude: 179}]);
ok(edge.x0 + edge.w <= 360 + 1e-9 && inside(edge, 70, 179), 'viewport clamps at the world edge and still contains the point');
// —— 点横跨全球 → 回退全世界（无从放大）
const world = [].concat(...[-170,-60,0,60,170].map(lon => [-60,0,60].map(lat => ({latitude: lat, longitude: lon}))));
ok(geomapViewport(world).w === 360, 'points spanning the globe -> whole world');
// —— 生产路径消费（drawGeoMap 的教训：纯函数可以是死代码）：卡片 draw() 必须取景并裁剪底图，
//    且底图定位必须走 transform——use 的 x/y 定位（取景后常为负）在个别引擎里不被应用，
//    底图落回原点而点留在放大位置（人工复核第 3 条截图：北美居左、点悬大西洋）。
ok(/geomapViewport\(placed\)/.test(src), 'geomap draw() derives its viewport from the placed rows');
ok(src.includes("el(s,'g',{'clip-path':'url(#geomap_plot)'})"),
   'base map is clipped via the wrapper g (it overflows the plot rect when zoomed)');
ok(src.includes("el(_mg,'use',{href:'#world_land'") &&
   src.includes('transform:`translate(${ox} ${oy}) scale(${sc})`'),
   'base map positions/scales via use transform (equivalent geometry on every engine)');
ok(!src.includes("'use',{href:'#world_land',x:ox"),
   'use must not position via x/y: negative x/y is dropped back to the origin by some engines');
console.log(FAIL.length ? 'FAIL\n' + FAIL.join('\n') : 'PASS');
process.exit(FAIL.length ? 1 : 0);
"""


@unittest.skipIf(shutil.which("node") is None, "node not installed here (runs on the server)")
class TestGeomapViewport(unittest.TestCase):
    def test_viewport_math_and_consumption(self):
        r = subprocess.run(["node", "-e", HARNESS, "--", str(JS)],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main()
