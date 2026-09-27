"""AN6：地图的坐标变换、BP→公元换算与年代筛选（跑 report_script.js 里的纯函数）。

这些逻辑住在构建器用的模板里，不能 import（脚本需要 DOM），所以用 node 抽出函数体再断言。
没有 node 时**跳过并说明**——按 §7，未执行不等于通过。
"""
import pathlib, shutil, subprocess, unittest

JS = pathlib.Path(__file__).resolve().parents[1] / "templates" / "report_script.js"

HARNESS = r"""
const fs = require('fs');
// node -e <script> -- <path>: argv[0]=node, argv[1]=path
const src = fs.readFileSync(process.argv[1], 'utf8');
const grab = (re) => { const m = src.match(re); if (!m) throw new Error('pattern not found: ' + re); return m[0]; };
const code = [
  grab(/function overlapsAge[\s\S]*?\n}/),
  grab(/function ageInRange[\s\S]*?\n}/),
  grab(/const GEO_VB[\s\S]*?Math\.abs\(b\) <= 180;\n};/),
].join('\n').replace(/\bconst /g, 'var ');   // const 在 eval 里是块作用域，取不出来
eval(code);
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };
// —— §7 公布的两条筛选用例
ok(overlapsAge(3000, 4000, 3500, 4500) === true, 'section 7: overlapping interval must match');
ok(overlapsAge(null, null, 3500, 4500) === false, 'section 7: unknown interval must not match');
ok(overlapsAge(4500, 5000, 3500, 4500) === true, 'touching at the boundary counts');
ok(overlapsAge(1000, 2000, 3500, 4500) === false, 'entirely outside must not match');
// —— 只有均值 / 无年代
ok(ageInRange({ date_mean_bp: 4000 }, 3500, 4500) === true, 'mean-only record uses a point');
ok(ageInRange({ date_min_bp: 3000, date_max_bp: 4000 }, 3500, 4500) === true, 'interval is preferred');
ok(ageInRange({}, 3500, 4500) === null, 'no date is neither matched nor excluded');
// —— 空值不能变成 (0,0)
ok(geoValid(null, null) === false, 'null coordinates are invalid (must not become 0,0)');
ok(geoValid('', '') === false, 'empty-string coordinates are invalid');
ok(geoValid(undefined, undefined) === false, 'undefined coordinates are invalid');
ok(geoValid(91, 0) === false && geoValid(0, 181) === false, 'out-of-range coordinates are invalid');
ok(geoValid(-90, 180) === true && geoValid('39.9', '116.4') === true, 'valid values pass (strings included)');
// —— 等距投影：四角与已知城市
ok(Math.abs(geoXY(90, -180).x) < 1e-9 && Math.abs(geoXY(-90, 180).y - 180) < 1e-9, 'world corners land in the viewBox');
const bj = geoXY(39.9, 116.4);
ok(bj.x > 290 && bj.y < 60, 'Beijing lands in the north-east quadrant');
ok(geoXY(0, 0).x === 180 && geoXY(0, 0).y === 90, 'equator/prime meridian maps to the centre');
// —— BP → 公元（1950 基准）：方向写反就会把"前 2050 年"说成"公元 2050 年"
const bp2ce = (v) => 1950 - v;
ok(bp2ce(2000) === -50, '2000 BP is 50 BCE, not 50 CE');
ok(bp2ce(4000) === -2050 && bp2ce(8000) === -6050, 'older dates are BCE');
ok(bp2ce(1950) === 0 && bp2ce(1000) === 950, 'the 1950 baseline and recent dates are CE');
if (FAIL.length) { console.error('FAILED:\n  ' + FAIL.join('\n  ')); process.exit(1); }
console.log('geo assertions passed');
"""


class TestGeoMap(unittest.TestCase):
    def test_report_geo_helpers(self):
        node = shutil.which("node")
        if not node:
            self.skipTest("node not available: the geo assertions were NOT run")
        r = subprocess.run([node, "-e", HARNESS, "--", str(JS)], capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, (r.stdout or "") + (r.stderr or ""))

    def test_base_map_matches_the_transform(self):
        """底图 path 必须落在 geoXY 产生的同一坐标系里，否则点会与轮廓错位。"""
        import re
        svg = pathlib.Path(__file__).resolve().parents[1] / "panel" / "world_land.svg"
        if not svg.exists():
            self.skipTest("panel/world_land.svg missing")
        text = svg.read_text(encoding="utf-8")
        self.assertIn("viewBox=\"0 0 360 180\"", text)
        self.assertIn("equirectangular", text, "投影规则要写在 metadata 里")
        self.assertIn("public domain", text, "许可证要写在 metadata 里")
        self.assertIn("sha256", text, "源文件校验和要写在 metadata 里")
        pts = re.findall(r"([0-9.]+) ([0-9.]+)", text.split('d="', 1)[1])
        xs = [float(a) for a, _ in pts]
        ys = [float(b) for _, b in pts]
        self.assertGreaterEqual(min(xs), 0.0)
        self.assertLessEqual(max(xs), 360.0)
        self.assertGreaterEqual(min(ys), 0.0)
        self.assertLessEqual(max(ys), 180.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
