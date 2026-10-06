"""人工复核回归：染色体画的片段颜色必须与图例色块同源（跑 report_script.js 里的真实表达式）。

用户在成品报告里看到"北方东亚的图层不见了"：bd3d3fc 把来源面板 A 的片段画成 c.bg——
即页面底色（浅色主题 #F7F2EB）——1458 个不透明米色矩形盖掉淡蓝底轨，图例却仍用实心
--ancn 深蓝色块代表北方东亚。数据完好（la_segments 里 NorthEA 1458 段俱在），纯渲染缺陷。
表达式住在模板里不能 import（需要 DOM），沿用 test_la_calib_card 的 node 抽取法；
没有 node 时**跳过并说明**——按 §7，未执行不等于通过。
"""
import pathlib, shutil, subprocess, unittest

TPL = pathlib.Path(__file__).resolve().parents[1] / "templates"
JS = TPL / "report_script.js"
BODY = TPL / "report_body.html"

HARNESS = r"""
const fs = require('fs');
const src = fs.readFileSync(process.argv[1], 'utf8');
const body = fs.readFileSync(process.argv[2], 'utf8');
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };
// 抓染色体画里片段的 fill 三元表达式（reveal('painting' 内唯一一处 g.anc 分支）
const m = src.match(/fill:\(_P1&&g\.anc===_P1\)\?a\.s:\(g\.anc===_P0\?\s*([A-Za-z_.]+)\s*:\s*([A-Za-z_.]+)\s*\)\}/);
if (!m) { console.log('FAIL\npainting segment fill expression not found'); process.exit(1); }
const pickP0 = m[1].trim(), pickOther = m[2].trim();
// 面板 A 的片段色必须是调色板主色 a.n（与图例 lg_n 色块 var(--ancn) 同源），
// 不能是 c.bg——页面底色当数据层，画在底轨上等于图层消失（本缺陷的正是形态）。
ok(pickP0 === 'a.n', `panel-A segment fill is ${pickP0}; must be the palette colour a.n, ` +
   `not the page background (the legend swatch lg_n is solid --ancn, so a c.bg layer is invisible)`);
ok(!/c\.bg/.test(m[0]), 'no branch of the segment fill may use c.bg (page background)');
ok(pickOther === 'c.fd', `non-source panels keep the neutral colour (got ${pickOther})`);
// AC() 把 a.n 绑定到 --ancn：这是"与图例同源"的另一半证据
ok(/n:V\('--ancn'\)/.test(src), "palette maps a.n to var(--ancn) (AC())");
// 图例侧的契约：lg_n 色块用实心 --ancn。片段色与色块都指 --ancn，二者才对得上。
ok(/<i style="background:var\(--ancn\)"><\/i><em data-i18n="lg_n"/.test(body),
   'legend swatch for lg_n (NorthEA) is solid var(--ancn) — the fill colour must match it');
console.log(FAIL.length ? 'FAIL\n' + FAIL.join('\n') : 'PASS');
process.exit(FAIL.length ? 1 : 0);
"""


@unittest.skipIf(shutil.which("node") is None, "node not installed here (runs on the server)")
class TestLAPaintingColors(unittest.TestCase):
    def test_panel_a_layer_is_visible_and_matches_legend(self):
        r = subprocess.run(["node", "-e", HARNESS, "--", str(JS), str(BODY)],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main()
