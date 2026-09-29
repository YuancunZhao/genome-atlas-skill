"""AN3/AN5 回归：局部祖源校准卡片的取数纯函数（跑 report_script.js 里的函数体）。

审查复现的问题：dayuCal 从 17 的 per_chrom **片段跨度比例**折算（例 0.9），再与 17b 的
**posterior** 校准均值/SD 比较——片段 0.9、posterior 目标 0.6、参照均值 0.5/SD 0.1 时页面显示
90%/+4 SD，同口径正解是 60%/+1 SD。修正后 KPI 与校准图都取 17b 写回的 target_north；
无校准时标签回退原始 LA 的来源面板，而不是变 '—'。函数住在模板里不能 import（需要 DOM），
用 node 抽出函数体再断言；没有 node 时**跳过并说明**——按 §7，未执行不等于通过。
"""
import pathlib, shutil, subprocess, unittest

JS = pathlib.Path(__file__).resolve().parents[1] / "templates" / "report_script.js"

HARNESS = r"""
const fs = require('fs');
const src = fs.readFileSync(process.argv[1], 'utf8');
const grab = (re) => { const m = src.match(re); if (!m) throw new Error('pattern not found: ' + re); return m[0]; };
const code = [
  grab(/function calibTarget[\s\S]*?\n}/),
  grab(/function laSdScore[\s\S]*?\n}/),
  grab(/function laPanelLabels[\s\S]*?\n  return \(rawPans\|\|\[\]\)\.filter[\s\S]*?\n}/),
].join('\n').replace(/\bconst /g, 'var ');
eval(code);
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };
// —— 审查的合成例：posterior 目标 0.6 vs 参照均值 0.5/SD 0.1 → +1 SD（不是片段 0.9 的 +4 SD）
ok(calibTarget([{population:'CHB',n:20,north_mean:0.5,north_sd:0.1,target_north:0.6}]) === 0.6,
   'the KPI value is 17b posterior target_north (0.6), not the 0.9 span proportion');
const z = laSdScore(calibTarget([{north_mean:0.5,north_sd:0.1,target_north:0.6}]), 0.5, 0.1);
ok(Math.abs(z - 1.0) < 1e-9, '0.6 vs mean 0.5 sd 0.1 is +1 SD, not the +4 SD the span proportion gave');
// —— target_north 缺失（17b 未运行）：null，不拿别的口径顶上
ok(calibTarget([]) === null, 'no calibration rows -> no KPI value');
ok(calibTarget([{north_mean:0.5,north_sd:0.1}]) === null,
   'rows without target_north (17b never attached it) -> null, never the span proportion');
ok(calibTarget(null) === null && calibTarget(undefined) === null, 'missing calibration is null');
// —— SD 计算的守卫：sd=0（组内同值）/均值缺失/目标缺失都不给 z
ok(laSdScore(0.6, 0.5, 0) === null, 'sd exactly 0 gives no z-score');
ok(laSdScore(0.6, null, 0.1) === null && laSdScore(null, 0.5, 0.1) === null,
   'missing target or mean gives no z-score');
ok(laSdScore(NaN, 0.5, 0.1) === null, 'non-finite target gives no z-score');
ok(Math.abs(laSdScore(0.3, 0.5, 0.1) + 2) < 1e-9, 'negative direction keeps its sign (-2 SD)');
// —— 标签回退：空 calibration_panels 不再变 '—'，回退到原始 LA 的来源面板
ok(JSON.stringify(laPanelLabels([], [{role:'source',id:'NorthEA'},{role:'source',id:'SouthEA'},{role:'control',id:'European'}]))
   === JSON.stringify(['NorthEA','SouthEA']),
   'empty calibration_panels falls back to the raw LA source panels, not dashes');
ok(JSON.stringify(laPanelLabels(['A','B'], [{role:'source',id:'NorthEA'}])) === JSON.stringify(['A','B']),
   'calibration panels take precedence when present (config la_labels)');
ok(JSON.stringify(laPanelLabels([], [])) === JSON.stringify([]), 'both empty -> empty (caller renders dashes)');
ok(JSON.stringify(laPanelLabels(undefined, undefined)) === JSON.stringify([]), 'missing fields degrade to empty');
console.log(FAIL.length ? 'FAIL\n' + FAIL.join('\n') : 'PASS');
process.exit(FAIL.length ? 1 : 0);
"""


@unittest.skipIf(shutil.which("node") is None, "node not installed here (runs on the server)")
class TestLACalibCard(unittest.TestCase):
    def test_pure_helpers(self):
        r = subprocess.run(["node", "-e", HARNESS, "--", str(JS)],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_helpers_are_extractable_from_the_shipped_template(self):
        src = JS.read_text(encoding="utf-8")
        for pat in ("function calibTarget", "function laSdScore", "function laPanelLabels"):
            self.assertIn(pat, src)
        # the card must consume calibTarget/laPanelLabels, not recompute a span proportion inline
        self.assertIn("const dayuCal=calibTarget(_CAL);", src)
        self.assertIn("const _PANELS=laPanelLabels(", src)
        self.assertNotIn("calLen", src, "the chromosome-length span folding is gone entirely")


if __name__ == "__main__":
    unittest.main(verbosity=2)
