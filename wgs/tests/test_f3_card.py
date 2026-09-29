"""W-T2 回归：f3 卡片的取数/判定纯函数（跑 report_script.js 里的函数体）。

审查复现的问题：卡片把"混合 f3 均为正"写死成本样本结论，且 _sg 对 null 值会打出
"+0.0"（Number(null)===0）。修正后的 28 会在零 SE/非有限时输出 null，模板的每个
数字和结论必须从行数据本身得出。这些函数住在模板里不能 import（脚本需要 DOM），
用 node 抽出函数体再断言；没有 node 时**跳过并说明**——按 §7，未执行不等于通过。
"""
import pathlib, shutil, subprocess, unittest

JS = pathlib.Path(__file__).resolve().parents[1] / "templates" / "report_script.js"

HARNESS = r"""
const fs = require('fs');
const src = fs.readFileSync(process.argv[1], 'utf8');
const grab = (re) => { const m = src.match(re); if (!m) throw new Error('pattern not found: ' + re); return m[0]; };
const code = [
  grab(/function f3Num[\s\S]*?\n}/),
  grab(/function f3Verdict[\s\S]*?\n}/),
  grab(/function f3Cell[\s\S]*?\n}/),
  grab(/function f3Note[\s\S]*?\n  return 'no-sig-neg';\n}/),
].join('\n').replace(/\bconst /g, 'var ');
eval(code);
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };
// —— 数字格式：null/NaN/Infinity 不得伪装成 0（旧 _sg 的 Number(null)===0 陷阱）
ok(f3Num(null, 1) === null, 'null must stay null, not print +0.0');
ok(f3Num(undefined, 1) === null, 'undefined must stay null');
ok(f3Num(NaN, 1) === null && f3Num(Infinity, 1) === null, 'non-finite must stay null');
ok(f3Num(0.5, 2) === '+0.50', 'positive keeps the sign marker');
ok(f3Num(-0.25, 2) === '-0.25' && f3Num(0, 1) === '0.0', 'negative and zero format correctly');
// —— 行判定：结论从行数据得出，不再写死"全为正"
ok(f3Verdict({f3:-0.01, se:0.001, z:-10}) === 'sig-neg', 'negative f3 with |Z|>3 is significant admixture evidence');
ok(f3Verdict({f3:0.0757, se:0.0004, z:184.3}) === 'ns', 'the real sample row (all-positive, huge Z) is non-significant, not auto-evidence');
ok(f3Verdict({f3:-0.01, se:0.001, z:-2}) === 'ns', 'negative but |Z|<3 is NOT evidence');
ok(f3Verdict({f3:0.0, se:0.0, z:null}) === 'unavailable', 'zero-SE row has no Z and no verdict');
ok(f3Verdict({f3:-0.01, se:0.0, z:null}) === 'unavailable', 'negative f3 without a Z must not claim significance');
ok(f3Verdict(null) === 'unavailable' && f3Verdict({}) === 'unavailable', 'missing rows/fields are unavailable');
// —— 单元格文案：三个数字齐全才给"f3 ±se (Z=z)"，否则显式不可用
ok(f3Cell({f3:null, se:0.1, z:1}) === '—' && f3Cell({f3:0.1, se:null, z:1}) === '—', 'any missing f3/se renders the unavailable marker');
ok(f3Cell({f3:0.0757, se:0.0004, z:184.3}) === '+0.0757 ±0.0004  (Z=+184.3)', 'complete row formats all three numbers');
ok(f3Cell({f3:0.1, se:0.02, z:null}) === '+0.1000 ±0.0200  (Z=—)', 'finite f3/se with null Z shows Z=— instead of +0.0');
// —— 整卡结论分类：不再永远说"均为正"
ok(f3Note([]) === 'none', 'no admixture rows -> no conclusion sentence');
ok(f3Note([{f3:null, se:null, z:null}]) === 'all-unavailable', 'only-unavailable rows -> explicit no-verdict note');
ok(f3Note([{f3:0.07, se:0.001, z:70}, {f3:null, se:null, z:null}]) === 'no-sig-neg', 'at least one available non-significant row -> non-detection note');
ok(f3Note([{f3:0.07, se:0.001, z:70}]) === 'no-sig-neg', 'all positive -> non-detection (wording derived, not assumed)');
ok(f3Note([{f3:0.07, se:0.001, z:70}, {f3:-0.01, se:0.002, z:-5}]) === 'has-sig-neg', 'any significant-negative row flips the card to the admixture-signal note');
console.log(FAIL.length ? 'FAIL\n' + FAIL.join('\n') : 'PASS');
process.exit(FAIL.length ? 1 : 0);
"""


@unittest.skipIf(shutil.which("node") is None, "node not installed here (runs on the server)")
class TestF3Card(unittest.TestCase):
    def test_pure_helpers(self):
        r = subprocess.run(["node", "-e", HARNESS, "--", str(JS)],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_helpers_are_extractable_from_the_shipped_template(self):
        # The harness regexes must keep matching the template as it ships; a renamed helper
        # fails here with 'pattern not found' rather than silently testing nothing.
        src = JS.read_text(encoding="utf-8")
        for pat in ("function f3Num", "function f3Verdict", "function f3Cell", "function f3Note"):
            self.assertIn(pat, src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
