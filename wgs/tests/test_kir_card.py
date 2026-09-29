"""M3/H2 回归：KIR 卡片的单倍型/配体文案纯函数（跑 report_script.js 里的函数体）。

审查复现的问题：卡片写死"单倍型：AA（两条 A 单倍型）"与示例样本的 C*03:04 /
C*01:02 / B*13:01 / B*54:01 配体句，与本样本真实的 Bx 及其 HLA 不一致；D.kir=[]
时也照印。修正后这些句子全部来自 24 的 kir_summary（经 30 的 D.kir_summary），
未运行是 unknown，不是 absent。函数住在模板里不能 import（需要 DOM），用 node
抽出函数体再断言；没有 node 时**跳过并说明**——按 §7，未执行不等于通过。
"""
import pathlib, shutil, subprocess, unittest

JS = pathlib.Path(__file__).resolve().parents[1] / "templates" / "report_script.js"

HARNESS = r"""
const fs = require('fs');
const src = fs.readFileSync(process.argv[1], 'utf8');
const grab = (re) => { const m = src.match(re); if (!m) throw new Error('pattern not found: ' + re); return m[0]; };
const code = [
  grab(/function kirRan[\s\S]*?\n}/),
  grab(/function kirHapLabel[\s\S]*?\n  return \(Z\?'单倍型：':'Haplotype: '\)\+hap;\n}/),
  grab(/function kirLigandLine[\s\S]*?\n  return seg\.join\(Z\?' · ':'; '\);\n}/),
].join('\n').replace(/\bconst /g, 'var ');
eval(code);
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };
// —— 是否运行：空表 = 未运行（unknown），不是全 absent
ok(kirRan([]) === false, 'empty kir.tsv means typing did not run');
ok(kirRan(undefined) === false && kirRan(null) === false, 'missing kir rows are not a run');
ok(kirRan([{gene:'KIR2DL1',present:true}]) === true, 'rows present means typing ran');
// —— 单倍型：来自 summary，两种真实取值都不写死
ok(kirHapLabel({haplotype:'Bx (at least one B haplotype)'}, true) === '单倍型：Bx（至少一条 B 单倍型）',
   'this sample is Bx: the card must print Bx, not the old hardcoded AA');
ok(kirHapLabel({haplotype:'AA (two A haplotypes)'}, false) === 'Haplotype: AA (two A haplotypes)',
   'an AA sample still prints AA -- the fix must not invert the claim either');
ok(kirHapLabel({haplotype:'unavailable (T1K was not run)'}, true) === '单倍型：未知（T1K 未运行）',
   'unavailable summary reads as unknown, never as a haplotype');
ok(kirHapLabel(null, true) === '单倍型：未知（KIR 分型未运行）' && kirHapLabel({}, false) === 'Haplotype: unknown (KIR typing not run)',
   'missing summary reads as unknown');
ok(kirHapLabel({haplotype:'Cx (future call)'}, true) === '单倍型：Cx (future call)',
   'an unrecognised value passes through verbatim instead of being coerced to AA/Bx');
// —— 配体行：来自本次运行的 C/B/A 分型组，不再引用示例等位基因
ok(kirLigandLine({C_ligands:'C1', B_epitopes:'Bw6/Bw6', A311_ligands:''}, false)
   === 'HLA-C ligand groups: C1; HLA-B epitopes: Bw6/Bw6; HLA-A3/A11 (KIR3DL2 ligands): —',
   'the real sample row renders C1 / Bw6/Bw6 with an explicit empty A3/A11');
ok(kirLigandLine({C_ligands:'C1/C2', B_epitopes:'Bw4/Bw6', A311_ligands:'A*11:01'}, true)
   .includes('C1/C2') && kirLigandLine({C_ligands:'C1/C2', B_epitopes:'Bw4/Bw6', A311_ligands:'A*11:01'}, true).includes('A*11:01'),
   'typed groups flow through in both languages');
ok(kirLigandLine(null, true) === 'HLA-C 配体组：— · HLA-B 表位：— · HLA-A3/A11（KIR3DL2 配体）：—',
   'no summary renders explicit placeholders, not an example genotype');
// —— 回归：示例样本的等位基因不再出现在任何输出里
['C*03:04','C*01:02','B*13:01','B*54:01'].forEach((al)=>{
  ok(kirHapLabel({haplotype:'Bx (at least one B haplotype)'}, true).indexOf(al) < 0
     && kirLigandLine({C_ligands:'C1',B_epitopes:'Bw6/Bw6'}, false).indexOf(al) < 0,
     'example allele '+al+' must not leak into the verdict lines');
});
console.log(FAIL.length ? 'FAIL\n' + FAIL.join('\n') : 'PASS');
process.exit(FAIL.length ? 1 : 0);
"""


@unittest.skipIf(shutil.which("node") is None, "node not installed here (runs on the server)")
class TestKirCard(unittest.TestCase):
    def test_pure_helpers(self):
        r = subprocess.run(["node", "-e", HARNESS, "--", str(JS)],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_helpers_are_extractable_from_the_shipped_template(self):
        src = JS.read_text(encoding="utf-8")
        for pat in ("function kirRan", "function kirHapLabel", "function kirLigandLine"):
            self.assertIn(pat, src)
        # the card must call the helpers, not print a verdict inline
        self.assertIn("kirHapLabel(D.kir_summary", src)
        self.assertIn("kirLigandLine(D.kir_summary", src)
        # the example sample's alleles must be gone entirely (the AA label itself stays as
        # the translation of a real possible value inside kirHapLabel -- that is the point)
        for gone in ("C*03:04", "C*01:02", "B*13:01", "B*54:01"):
            self.assertNotIn(gone, src, f"hardcoded example allele still in template: {gone}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
