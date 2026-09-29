"""AN6 回归：geomap 的 groups 归一化纯函数（跑 report_script.js 里的函数体）。

复审复现的问题：09b 的 groups 是 {modern:[...], ancient:[...]}，模板旧代码字典形状只取
.groups.ancient——现代主图（默认视图）按 kind='modern' 过滤后一个不剩，排名/前五名/列表
全空。修正后 analysisGroups 把两层的已算好结果拼起来（不重算排名，各行保留自己梯队的
rank），KPI 的 _ranked 与 geomap 卡片共用同一归一化。函数住在模板里不能 import（需要
DOM），用 node 抽出函数体再断言；没有 node 时**跳过并说明**——按 §7，未执行不等于通过。
"""
import pathlib, shutil, subprocess, unittest

JS = pathlib.Path(__file__).resolve().parents[1] / "templates" / "report_script.js"

HARNESS = r"""
const fs = require('fs');
const src = fs.readFileSync(process.argv[1], 'utf8');
const grab = (re) => { const m = src.match(re); if (!m) throw new Error('pattern not found: ' + re); return m[0]; };
const code = [
  grab(/function analysisGroups[\s\S]*?\n}/),
  grab(/const _ranked=[\s\S]*?;\n/),
].join('\n').replace(/\bconst /g, 'var ');
eval(code);
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };

// —— 09b 的字典形状：两层拼接，modern 在前（稳定并列顺序）
const M1={kind:'modern',label:'Han',rank:1,distance_mean:0.021,n:40,member_ids:['a']},
      M2={kind:'modern',label:'Japanese',rank:2,distance_mean:0.025,n:38,member_ids:['b']},
      A1={kind:'ancient',label:'China_Am',rank:1,distance_mean:0.019,n:6,member_ids:['c']},
      SMALL={kind:'modern',label:'Orogen',rank:null,small_group:true,n:1,distance_mean:0.02,member_ids:['d']};
const dict={analysis_id:'aadr',groups:{modern:[M1,M2,SMALL],ancient:[A1]}};
ok(analysisGroups(dict).length===4, 'dict shape concatenates both layers');
ok(analysisGroups(dict).slice(0,2)[0]===M1 && analysisGroups(dict).slice(0,2)[1]===M2,
   'modern layer comes first, ancient appended');
// —— 复审点名的那条：只有 modern 的字典（ancient 为空）旧代码返回 []，现代排名整图消失
ok(analysisGroups({groups:{modern:[M1],ancient:[]}}).length===1,
   'modern-only dict must yield the modern groups (old code returned [])');
ok(analysisGroups({groups:{modern:[],ancient:[A1]}}).length===1, 'ancient-only dict keeps ancient');
// —— 1000G（04b）的平铺数组原样通过
const flat=[M1,M2,A1];
ok(analysisGroups({groups:flat})===flat, 'array shape passes through untouched');
// —— 缺失形状是 []，不是异常
ok(analysisGroups(null).length===0 && analysisGroups({}).length===0
   && analysisGroups({groups:null}).length===0 && analysisGroups({groups:'x'}).length===0,
   'missing/malformed groups normalise to []');
// —— KPI 的 _ranked 与 geomap 卡片共用归一化：字典分析的现代排名非空（旧代码为空），
//     small_group 不进默认榜，两个梯队的 rank 各自保留
const ranked=_ranked(dict);
ok(ranked.length===3, 'ranked view drops the small_group row');
const modernView=ranked.filter(g=>String(g.kind||'')==='modern');
ok(modernView.length===2 && modernView[0].label==='Han' && modernView[0].rank===1,
   'modern layer view (the geomap default) is populated and ordered');
ok(ranked[0]===A1 && ranked[1]===M1 && ranked[0].rank===1 && ranked[1].rank===1,
   'cross-ladder rank ties break deterministically by distance_mean; computed ranks kept');
console.log(FAIL.length ? 'FAIL\n' + FAIL.join('\n') : 'PASS');
process.exit(FAIL.length ? 1 : 0);
"""


@unittest.skipIf(shutil.which("node") is None, "node not installed here (runs on the server)")
class TestAnalysisGroupsNormalisation(unittest.TestCase):
    def test_pure_functions_from_the_template(self):
        r = subprocess.run(["node", "-e", HARNESS, "--", str(JS)],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("PASS", r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
