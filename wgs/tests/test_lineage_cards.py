"""AN6-③ 回归：谱系卡"已知发现"观测成点、树时间线基准如实（跑 report_script.js 里的真实纯函数）。

§7 AN6 约束在此卡上的落点：观测用点（是**该支系的已发表记录**，不是本样本坐标）；formed/TMRCA
（树估计）与观测的报告年代是两套基准，不得混用或共用无说明的年份；没有带来源的迁移路线就不画
箭头。两个纯函数与卡片/测试同源：
  - lineageObsPoints(doc)：观测 → 统一地点词汇（latitude/longitude/locality/precision），缺坐标
    进 unlocated 不造点；relation/报告年代/publication 原样透传。
  - lineageTimelineNodes(path)：YFull formed 优先；formed 缺失（含 0）退 TMRCA 并记 basis；
    两者皆无的节点（如 phylotree 的 mt 节点）不进时间线——旧实现把 formed=0 当"现在"，K-M2335
    这类只有 TMRCA 的节点会被放到时间轴最右端。
沿用 node 抽取法；没有 node 时**跳过并说明**。另附源级断言：卡片必须消费这两个纯函数（drawGeoMap
的教训），发现地图的底图定位必须走 use transform（845ead2），且"无路线不画箭头"的诚实文案在。
"""

import pathlib, re, shutil, subprocess, unittest

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
  grab(/function lineageObsPoints[\s\S]*?\n}/),
  grab(/function lineagePendingNote[\s\S]*?\n}/),
  grab(/function lineageTimelineNodes[\s\S]*?\n}/),
].join('\n').replace(/\bconst /g, 'var ');
eval(code);
const FAIL = [];
const ok = (cond, msg) => { if (!cond) FAIL.push(msg); };

// —— 真实数据形状（本样本 lineage_history.json 的 Y 段）：白羊村/宗日古人点 + 北京现代个体，
//    外加一条无坐标观测。观测是支系的已发表记录，字段原样透传。
const doc = { history: { observations: [
  { record_id: 'BaiyangcunM1_W.SG', locality: 'Baiyangcun (Yunnan, West_YN)',
    coordinates: { latitude: 25.84, longitude: 100.59 }, precision: 'unknown', relation: 'exact',
    date_range: { mean: 3300, min: 3350, max: 3650 }, publication: 'WangFuScience2025' },
  { record_id: 'CSP048.AG', locality: 'Zongri (Qinghai)',
    coordinates: { latitude: 35.3, longitude: 100.4 }, precision: 'region', relation: 'exact',
    date_range: { mean: 4635 }, publication: 'WangFuSciAdv2023' },
  { record_id: 'NA18639.DG', locality: 'Beijing (Han)',
    coordinates: { latitude: 39.916667, longitude: 116.383333 }, relation: 'exact',
    date_range: { mean: 0 }, publication: 'The1000GenomesProjectConsortiumNature2015' },
  { record_id: 'NOGEO', locality: 'somewhere', coordinates: {}, relation: 'exact', publication: 'X' },
]}};
const op = lineageObsPoints(doc);
ok(op.points.length === 3, `3 located observation points (got ${op.points.length})`);
ok(op.unlocated.length === 1, `1 observation without coordinates -> unlocated (got ${op.unlocated.length})`);
ok(op.pending.length === 0, `no pending records in this doc (got ${op.pending.length})`);
const b = op.points[0];
ok(b.key === 'BaiyangcunM1_W.SG' && b.latitude === 25.84 && b.longitude === 100.59, 'identity + coordinates carry through');
ok(b.locality.indexOf('Baiyangcun') === 0, 'locality carried');
ok(b.precision === 'unknown', 'precision defaults to unknown, never invented as site');
ok(op.points[1].precision === 'region', 'region precision carried');
ok(b.relation === 'exact' && b.publication === 'WangFuScience2025', 'relation + publication carried');
ok(b.date_mean === 3300 && b.date_min === 3350 && b.date_max === 3650, 'reported date range carried');
ok(op.points[2].date_mean === 0, 'modern individual (mean 0) stays a point, labelled at render time');
ok(JSON.stringify(lineageObsPoints(null)) === '{"points":[],"unlocated":[],"pending":[]}', 'no doc -> empty, not an error');
ok(JSON.stringify(lineageObsPoints({ history: {} })) === '{"points":[],"unlocated":[],"pending":[]}', 'no observations -> empty');

// —— 复审 §3.2 P0-2b：待核对（pending_review）不是已发表观测。有坐标也不落点、不进
//    points/unlocated（nObs 与地图都不得把它算进去），单独进 pending 桶。
const doc2 = { history: { observations: [
  ...doc.history.observations,
  { record_id: 'PND1', locality: 'Zongri (Qinghai)', relation: 'pending_review',
    source_version: 'YFull12.03',
    pending_reason: 'label matches by name, but version equivalence ... is not proven',
    coordinates: { latitude: 35.3, longitude: 100.4 }, publication: 'AADRv66' },
  { record_id: 'PND2', locality: 'no coords', relation: 'pending_review', source_version: '',
    coordinates: {}, publication: 'AADRv66' },
]}};
const op2 = lineageObsPoints(doc2);
ok(op2.points.length === 3 && op2.unlocated.length === 1,
   `pending records never become points/unlocated even with coordinates (got ${op2.points.length}/${op2.unlocated.length})`);
ok(op2.pending.length === 2, `both pending records land in the pending bucket (got ${op2.pending.length})`);
ok(op2.pending[0].record_id === 'PND1' && op2.pending[1].record_id === 'PND2', 'pending records carry through');

// —— lineagePendingNote：空桶不渲染；两种语言都写出记录侧树版本与本项目树版本，
//    并说明"不计入已发表观测"。
ok(lineagePendingNote({tree_source:'YFull',tree_version:'14.06.0'}, op, true) === ''
   && lineagePendingNote({tree_source:'YFull',tree_version:'14.06.0'}, op, false) === '',
   'no pending records -> no note line');
const lin2 = {tree_source:'YFull', tree_version:'14.06.0'};
const zh = lineagePendingNote(lin2, op2, true), en = lineagePendingNote(lin2, op2, false);
ok(zh.indexOf('待核对 2 条') === 0, `zh note leads with the pending count (got "${zh.slice(0,12)}")`);
ok(zh.indexOf('YFull12.03') > 0 && zh.indexOf('YFull 14.06.0') > 0,
   'zh note names BOTH tree versions (record side and this project side)');
ok(zh.indexOf('不计入已发表观测') > 0, 'zh note says they are not counted as published observations');
ok(en.indexOf('2 record(s) pending review') === 0, 'en note leads with the pending count');
ok(en.indexOf('YFull12.03') > 0 && en.indexOf('YFull 14.06.0') > 0, 'en note names both tree versions');
ok(en.indexOf('not counted as published observations') > 0, 'en honesty phrasing survives');
// 无版本标注的记录（mt：AADR 不发布 mt 树版本）如实说"未标注"，不造一个
const mtOnly = lineageObsPoints({ history: { observations: [
  { record_id: 'M1', relation: 'pending_review', source_version: '', coordinates: {} }] } });
ok(lineagePendingNote({tree_source:'PhyloTree',tree_version:'rcrs@17.2'}, mtOnly, true).indexOf('未标注') > 0,
   'unlabelled source version is reported as such, never invented');
// 点集中在云南–青海–北京一带：取景应放大而非整幅世界（复用 geomapViewport 的同一份几何）。
const vp = geomapViewport(op.points);
ok(vp.w < 360 && vp.h < 180, `observation cluster zooms in (got ${vp.w}x${vp.h})`);

// —— 树时间线节点：formed 优先；formed=0（YFull 未给，如 K-M2335）退 TMRCA 并如实记基准；
//    两者皆无（mt 的 A13 只有 n_defining_sites）不进时间线。
const nodes = lineageTimelineNodes([
  { node: 'K-M2335', formed: 0, tmrca: 41500 },
  { node: 'N', formed: 36800, tmrca: 22000 },
  { node: 'A13', n_defining_sites: 22 },
  { node: 'N-CTS4714', formed: 6000, tmrca: 5700 },
]);
ok(nodes.length === 3, `nodes without any date are excluded (got ${nodes.length})`);
ok(nodes[0].node === 'K-M2335' && nodes[0].year === 41500 && nodes[0].basis === 'tmrca',
   'formed=0 falls back to TMRCA with basis recorded');
ok(nodes[1].basis === 'formed' && nodes[1].year === 36800, 'formed wins when present');
ok(nodes.every(n => n.basis === 'formed' || n.basis === 'tmrca'), 'every node declares its time basis');
ok(lineageTimelineNodes([]).length === 0 && lineageTimelineNodes(null).length === 0, 'empty path -> empty timeline');

// —— 生产路径消费（drawGeoMap 的教训：纯函数可以是死代码）与 §7 诚实约束。
const cardSrc = src.slice(src.indexOf('renderLineageCards'));
ok(/lineageObsPoints\(lin\)/.test(cardSrc), 'lineage card derives its points from lineageObsPoints');
ok(/lineagePendingNote\(lin,op,Z\)/.test(cardSrc), 'lineage card renders the pending-review note (P0-2b)');
ok(/lineageTimelineNodes\(path\)/.test(cardSrc), 'lineage timeline consumes lineageTimelineNodes');
ok(/obsMap\(op\.points,kind\)/.test(cardSrc), 'observation map is drawn from the placed points');
ok(/clip-lineage-\$\{kind\}/.test(src), 'mini map clips via a per-kind clipPath id');
ok(/href="#world_land"/.test(cardSrc) && /transform="translate\(\$\{ox/.test(cardSrc),
   'base map is positioned via use transform (negative x/y is dropped by some engines, 845ead2)');
ok(!/x="\$\{ox/.test(cardSrc), 'use must not position via interpolated x/y inside the lineage card');
ok(src.includes('no sourced route; distribution only'), 'no-route honesty line survives');
ok(/报告年代/.test(src) && /formed\/TMRCA/.test(src), 'reported dates and tree estimates are labelled as separate bases');

console.log(FAIL.length ? 'FAIL\n' + FAIL.join('\n') : 'PASS');
process.exit(FAIL.length ? 1 : 0);
"""


@unittest.skipIf(shutil.which("node") is None, "node not installed here (runs on the server)")
class TestLineageCards(unittest.TestCase):
    def test_obs_points_and_timeline_nodes(self):
        r = subprocess.run(["node", "-e", HARNESS, "--", str(JS)],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_source_wiring_static(self):
        """接线断言不依赖 node：纯函数定义存在、被卡片消费（防止提交只带测试没带实现）。"""
        src = JS.read_text(encoding="utf-8")
        self.assertRegex(src, r"function lineageObsPoints\(doc\)\{")
        self.assertRegex(src, r"function lineageTimelineNodes\(path\)\{")
        self.assertRegex(src, r"function lineagePendingNote\(lin,op,zh\)\{")
        self.assertIn("lineageObsPoints(lin)", src)
        self.assertIn("lineageTimelineNodes(path)", src)
        self.assertIn("lineagePendingNote(lin,op,Z)", src)


if __name__ == "__main__":
    unittest.main()
