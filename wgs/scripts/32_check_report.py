#!/usr/bin/env python
"""Static acceptance checks for the assembled WGS report (standard library only).

Usage: python3 32_check_report.py [report_data.json] [report.html]

Checks (P0_FIX_PLAN section 6, condensed):
  1. example-sample constants are gone from wgs/ (source scan, comments excluded)
  2. every D.<key> the template reads exists in the JSON (small allow-list for optional keys)
  3. the rendered HTML contains no NaN / undefined / null / unfilled {placeholder}
  4. key-shape invariants the figures depend on
  5. reference naming: chr-prefixed names must not appear in the key tables
  6. per-section status table (H1): status, reason code and bilingual reason are consistent
  7. the report JS runs to completion under node with zero failed figures

Exit code 0 = pass, 1 = at least one failure (details on stderr).
"""
import json, math, os, pathlib, re, shutil, subprocess, sys, tempfile

FAIL = []

# Minimal DOM stub so the report script can be executed outside a browser. It is deliberately
# dumb (every element is a fresh object) -- enough to catch a figure that throws at run time,
# which the static checks above cannot see.
DOM_STUB = r"""
const fakeEl = () => ({
  innerHTML: '', textContent: '', title: '',
  style: { setProperty(){}, cursor: '' },
  appendChild(){}, addEventListener(){}, removeEventListener(){},
  setAttribute(){}, getAttribute: () => null, removeAttribute(){},
  insertAdjacentHTML(){}, querySelectorAll: () => [], querySelector: () => null,
  classList: { add(){}, remove(){}, toggle(){} },
  getBBox: () => ({ x: 0, y: 0, width: 100, height: 10 }),
  getBoundingClientRect: () => ({ x: 0, y: 0, width: 900, height: 300 }),
  closest: () => null, contains: () => false, focus(){}, blur(){}, parentNode: null,
});
global.__boot = fakeEl();
global.window = { __renderErrors: [], matchMedia: null, addEventListener(){}, devicePixelRatio: 1,
                  location: { href: 'file://report.html' } };
global.document = {
  documentElement: { setAttribute(){}, style: { setProperty(){} } },
  body: fakeEl(),
  getElementById: (id) => (id === 'bootstate' ? global.__boot : fakeEl()),
  querySelectorAll: (sel) => (sel === '[data-i18n]' ? (global.__i18nEls || []) : []),
  querySelector: () => null,
  createElementNS: () => fakeEl(),
  createElement: () => fakeEl(),
  addEventListener(){},
  getElementsByTagName: () => [],
  title: '',
};
global.getComputedStyle = () => ({ getPropertyValue: () => '#123456' });
global.localStorage = { getItem: () => null, setItem(){} };
global.requestAnimationFrame = (cb) => 0;
global.navigator = { userAgent: 'node', language: process.env.PROBE_LANG || 'zh-CN' };
"""

NODE_RUNNER = """require(process.argv[2]);
const fs = require('fs');
// The page's copy lives in data-i18n attributes and is filled by renderAll through t(), which
// indexes UI[key][zh ? 0 : 1]. Stand in for those elements, or that loop is a no-op.
const html = fs.readFileSync(process.argv[4], 'utf8');
global.__i18nEls = [...new Set([...html.matchAll(/data-i18n="([^"]+)"/g)].map(m => m[1]))].map(k => ({
  __i18nKey: k, innerHTML: '', textContent: '',
  getAttribute: (n) => (n === 'data-i18n' ? k : null), setAttribute(){}, classList: { add(){}, remove(){} },
}));
console.log('I18N_KEYS=' + global.__i18nEls.length);
try { eval(fs.readFileSync(process.argv[3], 'utf8')); }
catch (e) { console.log('TOP-LEVEL ERROR: ' + e.message); process.exit(2); }
const errs = (global.window && global.window.__renderErrors) || [];
console.log('RENDER_ERRORS=' + JSON.stringify(errs));
console.log('BOOTSTATE=' + ((global.__boot && global.__boot.textContent) || ''));
"""

FAIL = []


def fail(msg):
    FAIL.append(msg)


def _code_of_js_line(line):
    line = re.sub(r"/\*.*?\*/", "", line)
    return line.split("//", 1)[0]


def check_constants(root):
    pat = re.compile(r"28\.88|4850\.\d|6\.52 ?亿|504 名|208 名|B\*13:01:01|O-F438|O-M122|6,919|2\.87 Gb")
    for p in sorted((root / "scripts").glob("*.py")):
        if p.name == "32_check_report.py":
            continue  # the checker carries the pattern list itself
        for i, l in enumerate(p.read_text(errors="ignore").splitlines(), 1):
            if pat.search(l.split("#", 1)[0]):
                fail(f"example constant {p.name}:{i}: {l.strip()[:90]}")
    for p in sorted((root / "templates").glob("*.js")):
        for i, l in enumerate(p.read_text(errors="ignore").splitlines(), 1):
            if pat.search(_code_of_js_line(l)):
                fail(f"example constant {p.name}:{i}: {l.strip()[:90]}")


def check_keys(root, D):
    js = (root / "templates/report_script.js").read_text(errors="ignore")
    keys = set(re.findall(r"\bD\.([A-Za-z_][A-Za-z0-9_]*)", js))
    optional = {"cnv", "sv_gene_dels"}  # template tolerates these when absent
    for k in sorted(keys - optional):
        if k not in D:
            fail(f"template reads D.{k} but report_data.json has no such key")


def check_payload_escaping(html_path):
    """载荷块里不能出现裸的闭合序列。

    31 把数据写成 JS 字符串字面量；如果数据里含 "</script"（例如某人的显示名），浏览器的解析器
    会在那里结束脚本块，剩下的载荷变成页面文本，图表全空 —— 而所有静态检查仍然通过。转义后
    文本里只有反斜杠 + 斜杠，数据本身不变。
    """
    if not html_path or not html_path.exists():
        return
    html = html_path.read_text(encoding="utf-8")
    blocks = re.findall(r"<script>window\.__c_\w+_\d+=(.*?)</script>", html, re.S)
    if not blocks:
        return          # 没有分块载荷的构建（例如只渲染片段）不适用
    bad = [b[:60] for b in blocks if "</" in b]
    if bad:
        fail(f"{len(bad)} payload chunk(s) contain a raw closing tag sequence, e.g. {bad[0]!r}")
    if "NaN" in "".join(blocks) or "Infinity" in "".join(blocks):
        fail("payload contains NaN/Infinity; it must be valid JSON (null for missing values)")


def check_html(html_path):
    if not html_path or not html_path.exists():
        return
    t = re.sub(r"<script.*?</script>", "", html_path.read_text(errors="ignore"), flags=re.S)
    for bad in ("NaN", "undefined", "null", "{ver}", "{cmp}"):
        if re.search(r"(?<![A-Za-z])" + re.escape(bad) + r"(?![A-Za-z])", t):
            fail(f"rendered HTML contains '{bad}'")


def check_print_and_narrow(html_path):
    """打印与窄屏（AN6）：交付的是单文件 HTML，读者常直接打印成 PDF。

    深色模式系统下打印若不强制亮底，会印出整页深蓝；过高的图不限制高度会只印出一半。
    这些都在生成产物里可静态断言，比一次性人工试打更可靠。
    """
    text = pathlib.Path(html_path).read_text(encoding="utf-8")
    for pat, why in ((r"@media\s*print", "a print media query"),
                     (r"break-inside:\s*avoid", "break-inside avoidance for figures and tables"),
                     (r"background:\s*#fff\s*!important", "forced light background when printing"),
                     (r"max-height:\s*\d+mm", "a max height so tall figures fit one page"),
                     (r"@media\s*\(max-width:560px\)", "a narrow-screen breakpoint")):
        if not re.search(pat, text):
            fail(f"report.html is missing {why} (pattern {pat})")


def check_shapes(D):
    if "chip_hotspots" in D and "alt_reads" not in D["chip_hotspots"]:
        fail("chip_hotspots.alt_reads missing (F.chipalt depends on it)")
    for k in ("la_segments",):
        for seg in (D.get(k) or [])[:50]:
            if isinstance(seg, dict) and "hap" not in seg:
                fail(f"{k} rows lack the 'hap' field the painting/circos figures need")
                break
    if D.get("roh_stats", {}).get("n_gt5") not in (None, 0) and not D.get("roh"):
        fail("roh_stats.n_gt5 > 0 without any roh rows")
    # M2: every assayed STR locus carries a class; one panel number is not a per-locus
    # "pathogenic" verdict, and loci without a rule must say so rather than pass as normal.
    str_classes = {"normal", "intermediate", "premutation", "full_mutation",
                   "below_threshold", "at_or_above_threshold", "no_local_rule"}
    for r in (D.get("str") or [])[:60]:
        if not isinstance(r, dict) or r.get("class") not in str_classes:
            fail(f"str row {r.get('locus') if isinstance(r, dict) else r!r}: missing/unknown class (got {r.get('class') if isinstance(r, dict) else None!r})")
            break
        if r.get("class") == "no_local_rule" and r.get("thr") is not None:
            fail(f"str row {r.get('locus')}: no_local_rule but thr is set")
    # M4: sv_gene_dels are whole-gene deletions (13's whole_gene_del column), one entry
    # per deleted gene; the event/gene counts must match the list they summarise.
    gd = D.get("sv_gene_dels") or []
    st = D.get("sv_gene_dels_stats") or {}
    if gd:
        if not st:
            fail("sv_gene_dels present without sv_gene_dels_stats (events/genes counts)")
        else:
            if st.get("genes") != len({e["gene"] for e in gd}):
                fail("sv_gene_dels_stats.genes does not match the unique genes in sv_gene_dels")
            if not (0 < st.get("events", 0) <= len(gd)):
                fail("sv_gene_dels_stats.events inconsistent with sv_gene_dels length")
        for e in gd[:20]:
            if not e.get("gene") or "," in e["gene"]:
                fail(f"sv_gene_dels entries must carry one gene name each, got {e.get('gene')!r}")
                break
    # M5: chip_hotspots must carry the non-zero per-hotspot detail (strongest first) and
    # its alt_reads total must match that list; the KPI states an observation, not a CHIP call.
    ch = D.get("chip_hotspots") or {}
    if ch.get("median_depth"):
        hs = ch.get("hotspots")
        if not isinstance(hs, list):
            fail("chip_hotspots lacks the per-hotspot list the KPI caption reads")
        else:
            if any(h.get("alt_reads", 0) <= 0 for h in hs):
                fail("chip_hotspots.hotspots lists entries with no alt reads")
            if sum(h.get("alt_reads", 0) for h in hs) != ch.get("alt_reads"):
                fail("chip_hotspots.alt_reads does not equal the sum of its per-hotspot list")
            if len(hs) > 1 and any(hs[i]["alt_reads"] < hs[i + 1]["alt_reads"] for i in range(len(hs) - 1)):
                fail("chip_hotspots.hotspots is not sorted strongest-first")
    # M5: an insufficient_reads telomere verdict must disclose the working 100000 gate.
    tel = next((s for s in (D.get("sections") or []) if s.get("id") == "telomere"), None)
    if tel and tel.get("code") == "insufficient_reads" and "100000" not in (tel.get("reason_zh") or "") + (tel.get("reason_en") or ""):
        fail("telomere insufficient_reads reason does not disclose the 100000 working gate")
    # H2: template contract. The mutation-spectrum figure colours bars via d.sub, the PRS
    # tooltip scales on r.pct_Han, the Human-Origins PCA indexes ho_modern rows positionally,
    # and the circos ring steps density counts by 5 Mb -- any drift silently blanks a figure.
    sub6 = {"C>A", "C>G", "C>T", "T>A", "T>C", "T>G"}
    for r in (D.get("spectrum") or [])[:96]:
        if r.get("sub") not in sub6:
            fail(f"spectrum row {r.get('ctx')}: sub missing or not a canonical substitution (got {r.get('sub')!r})")
            break
    for p in (D.get("prs") or [])[:5]:
        if "pct_Han" not in p:
            fail("prs rows lack pct_Han (template tooltip scales on it; step 12 calls it pct_sub)")
            break
    for p in (D.get("ho_modern") or [])[:5]:
        if not (isinstance(p, (list, tuple)) and len(p) == 3):
            fail("ho_modern rows must be [label, pc1, pc2] arrays -- the PCA indexes them positionally")
            break
    for d in (D.get("density") or [])[:25]:
        cl = D.get("chrlen", {}).get(d.get("chrom"), 0)
        if cl and not (abs(len(d.get("counts") or []) - (cl // 5_000_000 + 1)) <= 1):
            fail(f"density counts for {d.get('chrom')} are not 5 Mb bins ({len(d['counts'])} vs chrlen/{5_000_000:.0e})")
            break


    # W-T1/AN5: the affinity figure draws D.ho_affinity (modern + ancient on one distance axis) and the
    # three closest ancient groups with their members. These checks stay independent of any particular
    # reference panel: an earlier version demanded China_-prefixed labels, a north/south vocabulary and
    # a Han name_zh, which would fail for every non-Chinese panel and for the English rendering.
    aff = D.get("ho_affinity") or []
    # AN7: a build with the ancient-DNA projection disabled legitimately has no ho_affinity -- the
    # requirement applies only when the aadr section itself claims the analysis ran and is ok.
    # Unconditionally demanding it failed every non-AADR configuration as a false positive.
    _aadr = next((s for s in (D.get("sections") or []) if s.get("id") == "aadr"), None)
    if not aff:
        if _aadr and _aadr.get("status") == "ok":
            fail("ho_affinity missing or empty although the aadr section is ok (the ancestry figure has nothing to draw)")
    else:
        ds = [r.get("d") for r in aff]
        if ds != sorted(ds):
            fail("ho_affinity is not sorted by d")
        if not any(r.get("kind") == "modern" for r in aff):
            fail("ho_affinity carries no modern rows; the nearest-modern anchor would vanish")
        for r in aff:
            if r.get("kind") == "ancient" and int(r.get("n") or 0) < 2:
                fail(f"ho_affinity ancient row {r.get('label')} has n<2")
                break
            # region is a free-form grouping label from the location table (a colour key): it only has
            # to be a string. The vocabulary belongs to the reference panel, not to this checker.
            if r.get("region") is not None and not isinstance(r.get("region"), str):
                fail(f"ho_affinity {r.get('label')} has a non-string region {r.get('region')!r}")
                break
            # name_zh, when present, is whatever the location table holds. No script requirement: 7.1
            # says a record is not required to have a Han name, and an English report shows raw labels.
            nz = r.get("name_zh")
            if nz is not None and not isinstance(nz, str):
                fail(f"ho_affinity {r.get('label')} has a non-string name_zh {nz!r}")
                break
            if r.get("d") is not None and float(r.get("d") or 0) < 0:
                fail(f"ho_affinity {r.get('label')} has a negative distance")
                break
        strip = D.get("ho_affinity_strip") or []
        top3 = [x.get("label") for x in sorted([r for r in aff if r.get("kind") == "ancient"],
                                               key=lambda r: r["d"])[:3]]
        if [s2.get("label") for s2 in strip] != top3:
            fail(f"ho_affinity_strip {[s2.get('label') for s2 in strip]} != the 3 closest ancient groups {top3}")
        elif not all(s2.get("members") for s2 in strip):
            fail("ho_affinity_strip rows must carry their members")
        tg = D.get("ho_target_group")
        if tg and not (tg.get("members") and tg.get("mean_d") is not None):
            fail("ho_target_group must carry its members and its mean distance")
    # 28_f3_stats payload. Optional as a whole (null = step not run), but when present the
    # modern ranking must be sorted desc, the contrast list must reference exactly the
    # adjacent pairs of that ranking (the figure annotates them by position), and every
    # row's f3/se/z is either a finite number or an explicit null = "unavailable" (zero-SE /
    # non-finite; the card renders those as 不可用). NaN/Infinity parse as floats, so the
    # isfinite check is what actually keeps them out -- isinstance alone would wave them
    # through as numbers.
    f3 = D.get("f3")
    if f3:
        for side, rows_key in (("modern", "outgroup_f3"),):
            blk = f3.get(side) or {}
            rows = blk.get(rows_key) or []
            if not rows:
                continue
            for r in rows:
                for k in ("f3", "se", "z"):
                    v = r.get(k)
                    if v is not None and (not isinstance(v, (int, float)) or not math.isfinite(v)):
                        fail(f"f3 {side} row {r.get('set')}: {k} is neither a finite number nor null")
                if r.get("z") is not None and (r.get("f3") is None or r.get("se") is None):
                    fail(f"f3 {side} row {r.get('set')}: a Z without its f3/SE is not a statistic")
                if not (r.get("n") and r.get("label_zh") and r.get("label_en")):
                    fail(f"f3 {side} row {r.get('set')}: n and bilingual labels are required")
            vals = [r["f3"] for r in rows if isinstance(r.get("f3"), (int, float))]
            if len(vals) >= 2 and vals != sorted(vals, reverse=True):
                fail("f3 modern outgroup_f3 is not sorted descending")
            got = [(c.get("a"), c.get("b")) for c in (blk.get("contrasts") or [])]
            want = [(rows[i]["set"], rows[i + 1]["set"]) for i in range(len(rows) - 1)]
            if got != want:
                fail(f"f3 contrasts {got} do not match the adjacent ranking pairs {want}")
        for label, blk in (("modern", f3.get("modern") or {}), ("ancient", f3.get("ancient") or {})):
            for r in blk.get("admixture") or []:
                for k in ("f3", "se", "z"):
                    v = r.get(k)
                    if v is not None and (not isinstance(v, (int, float)) or not math.isfinite(v)):
                        fail(f"f3 {label} admixture row {r.get('a')}x{r.get('b')}: {k} is neither a finite number nor null")
                if not (r.get("a") and r.get("b")):
                    fail(f"f3 {label} admixture row lacks its source names")
    # KIR (step 24 via 30): gene rows and the haplotype/ligand summary must agree. Rows
    # without a summary made the old card invent a haplotype; a summary claiming AA/Bx with
    # no rows is stale output from another run. "unknown" and a genuine absent row are
    # different states and must not be swapped (M3/H2).
    kir_rows = D.get("kir") or []
    ks = D.get("kir_summary") or {}
    hap = ks.get("haplotype", "")
    if kir_rows and not hap:
        fail("kir rows present but kir_summary is missing: the haplotype/ligand line has no source")
    if not kir_rows and hap and not hap.lower().startswith("unavailable"):
        fail(f"kir_summary claims {hap!r} but the gene table is empty: typing did not run")
    if kir_rows and hap and hap.lower().startswith("unavailable"):
        fail("kir rows present but the summary says typing was unavailable")


def check_sections(D):
    """H1: every section carries a status; 'unavailable' must explain itself bilingually."""
    secs = D.get("sections")
    if not isinstance(secs, list) or not secs:
        fail("sections status table missing from report data")
        return
    ids = [s.get("id") for s in secs]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        fail(f"duplicate section ids: {sorted(dup)}")
    # AN7：词表必须等于 30 实际会写进 sections 的全部 reason code，缺一个就会把合法的
    # 结构化空态（如未运行可选模块的 missing_manifest）判成"未知码"。来源：
    #   30 的 _MISS_IN/_MISS_OUT/_MISS_STEP（no_input/missing_output/step_not_run）
    #   + 30 的 _LA_CODES（manifest 准入状态）+ run_all skipped() 的禁用理由
    #   + insufficient_reads（23 的血型低读段）。收紧或扩充词表时 30 与这里必须同步。
    known_codes = {"no_input", "missing_output", "step_not_run", "insufficient_reads",
                   "missing_manifest", "stale_result", "missing_result", "unreadable_result",
                   "disabled_by_config", "no_segments",
                   "local_ancestry_not_configured", "aadr_not_configured", "regional_axis_not_configured"}
    for s in secs:
        st = s.get("status")
        if st not in ("ok", "negative", "unavailable"):
            fail(f"section {s.get('id')}: bad status {st!r}")
        if st == "unavailable" and not (s.get("code") and s.get("reason_zh") and s.get("reason_en")):
            fail(f"section {s.get('id')}: unavailable without a code and bilingual reason")
        if st == "negative" and not s.get("evidence"):
            fail(f"section {s.get('id')}: negative without an evidence path")
        if s.get("code") and s["code"] not in known_codes:
            fail(f"section {s.get('id')}: unknown reason code {s['code']!r}")
        if not (s.get("name_zh") and s.get("name_en")):
            fail(f"section {s.get('id')}: missing bilingual display name (renderSections reads it)")


def check_sv_source(run_dir, D):
    """M4: every sv_gene_dels gene must come from 13's whole_gene_del column -- a gene the
    event merely overlaps is not a whole-gene deletion, however low the read depth is.
    run_dir is the outputs root of the report being checked (the directory holding its
    report_data.json), so a report built under a different work tree is checked against
    its own sources, not this repository's."""
    tsv = run_dir / "08_sv/sv_filtered.tsv"
    if not tsv.exists() or not (D.get("sv_gene_dels") or []):
        return
    import csv
    wgd = set()
    with open(tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for g in (row.get("whole_gene_del") or "").split(","):
                if g:
                    wgd.add(g)
    stray = {e["gene"] for e in D["sv_gene_dels"]} - wgd
    if stray:
        fail(f"sv_gene_dels names genes absent from sv_filtered whole_gene_del (overlap-only?): {sorted(stray)[:5]}")


def check_naming(run_dir):
    """Main chromosomes must use the reference naming (1..22,X,Y,MT); a chr-prefixed main
    contig silently breaks region queries everywhere downstream. Decoy/unplaced names in
    vendor style (chrUn_*, *_random, GL*) are tolerated. run_dir follows the report being
    checked, like check_sv_source."""
    import subprocess
    ok = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    bad_main = re.compile(r"^chr(?:[1-9][0-9]?|1[0-9]|2[0-2]|X|Y|M|MT)$")
    decoy = re.compile(r"(_random$|^chrUn|^GL0)", re.I)
    entries = []
    vcf = run_dir / "00_input/target.pass.vcf.gz"
    if vcf.exists():
        out = subprocess.run(["bcftools", "index", "-s", str(vcf)], capture_output=True, text=True).stdout
        entries += [("target.pass.vcf.gz", l.split("\t")[0]) for l in out.splitlines()]
    bed = run_dir / "00_input/callable.bed"
    if bed.exists():
        entries += [("callable.bed", l.split("\t")[0]) for l in open(bed)]
    seen = set()
    for src, c in entries:
        if not c or c in seen:
            continue
        seen.add(c)
        if bad_main.match(c):
            fail(f"chr-prefixed main contig '{c}' in {src} (region queries will silently miss)")
        elif c not in ok and not decoy.search(c):
            fail(f"unexpected chromosome name '{c}' in {src}")


def check_runtime(html_path):
    """Run the report JS under node with a minimal DOM stub; no figure may fail.

    Static checks cannot see a figure that throws at run time -- a circos crash on an
    X-chromosome gene deletion shipped past all of them. Nor can they see a failure that
    only happens in one language: the copy is indexed as UI[key][zh ? 0 : 1], so a malformed
    English entry throws inside t() and blanks the page in English alone. Run both.
    A missing node is a failure, not a skip: "跳过不可计通过" -- exiting 0 here let a
    machine without node hand over reports whose figures had never executed anywhere.
    """
    node = shutil.which("node")
    if not node:
        fail("node is not installed; the runtime figure check cannot run, and a skipped check must not count as passed")
        return
    if html_path is None or not html_path.exists():
        fail("report.html not found for the runtime check")
        return
    html = html_path.read_text(encoding="utf-8")
    # The page carries several <script> blocks (boot diagnostic, payload, logic, progress markers).
    # They share one global scope, so concatenating them reproduces what the browser executes.
    blocks = re.findall(r"<script>(.*?)</script>", html, re.S)
    if not blocks:
        fail("report.html has no inline <script> block")
        return
    script = "\n".join(blocks)
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td)
        (p / "stub.js").write_text(DOM_STUB, encoding="utf-8")
        (p / "runner.js").write_text(NODE_RUNNER, encoding="utf-8")
        (p / "report.js").write_text(script, encoding="utf-8")
        runs = {}
        for lang in ("zh-CN", "en-US"):
            try:
                r = subprocess.run([node, str(p / "runner.js"), str(p / "stub.js"), str(p / "report.js"),
                                    str(html_path)],
                                   capture_output=True, text=True, timeout=180,
                                   env={**os.environ, "PROBE_LANG": lang})
            except subprocess.TimeoutExpired:
                fail(f"the report JS did not finish within 180 s under node ({lang})")
                return
            runs[lang] = (r.stdout or "") + (r.stderr or "")
    for lang, out in runs.items():
        tag = " (zh)" if lang == "zh-CN" else " (en)"
        nk = re.search(r"I18N_KEYS=(\d+)", out)
        if not nk or int(nk.group(1)) == 0:
            fail(f"the runtime harness saw no data-i18n elements{tag}; the copy check would be vacuous")
        if "TOP-LEVEL ERROR" in out:
            line = next((l for l in out.splitlines() if "TOP-LEVEL ERROR" in l), out[-200:])
            fail(f"report JS throws at top level{tag}: {line}")
        elif "RENDER_ERRORS=[]" not in out:
            m = re.search(r"RENDER_ERRORS=(\[.*\])", out)
            fail(f"figures failed to render{tag}: {m.group(1) if m else out.strip()[-200:]}")
        elif "RENDER FAILED" in (re.search(r"BOOTSTATE=(.*)", out) or [None, ""])[1]:
            banner = re.search(r"BOOTSTATE=(.*)", out).group(1).strip()[:180]
            fail(f"the page banner reports a failed render{tag}: {banner}")


def main():
    root = pathlib.Path(__file__).resolve().parents[1]
    # The outputs root of the report being checked: report_data.json lives directly in it
    # (work/wgs). Defaulting to this repository's own tree kept only the no-argument case
    # honest -- a report passed from another work tree was still checked against this repo's
    # sources, so the source checks passed or failed on the wrong data (AN7).
    json_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else root / "work/wgs/report_data.json"
    run_dir = json_path.parent
    html_path = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else root / "work/report/report.html"

    check_constants(root)
    D = {}
    if json_path.exists():
        D = json.load(open(json_path))
        check_keys(root, D)
        check_shapes(D)
        check_sections(D)
        check_sv_source(run_dir, D)
    else:
        fail(f"{json_path} not found")
    check_html(html_path)
    check_payload_escaping(html_path)
    check_print_and_narrow(html_path)
    check_naming(run_dir)
    check_runtime(html_path)

    if FAIL:
        print("FAILED:", file=sys.stderr)
        for f in FAIL:
            print("  -", f, file=sys.stderr)
        return 1
    print("32_check_report: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
