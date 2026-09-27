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
import json, pathlib, re, shutil, subprocess, sys, tempfile

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
global.window = { __renderErrors: [], matchMedia: null, addEventListener(){}, devicePixelRatio: 1,
                  location: { href: 'file://report.html' } };
global.document = {
  documentElement: { setAttribute(){}, style: { setProperty(){} } },
  body: fakeEl(),
  getElementById: () => fakeEl(),
  querySelectorAll: () => [],
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
global.navigator = { userAgent: 'node' };
"""

NODE_RUNNER = """require(process.argv[2]);
const fs = require('fs');
try { eval(fs.readFileSync(process.argv[3], 'utf8')); }
catch (e) { console.log('TOP-LEVEL ERROR: ' + e.message); process.exit(2); }
const errs = (global.window && global.window.__renderErrors) || [];
console.log('RENDER_ERRORS=' + JSON.stringify(errs));
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


def check_html(html_path):
    if not html_path or not html_path.exists():
        return
    t = re.sub(r"<script.*?</script>", "", html_path.read_text(errors="ignore"), flags=re.S)
    for bad in ("NaN", "undefined", "null", "{ver}", "{cmp}"):
        if re.search(r"(?<![A-Za-z])" + re.escape(bad) + r"(?![A-Za-z])", t):
            fail(f"rendered HTML contains '{bad}'")


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
    known_codes = {"no_input", "missing_output", "step_not_run", "insufficient_reads"}
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


def check_sv_source(root, D):
    """M4: every sv_gene_dels gene must come from 13's whole_gene_del column -- a gene the
    event merely overlaps is not a whole-gene deletion, however low the read depth is."""
    tsv = root / "work/wgs/08_sv/sv_filtered.tsv"
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


def check_naming(root):
    """Main chromosomes must use the reference naming (1..22,X,Y,MT); a chr-prefixed main
    contig silently breaks region queries everywhere downstream. Decoy/unplaced names in
    vendor style (chrUn_*, *_random, GL*) are tolerated."""
    import subprocess
    ok = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    bad_main = re.compile(r"^chr(?:[1-9][0-9]?|1[0-9]|2[0-2]|X|Y|M|MT)$")
    decoy = re.compile(r"(_random$|^chrUn|^GL0)", re.I)
    run = root / "work/wgs"
    entries = []
    vcf = run / "00_input/target.pass.vcf.gz"
    if vcf.exists():
        out = subprocess.run(["bcftools", "index", "-s", str(vcf)], capture_output=True, text=True).stdout
        entries += [("target.pass.vcf.gz", l.split("\t")[0]) for l in out.splitlines()]
    bed = run / "00_input/callable.bed"
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
    X-chromosome gene deletion shipped past all of them. Skipped, with a note, when node
    is unavailable.
    """
    node = shutil.which("node")
    if not node:
        print("note: node was not found; the runtime figure check was skipped")
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
        try:
            r = subprocess.run([node, str(p / "runner.js"), str(p / "stub.js"), str(p / "report.js")],
                               capture_output=True, text=True, timeout=180)
        except subprocess.TimeoutExpired:
            fail("the report JS did not finish within 180 s under node")
            return
    out = (r.stdout or "") + (r.stderr or "")
    if "TOP-LEVEL ERROR" in out:
        line = next((l for l in out.splitlines() if "TOP-LEVEL ERROR" in l), out[-200:])
        fail(f"report JS throws at top level: {line}")
    elif "RENDER_ERRORS=[]" not in out:
        m = re.search(r"RENDER_ERRORS=(\[.*\])", out)
        fail(f"figures failed to render: {m.group(1) if m else out.strip()[-200:]}")


def main():
    root = pathlib.Path(__file__).resolve().parents[1]
    json_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else root / "work/wgs/report_data.json"
    html_path = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else root / "work/report/report.html"

    check_constants(root)
    D = {}
    if json_path.exists():
        D = json.load(open(json_path))
        check_keys(root, D)
        check_shapes(D)
        check_sections(D)
        check_sv_source(root, D)
    else:
        fail(f"{json_path} not found")
    check_html(html_path)
    check_naming(root)
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
