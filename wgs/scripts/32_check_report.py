#!/usr/bin/env python
"""Static acceptance checks for the assembled WGS report (standard library only).

Usage: python3 32_check_report.py [report_data.json] [report.html]

Checks (P0_FIX_PLAN section 6, condensed):
  1. example-sample constants are gone from wgs/ (source scan, comments excluded)
  2. every D.<key> the template reads exists in the JSON (small allow-list for optional keys)
  3. the rendered HTML contains no NaN / undefined / null / unfilled {placeholder}
  4. key-shape invariants the figures depend on
  5. reference naming: chr-prefixed names must not appear in the key tables

Exit code 0 = pass, 1 = at least one failure (details on stderr).
"""
import json, pathlib, re, sys

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


def main():
    root = pathlib.Path(__file__).resolve().parents[1]
    json_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else root / "work/wgs/report_data.json"
    html_path = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else None

    check_constants(root)
    D = {}
    if json_path.exists():
        D = json.load(open(json_path))
        check_keys(root, D)
        check_shapes(D)
    else:
        fail(f"{json_path} not found")
    check_html(html_path)
    check_naming(root)

    if FAIL:
        print("FAILED:", file=sys.stderr)
        for f in FAIL:
            print("  -", f, file=sys.stderr)
        return 1
    print("32_check_report: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
