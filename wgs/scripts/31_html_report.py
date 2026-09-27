#!/usr/bin/env python
"""Report v3 (whole-genome): bilingual single-file HTML. Numbers come from wgs/report_data_v3.json; text lives here."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import json, pathlib, re
P = P; S = pathlib.Path(__file__).resolve().parents[1]/"templates"
D = json.load(open(P/"wgs"/"report_data.json"))
OUT = REPORT/"report.html"; OUT.parent.mkdir(parents=True, exist_ok=True)

# ---- report copy: work/report_text.yaml when it exists, otherwise the bundled example (with a warning)
import yaml as _yaml
_cand = [P/"report_text.yaml", pathlib.Path(__file__).resolve().parents[1]/"example"/"report_text.yaml"]
_src = next(f for f in _cand if f.exists())
if _src != _cand[0]:
    print(f"WARNING: rendering with the bundled example copy at {_src}.\n"
          f"         Copy it to {_cand[0]} and rewrite it for this sample before delivering.", flush=True)
_TXT = _yaml.safe_load(open(_src, encoding="utf-8"))
# When the sample copy is missing and the bundled example is used, the delivered HTML must
# say so on its face -- a stdout warning is invisible to whoever opens the report file.
_EXAMPLE_BANNER = "" if _src == _cand[0] else (
    '<div style="margin:0 auto;max-width:900px;padding:10px 14px;border:1px solid #b45309;'
    'border-radius:8px;color:#b45309;background:#fffbeb;font:13px/1.6 system-ui">'
    '本报告的叙述文案来自<b>示例模板</b>，尚未针对该样本改写；数字由数据文件自动填充，'
    '但所有解释性文字在交付前必须逐段核对。 / This report&apos;s narrative text is the '
    '<b>bundled example copy</b>, not sample-specific wording; figures are data-driven but '
    'every explanatory paragraph must be reviewed before delivery.</div>')
_FILL = {"name_en": NAME_EN, "name_zh": NAME_ZH, "sample": SAMPLE,
         **{k: v for k, v in D.items() if isinstance(v, (str, int, float))}}
class _Keep(dict):
    def __missing__(self, k): return "{" + k + "}"
def _fmt(x):
    if isinstance(x, str):
        try: return x.format_map(_Keep(_FILL))
        except (IndexError, ValueError): return x
    if isinstance(x, list): return [_fmt(i) for i in x]
    if isinstance(x, dict): return {k: _fmt(v) for k, v in x.items()}
    return x
_TXT = _fmt(_TXT)
UI = _TXT["ui"]; FIND = _TXT["findings"]; PGX = _TXT["pgx"]; BLOOD = _TXT["blood"]; PRS_ZH = _TXT["prs_names"]






UI["sub_y"] = [s.replace("{ver}", D["versions"]["yfull"]) for s in UI["sub_y"]]

head = open(S/"report_head.html", encoding="utf-8").read()
head = re.sub(r"<title>.*?</title>", f"<title>{NAME_ZH}</title>", head, count=1)
body = open(S/"report_body.html", encoding="utf-8").read()

# --- build stamp + boot diagnostics. The report must work from file:// on a machine with no
# network, so a silent failure has to be impossible: the strip states whether the script ran,
# and if not, why not.
import datetime
_BUILD = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
_EARLY = """<script>
window.__earlyErrors = [];
window.addEventListener('error', function (e) { window.__earlyErrors.push((e.message || 'error') + ' @' + e.lineno + ':' + e.colno); });
document.addEventListener('DOMContentLoaded', function () {
  var b = document.getElementById('bootstate'); if (!b) return;
  if (b.textContent.indexOf('waiting') < 0) return;  // the report already reported itself -- do not overwrite it
  b.textContent = 'the report did not finish; last progress marker: ' + (window.__boot || 'none')
    + (window.__earlyErrors.length ? ' | errors: ' + window.__earlyErrors.join(' | ') : '');
  b.style.color = '#ffd0c8';
});
</script>
"""

# Progress markers: each one records how far the parser got, so a parse failure can be localised
# without a developer console.
_MARK = lambda s: '<script>window.__boot=' + json.dumps(s) + ';</script>\n'

# Last block: read every diagnostic back and show it in a fixed, unmissable banner whenever the
# report did not finish. The report can fail in ways that never reach the boot strip (a missing
# element, a failing innerHTML write), so this reports state instead of assuming it.
_SUMMARY = """<script>
(function () {
  var b = document.getElementById('bootstate');
  if (b && b.textContent.indexOf('rendered') >= 0) return;      // it worked -- stay out of the way
  var parts = [];
  parts.push('progress marker: ' + (window.__boot || 'none'));
  parts.push('all logic blocks reached: ' + (window.__codeDone || 'no'));
  parts.push('boot strip element: ' + (b ? 'found' : 'MISSING'));
  parts.push('document.title: ' + document.title);
  parts.push('early errors: ' + JSON.stringify(window.__earlyErrors || []));
  parts.push('render errors: ' + JSON.stringify(window.__renderErrors || []));
  parts.push('document.body: ' + (document.body ? 'found' : 'MISSING')
             + ' | svg count: ' + document.getElementsByTagName('svg').length);
  var e = document.getElementById('errorbox');
  if (e) { e.textContent = parts.join('\\n'); e.style.display = 'block'; }
})();
</script>
"""
_BOOT = ('<div id="boot" style="position:fixed;top:0;left:0;right:0;z-index:99999;background:#111;color:#eee;'
         'font:12px/1.7 ui-monospace,Menlo,monospace;padding:5px 10px;text-align:center">build ' + _BUILD +
         ' &middot; <span id="bootstate">HTML parsed; waiting for the report script&hellip;</span></div>'
         '<div id="errorbox" style="display:none;position:fixed;top:32px;left:0;right:0;z-index:99999;'
         'background:#7f1d1d;color:#fff;font:13px/1.7 ui-monospace,Menlo,monospace;padding:10px 14px;'
         'white-space:pre-wrap;word-break:break-all"></div>'
         '<div style="height:30px"></div>')
head = _EARLY + head
body = body + _BOOT
js = open(S/"report_script.js", encoding="utf-8").read()

# --- payload assembly. Safari refused to parse the inline data as one big literal (634 KB on a
# single line, then 669 KB spread over many lines). Each object is now emitted as a series of
# ~100 KB string chunks in small <script> blocks and rebuilt with JSON.parse, and every step
# records a progress marker, so a failure names the object it died on.
_CHUNK = 100_000
_DATA_OBJS = [("D", D), ("BLOODV3", BLOOD), ("UI", UI), ("FIND", FIND), ("PGX", PGX), ("PRS_EN", PRS_ZH)]
_PH = ("__DATA__", "__UI__", "__FIND__", "__PGX__", "__BLOOD__", "__PRS_EN__")
_lines = js.split("\n")
_decl = {i for i, l in enumerate(_lines) if any(p in l for p in _PH)}
if not _decl:
    raise SystemExit("31_html_report: no payload placeholder found in report_script.js")
js_code = "\n".join(l for i, l in enumerate(_lines) if i not in _decl)   # logic only; payload is built below

# On Safari the whole 54 KB logic block silently did nothing -- no error, no effect -- while the
# blocks after it ran fine. The logic is therefore emitted as a series of small blocks too, each
# with its own progress marker and error report, so a failure names the exact part.
_CODE_CAP = 12000
def _code_blocks():
    parts, buf, size = [], [], 0
    for lineno, ln in enumerate(js_code.split("\n"), 1):
        buf.append(ln)
        size += len(ln) + 1
        # Only break at a top-level statement boundary: a line with no leading indent that ends a
        # statement or a block. Breaking anywhere else split an argument list, which JavaScriptCore
        # rejected with 'Unexpected keyword catch'.
        if ln[:1] not in (" ", "\t") and ln.rstrip().endswith((";", "}")):
            # Record the source line of the last completed statement. Safari reports only
            # 'Script error. @0:0' for these blocks, so "how far did it get" is the only signal.
            buf.append("window.__boot='L%d';" % lineno)
            if size >= _CODE_CAP:
                parts.append("\n".join(buf)); buf, size = [], 0
    if buf:
        parts.append("\n".join(buf))
    out = []
    for i, part in enumerate(parts):
        # No try/catch per block: `try` introduces a block scope, so a `const` defined in one block
        # would be invisible in the next (JavaScriptCore: 'reveal is not defined'). The marker before
        # each block plus the window error listener is enough to localise a failure.
        out.append(_MARK("code-%d/%d" % (i + 1, len(parts))))
        out.append("<script>\n" + part + "\n</script>\n")
    return "".join(out)

def _payload_blocks():
    out = []
    for name, obj in _DATA_OBJS:
        s = json.dumps(obj, ensure_ascii=False)
        chunks = [s[i:i + _CHUNK] for i in range(0, len(s), _CHUNK)] or [""]
        out.append(_MARK(name + ":chunks"))
        for ci, ch in enumerate(chunks):
            out.append('<script>window.__c_%s_%d=%s;</script>\n'
                       % (name, ci, json.dumps(ch, ensure_ascii=False)))
        out.append('<script>window.__boot=%s;var %s=JSON.parse([%s].join(""));</script>\n'
                   % (json.dumps(name + ":parsed"), name,
                      ",".join("window.__c_%s_%d" % (name, i) for i in range(len(chunks)))))
    return "".join(out)
# The main script is wrapped so that a runtime failure reports itself: on a file:// page the
# window.onerror listener only sees "Script error." when the detail is suppressed, but a catch
# inside the same script always sees the real message. (Parse errors still only reach the console.)
OUT.write_text(head + _EXAMPLE_BANNER + body
               + _MARK("before-payload")
               + _payload_blocks()
               + _MARK("payload-loaded")
               + _code_blocks()
               + "<script>window.__codeDone='yes';</script>\n"
               + _SUMMARY, encoding="utf-8")
print("wrote", OUT, OUT.stat().st_size // 1024, "KB")
