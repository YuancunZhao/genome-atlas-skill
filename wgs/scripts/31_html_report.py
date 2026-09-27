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

# AN5（§7）：祖源与父母系的**事实**文案由数据生成，覆盖 yaml 里的手写版本。手写文案里常写死某个样本
# 的具体值（"最接近云南白羊村"、"N-M1845 之下进入分辨率极限"），换样本就成了错话；这里用实际的分析
# 结果重新表述，并打印被替代的键，方便对照。其余医学文案仍按原机制由 yaml 提供。
from ancestry_data import ancestry_copy as _ancestry_copy
_ANC_COPY = _ancestry_copy(D)
_ANC_REPLACED = sorted(k for k in _ANC_COPY if k in UI)
UI.update(_ANC_COPY)
if _ANC_REPLACED:
    print("31_html_report: ancestry copy generated from the data replaces: " + ", ".join(_ANC_REPLACED))






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
def _json_clean(obj):
    """递归把 NaN/Infinity 变成 null。

    json.dumps 默认把它们写成裸 NaN / Infinity —— 严格 JSON 解析器会直接拒绝，而浏览器的
    JSON.parse 也要求合法 JSON。这些值本来来自 pandas 的缺失，语义就是 null。
    """
    if isinstance(obj, float):
        return obj if obj == obj and obj not in (float("inf"), float("-inf")) else None
    if isinstance(obj, dict):
        return {k: _json_clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_clean(v) for v in obj]
    try:                                   # numpy 标量：转成 Python 原生类型
        import numpy as _np
        if isinstance(obj, _np.generic):
            return _json_clean(obj.item())
    except Exception:
        pass
    return obj


def _json_text(obj):
    """合法 JSON 文本，且不含会提前闭合 <script> 的序列。"""
    return json.dumps(_json_clean(obj), ensure_ascii=False, allow_nan=False)


def _js_string_literal(text):
    """把一段文本写成 JS 字符串字面量，并在斜杠前加反斜杠以断开闭合标签序列（JSON 与 JS 都
    把转义的斜杠当普通斜杠，所以数据不变，浏览器也不会提前结束脚本块）。"""
    return json.dumps(text, ensure_ascii=False).replace("</", "<\\/")


_DATA_OBJS = [("D", D), ("BLOODV3", BLOOD), ("UI", UI), ("FIND", FIND), ("PGX", PGX), ("PRS_EN", PRS_ZH)]
_PH = ("__DATA__", "__UI__", "__FIND__", "__PGX__", "__BLOOD__", "__PRS_EN__")
_lines = js.split("\n")
_decl = {i for i, l in enumerate(_lines) if any(p in l for p in _PH)}
if not _decl:
    raise SystemExit("31_html_report: no payload placeholder found in report_script.js")
js_code = "\n".join(l for i, l in enumerate(_lines) if i not in _decl)   # logic only; payload is built below

# Safari scopes a top-level `const` to its own <script> element, so splitting the logic broke every
# cross-block reference: the page reported "Can't find variable: reveal" and "renderAll", the figures
# stayed empty and document.title never changed, while Chrome, node and jsc all worked. The logic is
# therefore one block again. The payload stays chunked because it only ever assigns to window.*
# properties, which are shared. The per-statement markers keep the diagnostics.
def _code_blocks():
    out = []
    for lineno, ln in enumerate(js_code.split("\n"), 1):
        out.append(ln)
        # Record the source line of the last completed statement. Safari reports only
        # 'Script error. @0:0' for a failing inline script, so "how far did it get" is the signal.
        if ln[:1] not in (" ", "\t") and ln.rstrip().endswith((";", "}")):
            out.append("window.__boot='L%d';" % lineno)
    return "<script>\n" + "\n".join(out) + "\n</script>\n"

def _payload_blocks():
    out = []
    for name, obj in _DATA_OBJS:
        s = _json_text(obj)
        chunks = [s[i:i + _CHUNK] for i in range(0, len(s), _CHUNK)] or [""]
        out.append(_MARK(name + ":chunks"))
        for ci, ch in enumerate(chunks):
            out.append('<script>window.__c_%s_%d=%s;</script>\n'
                       % (name, ci, _js_string_literal(ch)))
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
