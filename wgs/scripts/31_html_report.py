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
  if (window.__earlyErrors.length) { b.textContent = 'SCRIPT ERROR: ' + window.__earlyErrors.join(' | '); b.style.color = '#ffd0c8'; }
  else if (b.textContent.indexOf('waiting') >= 0) { b.textContent = 'the report script never ran -- open the developer console (F12) for the reason'; b.style.color = '#ffd0c8'; }
});
</script>
"""
_BOOT = ('<div id="boot" style="position:fixed;top:0;left:0;right:0;z-index:99999;background:#111;color:#eee;'
         'font:12px/1.7 ui-monospace,Menlo,monospace;padding:5px 10px;text-align:center">build ' + _BUILD +
         ' &middot; <span id="bootstate">HTML parsed; waiting for the report script&hellip;</span></div>'
         '<div style="height:30px"></div>')
head = _EARLY + head
body = body + _BOOT
js = open(S/"report_script.js", encoding="utf-8").read()
js = (js.replace("__DATA__", json.dumps(D, ensure_ascii=False, separators=(",", ":"))).replace("__UI__", json.dumps(UI, ensure_ascii=False))
        .replace("__FIND__", json.dumps(FIND, ensure_ascii=False)).replace("__PGX__", json.dumps(PGX, ensure_ascii=False)).replace("__PRS_EN__", json.dumps(PRS_ZH, ensure_ascii=False)).replace("__BLOOD__", json.dumps(BLOOD, ensure_ascii=False)))
# The main script is wrapped so that a runtime failure reports itself: on a file:// page the
# window.onerror listener only sees "Script error." when the detail is suppressed, but a catch
# inside the same script always sees the real message. (Parse errors still only reach the console.)
OUT.write_text(head + body + "<script>\ntry{\n" + js
               + "\n}catch(e){var _bf=document.getElementById('bootstate');"
                 "if(_bf){_bf.textContent='MAIN SCRIPT FAILED: '+(e&&e.message?e.message:e);_bf.style.color='#ffd0c8';}"
                 "console.error('[main script]',e);throw e;}\n</script>\n", encoding="utf-8")
print("wrote", OUT, OUT.stat().st_size // 1024, "KB")
