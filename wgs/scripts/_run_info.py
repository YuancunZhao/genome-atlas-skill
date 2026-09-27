#!/usr/bin/env python
"""Record the facts of one run: run id, code revision, effective parameters, tool and reference versions.

H6 asks the pipeline to keep the run id, the code SHA (and whether the tree was dirty), the effective
parameters, the tool and reference versions, exit status and product association -- and for the report to
expose only what it needs. Everything except the report exposure lands here.

Why it matters: a report rendered months later cannot be interpreted without knowing which tree version
called the haplogroup, which reference build positions refer to, or whether the code that produced it had
uncommitted changes. "Installed version" is not "the version this run used", so versions are probed now
and stored next to the outputs rather than looked up later.

The file is written to work/wgs/run_info.json. Nothing here is printed into the report by default.
"""
import datetime, hashlib, json, os, pathlib, platform, shutil, socket, subprocess, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, BUILD, MIN_* ...
import wgsconfig as _c

# 记录哪些参数：影响结果的取样口径与阈值。跑完再想不起来当时设了什么，是最常见的复现障碍。
PARAM_KEYS = [
    "SAMPLE", "BUILD", "SEX",
    "MIN_DP", "MIN_MQ", "MIN_BQ", "MIN_VQ", "MPILEUP_MAX_DP", "MT_MAX_DP",
    "THREADS",
    "REGIONAL_ENABLED", "AADR_ENABLED", "LOCAL_ENABLED",
    "REF_SUPERPOP", "REF_SUBPOPS", "AXIS_POPS",
    "LA_LABELS", "AADR_MODERN", "AADR_ANCIENT_PREFIX", "AADR_ANNOTATION",
]
TOOL_KEYS = ["samtools", "bcftools", "tabix", "plink2", "mosdepth", "bwa", "fastp", "java", "node", "R"]


def _run(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        out = (r.stdout or r.stderr or "").strip().splitlines()
        return out[0].strip() if out else ""
    except (OSError, subprocess.SubprocessError):
        return ""


def _tool_versions():
    """Probe the tools this run will actually use (PATH-resolved), not what is merely installed somewhere."""
    out = {}
    for t in TOOL_KEYS:
        exe = shutil.which(t)
        if not exe:
            out[t] = None
            continue
        ver = ""
        for flag in (["--version"], ["-version"], ["version"]):
            ver = _run([exe, *flag])
            if ver:
                break
        out[t] = f"{ver} ({exe})" if ver else exe
    return out


def _code_revision(root):
    """Branch, commit and whether the working tree had uncommitted changes at the time of the run."""
    sha = _run(["git", "-C", str(root), "rev-parse", "HEAD"])
    branch = _run(["git", "-C", str(root), "rev-parse", "--abbrev-ref", "HEAD"])
    dirty = _run(["git", "-C", str(root), "status", "--porcelain"])
    return {"commit": sha or None, "branch": branch or None,
            "dirty": bool(dirty), "dirty_entries": len([l for l in dirty.splitlines() if l.strip()]) if dirty else 0}


def _reference_info():
    """Build plus a fingerprint of the reference index, so a later run can tell if the reference changed."""
    fai = pathlib.Path(str(FASTA) + ".fai")
    info = {"build": BUILD, "fasta": str(FASTA), "fai": str(fai) if fai.exists() else None}
    if fai.exists():
        h = hashlib.sha256()
        n = 0
        for line in open(fai, "rb"):
            h.update(line)
            n += 1
        info["fai_sha256"] = h.hexdigest()[:16]
        info["fai_sequences"] = n
        # 参考序列本身不重算（几十 GB）；索引指纹足以发现"参考换了"这件事。
    else:
        info["fai_sha256"] = None
        info["fai_sequences"] = None
    return info


def build(run_id=None, started=None):
    values = {k: getattr(_c, k, None) for k in PARAM_KEYS}
    return {
        "schema_version": 1,
        "run_id": run_id or datetime.datetime.now().strftime("%Y%m%d-%H%M%S"),
        "started_at": (started or datetime.datetime.now().astimezone()).isoformat(timespec="seconds"),
        "host": {"hostname": socket.gethostname(), "cpu": os.cpu_count(), "platform": platform.platform()},
        "code": _code_revision(_c.ROOT),
        "reference": _reference_info(),
        "parameters": values,
        "tools": _tool_versions(),
    }


def write(path=None, **kw):
    info = build(**kw)
    out = pathlib.Path(path) if path else pathlib.Path(_c.W) / "run_info.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(info, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    return out, info


if __name__ == "__main__":
    out, info = write()
    print(f"wrote {out}")
    print(f"  run {info['run_id']} · {info['code']['branch']}@{str(info['code']['commit'])[:8]}"
          f"{' (dirty)' if info['code']['dirty'] else ''} · {info['reference']['build']} · {info['host']['cpu']} cpu")
    missing = [k for k, v in info["tools"].items() if v is None]
    if missing:
        print(f"  tools not on PATH: {' '.join(missing)}")
