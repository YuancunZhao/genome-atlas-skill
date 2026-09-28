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
# 键名必须是 wgsconfig 真实导出的名字（SUPERPOP/SUBPOPS/AXIS）：曾用 REF_SUPERPOP/REF_SUBPOPS/
# AXIS_POPS 这三个不存在的名字，getattr 默认 None，每份 run_info.json 都在参数区记假值。
# build() 对此有硬守卫（名字不存在直接抛错），这里的测试同样逐键核验。
PARAM_KEYS = [
    "SAMPLE", "BUILD", "SEX",
    "MIN_DP", "MIN_MQ", "MIN_BQ", "MIN_VQ", "MPILEUP_MAX_DP", "MT_MAX_DP",
    "THREADS",
    "REGIONAL_ENABLED", "AADR_ENABLED", "LOCAL_ENABLED",
    "SUPERPOP", "SUBPOPS", "AXIS",
    "LA_LABELS", "AADR_MODERN", "AADR_ANCIENT_PREFIX", "AADR_ANNOTATION",
]
# 键 -> wgsconfig 里指名道姓该工具的属性。配置过的工具必须探**配置的那份**：PATH 上另装一份同名
# 工具时，记下的版本不是本次运行会调用的那个；配置指向的文件不存在时记 None，也不回落到 PATH
# 的同名者——那等于把别的二进制的版本安到本次运行头上。值为 None 的键按名字走 PATH（与步骤
# 脚本 `samtools`、`bcftools` 这类裸名调用方式一致）。
TOOL_KEYS = {
    "plink2": "PLINK2",
    "haplogrep3": "HAPLOGREP3",
    "java": "JAVA",
    "samtools": None, "bcftools": None, "tabix": None, "mosdepth": None,
    "bwa": None, "fastp": None, "node": None, "R": None,
}


def _run(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        out = (r.stdout or r.stderr or "").strip().splitlines()
        return out[0].strip() if out else ""
    except (OSError, subprocess.SubprocessError):
        return ""


def _tool_versions():
    """Probe the tools this run will actually use, not what is merely installed somewhere.

    Configured tools (wgsconfig names an explicit executable) are probed at that path and
    never fall back to a PATH binary of the same name; an absent configured tool is None.
    Name-resolved tools follow PATH, the way the step scripts invoke them."""
    out = {}
    for t, cfg_attr in TOOL_KEYS.items():
        candidate = str(getattr(_c, cfg_attr)) if cfg_attr else t
        exe = shutil.which(candidate)
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


def _input_fingerprints():
    """关键输入的大小与修改时间。整文件哈希对几十 GB 的 CRAM 不现实；size+mtime 足以发现"输入变了"。

    H6 要求样本/输入/参数/工具/参考变化后缓存失效。判定"变了"只需要一个稳定的指纹，不需要
    密码学强度——所以这里不读文件内容，读取成本与文件大小无关。
    """
    out = {}
    for name, val in (("reads", READS), ("fastq1", FASTQ1), ("fastq2", FASTQ2),
                      ("vendor_vcf", VENDOR_VCF), ("y_reads", Y_READS)):
        if not val:
            out[name] = None            # 空可选项显式记为 null，而不是缺键
            continue
        p = pathlib.Path(str(val))
        if p.exists():
            st = p.stat()
            out[name] = {"path": str(p), "size": st.st_size, "mtime": int(st.st_mtime)}
        else:
            out[name] = {"path": str(p), "missing": True}
    return out


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
    values = {}
    for k in PARAM_KEYS:
        # 硬守卫：键名打错时 getattr(_c, k, None) 会把参数记成静默的 None，每份记录都在撒谎。
        # 宁可当场抛错，也不要一份看起来完整的假配置。
        if not hasattr(_c, k):
            raise AttributeError(
                f"_run_info.PARAM_KEYS names {k!r}, which wgsconfig does not export; "
                "fix the key (wgsconfig exports e.g. SUPERPOP/SUBPOPS/AXIS), not the value")
        values[k] = getattr(_c, k)
    return {
        "schema_version": 1,
        "run_id": run_id or datetime.datetime.now().strftime("%Y%m%d-%H%M%S"),
        "started_at": (started or datetime.datetime.now().astimezone()).isoformat(timespec="seconds"),
        "host": {"hostname": socket.gethostname(), "cpu": os.cpu_count(), "platform": platform.platform()},
        "code": _code_revision(_c.ROOT),
        "reference": _reference_info(),
        "inputs": _input_fingerprints(),
        "parameters": values,
        "tools": _tool_versions(),
    }


def default_path():
    return pathlib.Path(_c.W) / "run_info.json"


def append_step(name, status, rc=None, seconds=None, products=None, note="", path=None):
    """记下一步的结束状态与产物，追加进 run_info.json 的 steps。

    H6：保留退出状态与产物关联、失败不复用半成品。做法是读-改-写并在失败时拒绝覆盖上一次的
    成功记录：一个被中断的步骤必须留下痕迹，而不是让下一轮以为它跑过。
    """
    out = pathlib.Path(path) if path else default_path()
    info = build()                      # 新文件也要带运行事实，不能只剩一个孤立的 steps
    if out.exists():
        try:
            info = json.loads(out.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            # 读不动就重建：宁可丢历史，也不要把坏文件当依据
            info = build()
        info.setdefault("steps", [])
    info.setdefault("steps", [])
    entry = {"name": str(name), "status": status}
    if rc is not None:
        entry["rc"] = int(rc)
    if seconds is not None:
        entry["seconds"] = round(float(seconds), 1)
    if products:
        entry["products"] = [str(x) for x in products]
    if note:
        entry["note"] = str(note)
    entry["at"] = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    info["steps"].append(entry)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(info, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    tmp.replace(out)          # 原子替换：读者不会看到半个 steps
    return out, entry


def write(path=None, **kw):
    info = build(**kw)
    out = pathlib.Path(path) if path else pathlib.Path(_c.W) / "run_info.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(info, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    return out, info


def _main(argv):
    import argparse
    ap = argparse.ArgumentParser(description="record run facts, or one step's outcome")
    ap.add_argument("--step", help="step name to record")
    ap.add_argument("--status", default="ok", choices=["ok", "failed", "skipped"])
    ap.add_argument("--rc", type=int, default=None)
    ap.add_argument("--seconds", type=float, default=None)
    ap.add_argument("--products", nargs="*", default=None)
    ap.add_argument("--note", default="")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    if a.step:
        out, entry = append_step(a.step, a.status, rc=a.rc, seconds=a.seconds,
                                 products=a.products, note=a.note, path=a.out)
        print(f"{a.status}: {entry['name']}" + (f" (rc={entry['rc']})" if "rc" in entry else ""))
        return 0
    return None


if __name__ == "__main__":
    _r = _main(sys.argv[1:])
    if _r is not None:
        sys.exit(_r)
    out, info = write()
    print(f"wrote {out}")
    print(f"  run {info['run_id']} · {info['code']['branch']}@{str(info['code']['commit'])[:8]}"
          f"{' (dirty)' if info['code']['dirty'] else ''} · {info['reference']['build']} · {info['host']['cpu']} cpu")
    missing = [k for k, v in info["tools"].items() if v is None]
    if missing:
        print(f"  tools not on PATH: {' '.join(missing)}")
