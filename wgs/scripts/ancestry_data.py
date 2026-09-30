"""Ancestry configuration and result bookkeeping (AN0).

Two jobs, both deliberately free of side effects at import time:

1. `read_options(cfg)` turns a raw configuration dict into validated switches and parameters.
   It must be a pure function -- no file reads, no directory creation, no dependency on
   `wgsconfig` (which creates the work tree when imported). "Explicitly configured" is decided
   from the keys the user actually wrote, never from a shared default: defaulting
   `ref_superpop` to EAS and then treating that default as a user choice would silently run
   regional East-Asian panels for every sample, which is exactly the failure this module exists
   to prevent.

2. Manifest read/write/match. Every analysis directory records what produced it; step 30 only
   consumes results whose manifest matches the current sample, configuration and references.
   A missing or mismatching manifest means "rebuild", never "reuse".

Nothing here knows a sample name, a haplogroup or a population count.
"""
import json, math, os, pathlib, random, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from gt_alleles import gt_alleles  # noqa: E402  复用它，不另写一个 GT 解析器

SCHEMA_VERSION = 1
STATES = ("ok", "disabled", "unavailable", "failed")
REQUIRED_MANIFEST_KEYS = ("sample_id",)
# 只支持 GRCh37。别名（hg19 / b37 / GRCh37.p13）归一化到规范名，避免"同一个 build 两种写法"
# 在配置校验与下游取参考时给出不同结论。新增 build 必须同时具备 FASTA/chain/1000G/AADR/注释，
# 否则位点会按错误的坐标系解释而没有任何一步失败。
SUPPORTED_BUILDS = ("GRCh37",)
BUILD_ALIASES = {"GRCh37": "GRCh37", "grch37": "GRCh37", "hg19": "GRCh37", "b37": "GRCh37",
                 "GRCh37.p13": "GRCh37"}


def normalize_build(value, default="GRCh37"):
    """把配置里的 build 写法归一化；不支持的取值抛出 ValueError（列出可接受的写法）。"""
    raw = (value if isinstance(value, str) else "") or default
    key = raw.strip()
    if key not in BUILD_ALIASES:
        raise ValueError(
            f"unsupported reference build {raw!r}: this pipeline's FASTA, chain files, 1000G panel, "
            f"AADR set and annotation GFF are GRCh37, so another build would be read in the wrong "
            f"coordinate system without any step failing. Accepted: GRCh37 (aliases "
            f"{sorted(set(BUILD_ALIASES) - {'GRCh37'})}).")
    return BUILD_ALIASES[key]

_LIST_KEYS = (
    "ref_subpops", "axis_pops", "local_ancestry_a", "local_ancestry_b", "local_ancestry_control",
    "local_ancestry_labels", "aadr_modern", "aadr_ancient_prefix",
    "local_ancestry_calibration_pops", "local_ancestry_calibration_chroms",
)


def _as_list(cfg, key, default):
    """配置里的列表键。字符串会被 list() 拆成字符，必须在入口挡住。"""
    if key not in cfg:
        return list(default)
    v = cfg[key]
    if v is None:
        return list(default)
    if isinstance(v, str):
        raise ValueError(
            f"{key} must be a YAML list (got the string {v!r}: it would be read as the characters "
            f"{list(v)!r}). Write it as [{v}] or as a - item block."
        )
    if not isinstance(v, (list, tuple)):
        raise ValueError(f"{key} must be a list, got {type(v).__name__}")
    return [str(x).strip() for x in v]


def _rate(cfg, key, default):
    """比例型门槛：取值 (0, 1]。"""
    v = cfg.get(key, default)
    if isinstance(v, str):
        try:
            v = float(v)
        except ValueError:
            raise ValueError(f"{key} must be a number in (0, 1], got {v!r}") from None
    if not isinstance(v, (int, float)) or isinstance(v, bool):
        raise ValueError(f"{key} must be a number in (0, 1], got {type(v).__name__}")
    if not (0 < float(v) <= 1):
        raise ValueError(f"{key} must be in (0, 1], got {v!r}")
    return float(v)


def _positive_int(cfg, key, default):
    v = cfg.get(key, default)
    if isinstance(v, str):
        try:
            v = int(v)
        except ValueError:
            raise ValueError(f"{key} must be a positive integer, got {v!r}") from None
    if isinstance(v, bool) or not isinstance(v, int) or v <= 0:
        raise ValueError(f"{key} must be a positive integer, got {v!r}")
    return int(v)


def _int(cfg, key, default):
    v = cfg.get(key, default)
    if isinstance(v, str):
        try:
            v = int(v)
        except ValueError:
            raise ValueError(f"{key} must be an integer, got {v!r}") from None
    if isinstance(v, bool) or not isinstance(v, int):
        raise ValueError(f"{key} must be an integer, got {v!r}")
    return int(v)


def _str(cfg, key, default=""):
    v = cfg.get(key, default)
    return "" if v is None else str(v)


def chromosome_lengths(candidates, fallback=None, log=print, build=None):
    """从 FASTA 索引（`.fai`）读染色体长度；读不到才用内置表，并**明确告警**。

    内置长度表只对 hg19/GRCh37 正确。换参考版本（例如 GRCh38）后，用 hg19 的长度去算百分比会
    静默给出偏大的数字——图上看着合理，含义全错。所以优先读索引，回退时开口说话。
    `.fai` 的键可能是 `1` 或 `chr1`，统一去掉 `chr` 前缀；副contig（含 `_`、随机/未定位序列）不收。

    build（规范名，如 "GRCh37"）给出后按**同 build**约束候选：元素必须是 (path, build_tag) 对，
    tag 归一化后不等于所需 build 的索引整条跳过并记 note——GRCh37 的任务读到 hg38 的长度是
    跨 build 错误，不是回退；宁可落到"没有同 build 索引"的告警。此时传裸路径（忘了声明 tag）
    直接抛 TypeError，逼调用方显式。不带 build 时维持旧语义（裸路径、先读到先得），供非生产
    调用与既有测试使用。回退表同样只在 GRCh37 任务上发出：其他 build 没有同 build 索引时
    返回空 dict（可见的缺），而不是发一张已知错误的表。
    """
    wanted = BUILD_ALIASES.get(build, build) if build else None
    for cand in candidates:
        if wanted is not None:
            if not (isinstance(cand, tuple) and len(cand) == 2):
                raise TypeError(
                    f"with build={build!r}, candidates must be (path, build_tag) pairs so an index "
                    f"of another build cannot slip in; got {cand!r}")
            path_s, tag = cand
            if BUILD_ALIASES.get(tag, tag) != wanted:
                if path_s:
                    log(f"skipped {path_s}: {tag} index, this task needs {wanted} (cross-build lengths are never a fallback)")
                continue
        else:
            path_s = cand
        if not path_s:
            continue          # 配置未给出路径时这里是空串；Path('') 会变成 '.' 并触发假的"读不了"告警
        try:
            import pathlib as _pl
            path = _pl.Path(path_s)
            if not path.exists():
                continue
            out = {}
            with open(path, encoding="utf-8") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 2 or not f[1].isdigit():
                        continue
                    name = f[0][3:] if f[0].startswith("chr") else f[0]
                    if any(ch in name for ch in "_.") or not name[:2].rstrip("XMY").isdigit() and name not in ("X", "Y", "M", "MT"):
                        continue
                    out[name] = int(f[1])
            if out:
                log(f"chromosome lengths read from {path} ({len(out)} sequences, build {wanted or 'unspecified'})")
                return out
        except OSError as e:
            log(f"warning: could not read {path_s} ({e})")
    if fallback:
        if wanted is None or wanted == "GRCh37":
            log("warning: no FASTA index found; falling back to the built-in hg19/GRCh37 length table, "
                "which is WRONG for any other reference build")
            return dict(fallback)
        log(f"warning: no FASTA index of build {wanted} found; the built-in table is hg19/GRCh37 and "
            f"would be wrong here, so no lengths are returned")
    return {}


def target_key(sample_id):
    """目标的内部唯一键。显示名（可重名、含空格与中文）不参与匹配，只有 sample_id 参与。"""
    return f"target:{sample_id}"


def read_options(cfg):
    """把配置翻译成开关与已校验参数；非法输入抛 ValueError。

    开关只看用户**写没写**、且写的是非空值——不回退到任何共享默认：
      regional_enabled  显式给出非空 ref_superpop
      aadr_enabled      显式给出非空 aadr_modern 或 aadr_ancient_prefix
      local_enabled     local_ancestry_a 与 local_ancestry_b 均非空、互不重叠、标签唯一
    """
    if not isinstance(cfg, dict):
        raise ValueError(f"configuration must be a mapping, got {type(cfg).__name__}")

    build = normalize_build(_str(cfg, "build", "GRCh37"))

    regional = _str(cfg, "ref_superpop", "")
    subpops = _as_list(cfg, "ref_subpops", [])
    axis = _as_list(cfg, "axis_pops", [])
    if axis and len(axis) != 2:
        raise ValueError(f"axis_pops needs exactly two populations, got {axis!r}")

    la_a = _as_list(cfg, "local_ancestry_a", [])
    la_b = _as_list(cfg, "local_ancestry_b", [])
    la_control = _as_list(cfg, "local_ancestry_control", [])
    if la_a and la_b:
        overlap = sorted(set(la_a) & set(la_b))
        if overlap:
            raise ValueError(f"local_ancestry_a and local_ancestry_b overlap on {overlap}")
        if len(set(la_a)) != len(la_a) or len(set(la_b)) != len(la_b):
            raise ValueError("local_ancestry_a/b contain duplicate populations")
        if set(la_a) & set(la_control) or set(la_b) & set(la_control):
            raise ValueError("the control panel overlaps a source panel; A/B must stay disjoint from it")
    labels_given = "local_ancestry_labels" in cfg
    labels = _as_list(cfg, "local_ancestry_labels", [])
    if labels_given:
        if len(labels) != 2:
            raise ValueError(f"local_ancestry_labels needs exactly two labels, got {labels!r}")
        if len(set(labels)) != len(labels) or any(not x for x in labels):
            raise ValueError(f"local_ancestry_labels must be two distinct non-empty labels, got {labels!r}")
    if la_a and la_b and not labels_given:
        labels = ["+".join(la_a), "+".join(la_b)]  # 由来源生成，不硬编码任何地理或成分语义

    modern = _as_list(cfg, "aadr_modern", [])
    ancient_prefix = _as_list(cfg, "aadr_ancient_prefix", [])
    if any(not p for p in ancient_prefix):
        raise ValueError("aadr_ancient_prefix must not contain an empty string (it would match everything)")

    pops = _as_list(cfg, "local_ancestry_calibration_pops", [])
    if not pops:
        pops = list(axis) if axis else ([la_a[0]] + [la_b[0]] if (la_a and la_b) else [])
    chroms = _as_list(cfg, "local_ancestry_calibration_chroms", ["1", "2", "6", "22"])
    if any(not c for c in chroms):
        raise ValueError("local_ancestry_calibration_chroms must not contain empty entries")

    return {
        "build": build,
        "sample_id": _str(cfg, "sample_id", "SAMPLE"),
        "name_zh": _str(cfg, "name_zh", "") or _str(cfg, "sample_id", "SAMPLE"),
        "name_en": _str(cfg, "name_en", "") or _str(cfg, "sample_id", "SAMPLE"),
        "threads": _str(cfg, "threads", "8"),
        "mem_gb": _str(cfg, "mem_gb", "30"),
        # ---- switches (explicit configuration only)
        "regional_enabled": bool(regional),
        "aadr_enabled": bool(modern or ancient_prefix),
        "local_enabled": bool(la_a and la_b),
        # ---- reference selection
        "ref_superpop": regional,
        "ref_subpops": subpops,
        "axis_pops": axis,
        "la_a": la_a,
        "la_b": la_b,
        "la_control": la_control,
        "la_labels": labels,
        "aadr_prefix": _str(cfg, "aadr_prefix", ""),
        "aadr_annotation": _str(cfg, "aadr_annotation", ""),
        "aadr_modern": modern,
        "aadr_ancient_prefix": ancient_prefix,
        # ---- thresholds
        "min_call_rate_modern": _rate(cfg, "ancestry_min_call_rate_modern", 0.95),
        "min_call_rate_target": _rate(cfg, "ancestry_min_call_rate_target", 0.95),
        "min_call_rate_ancient": _rate(cfg, "aadr_min_call_rate_ancient", 0.50),
        "min_projection_snps": _positive_int(cfg, "ancestry_min_projection_snps", 10000),
        "min_group_n": _positive_int(cfg, "ancestry_min_group_n", 2),
        # ---- calibration
        "calibration_pops": pops,
        "calibration_chroms": chroms,
        "calibration_n": _positive_int(cfg, "local_ancestry_calibration_n", 20),
        "calibration_seed": _int(cfg, "local_ancestry_calibration_seed", 1),
        # ---- lineage history (optional evidence file)
        "lineage_history_file": _str(cfg, "lineage_history_file", ""),
    }


# ---------------------------------------------------------------- manifest

def build_manifest(sample_id, analysis_id, state="ok", reason_code="", **fields):
    """一份 manifest 的最小骨架；调用方补齐 reference_release/parameters/... 等字段。"""
    if state not in STATES:
        raise ValueError(f"state must be one of {STATES}, got {state!r}")
    m = {"schema_version": SCHEMA_VERSION, "sample_id": str(sample_id), "analysis_id": str(analysis_id),
         "state": state, "reason_code": str(reason_code), "build": "GRCh37", "reference_release": "",
         "parameters": {}, "tool_versions": {}, "input_fingerprints": {}, "outputs": []}
    m.update(fields)
    return m


def disabled_manifest(sample_id, analysis_id, reason_code):
    """显式禁用也要留下记录：30 才能区分"没跑"与"跑了但被配置禁用"。"""
    return build_manifest(sample_id, analysis_id, state="disabled", reason_code=reason_code)


def write_manifest(path, manifest):
    """先写临时文件再原子替换：读者永远看不到半个 manifest。"""
    state = manifest.get("state")
    if state not in STATES:
        raise ValueError(f"state must be one of {STATES}, got {state!r} (a run that is not finished is not 'ok')")
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=1, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def clear_step_manifests(dir_path, *step_ids):
    """一个步骤本次**成功**运行后，撤下它自己的步骤级 manifest（manifest.<step-id>.json）。

    禁用→重新启用→成功之后，目录 manifest.json 已是本次结果，但旧的步骤级 disabled 记录还
    留在目录里，analysis_state 仍会返回 disabled（复审 AN0/AN5/H6 的生命周期缺口：成功必须
    结束对应步骤的非 ok 状态）。只删**自己**的记录——同目录其他步骤的失败/禁用证据原样保留，
    不粗暴清场。返回被撤下的文件名列表；记录本就不存在时是无操作。"""
    d = pathlib.Path(dir_path)
    removed = []
    for sid in step_ids:
        f = d / f"manifest.{sid}.json"
        if f.exists():
            f.unlink()
            removed.append(f.name)
    return removed
    return path


def read_manifest(path):
    path = pathlib.Path(path)
    if not path.exists():
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError):
        return None  # 损坏的 manifest 等同于没有：调用方必须重建


def manifest_matches(actual, expected):
    """这份产物能不能代表当前样本/配置/参考？

    规则：actual 必须存在、含 sample_id、且（若写了 state）state == "ok"；expected 里的每个键
    必须在 actual 中**存在且相等**。复审 AN0/AN5/H6：此前非 REQUIRED 键缺失会 `continue` 跳过——
    一份没有记录参考版本的 manifest 也能通过按 reference_release 的比对，等于没比对。任何不等、
    任何键缺失都返回 False——调用方据此重建，绝不复用旧结果。
    """
    if not isinstance(actual, dict) or not actual or not isinstance(expected, dict) or not expected:
        return False
    if actual.get("state") not in (None, "ok"):
        return False
    for k in REQUIRED_MANIFEST_KEYS:
        if k not in actual:
            return False
    for k, v in expected.items():
        if k not in actual:
            return False
        if actual[k] != v:
            return False
    return True


def _cli(argv=None):
    """最小入口，只给 run_all 用：把一个"按配置跳过"的分析写成 disabled manifest。

        python3 scripts/ancestry_data.py --disabled <analysis_id> --out <path> --sample <id> --reason <code>

    不是调度器，不扫目录、不推断状态：run_all 已经知道自己在跳过哪一步。
    """
    import argparse
    ap = argparse.ArgumentParser(description="ancestry helpers (AN0/AN3)")
    ap.add_argument("--pick-holdout", metavar="PSAM", help="select calibration holdouts from a .psam")
    ap.add_argument("--pops", default="", help="comma-separated populations to hold out")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--check-samples", metavar="VCF", help="verify a subset VCF carries exactly these samples")
    ap.add_argument("--expect", metavar="FILE", help="file with the expected sample names, one per line")
    ap.add_argument("--disabled", metavar="ANALYSIS_ID")
    ap.add_argument("--step-ok", metavar="ANALYSIS_ID",
                    help="write a state=ok step record (a successful optional step's receipt for its consumer)")
    ap.add_argument("--clear-step", metavar="ANALYSIS_ID",
                    help="remove this step's own manifest.<id>.json after a successful run (lifecycle)")
    ap.add_argument("--dir", help="directory holding the step-level manifest to clear")
    ap.add_argument("--out", help="output path (holdout list, or the manifest to write)")
    ap.add_argument("--sample")
    ap.add_argument("--reason", default="disabled_by_config")
    a = ap.parse_args(argv)
    if a.pick_holdout:
        rows = read_psam(a.pick_holdout)
        out = pick_holdout(rows, [p for p in a.pops.split(",") if p], a.n, a.seed)
        print(f"state={out['state']} reason={out['reason_code']} per_pop={out['per_pop']} "
              f"ids={len(out['ids'])} {out['detail']}")
        if a.out and out["state"] == "ok":
            pathlib.Path(a.out).write_text("\n".join(out["ids"]) + "\n", encoding="utf-8")
        elif a.out:
            pathlib.Path(a.out).write_text("", encoding="utf-8")
        return 0 if out["state"] == "ok" else 2
    if a.check_samples:
        got = set(read_vcf_samples(a.check_samples))
        want = {l.strip() for l in pathlib.Path(a.expect).read_text(encoding="utf-8").splitlines() if l.strip()}
        if got == want:
            print(f"sample set matches ({len(got)} names)")
            return 0
        print(f"sample set mismatch: missing={sorted(want - got)[:5]} unexpected={sorted(got - want)[:5]}",
              file=sys.stderr)
        return 3
    if a.clear_step:
        if not a.dir:
            ap.error("--clear-step requires --dir")
        removed = clear_step_manifests(a.dir, a.clear_step)
        print(f"cleared {', '.join(removed)}" if removed
              else f"no step manifest for {a.clear_step} (nothing to clear)")
        return 0
    if a.step_ok:
        if not a.out or not a.sample:
            ap.error("--step-ok requires --out and --sample")
        write_manifest(a.out, build_manifest(a.sample, a.step_ok, state="ok"))
        print(f"wrote {a.out} (state=ok)")
        return 0
    if not a.disabled or not a.out or not a.sample:
        ap.error("--disabled requires --out and --sample")
    m = disabled_manifest(a.sample, a.disabled, a.reason)
    write_manifest(a.out, m)
    print(f"wrote {a.out} (state=disabled, reason={a.reason})")
    return 0




# ─────────────────────────────────────────── AN1: metadata, dosage and locations

# The .anno header carries a paragraph of documentation in every column name, so columns are matched
# by prefix. This table is the only place where those raw names appear; nothing downstream may depend
# on them. Required columns fail loudly; optional ones become None rather than "" or 0.
_ANNO_ALIASES = (
    ("record_id", ("Genetic ID",)),
    ("master_id", ("Persistent Genetic ID",)),
    ("individual_id", ("Individual ID",)),
    ("source_population_id", ("Group ID",)),
    ("locality", ("Locality",)),
    ("political_entity", ("Political Entity",)),
    ("latitude", ("Latitude",)),
    ("longitude", ("Longitude",)),
    ("date_mean_bp", ("Date mean in BP",)),
    ("date_sd_bp", ("Date standard deviation in BP",)),
    ("date_raw", ("Full Date",)),
    ("date_basis", ("Method for Determining Date",)),
    ("genotype_representation", ("Suffices",)),
    ("data_type", ("Data type",)),
    ("y_hg_raw", ("Y haplogroup in terminal",)),
    ("y_hg_isogg_raw", ("Y haplogroup  in ISOGG",)),
    ("mt_hg_raw", ("mtDNA haplogroup",)),
    ("publication", ("Publication abbreviation",)),
    ("doi", ("doi for publication",)),
    ("assessment", ("ASSESSMENT WARNINGS",)),
)
_REQUIRED_ANNO = ("record_id",)
# 一个区间出现的位置（考古上下文区间、或 95.4% CI 的 calBCE 区间）。BP 基准固定 1950。
_DATE_RANGE = re.compile(r"(\d{3,5})\s*[-–]\s*(\d{3,5})\s*(?:cal)?\s*(BCE|BC|CE|BP)", re.I)


def _cell(value):
    """空单元格在 AADR 里写作 '..'；统一成 None。"""
    if value is None:
        return None
    v = str(value).strip()
    return None if v in ("", "..", "nan", "NA", "N/A") else v


def _coord(value, lo, hi):
    if value is None:
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x) or not (lo <= x <= hi):
        return None
    return x


def _number(value):
    if value is None:
        return None
    try:
        x = float(str(value).replace(",", ""))
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def _dates_from_raw(raw):
    """从 Full Date 里取可信区间，返回 (min_bp, max_bp)。BP 基准 1950。不可解析时 (None, None)。"""
    if not raw:
        return None, None
    m = _DATE_RANGE.search(raw)
    if not m:
        return None, None
    a, b = int(m.group(1)), int(m.group(2))
    axis = m.group(3).upper()
    if axis in ("BCE", "BC"):
        lo_bp, hi_bp = 1950 + min(a, b), 1950 + max(a, b)
    elif axis == "CE":
        # 复审 P1（AN1/AN6）：CE 年份此前落进"已经是 BP"分支，1000-1200 CE 被当成
        # 1000-1200 BP；正确换算是 1950-年 → 750-950 BP（CE 越晚 BP 越小）。
        lo_bp, hi_bp = 1950 - max(a, b), 1950 - min(a, b)
    else:  # 已经是 BP
        lo_bp, hi_bp = min(a, b), max(a, b)
    return lo_bp, hi_bp


def normalize_metadata(rows, dataset, release):
    """把 .anno 的原始行规范化成 7.3 的统一记录。

    只在这里适配原始列名。年代未知保持 None（不填 0 伪装现代）；坐标缺失或越界保持 None（不变成
    (0,0)，也不保留假位置）；来源没有给出的字段一律 None。
    """
    if not rows:
        return []
    cols = list(rows[0].keys())
    picked = {}
    for field, prefixes in _ANNO_ALIASES:
        for c in cols:
            if any(c.startswith(p) for p in prefixes):
                picked[field] = c
                break
    missing = [f for f in _REQUIRED_ANNO if f not in picked]
    if missing:
        raise ValueError(
            f"the .anno file has no column for {missing}; expected one starting with "
            f"{[p for f in missing for p in dict(_ANNO_ALIASES)[f]]}. Refusing to guess."
        )

    out = []
    for row in rows:
        get = lambda f: _cell(row.get(picked[f])) if f in picked else None  # noqa: E731
        lat = _coord(get("latitude"), -90.0, 90.0)
        lon = _coord(get("longitude"), -180.0, 180.0)
        raw = get("date_raw")
        mean, sd = _number(get("date_mean_bp")), _number(get("date_sd_bp"))
        basis = get("date_basis")
        lo, hi = _dates_from_raw(raw)
        if mean is None and lo is not None and hi is not None:
            mean = int(round((lo + hi) / 2))
            basis = "contextual_range_midpoint"   # 区间中点，不是把 SD 冒充 95% 区间
        out.append({
            "record_id": get("record_id"),
            "individual_id": get("individual_id") or get("master_id") or get("record_id"),
            "master_id": get("master_id") or get("record_id"),
            "source_version_id": get("record_id"),
            "genotype_representation": (get("genotype_representation") or "").upper() or None,
            "dataset": dataset, "reference_release": release,
            "kind": "unknown",                       # 由 assign_kind 依据面板选择决定，绝不按年代切
            "source_population_id": get("source_population_id"),
            "label": get("source_population_id"),
            "locality": get("locality"),
            "location_id": None,
            "latitude": lat, "longitude": lon,
            # 复审 P1（AN1/AN6）：AADR 只发布坐标本身，不发布精度等级——"有坐标"冒充 site 让
            # 模板的"未知精度"防护失效。site/region 只能由人工地点覆盖表（30 的细化层）给出；
            # 这里如实记 unknown，坐标来源单独用 location_source 表达。
            "location_precision": "unknown",
            "location_source": "anno" if (lat is not None and lon is not None) else None,
            "date_mean_bp": mean,
            "date_sd_bp": sd,
            "date_min_bp": lo, "date_max_bp": hi,
            "date_basis": basis, "date_raw": raw,
            "publication": get("publication"),
            "y_hg_raw": get("y_hg_raw"), "mt_hg_raw": get("mt_hg_raw"),
            "hg_source_tree": "YFull12.03" if get("y_hg_raw") else None,
            "hg_call_source": "aadr_automatic" if get("y_hg_raw") else None,
            "hg_qc": get("assessment"),
            "pc1": None, "pc2": None, "pc3": None, "pc4": None,
            "n_reference_snps": None, "n_called_snps": None, "call_rate": None,
            "eligible": None, "exclusion_reason": None,
            "distance_to_target": None, "group_id": None,
        })
    return out


def location_key_of(record):
    """(dataset, population, locality) —— 同一群体的不同遗址必须是不同地点，绝不共用一个坐标。"""
    return "{}:{}:{}".format(record.get("dataset") or "", record.get("source_population_id") or "",
                             (record.get("locality") or "").strip().lower() or "unknown")


def _representation_rank(rep):
    """来源推荐的表示优先：SG/DG（高覆盖二倍体）> HO > AG/TW/其它。"""
    rep = (rep or "").upper()
    if rep in ("SG", "DG"):
        return 3
    if rep == "HO":
        return 2
    if rep in ("AG", "TW", "WGC", "AA", "EC"):
        return 1
    return 0


def dedupe_by_master_id(records):
    """同一个人的多种技术表示只算一个人；丢弃者带原因保留，供详情使用。"""
    groups = {}
    for r in records:
        groups.setdefault(r.get("master_id") or r.get("record_id"), []).append(r)
    kept, dropped = [], []
    for master, rows in groups.items():
        if len(rows) == 1:
            kept.append(rows[0])
            continue
        ranked = sorted(rows, key=lambda r: (-_representation_rank(r.get("genotype_representation")),
                                             str(r.get("record_id"))))
        kept.append(ranked[0])
        for r in ranked[1:]:
            dropped.append({"record_id": r.get("record_id"), "master_id": master,
                            "reason_code": "duplicate_representation",
                            "detail": f"{r.get('genotype_representation')} not used; "
                                      f"{ranked[0].get('record_id')} ({ranked[0].get('genotype_representation')}) kept"})
    kept.sort(key=lambda r: str(r.get("record_id")))
    dropped.sort(key=lambda r: str(r.get("record_id")))
    return {"kept": kept, "dropped": dropped}


def assign_kind(records, ancient_prefixes, modern_groups):
    """kind 由面板选择决定：现代/古代名单来自配置，不用年代阈值切片。

    两边都命中说明配置重叠，直接拒绝——否则同一条记录会同时进两个 PCA 集合。
    """
    modern = [str(x) for x in (modern_groups or [])]
    prefixes = [str(x) for x in (ancient_prefixes or [])]
    for p in prefixes:
        clash = [m for m in modern if m.startswith(p)]
        if clash:
            raise ValueError(f"aadr_modern and aadr_ancient_prefix overlap: {clash} would match prefix {p!r}")
    for r in records:
        pop = str(r.get("source_population_id") or "")
        if pop in modern:
            r["kind"] = "modern"
        elif any(pop.startswith(p) for p in prefixes):
            r["kind"] = "ancient"
        else:
            r["kind"] = "unknown"      # 不进现代/古代 PCA 集合
    return records


def unmatched_panel_entries(records, modern_groups, ancient_prefixes):
    """配置里点名了、但一条记录都匹配不到的条目（modern 精确匹配，ancient 前缀匹配）。

    匹配语义与 assign_kind 完全一致。一个拼错的群体名会静默掏空那一侧的 PCA 集合，
    分析表面上照常出结果——这里把没命中的配置点名交回去，让 08 打印并写进 manifest。
    """
    modern = [str(x) for x in (modern_groups or [])]
    prefixes = [str(x) for x in (ancient_prefixes or [])]
    pops = {str(r.get("source_population_id") or "") for r in records}
    return {"modern": [m for m in modern if m not in pops],
            "ancient_prefix": [p for p in prefixes if not any(s.startswith(p) for s in pops)]}


LOCATION_COLUMNS = ("dataset", "source_id", "location_id", "label_zh", "label_en", "locality",
                    "latitude", "longitude", "precision", "source_url", "note")


def load_locations(path):
    """读人工地点覆盖表，返回 {location_id: row}。缺列直接报错，不静默降级。

    表只用于可核查的纠错与翻译：没有中文名就留空，显示时回退到原始名称（不要求中文）。
    """
    path = pathlib.Path(path)
    header, rows = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if header is None:
            header = [p.strip() for p in parts]
            missing = [c for c in LOCATION_COLUMNS if c not in header]
            if missing:
                raise ValueError(f"{path.name} is missing columns {missing}; expected {list(LOCATION_COLUMNS)}")
            continue
        rows.append(dict(zip(header, parts)))
    out = {}
    for r in rows:
        lid = (r.get("location_id") or "").strip()
        if not lid:
            continue
        lat = _coord(r.get("latitude"), -90.0, 90.0)
        lon = _coord(r.get("longitude"), -180.0, 180.0)
        out[lid] = {"dataset": (r.get("dataset") or "").strip(), "source_id": (r.get("source_id") or "").strip(),
                    "location_id": lid, "label_zh": (r.get("label_zh") or "").strip(),
                    "label_en": (r.get("label_en") or "").strip(), "locality": (r.get("locality") or "").strip(),
                    "latitude": lat, "longitude": lon,
                    # 复审 P1（AN1/AN6）：人工表留空的精度不得因"有坐标"默认成 site——本表头
                    # 约定精度一律 region；留空=未核定，如实 unknown。
                    "precision": (r.get("precision") or "").strip() or "unknown",
                    "source_url": (r.get("source_url") or "").strip(), "note": (r.get("note") or "").strip()}
    return out


# ─────────────────────────────────────────── allele dosage

def a1_dosage_with_reason(gt, ref, alt, a1):
    """返回 (dosage, reason)。dosage 是 a1 在二倍体二等位 GT 里的拷贝数，无法确定时为 None。

    只接受 a1 与 REF/ALT 明确相同的记录；缺失、三倍性、越界索引、多等位一律记 None 并给出原因。
    不猜链翻转：无法解释的等位基因不会退化成"参考纯合"。
    """
    alts = [a for a in str(alt).split(",") if a not in ("", ".")]
    if len(alts) != 1 or ref is None:
        return None, "not_biallelic"
    alleles = gt_alleles(gt, ref, alts[0])
    if alleles is None:
        return None, "all_missing"
    if len(alleles) != 2:
        return None, "not_diploid"
    if any(a is None for a in alleles):
        return None, "partial_missing"
    if any(a not in (ref, alts[0]) for a in alleles):
        return None, "allele_not_observable"
    if a1 == ref:
        return sum(1 for a in alleles if a == ref), ""
    if a1 == alts[0]:
        return sum(1 for a in alleles if a == alts[0]), ""
    return None, "a1_not_ref_or_alt"


def a1_dosage(gt, ref, alt, a1):
    return a1_dosage_with_reason(gt, ref, alt, a1)[0]


# ─────────────────────────────────────────── AN2: one eligible set, one grouping

def coverage_from_smiss(imiss):
    """plink2 .smiss 表 → n_called_snps / call_rate（09b 与 04b 共用的唯一口径）。

    formats#smiss（并对照真实文件核过：OBS_CT == prune.prune.in 行数）：OBS_CT 是**分母**（总位
    点数），MISSING_CT 才是缺失数。已调用数是两者之差；把 OBS_CT 当成已调用数会把缺失位点计入
    n_called_snps，位点门槛随之失真。F_MISS = MISSING_CT/OBS_CT，故 call_rate = 1 - F_MISS 同口径。
    输入是已按 IID 重命名好 iid 列（或保留 IID）的 DataFrame；返回同一张表附上两列。
    """
    imiss = imiss.copy()
    if "IID" in imiss.columns and "iid" not in imiss.columns:
        imiss = imiss.rename(columns={"IID": "iid"})
    imiss["n_called_snps"] = imiss["OBS_CT"].astype(int) - imiss["MISSING_CT"].astype(int)
    imiss["call_rate"] = 1.0 - imiss["F_MISS"].astype(float)
    return imiss


def _number_or_none(value):
    """pandas/numpy 的 NaN 与 inf 都算缺失：它们不能拿去比较，也不该被 int() 撞出异常。"""
    if value is None:
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(x) or math.isinf(x)) else x


def eligible_records(records, kind, min_rate=None, min_snps=None):
    """默认榜单、PCA 与地图共用的唯一合格集合。

    就地标注 `eligible` / `exclusion_reason`（被排除的记录留在原地供详情使用），返回合格记录并
    保持输入顺序。标注**只作用于 kind 与本次请求相同的记录**：09b 拿同一份列表连筛 modern 再筛
    ancient，若每次都改写全部记录，第二遍会把第一遍判合格的 modern 全部改成 other_kind——两种
    kind 的结论因此互不覆盖，每条记录最终带的是它自己那一类的判定；kind 之外的记录（含 target
    与未知 kind）不带任何合格标注，也不混进任何一方的榜单。门槛只对调用方显式给出的维度生效
    （None = 该维度不设门槛）。
    """
    out = []
    for r in records:
        if r.get("kind") != kind:
            continue
        reason = None
        if _number_or_none(r.get("distance_to_target")) is None:
            reason = "no_distance"
        elif min_rate is not None:
            rate = _number_or_none(r.get("call_rate"))
            if rate is None or rate < float(min_rate):
                reason = "low_call_rate"
        if reason is None and min_snps is not None:
            n = _number_or_none(r.get("n_called_snps"))
            if n is None or n < float(min_snps):
                reason = "insufficient_sites"
        r["eligible"] = reason is None
        r["exclusion_reason"] = reason
        if reason is None:
            out.append(r)
    return out


def group_summaries(records, min_group_n):
    """按 group_id 汇总：先有每个个体的距离，再取组均值。

    不用"群体质心 vs 目标"另算一套距离——那与个体距离不是同一口径，会让两个榜单给出不同的
    近邻顺序。`rank` 只给达到 min_group_n 的群体；样本不足的组仍然可查（small_group=True，
    rank=None），但不混进默认前几名。平手按稳定的 group_id 排序。
    """
    buckets = {}
    for r in records:
        gid = r.get("group_id")
        if gid is None or r.get("distance_to_target") is None:
            continue
        buckets.setdefault(str(gid), []).append(r)
    out = []
    for gid, rows in buckets.items():
        ds = [float(r["distance_to_target"]) for r in rows]
        # 成员键：AADR 记录有 record_id；1000G 记录没有，只有 sscore 的 iid。直接 str(record_id)
        # 会把每个成员都写成 "None" 字面量——下游按 member_ids 关联记录时全体失配。
        _mem = lambda r: str(r.get("record_id") if r.get("record_id") is not None else r.get("iid"))
        rows_sorted = sorted(rows, key=_mem)
        first = rows_sorted[0]
        dmin = [r.get("date_min_bp") for r in rows if r.get("date_min_bp") is not None]
        dmax = [r.get("date_max_bp") for r in rows if r.get("date_max_bp") is not None]
        out.append({
            "group_id": gid,
            "label": first.get("label") or gid,
            "kind": first.get("kind"),
            "location_id": first.get("location_id"),
            "n": len(rows),
            "member_ids": [_mem(r) for r in rows_sorted],
            "rank": None,
            "small_group": len(rows) < int(min_group_n),
            "distance_mean": sum(ds) / len(ds),
            "distance_min": min(ds),
            "distance_max": max(ds),
            "date_min_bp": min(dmin) if dmin else None,
            "date_max_bp": max(dmax) if dmax else None,
        })
    ranked = sorted((g for g in out if not g["small_group"]),
                    key=lambda g: (g["distance_mean"], g["group_id"]))
    for i, g in enumerate(ranked, 1):
        g["rank"] = i
    out.sort(key=lambda g: (g["small_group"], g["distance_mean"], g["group_id"]))
    return out


# ─────────────────────────────────────────── AN3: holdout selection

def pick_holdout(rows, pops, n, seed):
    """按固定 seed 无放回抽取校准 holdout。

    每个配置群体取 `min(n, floor(群体人数/5))` 人：留下 1/5 在外面，剩下的仍是可用的参考面板。
    抽不满 2 人、或抽取后该群体剩余参考不足 2 人 → 返回 unavailable，而不是拿一两个人硬估 SD；
    未配置的群体不参与；配置里写了参考面板没有的群体则明确报错。
    """
    by_pop = {}
    for r in rows:
        pop, iid = str(r.get("population") or ""), str(r.get("iid") or "")
        if pop and iid:
            by_pop.setdefault(pop, []).append(iid)
    wanted = [str(p) for p in (pops or [])]
    if not wanted:
        return {"state": "unavailable", "reason_code": "no_populations_configured",
                "detail": "local_ancestry_calibration_pops is empty and no source panel is configured",
                "ids": [], "per_pop": {}}
    missing = [p for p in wanted if p not in by_pop]
    if missing:
        return {"state": "unavailable", "reason_code": "unknown_population",
                "detail": f"not present in the reference panel: {missing}", "ids": [], "per_pop": {}}
    rng = random.Random(int(seed))
    ids, per_pop = [], {}
    for p in wanted:
        pool = sorted(set(by_pop[p]))          # 排序后抽样只依赖 iid 集合与 seed，与输入顺序无关
        want = min(int(n), len(pool) // 5)
        if want < 2:
            return {"state": "unavailable", "reason_code": "insufficient_holdout",
                    "detail": f"{p}: {len(pool)} individuals give {want} holdouts, need >=2",
                    "ids": [], "per_pop": {}}
        if len(pool) - want < 2:
            return {"state": "unavailable", "reason_code": "insufficient_reference",
                    "detail": f"{p}: {len(pool) - want} individuals would remain in the reference, need >=2",
                    "ids": [], "per_pop": {}}
        ids += rng.sample(pool, want)
        per_pop[p] = want
    return {"state": "ok", "reason_code": "", "detail": "", "ids": ids, "per_pop": per_pop}


def read_psam(path, population_column="Population"):
    """读 1000G .psam（#IID 与 Population 两列）。只做解析，不做选择。"""
    rows, header = [], None
    for line in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if header is None:
            header = [p.lstrip("#") for p in parts]
            continue
        rec = dict(zip(header, parts))
        rows.append({"iid": rec.get("IID", ""), "population": rec.get(population_column, "")})
    return rows


def read_vcf_samples(path):
    """读 VCF/BCF 的样本名列表（只解析 #CHROM 头；bgzip 由调用方负责解压）。"""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                return cols[9:]
            break
    return []


if __name__ == "__main__":
    raise SystemExit(_cli())


# ─────────────────────────────────────────── AN5: 事实文案由数据生成

def _n(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def _fmt_d(v, digits=4):
    try:
        return f"{float(v):.{digits}f}"
    except (TypeError, ValueError):
        return "—"


# 1000G 面板的群体标签是拉丁转写；zh 文案按面板命名规则转写，不认识的标签原样返回（宁缺毋滥）。
_POP_ZH = {"Chongqing": "重庆", "Sichuan": "四川", "Hubei": "湖北", "Beijing": "北京",
           "Southern": "南方", "Northern": "北方", "Fujian": "福建"}

def _pop_zh(label):
    """Han_Chongqing → 重庆汉族；China_* 的古代群体已有 name_zh，由调用方使用。"""
    if not isinstance(label, str):
        return str(label or "")
    if label.startswith("Han_") and len(label) > 4:
        place = label[4:]
        return f"{_POP_ZH.get(place, place)}汉族"
    return label

def _pop_en(label):
    """Han_Chongqing → Han from Chongqing；其余原样。"""
    if isinstance(label, str) and label.startswith("Han_") and len(label) > 4:
        return f"Han from {label[4:].replace('_', ' ')}"
    return label

def _iid_readable(iid):
    """BaiyangcunM13.SG → Baiyangcun M13：去掉 .SG/.AG/.DG 等技术后缀，站点名与编号间加空格。"""
    s = str(iid).split(".")[0]
    import re
    m = re.match(r"^([A-Za-z][a-z]+?)([A-Z]\d+)$", s)
    return f"{m.group(1)} {m.group(2)}" if m else s

def _grp_name_zh(anc_rows, group):
    """按 group 标签找古代群体的 name_zh（如 China_MLBA → 中国 中晚期青铜）；没有就用原标签。"""
    for r in anc_rows:
        if r.get("label") == group and r.get("name_zh"):
            return r["name_zh"]
    return group or ""


def ancestry_copy(D):
    """由数据生成祖源与父母系的**事实**文案，返回 {key: [zh, en]}。

    这些键的手写文案会被这里的值覆盖（构建器会打印被替代的键），因为手写的版本里往往写死了某个
    样本的具体值——「最接近云南白羊村」「N-M1845 之下进入分辨率极限」——换成别的样本就成了错话。
    措辞只声明数据支持的内容：距离在"所选参考面板中"的排名、明确的样本数、实际的单倍群与树版本；
    不把"最接近"读成族群身份，也不在没有证据时给地理结论。
    """
    out = {}
    aff = list(D.get("ho_affinity") or [])
    mods = [r for r in aff if r.get("kind") == "modern"]
    ancs = [r for r in aff if r.get("kind") == "ancient"]
    near_i = (D.get("ho_near_individual") or [None])[0]
    acc = D.get("ho_accession") or {}
    n_rec = _n((acc.get("counts") or {}).get("selected")) if isinstance(acc, dict) else None
    n_elig = _n((acc.get("counts") or {}).get("eligible")) if isinstance(acc, dict) else None
    n_coord = _n((acc.get("counts") or {}).get("mapped")) if isinstance(acc, dict) else None

    # --- 祖源标题（W3 A1 口径）：最近的现代**群体**与最近的古代**基因组**分开说，不混用统计量
    if aff:
        m = mods[0] if mods else None
        a = ancs[0] if ancs else None
        zh_parts, en_parts = [], []
        if m:
            zh_parts.append(f"现代人群里最近的是{_pop_zh(m['label'])}")
            en_parts.append(f"the nearest present-day group is {_pop_en(m['label'])}")
        if near_i:
            zh_parts.append(f"古代基因组里距离最近的是 {_iid_readable(near_i['iid'])}"
                            f"（{_grp_name_zh(ancs, near_i.get('group'))}）")
            en_parts.append(f"the nearest ancient genome is {_iid_readable(near_i['iid'])} ({near_i.get('group')})")
        elif a:
            zh_parts.append(f"最近的古代群体是 {_pop_zh(a['label'])}（n={a['n']}，d={_fmt_d(a['d'])}）")
            en_parts.append(f"nearest ancient group is {a['label']} (n={a['n']}, d={_fmt_d(a['d'])})")
        out["c_anc"] = ["在所选参考面板中，" + "；".join(zh_parts),
                        "Within the selected reference panel: " + "; ".join(en_parts)]
    else:
        out["c_anc"] = ["本版没有可用的古 DNA 投影结果",
                        "No usable ancient-DNA projection in this build"]

    # --- 祖源副标题：面板规模与坐标覆盖率（明确 n，不写"130 万位点"这类只有某次运行才知道的数）
    bits_zh, bits_en = [], []
    if n_rec:
        bits_zh.append(f"{n_rec} 条参考记录"); bits_en.append(f"{n_rec} reference records")
    if n_elig:
        bits_zh.append(f"合格 {n_elig}"); bits_en.append(f"{n_elig} eligible")
    if n_coord is not None and n_rec:
        bits_zh.append(f"有坐标 {n_coord}"); bits_en.append(f"{n_coord} with coordinates")
    if mods or ancs:
        bits_zh.append(f"群体 {len(mods)} 现代 / {len(ancs)} 古代")
        bits_en.append(f"{len(mods)} modern / {len(ancs)} ancient groups")
    if bits_zh:
        out["sub_anc"] = [" · ".join(bits_zh), " · ".join(bits_en)]

    # --- 祖源脚注（W3 A2 口径）：最近的现代人群（前几个）、最近古代个体（带 call_rate）、
    #     个体噪声与群体均值不同量级的提示、不做亲缘/族群推断的边界声明
    if aff and (mods or ancs):
        top_m = mods[:3]
        zh_m, en_m = "", ""
        if top_m:
            places = [r["label"][4:] for r in top_m if str(r["label"]).startswith("Han_")]
            zh_places = [_POP_ZH.get(p, p) for p in places]
            if len(places) == len(top_m):          # 前几名全是 Han_* 时合并成「重庆、四川的汉族样本」
                zh_m = "、".join(zh_places) + "的汉族样本"
                en_m = "Han samples from " + ", ".join(p.replace("_", " ") for p in places)
            else:
                zh_m = "、".join(_pop_zh(r["label"]) for r in top_m)
                en_m = ", ".join(_pop_en(r["label"]) for r in top_m)
        zh_bits, en_bits = [], []
        if zh_m:
            zh_bits.append(f"最近的现代人群是{zh_m}")
            en_bits.append(f"The nearest present-day groups are {en_m}")
        if near_i:
            _cr = None
            for mem in (D.get("ho_target_group") or {}).get("members") or []:
                if mem.get("iid") == near_i.get("iid") and _n(mem.get("call_rate")) is not None:
                    _cr = float(mem["call_rate"])
                    break
            zh_bits.append(f"距离最近的古代基因组是 {_iid_readable(near_i['iid'])}"
                           f"（d={_fmt_d(near_i['d'], 3)}，call_rate {f'{_cr:.2f}' if _cr is not None else '—'}）")
            en_bits.append(f"The nearest ancient genome is {_iid_readable(near_i['iid'])} "
                           f"(d={_fmt_d(near_i['d'], 3)}, call_rate {f'{_cr:.2f}' if _cr is not None else '—'})")
            _gz = _grp_name_zh(ancs, near_i.get("group"))
            _gd = next((r["d"] for r in ancs if r.get("label") == near_i.get("group")), None)
            zh_bits.append(f"单个基因组的距离噪声大，且与其所属群体均值（{_gz}，{_fmt_d(_gd, 3) if _gd is not None else '—'}）"
                           f"不在同一量级")
            en_bits.append(f"A single genome's distance is noisy and is not on the same scale as its group's "
                           f"mean ({near_i.get('group')}, {_fmt_d(_gd, 3) if _gd is not None else '—'})")
        zh_bits.append("这里只描述“与哪些人群最接近”，不由此做亲缘或族群推断")
        en_bits.append("This section states only which groups lie closest and draws no kinship or ethnic inference")
        out["n_anc"] = ["；".join(zh_bits) + "。", ". ".join(en_bits) + "."]

    # --- 父系正文：从实际路径生成，不写死任何支系名
    ypath = list(D.get("ypath") or [])
    yt, yc = D.get("y_terminal"), D.get("y_conservative")
    unc = list(D.get("y_uncertain") or [])
    if ypath and yt:
        solid = [p for p in ypath if _n(p.get("der")) is not None and int(p["der"]) >= 5]
        trunk = f"{solid[0]['snp']}→{solid[-1]['snp']}" if solid else ypath[0].get("snp")
        zh = (f"Y 染色体只从父亲传给儿子，记录的是父系一条线。本次判定沿 {D.get('y_tree') or 'YFull'} 树逐级下行，"
              f"共 {len(ypath)} 级；支持充分的节点从 {trunk}。")
        en = (f"The Y chromosome passes only from father to son, so it records one paternal line. This call "
              f"walks {len(ypath)} levels of the {D.get('y_tree') or 'YFull'} tree; the solidly supported nodes run {trunk}.")
        if unc:
            zh += f"末端 {len(unc)} 级的支持位点不足（{'、'.join(str(u) for u in unc)}），保守回退到 {yc}。"
            en += (f" The last {len(unc)} level(s) rest on too few sites ({', '.join(str(u) for u in unc)}); "
                   f"the conservative fallback is {yc}.")
        elif yc and yc != yt:
            zh += f"保守回退为 {yc}。"; en += f" Conservative fallback: {yc}."
        out["n_y"] = [zh, en]
        out["c_y"] = [f"父系：Y 染色体属于 {yt}" + (f"（保守回退 {yc}）" if yc and yc != yt else ""),
                      f"Father's line: {yt}" + (f" (conservative fallback {yc})" if yc and yc != yt else "")]
    elif D.get("y_state") == "unavailable":
        out["c_y"] = ["父系：本次交付没有得到可用的 Y 判定（原因见状态表）",
                      "Father's line: no usable Y call in this delivery (see the status table)"]
    return out


# ─────────────────────────────────────────── AN5: 载荷卫生与分析状态（供 30/31 共用）

def json_clean(obj):
    """递归把 NaN/Infinity 变成 null，并把 numpy 标量转成原生类型。

    json.dumps 默认把 NaN 写成裸 NaN —— 不是合法 JSON，浏览器的 JSON.parse 会拒绝；而这些值来自
    pandas 的缺失，语义本来就是 null。
    """
    if isinstance(obj, float):
        return obj if (obj == obj and obj not in (float("inf"), float("-inf"))) else None
    if isinstance(obj, dict):
        return {k: json_clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_clean(v) for v in obj]
    try:
        import numpy as _np
        if isinstance(obj, _np.generic):
            return json_clean(obj.item())
    except Exception:
        pass
    return obj


def json_text(obj):
    """合法 JSON 文本（allow_nan=False / NaN 已转 null）。"""
    return json.dumps(json_clean(obj), ensure_ascii=False, allow_nan=False)


def js_string_literal(text):
    """写成 JS 字符串字面量，并断开闭合标签序列。

    数据里若含 "</script"（例如某个显示名），浏览器会在那里结束脚本块，剩下的载荷变成页面文字、
    图表全空 —— 而所有静态检查仍然通过。把斜杠转义后 JSON 与 JS 都把它当普通斜杠，数据不变、
    解析器不再提前收尾。
    """
    return json.dumps(text, ensure_ascii=False).replace("</", "<\\/")


def analysis_state(dir_path, expected=None, names=("summary.json", "local_ancestry.json"),
                   expected_parameters=None):
    """读某分析目录的 manifest 与结构化结果，返回 (state, reason_code, doc)。

    缺 manifest → missing_manifest（旧结果必须重建，不能当current用）；指纹不符 → stale_result；
    结果读不出来 → unreadable_result；没有结果文件 → missing_result。**不**从"文件在不在"推断，
    也不读旧路径猜结果。

    expected_parameters：把 manifest.parameters 里的**配置绑定键**逐一比对（复审 AN0/AN5/H6：
    30 此前只比 sample_id，同一样本换门槛/换参考/换 prune 集后旧结果照样进场）。调用方传它**当前
    有效配置**能推导出的键；任一键缺失或不等 → stale_result。数据依赖键（如 missing_chroms）不
    属于此列，不要传。
    """
    d = pathlib.Path(dir_path)
    # 步骤级 manifest（manifest.<step-id>.json）由 run_all 在"该步被配置禁用或失败"时写。它**优先于**
    # 目录里的 manifest.json：那一步没有产出，目录里若有结果只能来自更早的运行，用它等于让禁用状态
    # 失效（复审 AN0+AN5+H6）。只有所有步骤级记录都是 ok，才回落到目录级 manifest。
    # 第一版把这段放在"找不到 manifest.json 时"才走，实测两种位置都返回 ok——因为真实 manifest 存在，
    # glob 根本不会被走到。
    man = None
    _step_non_ok = False
    for extra in sorted(d.glob("manifest.*.json")) if d.exists() else []:
        cand = read_manifest(extra)
        if cand and str(cand.get("state") or "ok") not in ("ok",):
            man = cand
            _step_non_ok = True
            break
        if cand and man is None:
            man = cand
    if not man:
        man = read_manifest(d / "manifest.json")
    if not man:
        return "unavailable", "missing_manifest", None
    # 逐键比对标识（sample_id / analysis_id / 参考版本…），但**不**把 state 当准入条件：state 是这份
    # manifest 要如实报告的结论，disabled/failed 必须原样传出去，不能一律折叠成 unavailable。
    # （manifest_matches 另有一条 "state 必须是 ok" 的规则，适用于"能不能复用结果"，不是这里。）
    for k, v in (expected or {}).items():
        if v in (None, ""):
            continue
        if str(man.get(k)) != str(v):
            return "unavailable", "stale_result", None
    # 配置绑定参数：manifest.parameters 里缺键或值不符同样是 stale_result——"没记录"不等于"一致"。
    _params = man.get("parameters") or {}
    for k, v in (expected_parameters or {}).items():
        if k not in _params or str(_params[k]) != str(v):
            return "unavailable", "stale_result", None
    # 步骤级禁用/失败记录压过目录 manifest 时，目录里的结果文件属于**上一次**运行——状态照实
    # 报告，但结果不能作为 doc 交出（与目录自身 manifest.json 标 disabled、由该步写入原因文件的
    # 情形不同：那种结果属于本次运行，读出来是安全的）。
    if _step_non_ok:
        return str(man.get("state") or "ok"), str(man.get("reason_code") or ""), None
    for name in names:
        f = d / name
        if f.exists():
            try:
                return str(man.get("state") or "ok"), str(man.get("reason_code") or ""), json.loads(f.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                return "unavailable", "unreadable_result", None
    return "unavailable", "missing_result", None
