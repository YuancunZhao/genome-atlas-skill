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
SUPPORTED_BUILDS = ("GRCh37",)

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

    build = _str(cfg, "build", "GRCh37") or "GRCh37"
    if build not in SUPPORTED_BUILDS:
        raise ValueError(
            f"build {build!r} is not supported (this pipeline is GRCh37 only); refusing before any work"
        )

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

    规则：actual 必须存在、含 sample_id、且（若写了 state）state == "ok"；expected 里的每个
    键只要在 actual 中出现就必须相等。任何不等、任何关键键缺失都返回 False——调用方据此重建，
    绝不复用旧结果。
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
            if k in REQUIRED_MANIFEST_KEYS:
                return False
            continue
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
            "location_precision": "site" if (lat is not None and lon is not None) else "unknown",
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
                    "precision": (r.get("precision") or "").strip() or ("site" if lat is not None else "unknown"),
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
    保持输入顺序。`kind` 之外的一律排除，所以现代/古代/未分类不会互相混进对方的榜单；门槛只对
    调用方显式给出的维度生效（None = 该维度不设门槛）。
    """
    out = []
    for r in records:
        reason = None
        if r.get("kind") != kind:
            reason = "other_kind"
        elif _number_or_none(r.get("distance_to_target")) is None:
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
        rows_sorted = sorted(rows, key=lambda r: str(r.get("record_id")))
        first = rows_sorted[0]
        dmin = [r.get("date_min_bp") for r in rows if r.get("date_min_bp") is not None]
        dmax = [r.get("date_max_bp") for r in rows if r.get("date_max_bp") is not None]
        out.append({
            "group_id": gid,
            "label": first.get("label") or gid,
            "kind": first.get("kind"),
            "location_id": first.get("location_id"),
            "n": len(rows),
            "member_ids": [str(r.get("record_id")) for r in rows_sorted],
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
