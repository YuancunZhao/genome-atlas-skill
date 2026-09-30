#!/usr/bin/env python
"""Paternal/maternal lineage results and their historical evidence (AN4).

Three jobs, all shaped by the same rule: never turn missing evidence into a finding.

* `match_lineage` answers what a stored observation is relative to the sample's own branch, and it
  refuses to guess: without tree edges only a trusted, same-version canonical id can be `exact`,
  everything else is `unresolved`. String prefixes are never used -- "A13".startswith("A1") is exactly
  the kind of shortcut that puts one branch's records on another branch.
* `lineage_summary` normalises what 05/06 wrote. The deepest node in the tree is *not* automatically the
  proven terminal: the caller's call, the conservative fallback and the nodes that need review are kept
  apart, and an empty path is an unavailable result rather than a crash (the report used to read
  `path[-1]`).
* `lineage_observations` / `history_view` assemble the historical view from the normalised metadata. It
  deliberately does *not* apply the autosomal call-rate gate: a 0.49-coverage ancient genome with a
  published mt label is evidence about that lineage, while the autosomal gate is about PCA. Routes are
  only ever shown when a versioned, sourced route exists; otherwise the view says it is a distribution.

Y and mt names are namespaced (`y:`, `mt:`) because the same label can exist in both trees.
"""
import json, os, pathlib, sys


# AADR 用这些字符串表示"这一行没有单倍群标签"（女性和性别未知的个体）。它们不是支系名，必须当缺失
# 处理，否则 "n/a (female)" 会被当成一个群体名去和树比较。
MISSING_HG_PLACEHOLDERS = ("", ".", "..", "nan", "n/a", "n/a (female)", "n/a (sex unknown)",
                           "n/a (sex)", "unknown", "na")


def is_missing_hg(value):
    return str(value or "").strip().lower() in MISSING_HG_PLACEHOLDERS


def _ns(node):
    s = str(node or "")
    return s.split(":", 1)[0] if ":" in s else ""


def _ancestors(node, edges):
    """(ancestors, in_tree) —— 沿 edges 向上走；环或超长链会被截断，不无限循环。"""
    out, cur, guard = [], node, 0
    while cur is not None and cur in edges and guard < 10000:
        cur = edges[cur]
        guard += 1
        if cur is not None:
            out.append(cur)
    return out, (node in edges)


def match_lineage(query, record, parents, same_tree=True):
    """记录相对查询支系的关系：exact | descendant | ancestor | unrelated | unresolved。

    `parents` 是有来源的树边（子 → 父）。缺边时不猜：只有同版本且 ID 完全相同的记录可以 exact，
    其余一律 unresolved。查询与记录必须属于同一命名空间（Y 与 mt 可以重名）。
    """
    q, r = str(query or ""), str(record or "")
    if not q or not r or _ns(q) != _ns(r):
        return "unresolved"
    if q == r:
        return "exact" if same_tree else "unresolved"
    edges = parents or {}
    aq, q_in = _ancestors(q, edges)
    ar, r_in = _ancestors(r, edges)
    if r in aq:
        return "ancestor"          # 记录是查询支系的祖先（大支系背景）
    if q in ar:
        return "descendant"        # 记录在查询支系下游
    if q_in and r_in:
        return "unrelated"
    return "unresolved"


def load_tree_nodes(tree_path, kind):
    """从 YFull 的 current_tree.json 收集规范节点 ID（带命名空间），用于判断"这个标签在当前树里吗"。"""
    try:
        tree = json.loads(pathlib.Path(tree_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    out = set()

    def walk(n):
        nid = str(n.get("id") or "")
        if nid:
            out.add(f"{kind}:{nid}")
        for c in n.get("children", []):
            walk(c)
    walk(tree)
    return out


def canonicalize(node, kind, history, known_nodes=None):
    """把历史标签规范化到当前树版本的节点 ID，返回 (canonical, note)。

    note 说明做了什么：missing（本来就没有标签）/ alias（按别名表改写）/ none（树里有这个节点）/
    version_mismatch（树里没有——多半是标签来自另一个树版本，例如 AADR 的 YFull 12.03）/
    unverified（没拿到树节点集合，无法判断）。**绝不**假设"表里没有就等于无关支系。
    """
    raw = str(node or "").strip()
    if is_missing_hg(raw):
        return "", "missing"
    key = raw if ":" in raw else f"{kind}:{raw}"
    aliases = (history or {}).get("aliases") or {}
    if key in aliases:
        return str(aliases[key]), "alias"
    if known_nodes is None:
        return key, "unverified"
    if key in known_nodes:
        return key, "none"
    return key, "version_mismatch"


# ---------------------------------------------------------------- 05/06 的结果规范化

def unavailable_lineage(kind, reason_code, detail=""):
    """不可用也要有明确原因（D1：未测量不是阴性）。"""
    return {"kind": str(kind), "state": "unavailable", "reason_code": str(reason_code), "detail": str(detail),
            "reported_hg": None, "conservative_hg": None, "tree_source": "", "tree_version": "",
            "call_quality": "none", "supported_path": [], "uncertain_nodes": [], "conflicts": [],
            "route_review": "none"}


def lineage_summary(result):
    """把 05/06 的原始结果折成报告层要的形状；空路径、缺判定都变成明确的不可用状态。"""
    r = dict(result or {})
    path = list(r.get("supported_path") or [])
    conflicts = list(r.get("conflicts") or [])
    reported = r.get("reported_hg")
    state = str(r.get("state") or ("ok" if reported else "unavailable"))
    reason = str(r.get("reason_code") or "")
    if state == "ok" and not path:
        # 有判定但没有可展示的支持路径：仍然可用（例如调用来自复核记录），但要说清没有路径
        if not reported:
            state, reason = "unavailable", reason or "no_supported_path"
    if state == "ok" and not reported:
        state, reason = "unavailable", reason or "no_haplogroup_call"
    return {
        "kind": str(r.get("kind") or ""),
        "state": state, "reason_code": reason,
        "reported_hg": reported,
        "conservative_hg": r.get("conservative_hg"),
        "terminal": (str(path[-1].get("node")) if path and path[-1].get("node") else None),
        "call_quality": str(r.get("call_quality") or "none"),
        "tree_source": str(r.get("tree_source") or ""), "tree_version": str(r.get("tree_version") or ""),
        "supported_path": path,
        "uncertain_nodes": list(r.get("uncertain_nodes") or []),
        "conflicts": conflicts,
        # 有未解决的冲突时，严格匹配视图要排除这条记录，而不是挑一个说法继续用
        "strict_match_allowed": not conflicts,
        "route_review": str(r.get("route_review") or "none"),
    }


# ---------------------------------------------------------------- 历史观测与路线

def _dedupe_by_master(rows):
    """同一个人可能有多条技术表示（.SG/.AG/…）；历史视图里只算一次，优先高覆盖表示。"""
    def rank(r):
        rep = str(r.get("genotype_representation") or r.get("record_id", "")).upper()
        for i, tok in enumerate(("SG", "DG", "HO", "AG", "TW")):
            if rep.endswith(tok) or rep == tok:
                return i
        return 99
    best = {}
    for r in rows:
        key = str(r.get("master_id") or r.get("individual_id") or r.get("record_id") or "")
        if not key:
            continue
        if key not in best or rank(r) < rank(best[key]):
            best[key] = r
    return [best[k] for k in sorted(best)]


def lineage_observations(rows, query, parents, tree_kind, same_tree=True, history=None,
                         known_nodes=None, stats=None):
    """筛出与该支系相关的历史记录（含祖先/后代/未定），并带上坐标、年代与来源。

    这里**不用**常染色体 call_rate 门槛：来源对 Y/mt 的可用性决定记录是否可用，覆盖率低只影响
    常染色体 PCA。没有相应单倍群标签的记录直接跳过——不猜。
    """
    kind = str(tree_kind)
    field = "y_hg_raw" if kind == "y" else "mt_hg_raw"
    field_tree = "hg_source_tree"
    out = []
    for r in _dedupe_by_master([x for x in (rows or []) if not is_missing_hg(x.get(field))]):
        node, note = canonicalize(r.get(field), kind, history, known_nodes)
        if stats is not None and note != "missing":
            stats[note] = stats.get(note, 0) + 1     # 进了比较就记账，包括版本不匹配的那些
        rel = match_lineage(query, node, parents, same_tree=bool(same_tree))
        if note == "version_mismatch":
            # 复审 AN4：标签不在当前树节点集（多半来自另一个树版本，如 AADR 的 YFull 12.03）时，
            # 与查询支系字符串相等不等于版本等价——此前 version_mismatch 只记计数，same_tree=True
            # 照样可判 exact。未验证版本等价的记录不进默认视图，宁可 unresolved。
            rel = "unresolved"
        lat, lon = r.get("latitude"), r.get("longitude")
        out.append({
            "record_id": str(r.get("record_id") or ""),
            "node_id": node,
            "relation": rel,
            "locality": r.get("locality"),
            "coordinates": ({"latitude": lat, "longitude": lon} if (lat is not None and lon is not None) else None),
            "precision": str(r.get("location_precision") or "unknown"),
            "date_range": {"mean": r.get("date_mean_bp"), "min": r.get("date_min_bp"),
                           "max": r.get("date_max_bp")},
            "date_basis": r.get("date_basis"),
            "call_source": str(r.get(field_tree) or r.get("hg_call_source") or ""),
            "publication": r.get("publication"),
        })
    # 只保留能与查询支系建立关系的记录：exact 与后代默认显示，祖先作为背景层（视图决定怎么展开）。
    # unrelated 与 unresolved 不进默认视图——缺树边时无法确认，宁可说"未定"也不把别的支系算进来。
    out = [o for o in out if o["relation"] in ("exact", "descendant", "ancestor")]
    order = {"exact": 0, "descendant": 1, "ancestor": 2}
    out.sort(key=lambda o: (order.get(o["relation"], 9), o["record_id"]))
    return out


def conservative_from_path(path, solid=5, tail=4):
    """没有复核记录时的通用保守落点。

    判据不是"某个节点的支持数够不够"，而是"它是否位于末端一条弱链的末尾"。本项目的真实路径是
    12 / 1 / 0 / 1 / 5：末端那个 5 恰好达到阈值，但它上面三级分别是 1、0、1——它们是同一段分辨率
    极限，5 只是这段弱链的末尾，不是独立证据。只看单节点阈值会把它当成可靠的末端（这正是第一版
    规则错的地方）。所以：末端 tail 级里若有支持不足的节点，就从最靠上的那个再往上退一级。

    path 每项至少含 (node, der, ...)。
    """
    if not path:
        return None
    seg = list(path[-tail:]) if tail and tail > 0 else list(path)
    weak = [i for i, p in enumerate(seg) if int(p[1]) < int(solid)]
    if not weak:
        return str(path[-1][0])
    first_weak = seg[weak[0]]
    pos = list(path).index(first_weak)
    return str(path[pos - 1][0]) if pos > 0 else None


def reviewed_call(history, kind, reported_hg, sample_id=None, log=None):
    """带理由的复核记录（panel/lineage_history.json 的 reviewed_calls）优先于任何自动规则。

    §7 要求 conservative_hg 来自实际证据或带理由的复核记录，而不是"树里最深节点"或某个阈值。
    键形如 `y:N-CTS4714`。

    **样本限定（复审 AN4）**：键只写支系，所以一个样本的人工复核结论会被**任何**调出同一支系的
    新样本继承——那是把一个样本的判断当成通用规则。因此记录里必须写明 sample_id，且只在**当前样本
    与之相符**时才生效；没有 sample_id 的旧记录一律不套用（并说明原因），由调用方回退到通用规则。
    """
    hist = history or {}
    node = str(reported_hg or "")
    if not node:
        return None, None
    key = node if ":" in node else f"{kind}:{node}"
    rec = (hist.get("reviewed_calls") or {}).get(key)
    if not isinstance(rec, dict) or not rec.get("conservative_hg"):
        return None, None
    owner = rec.get("sample_id")
    if not owner:
        # 旧记录没有样本标记：不得当成通用规则套给任何样本。
        if log:
            log(f"reviewed call for {key} carries no sample_id and is NOT applied; "
                f"it is one sample's review, not a rule (add sample_id to use it)")
        return None, None
    if sample_id is None or str(owner) != str(sample_id):
        if log:
            log(f"reviewed call for {key} belongs to {owner}, not to this sample; "
                f"falling back to the generic rule")
        return None, None
    return str(rec["conservative_hg"]), str(rec.get("reason") or "reviewed by hand")


HISTORY_STATES = ("ok", "distribution_only", "unavailable")


def load_history(doc):
    """读 panel/lineage_history.json：版本、来源、树边/别名与经核对的路线。"""
    d = dict(doc or {})
    parents = {str(k): (str(v) if v is not None else None) for k, v in (d.get("parents") or {}).items()}
    routes = []
    for r in (d.get("routes") or []):
        if not isinstance(r, dict):
            continue
        if not r.get("route_id") or not r.get("applicable_node") or not r.get("tree_version"):
            continue
        if not (r.get("waypoints") or []) or not (r.get("sources") or []):
            continue          # 没有航点或没有出处的"路线"不是路线
        routes.append(r)
    return {"schema_version": int(d.get("schema_version") or 1),
            "tree_source": str(d.get("tree_source") or ""),
            "tree_version": str(d.get("tree_version") or ""),
            "parents": parents,
            "aliases": {str(k): str(v) for k, v in (d.get("aliases") or {}).items()},
            "reviewed_calls": {str(k): dict(v) for k, v in (d.get("reviewed_calls") or {}).items()
                               if isinstance(v, dict)},
            "routes": routes,
            "sources": list(d.get("sources") or [])}


def history_view(query, history, observations, parents, same_tree=True):
    """历史视图：观测 + 与查询支系匹配的路线；没有路线时明说是分布视图，不造故事。"""
    hist = history or {}
    walk = dict(parents or hist.get("parents") or {})
    matched, seen = [], set()
    for r in (hist.get("routes") or []):
        rid = str(r.get("route_id"))
        if rid in seen:
            continue
        if match_lineage(query, r.get("applicable_node"), walk, same_tree=same_tree) in ("exact", "descendant"):
            matched.append(r)
            seen.add(rid)
    obs = list(observations or [])
    if matched:
        state, reason = "ok", ""
    elif obs:
        state, reason = "distribution_only", "no_sourced_route_for_this_branch"
    else:
        # 没有记录不等于该支系历史上不存在（7.1）：这是分布视图的常态，不是"不可用"。
        state, reason = "distribution_only", "no_records_in_this_dataset"
    return {"history_state": state, "history_reason_code": reason,
            "routes": matched, "observations": obs,
            "tree_source": hist.get("tree_source", ""), "tree_version": hist.get("tree_version", ""),
            "sources": list(hist.get("sources") or [])}


# ---------------------------------------------------------------- CLI（可单独运行，AN4 要求）

def _cli(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="lineage history helper (AN4)")
    ap.add_argument("--history", default="", help="panel/lineage_history.json")
    ap.add_argument("--yard", default="", help="03_haplo (where y_result.json / mt_result.json live)")
    ap.add_argument("--out", required=True, help="where to write lineage_history.json")
    ap.add_argument("--sample", default="")
    ap.add_argument("--ytree", default="", help="YFull current_tree.json, for node-membership checks")
    ap.add_argument("--rows", default="", help="JSON with the AADR records (11_aadr/summary.json); "
                                              "their y_hg_raw/mt_hg_raw are what the observations are built from")
    a = ap.parse_args(argv)

    hist = load_history(json.loads(pathlib.Path(a.history).read_text(encoding="utf-8"))
                        if a.history and pathlib.Path(a.history).exists() else {})
    parents = hist["parents"]
    nodes = load_tree_nodes(a.ytree, "y") if a.ytree else None
    stats = {}
    if nodes:
        print(f"tree nodes loaded: {len(nodes)}")
    # 复审 AN4-P1：observations 此前恒为空列表（history_view(key, hist, [], parents)），于是报告里的
    # "已发表发现记录 / 迁移路线"永远是空占位——不是"没有记录"，而是**从来没查过**。这里读入 AADR 的
    # 记录（含 y_hg_raw/mt_hg_raw），交给 lineage_observations 去筛。查不到就如实报 0，并把原因带上。
    rows = []
    # 缺省指向 11_aadr 的结果：单独运行本步时也不必记住路径，而它正是观测的来源。
    _rows_path = a.rows or str(pathlib.Path(a.yard).parent / "11_aadr" / "summary.json")
    if _rows_path:
        _rp = pathlib.Path(_rows_path)
        if _rp.exists():
            try:
                _doc = json.loads(_rp.read_text(encoding="utf-8"))
                rows = _doc.get("records") or []
                print(f"observation source: {_rp.name}, {len(rows)} record(s)")
            except (json.JSONDecodeError, OSError) as e:
                print(f"WARNING: could not read {_rp}: {e}; observations will be empty", file=sys.stderr)
        else:
            print(f"WARNING: {_rp} not found; observations will be empty (has 09b run?)", file=sys.stderr)
    out = {"y": None, "mt": None}
    for kind in ("y", "mt"):
        p = pathlib.Path(a.yard) / f"{kind}_result.json"
        if not p.exists():
            out[kind] = lineage_summary(unavailable_lineage(kind, "missing_result", f"{p.name} not found"))
            continue
        try:
            raw = json.loads(p.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as e:
            out[kind] = lineage_summary(unavailable_lineage(kind, "unreadable_result", str(e)))
            continue
        s = lineage_summary(raw)
        node = s.get("terminal") or s.get("reported_hg")
        key = node if (node and ":" in str(node)) else (f"{kind}:{node}" if node else "")
        # 复审 AN4：树节点集只对 Y 有效（--ytree 是 YFull 树）；mt 没有 correspond 的树文件，
        # 传 Y 的节点集会把每个 mt 标签都误判 version_mismatch。mt 保持 unverified。
        _obs = lineage_observations(rows, key, parents, kind, history=hist,
                                    known_nodes=(nodes if kind == "y" else None), stats=stats) \
            if (rows and key) else []
        _hist = history_view(key, hist, _obs, parents)
        if stats:
            _hist["label_notes"] = dict(stats)
        out[kind] = dict(s, history=_hist)
    _out_path = pathlib.Path(a.out)
    _out_path.parent.mkdir(parents=True, exist_ok=True)
    _out_path.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    # 复审 H6/AN5 残留：30 此前按"文件存在"读 lineage_history.json，绕过任何状态。写一份
    # manifest（analysis_id=09d-lineage-history）让 30 走 analysis_state 准入；parameters 只记
    # 实际输入的 provenance，不用路径字符串冒充内容指纹。
    try:
        import ancestry_data as _adm
        _adm.write_manifest(_out_path.parent / "manifest.json", _adm.build_manifest(
            a.sample or "", "09d-lineage-history", state="ok",
            parameters={"history": a.history or "", "rows": _rows_path or "", "ytree": a.ytree or ""},
            outputs=[_out_path.name]))
    except Exception as _e:      # manifest 写失败不阻断结果文件，但要留下痕迹
        print(f"WARNING: could not write the 09d manifest: {_e}", file=sys.stderr)
    print(f"wrote {a.out} (y={out['y']['state']}, mt={out['mt']['state']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
