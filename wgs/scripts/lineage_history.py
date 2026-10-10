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
* `observation_source` decides where the records come from and records whether the query happened at
  all: the full 08-normalised metadata first, the configured .anno when AADR PCA is off, the 09b
  summary (a screened subset) only as an explicit fallback. A query that never ran is reported as
  "history source not available" -- never as "the dataset has no records".

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

def _dedupe_by_master(rows, field):
    """同一个人的多种技术表示在历史视图里只算一次，按固定技术偏好挑一份（SG/DG/HO/AG/TW）。

    field 是该 kind 的标签字段（"y_hg_raw"/"mt_hg_raw"）。复审 §3.2 P1 AN1：挑选不再悄悄
    抹掉同 Master 的其它表示。没有该 kind 标签的表示不参与挑选（缺标签不是另一种支系说法），
    但保留在证据里；返回的每条记录带 `_representations`（该人**全部**原始表示：
    record_id/该 kind 原始标签/QC）、`_label_conflicts`（其它表示的非空标签与所选表示不同）、
    `_qc_conflicts`（其它表示的非空 QC 结论与所选不同，含所选无 QC 的情况）——同一人两条
    表示给出不同支系标签或质量结论时，证据必须跟着观测走，不能在排序里无声解决。
    """
    def rank(r):
        rep = str(r.get("genotype_representation") or r.get("record_id", "")).upper()
        for i, tok in enumerate(("SG", "DG", "HO", "AG", "TW")):
            if rep.endswith(tok) or rep == tok:
                return i
        return 99
    def label_of(r):
        return None if is_missing_hg(r.get(field)) else str(r.get(field))
    def qc_of(r):
        return str(r.get("hg_qc") or "").strip() or None
    groups = {}
    for r in rows:
        key = str(r.get("master_id") or r.get("individual_id") or r.get("record_id") or "")
        if not key:
            continue
        groups.setdefault(key, []).append(r)
    kept = []
    for key in sorted(groups):
        grp = groups[key]
        candidates = [r for r in grp if label_of(r) is not None]
        if not candidates:
            continue                      # 这个人没有该 kind 的标签：不生成观测
        best = min(candidates, key=rank)
        best = dict(best)                 # 证据字段不写回调用方的行（测试夹具是共享 dict）
        kept_label, kept_qc = label_of(best), qc_of(best)
        best["_representations"] = [{"record_id": str(r.get("record_id") or ""),
                                     "label": label_of(r), "qc": qc_of(r)} for r in grp]
        best["_label_conflicts"] = [
            {"record_id": str(r.get("record_id") or ""), "label": label_of(r),
             "kept_label": kept_label}
            for r in grp if label_of(r) is not None and label_of(r) != kept_label]
        best["_qc_conflicts"] = [
            {"record_id": str(r.get("record_id") or ""), "qc": qc_of(r), "kept_qc": kept_qc}
            for r in grp if qc_of(r) is not None and qc_of(r) != kept_qc]
        kept.append(best)
    return kept


def _record_source_version(record, kind):
    """记录标签的来源树版本，按 Y/mt 分开取（复审 §3.2 P0-2b）。

    AADR 只为 Y 标签发布来源（.anno 的 "based on Y-full 12.03" → hg_source_tree）；mt 标签没有
    对应的树版本字段。mt 的来源版本如实为空——绝不把 Y 的 YFull12.03 安到 mt 头上冒充来源。
    返回 (source_version, call_source)：前者进版本准入，后者只是"标签怎么来的"。
    """
    if kind == "y":
        return str(record.get("hg_source_tree") or ""), \
            str(record.get("hg_source_tree") or record.get("hg_call_source") or "")
    return "", str(record.get("mt_call_source") or record.get("hg_call_source") or "")


def _version_proven(source_version, tree_source, tree_version, maps):
    """记录来源版本与结果树版本是否**已证明**等价（复审 §3.2 P0-2b）。

    只有两种情况算证明：两边的版本字面一致（AADR 写 "YFull12.03"，结果写 tree_source="YFull"
    + tree_version="12.03"，拼起来正好相同），或面板里有一条**有来源、明确起止版本**的映射
    （version_maps 的 from/to/source 齐全）。同名节点存在本身不是版本等价——YFull 12.03 的
    N-CTS4714 与 14.06.0 的 N-CTS4714 只是同名。
    """
    s, v = str(source_version or "").strip(), str(tree_version or "").strip()
    if not s or not v:
        return False
    name = str(tree_source or "").strip()
    target = f"{name}{v}" if name else v
    if s == v or s == target:
        return True
    return any(str(m.get("from") or "") == s and str(m.get("to") or "") == target
               and str(m.get("source") or "").strip()
               for m in (maps or []))


def lineage_observations(rows, query, parents, tree_kind, same_tree=True, history=None,
                         known_nodes=None, stats=None, result_tree_source="",
                         result_tree_version="", strict=True):
    """筛出与该支系相关的历史记录（含祖先/后代/未定），并带上坐标、年代与来源。

    这里**不用**常染色体 call_rate 门槛：来源对 Y/mt 的可用性决定记录是否可用，覆盖率低只影响
    常染色体 PCA。没有相应单倍群标签的记录直接跳过——不猜。

    版本准入（复审 §3.2 P0-2b）：exact/父子关系只在记录来源版本与**结果**的
    tree_source/tree_version 已证明等价时成立（字面一致，或有来源、明确起止版本的 version_maps
    映射）。未证明时，与查询同名的记录进"待核对"（relation=pending_review，带原因），其余
    unresolved——同名节点存在不等于版本等价，宁可待核对也不冒充已证实的观测。
    strict=False（本样本判定有未解决冲突）时关系一律不作数：strict_match_allowed 要真正约束
    查询，不是摆设。
    """
    kind = str(tree_kind)
    field = "y_hg_raw" if kind == "y" else "mt_hg_raw"
    hist = history or {}
    maps = hist.get("version_maps") or []
    out = []
    # 复审 §3.2 P1 AN1：去重交给 _dedupe_by_master（带 field），同 Master 的其它表示连同
    # 标签/QC 冲突证据挂在所选表示上；没有该 kind 标签的表示不生成观测。
    for r in _dedupe_by_master(rows or [], field):
        node, note = canonicalize(r.get(field), kind, history, known_nodes)
        if stats is not None and note != "missing":
            stats[note] = stats.get(note, 0) + 1     # 进了比较就记账，包括版本不匹配的那些
        src_version, call_src = _record_source_version(r, kind)
        proven = _version_proven(src_version, result_tree_source, result_tree_version, maps)
        rel = match_lineage(query, node, parents, same_tree=bool(same_tree) and proven)
        pending_reason = ""
        if note == "version_mismatch":
            # 复审 AN4：标签不在当前树节点集（多半来自另一个树版本，如 AADR 的 YFull 12.03）时，
            # 与查询支系字符串相等不等于版本等价——未验证版本等价的记录不进默认视图，宁可
            # unresolved。
            rel = "unresolved"
        elif not proven:
            # 复审 §3.2 P0-2b：版本等价未证明。同名记录不是 exact，是"待核对"——留给人对着
            # 两个版本的树核一遍，再经面板 version_maps（带来源）收录。
            if node and node == str(query or ""):
                rel = "pending_review"
                got = src_version or "an unlabelled tree"
                want = (f"{result_tree_source}{result_tree_version}" if result_tree_source
                        else str(result_tree_version or "an unknown tree"))
                pending_reason = (f"label matches by name, but version equivalence between its "
                                  f"source tree ({got}) and this result's tree ({want}) "
                                  f"is not proven")
            else:
                rel = "unresolved"
        if (not strict) and rel in ("exact", "descendant", "ancestor"):
            # strict_match_allowed=False：查询支系自身有未解决冲突，严格视图不冒充任何关系
            if stats is not None:
                stats["strict_match_blocked"] = stats.get("strict_match_blocked", 0) + 1
            rel = "unresolved"
        if rel == "pending_review" and stats is not None:
            stats["pending_review"] = stats.get("pending_review", 0) + 1
        lat, lon = r.get("latitude"), r.get("longitude")
        out.append({
            "record_id": str(r.get("record_id") or ""),
            "node_id": node,
            "relation": rel,
            "pending_reason": pending_reason,
            "source_version": src_version,
            "locality": r.get("locality"),
            "coordinates": ({"latitude": lat, "longitude": lon} if (lat is not None and lon is not None) else None),
            "precision": str(r.get("location_precision") or "unknown"),
            "date_range": {"mean": r.get("date_mean_bp"), "min": r.get("date_min_bp"),
                           "max": r.get("date_max_bp")},
            "date_basis": r.get("date_basis"),
            # 复审 §3.2 P1 AN6：原始日期串一并透传——模板呈现"报告年代"时依据（date_basis）
            # 与原始写法（date_raw）都在，不只是换算后的数字。
            "date_raw": r.get("date_raw"),
            "call_source": call_src,
            "publication": r.get("publication"),
            # 复审 §3.2 P1 AN1：原始表示与同 Master 表示间的标签/QC 冲突跟着观测走（证据
            # 保留在交付数据里；是否展示由视图决定，绝不在去重时丢弃）。
            "representations": list(r.get("_representations") or []),
            "label_conflicts": list(r.get("_label_conflicts") or []),
            "qc_conflicts": list(r.get("_qc_conflicts") or []),
        })
    # 只保留能与查询支系建立关系的记录：exact 与后代默认显示，祖先作为背景层（视图决定怎么展开）。
    # 待核对（pending_review）保留在载荷里交给模板诚实展示，但不计入已发表观测。
    # unrelated 与 unresolved 不进默认视图——缺树边时无法确认，宁可说"未定"也不把别的支系算进来。
    out = [o for o in out if o["relation"] in ("exact", "descendant", "ancestor", "pending_review")]
    order = {"exact": 0, "descendant": 1, "ancestor": 2, "pending_review": 3}
    out.sort(key=lambda o: (order.get(o["relation"], 9), o["record_id"]))
    return out


def conservative_from_path(path, solid, tail):
    """没有复核记录时的通用保守落点。

    判据不是"某个节点的支持数够不够"，而是"它是否位于末端一条弱链的末尾"。本项目的真实路径是
    12 / 1 / 0 / 1 / 5：末端那个 5 恰好达到阈值，但它上面三级分别是 1、0、1——它们是同一段分辨率
    极限，5 只是这段弱链的末尾，不是独立证据。只看单节点阈值会把它当成可靠的末端（这正是第一版
    规则错的地方）。所以：末端 tail 级里若有支持不足的节点，就从最靠上的那个再往上退一级。

    solid/tail 是**每次运行自己的配置**（config 的 lineage_solid_min / lineage_tail_levels），
    不是函数默认值：此前硬编码的 5/4 是从这个样本的真实路径归纳出来的，把它当通用默认等于把
    一个样本的分辨率极限写进所有样本（AN4-a）。未配置（None）时本规则不生效——没有复核记录
    也没有配置阈值，就没有 conservative_hg，不硬造一个。

    path 每项至少含 (node, der, ...)。
    """
    if not path or solid is None or tail is None:
        return None
    seg = list(path[-tail:]) if tail and tail > 0 else list(path)
    weak = [i for i, p in enumerate(seg) if int(p[1]) < int(solid)]
    if not weak:
        return str(path[-1][0])
    first_weak = seg[weak[0]]
    pos = list(path).index(first_weak)
    return str(path[pos - 1][0]) if pos > 0 else None


def load_review(path):
    """读本样本的人工复核文件（缺省 03_haplo/lineage_review.json，或 config lineage_review_file）。

    复核记录是**样本私有**的（它核对的是这个样本的末端 SNP 表和当时的树版本），所以它住在 work
    目录、不进共享 panel——panel/lineage_history.json 的 reviewed_calls 只保留空表。文件可以不
    存在：没有复核就没有 conservative_hg 的复核来源，仅此而已，不报错不阻断。
    """
    try:
        doc = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"reviewed_calls": {}}
    recs = {str(k): dict(v) for k, v in ((doc or {}).get("reviewed_calls") or {}).items()
            if isinstance(v, dict)}
    return {"reviewed_calls": recs}


def reviewed_call(review, kind, reported_hg, sample_id=None, tree_version=None, evidence_sha=None,
                  log=None):
    """带理由的复核记录（reviewed_calls）优先于任何自动规则——但只在它核对过的东西与本次运行一致时。

    §7 要求 conservative_hg 来自实际证据或带理由的复核记录，而不是"树里最深节点"或某个阈值。
    键形如 `y:N-CTS4714`。

    **样本限定（复审 AN4）**：键只写支系，所以一个样本的人工复核结论会被**任何**调出同一支系的
    新样本继承——那是把一个样本的判断当成通用规则。因此记录里必须写明 sample_id，且只在**当前样本
    与之相符**时才生效；没有 sample_id 的旧记录一律不套用（并说明原因），由调用方回退到通用规则。

    **证据/树版本绑定（AN4-a）**：人工复核是对着**这个样本的末端 SNP 表**（evidence_sha256）和
    **当时的项目树版本**（tree_version）做的。记录必须带这两个绑定字段，且调用方必须传本次运行的
    对应值；绑定缺失或对不上都不套用——对着旧证据/旧树做的结论不是本次运行的结论，宁可没有
    conservative_hg。
    """
    recs = (review or {}).get("reviewed_calls") or {}
    node = str(reported_hg or "")
    if not node:
        return None, None
    key = node if ":" in node else f"{kind}:{node}"
    rec = recs.get(key)
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
    rtv = str(rec.get("tree_version") or "")
    if not tree_version or not rtv or rtv != str(tree_version):
        if log:
            log(f"reviewed call for {key} was recorded under tree {rtv or '<unrecorded>'}, "
                f"this run classifies on {tree_version or '<unknown>'}; not applied")
        return None, None
    rsha = str(rec.get("evidence_sha256") or "")
    if not evidence_sha or not rsha or rsha != str(evidence_sha):
        if log:
            log(f"reviewed call for {key} is bound to different caller evidence "
                f"(recorded ...{rsha[-8:] if rsha else 'none'}, this run ...{str(evidence_sha)[-8:] if evidence_sha else 'none'}); not applied")
        return None, None
    return str(rec["conservative_hg"]), str(rec.get("reason") or "reviewed by hand")


# ---------------------------------------------------------------- 观测来源（AN4：完整元数据独立查询）

# reference_metadata.tsv 是 pandas 写的：数值列在 TSV 里都成了字符串。历史视图要拿这些字段做
# 坐标/年代（模板还要做投影运算），在这里转回数值；转不动的按缺失处理，不猜。
_NUMERIC_FIELDS = ("latitude", "longitude", "date_mean_bp", "date_min_bp", "date_max_bp")


def _coerce_numeric(row):
    for k in _NUMERIC_FIELDS:
        v = row.get(k)
        if v is None or str(v).strip() == "":
            row[k] = None
            continue
        try:
            f = float(v)
            row[k] = int(f) if f.is_integer() else f
        except (TypeError, ValueError):
            row[k] = None
    return row


def read_metadata_tsv(path):
    """读 08/09b 规范化后写出的 reference_metadata.tsv（**完整**元数据，未过 08/09 的 PCA 筛选）。"""
    import csv
    with open(path, newline="", encoding="utf-8") as fh:
        return [_coerce_numeric(dict(r)) for r in csv.DictReader(fh, delimiter="\t")]


def read_anno_records(path):
    """AADR PCA 关闭（08/09b 没跑）时直接读已配置的 .anno——历史查询不随 PCA 禁用而失效。

    复用 ancestry_data.normalize_metadata：08/09b 用同一函数写 reference_metadata.tsv，这里用
    同一函数读 .anno，字段名不会漂。.anno 缺必需列时 normalize_metadata 抛 ValueError，由
    调用方按"来源不可读"处理（如实报错，不静默变成空数据集）。
    """
    import csv
    import ancestry_data as ad
    with open(path, newline="", encoding="utf-8") as fh:
        raw = [{(k or "").strip(): v for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]
    return ad.normalize_metadata(raw, dataset="AADR", release="from-anno")


def _metadata_annotation_sha(metadata_path):
    """08 在 metadata 同目录 manifest.json 里记录的注释内容指纹；读不到/没记录 = 无法证明。"""
    import ancestry_data as ad
    m = ad.read_manifest(pathlib.Path(metadata_path).parent / "manifest.json")
    if not m or m.get("analysis_id") != "08-aadr-extract" or m.get("state") != "ok":
        return ""
    return str((m.get("parameters") or {}).get("annotation_sha") or "")


def observation_source(metadata="", anno="", rows="", yard=""):
    """决定观测记录从哪里来，并如实记录"到底查没查、查的是什么"。

    优先级：--metadata（08 的完整元数据）> --anno（已配置注释文件，AADR 禁用时也能查）>
    --rows（09b 的 summary——那一份经过 08/09 的 PCA 筛选，只剩面板里的记录；保留作显式回退，
    并在 query.note 里说明它是筛选子集）。都没给时用缺省位置：yard 旁的 11_aadr/
    reference_metadata.tsv，再退 11_aadr/summary.json。**都没有 = 没查过**，state=not_available：
    没查过不等于数据集没有记录，这个状态会一路带到 history_state，报告按"历史资料不可用"呈现。

    缓存身份（复审 §3.2 P1 AN1）：metadata 层是缓存，**不是无条件优先**。同时配置了 .anno 时，
    先核对 08 manifest 记录的 annotation_sha 与当前注释的内容指纹：对不上（或没记录——修复前
    的旧 manifest）就证明不了这份缓存出自当前配置的注释，直接规范化当前 .anno，note 说明原因。
    指纹对上（注释换回同一份，A→B→A 的最后一跳）才用缓存，省一次全量规范化。

    返回 (rows, query)，query = {state: ok|unreadable|not_available, source_kind, source, n_rows, note}。
    """
    def _q(state, kind, path, rows_list, note=""):
        return rows_list, {"state": state, "source_kind": kind, "source": str(path or ""),
                           "n_rows": len(rows_list), "note": str(note)}

    _anno_p = str(anno or "")
    _anno_ok = bool(_anno_p) and pathlib.Path(_anno_p).is_file()

    def _metadata_tier(p):
        """metadata 层：身份可证（或无注释可比对）才读；否则改查当前 .anno。"""
        import ancestry_data as _ad
        if _anno_ok:
            recorded, current = _metadata_annotation_sha(p), _ad.file_sha(_anno_p)
            if recorded != current:
                note = (f"existing metadata not bound to the current annotation "
                        f"(recorded sha {recorded or 'none'} != {current}); "
                        f"normalised the configured .anno directly")
                try:
                    return _q("ok", "anno", _anno_p, read_anno_records(_anno_p), note)
                except (OSError, ValueError) as e:
                    return _q("unreadable", "anno", _anno_p, [], f"{type(e).__name__}: {e}")
        try:
            return _q("ok", "reference_metadata", p, read_metadata_tsv(p))
        except (OSError, ValueError) as e:
            return _q("unreadable", "reference_metadata", p, [], f"{type(e).__name__}: {e}")

    _SCREENED = "09b summary rows: screened by the 08/09 PCA filters, not the full dataset"
    for path, kind, note in ((metadata, "reference_metadata", ""), (anno, "anno", "")):
        p = str(path or "")
        if not p or not pathlib.Path(p).exists():
            continue
        if kind == "reference_metadata":
            return _metadata_tier(p)
        try:
            return _q("ok", kind, p, read_anno_records(p), note)
        except (OSError, ValueError) as e:
            return _q("unreadable", kind, p, [], f"{type(e).__name__}: {e}")
    if rows:
        p = str(rows)
        if pathlib.Path(p).exists():
            try:
                doc = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
                return _q("ok", "summary", p, list(doc.get("records") or []), _SCREENED)
            except (OSError, json.JSONDecodeError) as e:
                return _q("unreadable", "summary", p, [], f"{type(e).__name__}: {e}")
    _meta_default = str(pathlib.Path(yard).parent / "11_aadr" / "reference_metadata.tsv")
    if pathlib.Path(_meta_default).exists():
        return _metadata_tier(_meta_default)
    _sum_default = str(pathlib.Path(yard).parent / "11_aadr" / "summary.json")
    if pathlib.Path(_sum_default).exists():
        try:
            doc = json.loads(pathlib.Path(_sum_default).read_text(encoding="utf-8"))
            return _q("ok", "summary", _sum_default, list(doc.get("records") or []), _SCREENED)
        except (OSError, json.JSONDecodeError) as e:
            return _q("unreadable", "summary", _sum_default, [], f"{type(e).__name__}: {e}")
    return _q("not_available", "none", "", [],
              "no reference_metadata.tsv / configured .anno / 09b summary was available to query")


HISTORY_STATES = ("ok", "distribution_only", "unavailable")


def load_history(doc):
    """读 panel/lineage_history.json：版本、来源、树边/别名、版本映射与经核对的路线。"""
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
    # 版本映射（复审 §3.2 P0-2b）：把另一个树版本的标签等价为当前版本。from/to/source 三者
    # 齐全才算证据——没有来源的映射不是映射，是愿望。
    maps = []
    for m in (d.get("version_maps") or []):
        if not isinstance(m, dict):
            continue
        if not str(m.get("from") or "").strip() or not str(m.get("to") or "").strip():
            continue
        if not str(m.get("source") or "").strip():
            continue
        maps.append({"from": str(m["from"]), "to": str(m["to"]), "source": str(m["source"])})
    return {"schema_version": int(d.get("schema_version") or 1),
            "tree_source": str(d.get("tree_source") or ""),
            "tree_version": str(d.get("tree_version") or ""),
            "parents": parents,
            "aliases": {str(k): str(v) for k, v in (d.get("aliases") or {}).items()},
            "version_maps": maps,
            "reviewed_calls": {str(k): dict(v) for k, v in (d.get("reviewed_calls") or {}).items()
                               if isinstance(v, dict)},
            "routes": routes,
            "sources": list(d.get("sources") or [])}


def history_view(query, history, observations, parents, same_tree=True, source=None):
    """历史视图：观测 + 与查询支系匹配的路线；没有路线时明说是分布视图，不造故事。

    `source` 是 observation_source 的查询记录。**没查过不等于数据集没有记录**（7.1）：来源缺失/
    不可读时视图整体 unavailable 并带原因；只有真正查过而一无所获，才可以说 no_records_in_this_dataset。
    """
    hist = history or {}
    base = {"tree_source": hist.get("tree_source", ""), "tree_version": hist.get("tree_version", ""),
            "sources": list(hist.get("sources") or []), "query": (dict(source) if source else None)}
    if source is not None and str(source.get("state")) != "ok":
        # 查询从未发生（not_available）或读了但读不出来（unreadable）：历史资料不可用，而不是空数据集。
        rc = "history_source_unreadable" if source.get("state") == "unreadable" else "history_source_not_available"
        return dict(base, history_state="unavailable", history_reason_code=rc, routes=[], observations=[])
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
        # 查过了、该支系一条相关记录都没有：这是分布视图的常态（7.1），不是"不可用"。
        state, reason = "distribution_only", "no_records_in_this_dataset"
    return dict(base, history_state=state, history_reason_code=reason, routes=matched, observations=obs)


# ---------------------------------------------------------------- CLI（可单独运行，AN4 要求）

def _cli(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="lineage history helper (AN4)")
    ap.add_argument("--history", default="", help="panel/lineage_history.json（LINEAGE_HISTORY_FILE）")
    ap.add_argument("--yard", default="", help="03_haplo (where y_result.json / mt_result.json live)")
    ap.add_argument("--out", required=True, help="where to write lineage_history.json")
    ap.add_argument("--sample", default="")
    ap.add_argument("--ytree", default="", help="YFull current_tree.json, for node-membership checks")
    ap.add_argument("--metadata", default="",
                    help="11_aadr/reference_metadata.tsv (08's FULL normalised metadata); preferred source")
    ap.add_argument("--anno", default="",
                    help="configured AADR .anno, normalised in-place when AADR PCA (08/09b) never ran")
    ap.add_argument("--rows", default="",
                    help="09b summary.json (a SCREENED subset of the metadata); explicit fallback only")
    a = ap.parse_args(argv)

    # 失效先于计算（复审 §3.2 P0 失败生命周期）：--history JSON 损坏、树读取失败、查询/合并/
    # 写文件中途崩溃时，旧 ok manifest 与旧 lineage_history.json 不得继续充当本次结果。写不出
    # 失效记录就不开工——否则旧 ok 恰好在最需要失效的路径上存活。
    import ancestry_data as _adm
    _adm.begin_run_manifest(pathlib.Path(a.out).parent / "manifest.json",
                            a.sample or "", "09d-lineage-history")

    hist = load_history(json.loads(pathlib.Path(a.history).read_text(encoding="utf-8"))
                        if a.history and pathlib.Path(a.history).exists() else {})
    parents = hist["parents"]
    nodes = load_tree_nodes(a.ytree, "y") if a.ytree else None
    stats = {}
    if nodes:
        print(f"tree nodes loaded: {len(nodes)}")
    # 复审 AN4-P1：observations 此前恒为空列表，于是报告里的"已发表发现记录 / 迁移路线"永远是空
    # 占位——不是"没有记录"，而是**从来没查过**。AN4 查询扩展：观测默认来自 08 的完整元数据
    # （reference_metadata.tsv，未过 08/09 的 PCA 筛选），AADR 禁用时退到已配置 .anno；09b 的
    # summary 只是显式回退（它是筛选子集）。没查过/读不出来都记入 query 状态，报告按"历史资料
    # 不可用"呈现——绝不写成"数据集没有记录"。
    rows, _query = observation_source(metadata=a.metadata, anno=a.anno, rows=a.rows, yard=a.yard)
    if _query["state"] == "ok":
        print(f"observation source: {_query['source_kind']} ({_query['source']}), "
              f"{_query['n_rows']} record(s)" + (f"; note: {_query['note']}" if _query["note"] else ""))
    else:
        print(f"WARNING: observation source {_query['state']} ({_query['note']}); "
              "history views will be delivered as unavailable, not as an empty dataset", file=sys.stderr)
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
        # 复审 §3.2 P0-2b：把结果的 tree_source/tree_version 与 strict_match_allowed 一并传进
        # 观测比较——exact/父子关系只在版本等价已证明时成立；本样本判定有冲突时严格视图不作数。
        _obs = lineage_observations(rows, key, parents, kind, history=hist,
                                    known_nodes=(nodes if kind == "y" else None), stats=stats,
                                    result_tree_source=s.get("tree_source") or "",
                                    result_tree_version=s.get("tree_version") or "",
                                    strict=bool(s.get("strict_match_allowed", True))) \
            if (rows and key) else []
        _hist = history_view(key, hist, _obs, parents, source=_query)
        if stats:
            _hist["label_notes"] = dict(stats)
        out[kind] = dict(s, history=_hist)
    _out_path = pathlib.Path(a.out)
    _out_path.parent.mkdir(parents=True, exist_ok=True)
    _out_path.write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    # 复审 H6/AN5 残留：30 此前按"文件存在"读 lineage_history.json，绕过任何状态。写一份
    # manifest（analysis_id=09d-lineage-history）让 30 走 analysis_state 准入；parameters 记
    # 观测来源的实际 provenance（层级/行数/查询状态），不用路径字符串冒充内容指纹。
    try:
        _adm.write_manifest(_out_path.parent / "manifest.json", _adm.build_manifest(
            a.sample or "", "09d-lineage-history", state="ok",
            parameters={
                # 复审 §3.2 P0 指纹绑定：路径字符串不充当指纹——记 history/ytree 的**内容**
                # sha（缺失/传目录时为 ""，30 比对侧同一算法）。观测来源记 provenance（层级/
                # 行数/查询状态），30 不把它们当比对键。
                "history_sha": _adm.file_sha(a.history),
                "ytree_sha": _adm.file_sha(a.ytree),
                "source_kind": _query["source_kind"], "source": _query["source"],
                "n_rows": _query["n_rows"], "query_state": _query["state"]},
            outputs=[_out_path.name]))
    except Exception as _e:      # manifest 写失败不阻断结果文件，但要留下痕迹
        print(f"WARNING: could not write the 09d manifest: {_e}", file=sys.stderr)
    print(f"wrote {a.out} (y={out['y']['state']}, mt={out['mt']['state']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
