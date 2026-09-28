#!/usr/bin/env python
"""Infer the sex-chromosome constitution from the data, instead of assuming it from the configuration.

H4: "degrade by evidence and applicability, do not infer chromosome constitution from configuration". The
configuration says what the sample was declared as; it cannot say whether this CRAM actually has Y reads,
whether the X was sequenced deeply enough to see two copies, or whether the library lost the Y. A pipeline
that trusts `sex: male` and computes Y-based statistics on a female sample produces numbers that look
ordinary and mean nothing.

Evidence used here is a read-count ratio from the alignment index, which is cheap and needs no extra
passes: the ratio of Y-derived to autosomal-derived reads. Thresholds are wide on purpose -- this decides
whether a downstream analysis is *applicable*, not a clinical karyotype. Output is a record with the
evidence and the conclusion, so the report can say why something was skipped or degraded.
"""
import json, pathlib, subprocess, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- CRAM, SAMPLE, WGS, THREADS ...
import wgsconfig as _c

# 判定阈值。宽区间是刻意的：这里决定"下游分析适不适用"，不是核型诊断。
Y_PRESENT_MIN = 0.05      # Y_reads / autosome_reads 高于此值 → 有 Y 证据
Y_ABSENT_MAX = 0.02       # 低于此值 → 无 Y 证据
X_SINGLE_MAX = 0.75       # X 相对常染色体的读取比低于此值 → 单条 X 的可能


def _depth_means(summary=None):
    """每条染色体的平均深度，来自 01_qc 的 mosdepth summary。返回 (means, reason)。

    原先用 `samtools idxstats` 数读段：在这份 47 GB 的 CRAM 上它要全文件解码，跑满 100 秒仍未回来
    （CPU 99%），一个"判定分析适不适用"的小步骤不该有这种代价。深度汇总已经在 00/01 步算过，读它
    既快又更贴近问题本身——决定 X/Y 分析适不适用的是**测到多深**，不是比对索引里有多少条记录。
    """
    path = pathlib.Path(summary) if summary else pathlib.Path(_c.W) / "01_qc" / "depth.mosdepth.summary.txt"
    if not path.exists():
        return None, "no_depth_summary"
    means = {}
    try:
        with open(path, encoding="utf-8") as fh:
            head = fh.readline().split("\t")
            try:
                ci, mi = head.index("chrom"), head.index("mean")
            except ValueError:
                return None, "depth_summary_unexpected_columns"
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(ci, mi):
                    continue
                name = f[ci]
                if name.endswith("_region"):     # mosdepth 会为 regions 再写一份，丢掉
                    continue
                name = name[3:] if name.startswith("chr") else name
                try:
                    means[name] = float(f[mi])
                except ValueError:
                    continue
    except OSError as e:
        return None, f"depth_summary_unreadable:{type(e).__name__}"
    return (means, "") if means else (None, "depth_summary_empty")


def classify_depths(means, declared=None, reason=""):
    """由深度汇总判断性染色体构成。返回证据记录；判不了就说为什么判不了。"""
    if not means:
        return {"state": "unavailable", "reason_code": reason or "no_depth_summary",
                "declared": declared, "evidence": None}
    autos = [v for k, v in means.items() if k.isdigit() and 1 <= int(k) <= 22 and v > 0]
    if not autos:
        return {"state": "unavailable", "reason_code": "no_autosomal_depth",
                "declared": declared, "evidence": {"chromosomes": len(means)}}
    auto = sum(autos) / len(autos)
    y = means.get("Y")
    x = means.get("X")
    if y is None or x is None:
        return {"state": "unavailable", "reason_code": "no_sex_chromosomes_in_summary",
                "declared": declared, "evidence": {"autosome_mean": round(auto, 3)}}
    yr, xr = y / auto, x / auto
    # 男性 Y 的相对深度约 0.5（可比对区约占一半），女性接近 0；阈值取宽区间，只用于判定"适不适用"。
    if yr >= Y_PRESENT_MIN:
        inferred, why = "has_y", f"Y/autosome mean depth {yr:.4f} >= {Y_PRESENT_MIN}"
    elif yr <= Y_ABSENT_MAX:
        inferred, why = "no_y", f"Y/autosome mean depth {yr:.4f} <= {Y_ABSENT_MAX}"
    else:
        inferred, why = "ambiguous", f"Y/autosome mean depth {yr:.4f} lies between the thresholds"
    agree = None
    if declared:
        d = "has_y" if str(declared).lower().startswith("m") else "no_y"
        agree = (d == inferred)
    return {"state": "ok", "reason_code": "", "declared": declared,
            "inferred": inferred, "agrees_with_declared": agree, "why": why,
            "evidence": {"autosome_mean_depth": round(auto, 3), "x_mean_depth": round(x, 3),
                         "y_mean_depth": round(y, 3), "y_ratio": round(yr, 5), "x_ratio": round(xr, 5),
                         "chromosomes": len(means)}}


def classify(counts, declared=None, reason=""):
    """把 idxstats 的计数变成证据记录。缺 X/Y 或总数太低时说明原因，而不是给一个默认结论。"""
    if not counts:
        return {"state": "unavailable", "reason_code": reason or "no_index_stats",
                "declared": declared, "evidence": None}
    # 只按名字取，不假定 chrom 命名：任何以 Y/MT 结尾的键都可能出现（chrY / Y）
    def _get(name, alt):
        for k in (name, alt, "chr" + name):
            if k in counts:
                return counts[k]
        return 0
    auto = sum(v for k, v in counts.items()
               if k.replace("chr", "").isdigit() and 1 <= int(k.replace("chr", "")) <= 22)
    x = _get("X", "X")
    y = _get("Y", "Y")
    if auto <= 0:
        return {"state": "unavailable", "reason_code": "no_autosomal_reads",
                "declared": declared, "evidence": {"x": x, "y": y, "autosome": auto}}
    yr = y / auto
    xr = x / auto
    if yr >= Y_PRESENT_MIN:
        inferred, why = "has_y", f"Y/autosome read ratio {yr:.4f} >= {Y_PRESENT_MIN}"
    elif yr <= Y_ABSENT_MAX:
        inferred, why = "no_y", f"Y/autosome read ratio {yr:.4f} <= {Y_ABSENT_MAX}"
    else:
        inferred, why = "ambiguous", f"Y/autosome read ratio {yr:.4f} lies between the thresholds"
    agree = None
    if declared:
        d = "has_y" if str(declared).lower().startswith("m") else "no_y"
        agree = (d == inferred)
    return {"state": "ok", "reason_code": "", "declared": declared,
            "inferred": inferred, "agrees_with_declared": agree, "why": why,
            "evidence": {"x": x, "y": y, "autosome": auto,
                         "y_ratio": round(yr, 5), "x_ratio": round(xr, 5)}}


def check(summary=None, declared=None):
    """判定本样本的性染色体构成。cram 参数已移除：读深度汇总而不是重扫比对文件。"""
    declared = declared if declared is not None else _c.SEX
    means, reason = _depth_means(summary)
    out = classify_depths(means, declared, reason)
    out["sample_id"] = _c.SAMPLE
    out["evidence_source"] = str(summary) if summary else str(pathlib.Path(_c.W) / "01_qc" / "depth.mosdepth.summary.txt")
    return out


if __name__ == "__main__":
    rec = check()
    path = pathlib.Path(_c.W) / "sex_evidence.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {path}")
    if rec["state"] != "ok":
        print(f"  sex-chromosome evidence unavailable: {rec['reason_code']}")
        sys.exit(0)
    print(f"  declared={rec['declared']} inferred={rec.get('inferred')} "
          f"agrees={rec['agrees_with_declared']} · {rec['why']}")
    if rec["agrees_with_declared"] is False:
        # 不自行改结论：把冲突说出来，由下游按适用性降级，并让报告能解释为什么。
        print("  WARNING: the data disagrees with the declared sex; downstream Y-based results must be "
              "degraded rather than reported as ordinary findings", file=sys.stderr)
