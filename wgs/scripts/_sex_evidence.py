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
# 工作门槛（H4 复审补）：相对比可以被低深度的噪声推高——常染色体只有 2× 的文库里，
# Y 0.15× 也能算出 0.075 的"达标"比值。绝对深度不足的"有 Y"不断言，按证据不足处理。
Y_MIN_ABS_DEPTH = 1.0     # Y 平均深度的绝对下限，低于它的 has_y 降级为 ambiguous
# X 互证（H4 复审补）：Y 存在性不是 X 倍性的充分证据。有 Y 而 X 深度像两条、或无 Y 而
# X 深度像一条，都是性染色体证据自相矛盾——记录两边证据并标记 x_conflict，由 evidence_sex
# 回退声明值，绝不把 Y 一条证据当核型结论发布。
X_DIPLOID_MIN = 0.90      # X/常染色体 高于此值 → X 看起来是两条


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
    # 绝对深度下限：低深度文库的噪声足以把比值推过线（2× 文库的 Y 0.15× 也能算出 0.075 的
    # "达标"比值）。Y 本身太浅的"有 Y"撑不起任何下游结论，降级为 ambiguous 走声明回退。
    if inferred == "has_y" and y < Y_MIN_ABS_DEPTH:
        inferred = "ambiguous"
        why += (f"; Y mean depth {y:.3f} is below the absolute floor {Y_MIN_ABS_DEPTH}"
                f" (the ratio can be inflated by noise in a shallow library)")
    # X 互证：Y 存在性不是 X 倍性的充分证据。矛盾时保留两侧证据并标记 x_conflict，
    # 结论交给 evidence_sex 回退声明——这里不是核型判定器。
    x_conflict = (inferred == "has_y" and xr >= X_DIPLOID_MIN) or \
                 (inferred == "no_y" and xr <= X_SINGLE_MAX)
    if x_conflict:
        side = ("X/autosome %.4f looks diploid despite Y evidence" % xr) if inferred == "has_y" \
            else ("X/autosome %.4f looks haploid without Y evidence" % xr)
        why += f"; conflicting X evidence: {side}"
    agree = None
    if declared:
        d = "has_y" if str(declared).lower().startswith("m") else "no_y"
        agree = (d == inferred)
    return {"state": "ok", "reason_code": "", "declared": declared,
            "inferred": inferred, "agrees_with_declared": agree, "why": why,
            "x_conflict": x_conflict,
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


def evidence_sex(declared=None, summary=None, log=print):
    """本样本应按哪种性染色体构成处理，以及依据记录。

    有证据就用证据；证据不可用时回退到声明值，但把回退原因写进记录——回退是可以的（否则一个缺
    深度汇总的环境会连带停掉 Y 相关分析），沉默地回退不行。与声明冲突时告警一次。
    """
    rec = check(summary=summary, declared=declared)
    inf = rec.get("inferred")
    if rec.get("x_conflict"):
        # X 与 Y 证据自相矛盾（有 Y 而 X 像两条 / 无 Y 而 X 像一条）：单条染色体的读数不
        # 足以定倍性，回退声明值并说明。这挡住的是"拿 Y 一条证据当核型结论发布"。
        sex = "male" if str(declared if declared is not None else _c.SEX).lower().startswith("m") else "female"
        log(f"warning: X and Y evidence conflict ({rec.get('why')}); not deciding ploidy from a "
            f"single chromosome -- falling back to the declared sex ({sex})")
    elif inf == "has_y":
        sex = "male"
    elif inf == "no_y":
        sex = "female"
    else:
        sex = "male" if str(declared if declared is not None else _c.SEX).lower().startswith("m") else "female"
        if rec.get("state") == "ok":
            log(f"warning: sex-chromosome evidence is not decisive ({rec.get('why')}); "
                f"falling back to the declared sex ({sex})")
        else:
            log(f"warning: sex-chromosome evidence unavailable ({rec.get('reason_code')}); "
                f"falling back to the declared sex ({sex})")
    if rec.get("agrees_with_declared") is False:
        log(f"WARNING: the data disagrees with the declared sex: {rec.get('why')}; "
            f"downstream sex-specific results must be degraded")
    rec["used_sex"] = sex
    return sex, rec


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
    if rec.get("x_conflict"):
        print("  NOTE: X and Y evidence conflict; evidence_sex falls back to the declared sex "
              "rather than publish a ploidy call from one chromosome")
    if rec["agrees_with_declared"] is False:
        # 不自行改结论：把冲突说出来，由下游按适用性降级，并让报告能解释为什么。
        print("  WARNING: the data disagrees with the declared sex; downstream Y-based results must be "
              "degraded rather than reported as ordinary findings", file=sys.stderr)
