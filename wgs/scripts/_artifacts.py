#!/usr/bin/env python
"""Several products in this pipeline are consumed but not produced by any script in the repository (H5).

Their producers were run by hand on the server, so the chain cannot be re-run end to end and a missing
product is discovered late -- as an empty table or a quietly degraded result rather than as a failure.
This module names each such product, what should produce it, and which steps consume it, and can check
that the ones a given step needs are actually present.

It deliberately does not invent the producer commands. Where a producer is unknown it says so, because a
guessed Delly/ExpansionHunter invocation that differs from the one used for the delivered report would be
worse than an honest gap: the numbers would change without anyone being told.
"""
import json, pathlib, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import wgsconfig as _c

W = pathlib.Path(_c.W)

# product -> (relative path, producer, consumers, note)
ARTIFACTS = {
    "delly_sv": {
        "path": "08_sv/delly",
        "producer": None,
        "produced_by_hand": True,
        "consumers": ["13_sv_filter"],
        "note": "Delly per-sample SV calls; 13 filters them by depth agreement",
    },
    "expansionhunter": {
        "path": "08_sv/eh",
        "producer": None,
        "produced_by_hand": True,
        "consumers": ["13_sv_filter", "30_build_report_data"],
        "note": "repeat expansions (incl. STR loci reported in section 05)",
    },
    "smn": {
        "path": "08_sv/smn",
        "producer": None,
        "produced_by_hand": True,
        "consumers": ["30_build_report_data"],
        "note": "SMN1/SMN2 copy number; the M1 report conclusion depends on it",
    },
    "roh": {
        "path": "09_misc/plink_roh",
        "producer": None,
        "produced_by_hand": True,
        "consumers": ["30_build_report_data"],
        "note": "runs of homozygosity",
    },
    "cyrius": {
        "path": "06_pgx/cyrius",
        "producer": None,
        "produced_by_hand": True,
        "consumers": ["11_pgx_extra"],
        "note": "CYP2D6 star alleles (PharmCAT does not resolve 2D6 on its own)",
    },
    "t1k": {
        "path": "06_pgx/t1k",
        "producer": None,
        "produced_by_hand": True,
        "consumers": ["11_pgx_extra"],
        "note": "HLA/KIR typing used by the immunogenetics section",
    },
    "pharmcat_summary": {
        "path": "06_pgx/pharmcat",
        "producer": None,
        "produced_by_hand": True,
        "consumers": ["11_pgx_extra", "30_build_report_data"],
        "note": "PharmCAT report and diplotypes",
    },
}


def check(keys=None, root=None, verbose=False):
    """Return (ok, report). Missing products are reported, not treated as empty results."""
    root = pathlib.Path(root) if root else W
    keys = list(keys) if keys else list(ARTIFACTS)
    rows, missing = [], []
    for k in keys:
        spec = ARTIFACTS.get(k)
        if not spec:
            rows.append({"key": k, "state": "unknown_product"})
            missing.append(k)
            continue
        p = root / spec["path"]
        exists = p.exists() and any(p.iterdir()) if p.is_dir() else p.exists()
        row = {"key": k, "path": str(p), "state": "ok" if exists else "missing",
               "produced_by_hand": spec["produced_by_hand"], "note": spec["note"]}
        if not exists:
            missing.append(k)
        rows.append(row)
        if verbose:
            print(f"  [{'ok' if exists else 'MISSING'}] {k:18s} {spec['path']}"
                  f"{'  (no producer in the repository)' if spec['produced_by_hand'] else ''}")
    return (not missing), {"checked": len(keys), "missing": missing, "rows": rows}


def require(keys, root=None):
    """Assert the named products exist; raise with a message that names the consequence and the gap."""
    ok, rep = check(keys, root=root)
    if not ok:
        lines = [f"missing product(s): {', '.join(rep['missing'])}"]
        for r in rep["rows"]:
            if r["state"] != "ok":
                spec = ARTIFACTS.get(r["key"], {})
                lines.append(f"  - {r['key']} expected at {spec.get('path')}: {spec.get('note', '')}")
                if spec.get("produced_by_hand"):
                    lines.append("    no script in this repository produces it (run by hand on the server)")
        raise FileNotFoundError("\n".join(lines))
    return rep


if __name__ == "__main__":
    keys = sys.argv[1:] or None
    ok, rep = check(keys, verbose=True)
    out = pathlib.Path(_c.W) / "artifacts.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    if not ok:
        print(f"  {len(rep['missing'])} product(s) missing: {' '.join(rep['missing'])}")
        sys.exit(1)
