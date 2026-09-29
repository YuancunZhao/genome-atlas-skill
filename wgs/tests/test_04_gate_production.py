"""P1（复审 AN0/AN2/AN5）验收：04 只在显式配置 ref_superpop 时训练区域空间。

复审原话："04 仍无 REGIONAL_ENABLED 判断，未配置区域时仍按共享 SUPERPOP 默认训练。"本文件
对**生产 04_ancestry_pca.sh** 断言：假 plink2 替身按 --out 落盘常规产物，脚本在
REGIONAL_ENABLED=0（配置未给 ref_superpop）时只跑 global——不创建 regional.*/prune.regional
任何文件；显式给出 ref_superpop 时区域空间照常训练。替身写带 SuperPop/Population 列的
sscore 头，与真实 plink2 输出同构（服务器实测）。
"""
import json, os, pathlib, subprocess, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "04_ancestry_pca.sh"

# 假 plink2：--indep-pairwise 写剪枝集，--score 写 sscore，其余按 --out 补齐常用副产物。
_STUB = r"""#!/bin/bash
prev=""; out=""
for a in "$@"; do [ "$prev" = "--out" ] && out="$a"; prev="$a"; done
case "$*" in
  *--indep-pairwise*) printf 'rs1\nrs2\n' > "$out.prune.in"; : > "$out.prune.out"; exit 0;;
esac
case "$out" in
  *.proj) printf '#FID\tIID\tSuperPop\tPopulation\tPC1_AVG\tPC2_AVG\tPC3_AVG\tPC4_AVG\nTEST\tTEST\tEAS\tCHB\t0.10\t0.01\t0.00\t0.00\n' > "$out.sscore"; exit 0;;
  *) for e in pgen psam pvar afreq "eigenvec.allele" log; do : > "$out.$e"; done; exit 0;;
esac
"""


def _setup(td, **over):
    """临时工程 + 假 plink2 + 两行 psam（1 EAS / 1 EUR）。返回 (config, 04 输出目录)。"""
    work = pathlib.Path(td) / "work"
    bin_dir = pathlib.Path(td) / "bin"
    bin_dir.mkdir(parents=True)
    stub = bin_dir / "plink2"
    stub.write_text(_STUB, encoding="utf-8")
    stub.chmod(0o755)
    ref = pathlib.Path(td) / "panel" / "all_phase3"
    ref.parent.mkdir(parents=True)
    (ref.parent / "all_phase3.psam").write_text(
        "#FID\tIID\tSEX\tSuperPop\tPopulation\nFAM1\tHAN1\t1\tEAS\tCHB\nFAM2\tEUR1\t1\tEUR\tGBR\n",
        encoding="utf-8")
    cfg = {"sample_id": "TESTSAMPLE", "work_dir": str(work), "plink2": str(stub),
           "kg_pfile": str(ref), "threads": 1, "mem_gb": 1}
    cfg.update(over)
    p = pathlib.Path(td) / "config.yaml"
    p.write_text("".join(f"{k}: {json.dumps(v)}\n" for k, v in cfg.items()), encoding="utf-8")
    return p, work / "wgs" / "04_ancestry"


def _run_04(cfg):
    return subprocess.run(["bash", str(SCRIPT)], capture_output=True, text=True,
                         env={**os.environ, "WGS_CONFIG": str(cfg)}, timeout=120)


class TestRegionalGate(unittest.TestCase):

    def test_no_ref_superpop_runs_global_only(self):
        with tempfile.TemporaryDirectory() as td:
            cfg, out = _setup(td)                       # 未给 ref_superpop
            r = _run_04(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("global PCA only", r.stdout)
            self.assertTrue((out / "kg.proj.sscore").exists(), "global 空间照常训练")
            self.assertFalse(list(out.glob("regional.*")),
                             "未配置区域时不得产出任何 regional.* 文件")
            self.assertFalse((out / "prune.regional.prune.in").exists())
            self.assertFalse((out / "regional.ids").exists())

    def test_explicit_ref_superpop_trains_the_regional_space(self):
        with tempfile.TemporaryDirectory() as td:
            cfg, out = _setup(td, **{"ref_superpop": "EAS"})
            r = _run_04(cfg)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("EAS pruned SNPs", r.stdout)   # wc -l 输出有列宽填充，不断言数字对齐
            self.assertTrue((out / "regional.proj.sscore").exists())
            self.assertTrue((out / "prune.regional.prune.in").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
