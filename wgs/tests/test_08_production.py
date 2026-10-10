"""AN1 验收：对着**生产 08_aadr_extract.py** 子进程断言，不对着复制品。

复审指出：bcftools 非零退出只打印后继续、TGENO 只断言总尺寸——一个先吐合法 stdout 再
exit 1 的 bcftools 替身、或一个等长损坏头的 .geno，都能让 08 照常写 manifest.state=ok。
这里用真实脚本 + 合成面板 + bcftools 替身复现这两条，断言失败落 failed manifest、非零退
出、不产半份面板；附 happy-path 对照证明替身真的走通了生产路径（09b 的生产测试不执行
08，覆盖不到提取链）。validate_tgeno_header 另有纯函数级断言——头格式取自服务器实测的
真实 v66 文件：b"TGENO   27594  584131 <hex> <hex>" + NUL 填充到 48 字节。

生产脚本 import numpy/pandas，config 里的 aadr_modern 需要真正的列表（yaml 解析）——
缺任一依赖时子进程用例整文件跳过（本地套件跳过、服务器套件执行，与 test_09b_production
同一约定）；纯函数用例无此限制。
"""
import importlib.util, json, os, pathlib, re, subprocess, sys, tempfile, unittest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "08_aadr_extract.py"
_SKIP = ""
for _mod in ("numpy", "pandas", "yaml"):
    if importlib.util.find_spec(_mod) is None:
        _SKIP = f"{_mod} not available (the production script or its config needs it)"
        break

# 合成面板：两个 Han 现代个体 + 一个未配置个体；6 个 1 号染色体双等位位点。
N_SNP = 6
INDS = [("HAN1", "M", "Han"), ("HAN2", "F", "Han"), ("UNK1", "F", "Other")]
SITES = [(f"rs{i + 1}", "1", str(1000 + i), "A", "C") for i in range(N_SNP)]


def _write_config(td):
    cfg = {
        "sample_id": "TESTSAMPLE",
        "work_dir": str(pathlib.Path(td) / "work"),
        "aadr_prefix": str(pathlib.Path(td) / "work/data/ref/aadr/panel"),
        "aadr_modern": ["Han"],
    }
    p = pathlib.Path(td) / "config.yaml"
    p.write_text("".join(f"{k}: {json.dumps(v)}\n" for k, v in cfg.items()), encoding="utf-8")
    return p


def _write_panel(td, header=None, n_rows=None):
    """合成 .ind/.snp/.anno/.geno。header=None 用合法 TGENO 头；n_rows=None 与 .ind 行数一致
    （头-尺寸自洽）。注入损坏用这两个参数。"""
    pref = pathlib.Path(td) / "work/data/ref/aadr/panel"
    pref.parent.mkdir(parents=True, exist_ok=True)
    (pref.parent / "panel.ind").write_text(
        "".join(f"{i}\t{s}\t{p}\n" for i, s, p in INDS), encoding="utf-8")
    (pref.parent / "panel.snp").write_text(
        "".join(f"{r}\t{c}\t0.0\t{pos}\t{a1}\t{a2}\n" for r, c, pos, a1, a2 in SITES), encoding="utf-8")
    (pref.parent / "panel.anno").write_text(
        "Genetic ID\tPersistent Genetic ID\tGroup ID\tDate mean in BP\n"
        + "".join(f"{iid}\t{iid}\t{grp}\t\n" for iid, _, grp in INDS), encoding="utf-8")
    rlen = (N_SNP + 3) // 4
    n_rows = len(INDS) if n_rows is None else n_rows
    if header is None:
        header = f"TGENO   {len(INDS)}  {N_SNP} 00000000 00000000"
    head = header.encode()[:48].ljust(48, b"\x00")
    (pref.parent / "panel.geno").write_bytes(head + bytes(rlen * n_rows))
    return pref.parent / "panel.geno"


def _write_bcftools(td, mode):
    """bcftools 替身：ok=按 -R 位点表出全量行；fail-exit1=先吐一条完全合法的行再失败；
    empty=exit 0 但一行不出。"""
    stub = pathlib.Path(td) / "bin" / "bcftools"
    stub.parent.mkdir(parents=True, exist_ok=True)
    stub.write_text("#!/bin/bash\n"
                    'sites=""; prev=""\n'
                    'for a in "$@"; do [ "$prev" = "-R" ] && sites="$a"; prev="$a"; done\n'
                    f'mode="{mode}"\n'
                    'if [ "$mode" = "fail-exit1" ]; then\n'
                    '  printf \'1\\t1000\\tA\\tC\\t0|1\\n\'\n'
                    '  echo "stub bcftools: simulated failure after valid output" >&2\n'
                    "  exit 1\n"
                    "fi\n"
                    'if [ "$mode" = "empty" ]; then exit 0; fi\n'
                    'awk \'BEGIN{OFS="\\t"} {print $1, $2, "A", "C", "0|1"}\' "$sites"\n'
                    "exit 0\n", encoding="utf-8")
    stub.chmod(0o755)
    return stub.parent


def _run_08(cfg, bin_dir):
    return subprocess.run(
        [sys.executable, str(SCRIPT)], capture_output=True, text=True,
        env={**os.environ, "WGS_CONFIG": str(cfg),
             "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"},
        timeout=180)


def _manifest(td):
    return json.loads(
        (pathlib.Path(td) / "work/wgs/11_aadr/manifest.json").read_text(encoding="utf-8"))


def _seed_stale_ok(td):
    """上一次运行留下的 ok manifest——本次失败必须原子覆盖它，不能让旧成功状态可误读。"""
    w = pathlib.Path(td) / "work/wgs/11_aadr"
    w.mkdir(parents=True, exist_ok=True)
    (w / "manifest.json").write_text(
        json.dumps({"schema_version": "1", "sample_id": "TESTSAMPLE",
                    "analysis_id": "08-aadr-extract", "state": "ok", "stale": True}),
        encoding="utf-8")


@unittest.skipIf(_SKIP, _SKIP)
class TestBcftoolsFailureIsAFailure(unittest.TestCase):
    """复审 AN1：bcftools 非零退出只打印后继续，08 仍写 manifest.state=ok。"""

    def test_valid_stdout_then_exit1_fails_and_invalidates_stale_ok(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_panel(td)
            bin_dir = _write_bcftools(td, "fail-exit1")
            _seed_stale_ok(td)
            r = _run_08(cfg, bin_dir)
            self.assertNotEqual(r.returncode, 0, "先吐合法行再 exit 1 的 bcftools 必须让 08 失败")
            self.assertIn("bcftools_query_failed", r.stderr)
            m = _manifest(td)
            self.assertEqual(m["state"], "failed")
            self.assertEqual(m["reason_code"], "bcftools_query_failed")
            self.assertEqual(m["sample_id"], "TESTSAMPLE")
            self.assertFalse((pathlib.Path(td) / "work/wgs/11_aadr/aadr.bed").exists(),
                             "失败不得产出半份面板")
            self.assertNotIn("stale", m, "失败 manifest 必须覆盖上次的 ok，不是并列留下")

    def test_exit0_but_no_rows_is_a_failure(self):
        """exit 0 + 空 stdout：输入坏了（错的 VCF/build），不是“目标恰好全缺失”。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_panel(td)
            bin_dir = _write_bcftools(td, "empty")
            r = _run_08(cfg, bin_dir)
            self.assertNotEqual(r.returncode, 0)
            m = _manifest(td)
            self.assertEqual(m["state"], "failed")
            self.assertEqual(m["reason_code"], "bcftools_query_empty")

    def test_happy_path_writes_ok_manifest_and_panel(self):
        """对照：替身按位点表出全量行时，真实 08 走通整条链（头校验→提取→bed/manifest）。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_panel(td)
            bin_dir = _write_bcftools(td, "ok")
            r = _run_08(cfg, bin_dir)
            self.assertEqual(r.returncode, 0, r.stderr)
            m = _manifest(td)
            self.assertEqual(m["state"], "ok")
            fam = (pathlib.Path(td) / "work/wgs/11_aadr/aadr.fam").read_text(encoding="utf-8")
            iids = [line.split()[1] for line in fam.strip().splitlines()]
            self.assertEqual(iids, ["HAN1", "HAN2", "TESTSAMPLE"], "面板 2 人 + 稳定 sample_id 目标")
            samples = (pathlib.Path(td) / "work/wgs/11_aadr/samples.tsv").read_text(encoding="utf-8")
            self.assertIn("TESTSAMPLE\tTESTSAMPLE\ttarget", samples)


@unittest.skipIf(_SKIP, _SKIP)
class TestTgenoHeaderValidated(unittest.TestCase):
    """复审 AN1：TGENO 只断言总尺寸——等长损坏头、与 .ind/.snp 不符的头维度都放行。"""

    def test_corrupt_equal_length_header_fails(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_panel(td, header="XGENO   2  6 00000000 00000000")  # 尺寸完全一致，魔数损坏
            bin_dir = _write_bcftools(td, "ok")
            r = _run_08(cfg, bin_dir)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("tgeno_header_invalid", r.stderr)
            self.assertEqual(_manifest(td)["reason_code"], "tgeno_header_invalid")

    def test_header_dimensions_must_match_ind_rows(self):
        """头自称 5 人、文件 3 行（与 .ind 一致，尺寸断言通过）：只有头维度校验抓得到——
        这是纯尺寸断言的盲区。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_panel(td, header="TGENO   5  6 00000000 00000000")
            bin_dir = _write_bcftools(td, "ok")
            r = _run_08(cfg, bin_dir)
            self.assertNotEqual(r.returncode, 0)
            m = _manifest(td)
            self.assertEqual(m["reason_code"], "tgeno_header_invalid")
            self.assertIn("header n_ind 5 != .ind rows 3", m["parameters"]["detail"])

    def test_truncated_geno_fails_cleanly(self):
        """截断一行的 .geno：旧 assert 裸栈崩掉且不留状态；现在带原因落 failed manifest。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            geno = _write_panel(td)
            rlen = (N_SNP + 3) // 4
            geno.write_bytes(geno.read_bytes()[: -rlen])
            bin_dir = _write_bcftools(td, "ok")
            r = _run_08(cfg, bin_dir)
            self.assertNotEqual(r.returncode, 0)
            m = _manifest(td)
            self.assertEqual(m["reason_code"], "tgeno_header_invalid")
            self.assertIn("size", m["parameters"]["detail"])


def _write_dup_panel(td, ind_ids):
    """同一 Master（Persistent Genetic ID=P）的两种表示：P.SG（Suffices=SG，偏好最高）与
    P.HO（Suffices=HO）。ind_ids 决定 .ind 里真有谁——SG 不在 .ind 时，HO 是唯一可用表示。"""
    pref = pathlib.Path(td) / "work/data/ref/aadr/panel"
    pref.parent.mkdir(parents=True, exist_ok=True)
    (pref.parent / "panel.ind").write_text(
        "".join(f"{i}\tM\tHan\n" for i in ind_ids), encoding="utf-8")
    (pref.parent / "panel.snp").write_text(
        "".join(f"{r}\t{c}\t0.0\t{pos}\t{a1}\t{a2}\n" for r, c, pos, a1, a2 in SITES), encoding="utf-8")
    (pref.parent / "panel.anno").write_text(
        "Genetic ID\tPersistent Genetic ID\tGroup ID\tDate mean in BP\tSuffices\n"
        "P.SG\tP\tHan\t\tSG\n"
        "P.HO\tP\tHan\t\tHO\n", encoding="utf-8")
    rlen = (N_SNP + 3) // 4
    head = f"TGENO   {len(ind_ids)}  {N_SNP} 00000000 00000000".encode()[:48].ljust(48, b"\x00")
    (pref.parent / "panel.geno").write_bytes(head + bytes(rlen * len(ind_ids)))
    return pref.parent


@unittest.skipIf(_SKIP, _SKIP)
class TestDedupStaysInsideTheAvailablePanel(unittest.TestCase):
    """复审 §3.2 P1 AN1：08 先对全 .anno 按固定表示偏好去重、再与 .ind 相交。

    同一 Master 的 SG 不在 .ind、HO 在 .ind 时，SG 被偏好留下、HO 被当作重复丢掉——
    唯一有基因型的表示从面板里无声消失。去重只能在**可用**候选内选择；且
    reference_metadata.tsv 是历史查询的输入，必须保留全部原始技术表示，不能再写成
    去重后的子集（09d 的"同一个人只算一次"由消费侧自己做，并负责暴露表示间冲突）。"""

    def test_only_available_representation_is_not_lost(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_dup_panel(td, ind_ids=["P.HO"])
            bin_dir = _write_bcftools(td, "ok")
            r = _run_08(cfg, bin_dir)
            self.assertEqual(r.returncode, 0, r.stderr)
            fam = (pathlib.Path(td) / "work/wgs/11_aadr/aadr.fam").read_text(encoding="utf-8")
            iids = [line.split()[1] for line in fam.strip().splitlines()]
            self.assertIn("P.HO", iids, "SG 不在 .ind 时，唯一可用的 HO 不得在去重时被丢掉")
            # 历史视图的输入保留两种原始表示（谁被提取、谁没被提取是提取侧的决定）
            meta = (pathlib.Path(td) / "work/wgs/11_aadr/reference_metadata.tsv").read_text(
                encoding="utf-8")
            self.assertEqual({row.split("\t")[0] for row in meta.strip().splitlines()[1:]},
                             {"P.SG", "P.HO"},
                             "reference_metadata.tsv 必须带全部原始表示，不是去重后的子集")

    def test_ok_manifest_records_the_annotation_content_sha(self):
        """复审 §3.2 P1 AN1（缓存绑定）：ok manifest 必须记录实际规范化的 .anno 内容指纹。
        09d 拿它与当前配置注释比对；没这份指纹，旧 metadata 换了注释也照样被当成当前查询。"""
        import hashlib
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_panel(td)
            bin_dir = _write_bcftools(td, "ok")
            r = _run_08(cfg, bin_dir)
            self.assertEqual(r.returncode, 0, r.stderr)
            anno = pathlib.Path(td) / "work/data/ref/aadr/panel.anno"
            want = hashlib.sha256(anno.read_bytes()).hexdigest()[:12]
            m = _manifest(td)
            self.assertEqual(m["parameters"].get("annotation_sha"), want,
                             "记录的是实际读的 .anno 的内容指纹（file_sha 同款 12 位）")

    def test_dedup_still_applies_when_both_representations_are_available(self):
        """两种表示都在 .ind 时仍去重（一人一份基因型）：修可用性排序不得顺手取消去重。"""
        with tempfile.TemporaryDirectory() as td:
            cfg = _write_config(td)
            _write_dup_panel(td, ind_ids=["P.SG", "P.HO"])
            bin_dir = _write_bcftools(td, "ok")
            r = _run_08(cfg, bin_dir)
            self.assertEqual(r.returncode, 0, r.stderr)
            fam = (pathlib.Path(td) / "work/wgs/11_aadr/aadr.fam").read_text(encoding="utf-8")
            iids = [line.split()[1] for line in fam.strip().splitlines()]
            self.assertEqual(iids, ["P.SG", "TESTSAMPLE"], "可用候选内仍按偏好取一份（SG 优于 HO）")
            meta = (pathlib.Path(td) / "work/wgs/11_aadr/reference_metadata.tsv").read_text(
                encoding="utf-8")
            self.assertEqual({row.split("\t")[0] for row in meta.strip().splitlines()[1:]},
                             {"P.SG", "P.HO"})


def _validate_fn():
    """从生产脚本里提取 validate_tgeno_header 本体做纯函数断言（08 是顶层脚本不能 import）。"""
    src = SCRIPT.read_text(encoding="utf-8")
    m = re.search(r"^def validate_tgeno_header\(.*?\n    return None$", src, re.M | re.S)
    if not m:
        raise AssertionError("validate_tgeno_header not found in the production script")
    ns = {"os": os}
    exec(m.group(0), ns)  # noqa: S102 -- 断言对象就是生产函数本体
    return ns["validate_tgeno_header"]


class TestValidateTgenoHeaderPure(unittest.TestCase):
    """头格式来自服务器实测的真实 v66 .geno 前 48 字节。"""

    def _file(self, header, n_ind, n_snp, extra_rows=0):
        rlen = (n_snp + 3) // 4
        head = header.encode()[:48].ljust(48, b"\x00")
        with tempfile.TemporaryDirectory() as td:
            p = pathlib.Path(td) / "p.geno"
            p.write_bytes(head + bytes(rlen * (n_ind + extra_rows)))
            yield str(p)

    def test_real_v66_shaped_header_passes(self):
        fn = _validate_fn()
        for p in self._file("TGENO   27594  584131 8c17d6d1 80974215", 27594, 584131):
            self.assertIsNone(fn(p, 27594, 584131))

    def test_rejects_bad_magic(self):
        fn = _validate_fn()
        for p in self._file("XGENO   2  6 0 0", 2, 6):
            self.assertIn("magic", fn(p, 2, 6))

    def test_rejects_non_numeric_dimensions(self):
        fn = _validate_fn()
        for p in self._file("TGENO  two six 0 0", 2, 6):
            self.assertIn("non-numeric", fn(p, 2, 6))

    def test_rejects_header_without_dimension_fields(self):
        fn = _validate_fn()
        for p in self._file("TGENO", 2, 6):
            self.assertIn("malformed", fn(p, 2, 6))

    def test_rejects_truncated_header(self):
        fn = _validate_fn()
        with tempfile.TemporaryDirectory() as td:
            p = pathlib.Path(td) / "short.geno"
            p.write_bytes(b"TGENO   2  6 0")
            self.assertIn("truncated", fn(str(p), 2, 6))

    def test_rejects_missing_file(self):
        fn = _validate_fn()
        self.assertIn("unreadable", fn("/nonexistent/x.geno", 2, 6))

    def test_rejects_dimension_mismatch_even_when_size_is_self_consistent(self):
        """头自称 5 人且文件真是 5 行（头-尺寸自洽）：与 .ind 的 2 行不符照样拒。"""
        fn = _validate_fn()
        for p in self._file("TGENO   5  6 0 0", 2, 6):
            self.assertIn("header n_ind 5 != .ind rows 2", fn(p, 2, 6))
        for p in self._file("TGENO   2  9 0 0", 2, 6):
            self.assertIn("header n_snp 9 != .snp rows 6", fn(p, 2, 6))

    def test_rejects_size_mismatch(self):
        fn = _validate_fn()
        for p in self._file("TGENO   2  6 0 0", 2, 6, extra_rows=1):
            self.assertIn("size", fn(p, 2, 6))


if __name__ == "__main__":
    unittest.main(verbosity=2)
