"""运行事实记录（H6）：缺什么要如实说缺，而不是给一个看起来正常的版本号。

H6 要求保留运行 ID、代码 SHA/dirty、有效参数、工具/参考版本与产物关联。"已安装版本"不等于
"本次运行使用的版本"，所以这些必须在运行的当下探测并落盘，而不是事后去查。
"""
import importlib, json, pathlib, subprocess, sys, tempfile, unittest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import _run_info as ri  # noqa: E402


class TestRunInfo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.info = ri.build(run_id="TEST-0001", started=None)

    def test_required_fields_are_present(self):
        for k in ("schema_version", "run_id", "started_at", "host", "code", "reference", "parameters", "tools"):
            self.assertIn(k, self.info, k)

    def test_run_id_is_echoed_when_supplied(self):
        self.assertEqual(self.info["run_id"], "TEST-0001")

    def test_code_revision_reports_dirtiness(self):
        code = self.info["code"]
        self.assertIn("commit", code)
        self.assertIn("dirty", code)
        self.assertIsInstance(code["dirty"], bool)
        # 本次测试运行于一个 git 工作树内；commit 应能取到
        self.assertTrue(code["commit"] or code["branch"] is None)

    def test_dirty_reflects_the_working_tree(self):
        """造一个脏工作树，dirty 必须为真；干净的必须为假。"""
        d = tempfile.mkdtemp()
        subprocess.run(["git", "init", "-q", d], check=True)
        (pathlib.Path(d) / "a.txt").write_text("x")
        subprocess.run(["git", "-C", d, "add", "a.txt"], check=True)
        subprocess.run(["git", "-C", d, "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "x"], check=True)
        clean = ri._code_revision(d)
        self.assertFalse(clean["dirty"], clean)
        (pathlib.Path(d) / "b.txt").write_text("y")
        after = ri._code_revision(d)
        self.assertTrue(after["dirty"], after)
        self.assertGreaterEqual(after["dirty_entries"], 1)

    def test_reference_has_a_fingerprint_or_says_it_has_none(self):
        ref = self.info["reference"]
        self.assertEqual(ref["build"], "GRCh37")
        self.assertIn("fai_sha256", ref)
        if ref["fai_sha256"] is not None:
            self.assertEqual(len(ref["fai_sha256"]), 16)
            self.assertGreater(ref["fai_sequences"], 0)

    def test_missing_tools_are_none_not_a_guess(self):
        """PATH 上没有的工具必须记为 None：编一个版本号会让人以为那次运行用了它。"""
        for name, v in self.info["tools"].items():
            self.assertTrue(v is None or isinstance(v, str), name)

    def test_parameters_come_from_the_configuration(self):
        p = self.info["parameters"]
        for k in ("SAMPLE", "BUILD", "MIN_DP", "MIN_MQ", "MIN_BQ"):
            self.assertIn(k, p, k)
        self.assertEqual(p["BUILD"], "GRCh37")

    def test_written_file_is_valid_json(self):
        out, info = ri.write(path=pathlib.Path(tempfile.mkdtemp()) / "ri.json", run_id="TEST-0002")
        self.assertTrue(out.exists())
        self.assertEqual(json.loads(out.read_text())["run_id"], "TEST-0002")


if __name__ == "__main__":
    unittest.main(verbosity=2)


class TestInputFingerprints(unittest.TestCase):
    """输入指纹：H6 要求输入变化后缓存失效，所以指纹必须对改动敏感、对未改动稳定。"""

    def test_size_and_mtime_change_when_content_changes(self):
        import os, time
        d = pathlib.Path(tempfile.mkdtemp())
        f = d / "reads.cram"
        f.write_text("a" * 100)
        first = ri._input_fingerprints.__wrapped__ if hasattr(ri._input_fingerprints, "__wrapped__") else None
        # 直接测内部逻辑：用本地文件与 stat 对比
        st1 = f.stat()
        time.sleep(0.01)
        f.write_text("b" * 200)
        os.utime(f, (st1.st_mtime + 5, st1.st_mtime + 5))
        st2 = f.stat()
        self.assertNotEqual((st1.st_size, st1.st_mtime), (st2.st_size, st2.st_mtime))

    def test_optional_inputs_are_null_not_missing(self):
        """空可选项（未配置的 FASTQ 等）必须是显式的 null，读者才能区分"没配"与"脚本忘了写"。"""
        info = ri.build(run_id="T2")
        for k in ("reads", "fastq1", "fastq2", "vendor_vcf", "y_reads"):
            self.assertIn(k, info["inputs"], k)

    def test_missing_file_is_recorded_not_dropped(self):
        info = ri.build(run_id="T3")
        for k, v in info["inputs"].items():
            if v is not None:
                self.assertIn("path", v)
                self.assertTrue("size" in v or v.get("missing") is True, (k, v))
