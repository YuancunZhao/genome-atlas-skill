"""旧缺陷 #8a/#8b 回归：setup/00_install_tools.sh 的下载校验与 haplogrep3 树安装。

#8a 旧况：plink2 与 haplogrep3 走裸 curl 到 /tmp 再解压（截断下载会被后续步骤当真用），且
haplogrep3 的 URL 是 `releases/latest/download/` + 固定 3.2.2 文件名的矛盾组合——上游 tag
是 v 前缀（v3.2.2），今天 latest 已是 v3.3.2 且资产改 .tar.gz，旧 URL 直接 404。
#8b 旧况：安装脚本完全不装 phylotree 树。06_mtdna.py 的 `classify --tree phylotree-rcrs@17.2`
依赖 haplogrep3 自己的树目录（jar 旁 trees/，`--tree` 只接受 tree ID、显式路径被拒，服务器
2026-10-07 实测），全新机器上 06 会在装完二进制后仍然失败。

两层验证：静态断言锁住生产脚本的 URL/路由/树安装与正向校验；行为层把生产 `_fetch_tool`
函数本体从脚本文本中抽出来真跑（curl/stat 桩 + 真 unzip），证明 .part/大小下限/zip 完整性
检查真实生效——旧脚本在这两层都有断言失败。
"""
import pathlib, re, subprocess, tempfile, unittest, zipfile

REPO = pathlib.Path(__file__).resolve().parents[1]
INSTALL = REPO / "setup" / "00_install_tools.sh"
MTDNA = REPO / "scripts" / "06_mtdna.py"

_FETCH_FN_RE = re.compile(r"(?m)^_fetch_tool\(\)\{\n.*?\n\}\n", re.S)


def _src():
    return INSTALL.read_text(encoding="utf-8")


class TestInstallScriptSource(unittest.TestCase):
    """静态断言：旧脚本（/latest/download + 裸 curl + 无树安装）在此全部失败。"""

    def setUp(self):
        self.src = _src()

    def test_haplogrep3_url_is_pinned_to_a_v_tag(self):
        """上游 tag 是 v 前缀；`releases/latest/download/` + 固定文件名已 404（实测 latest 为
        v3.3.2 且资产为 .tar.gz）。"""
        self.assertNotIn("releases/latest/download", self.src)
        self.assertIn("HV=${HAPLOGREP3_VERSION:-3.2.2}", self.src)
        self.assertIn('releases/download/v$HV/haplogrep3-$HV-linux.zip', self.src)

    def test_bare_curl_downloads_are_gone(self):
        """plink2/haplogrep3 不再裸 curl 到 /tmp：必须过 _fetch_tool（.part + 大小下限 + zip 校验）。"""
        self.assertNotIn("curl -fsSL --retry 3 -o /tmp/p2.zip", self.src)
        self.assertNotIn("curl -fsSL --retry 3 -o /tmp/h3.zip", self.src)
        self.assertRegex(self.src, r"_fetch_tool https://s3\.amazonaws\.com/plink2-assets/")
        self.assertRegex(self.src, r'_fetch_tool "https://github\.com/genepi/haplogrep3/')

    def test_tree_install_with_positive_check(self):
        """#8b：fresh 安装必须装 phylotree-rcrs@17.2 并正向校验 tree.yaml/rcrs.fasta 落地。"""
        self.assertIn("install-tree phylotree-rcrs@17.2", self.src)
        self.assertIn('H3TREE="$TOOLS/trees/phylotree-rcrs/17.2"', self.src)
        self.assertIn('[ -f "$H3TREE/tree.yaml" ]', self.src)
        self.assertIn('[ -s "$H3TREE/rcrs.fasta" ]', self.src)

    def test_mtdna_keeps_tool_id_tree_form(self):
        """--tree 只接受 tree ID（显式路径被 3.3.2 拒绝，2026-10-07 实测）：06 必须保持 ID 形式，
        树的本地化由 jar 旁 trees/ 目录 + install-tree 承担。把它改成路径会弄坏 classify。"""
        self.assertIn("phylotree-rcrs@17.2", MTDNA.read_text(encoding="utf-8"))

    def test_shell_syntax(self):
        for sh in (INSTALL, REPO / "scripts" / "01_normalize.sh"):
            r = subprocess.run(["bash", "-n", str(sh)], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"{sh}: {r.stderr}")


class TestFetchToolBehaviour(unittest.TestCase):
    """行为层：抽取生产 `_fetch_tool` 本体，配 curl/stat 桩与真 unzip 真跑。

    stat 桩把 GNU `-c %s` 翻译成本机 stat（macOS 无 -c），保证本地与服务器同一路径。
    """

    def _harness(self, td):
        m = _FETCH_FN_RE.search(_src())
        self.assertIsNotNone(m, "生产脚本里找不到 _fetch_tool 函数本体")
        sh = pathlib.Path(td) / "fetch.sh"
        sh.write_text(
            "#!/bin/bash\nset -uo pipefail\n" + m.group(0) +
            '\n_fetch_tool "$@"\nrc=$?\necho "FETCH_RC=$rc"\nexit 0\n', encoding="utf-8")
        bin_ = pathlib.Path(td) / "bin"
        bin_.mkdir(exist_ok=True)
        (bin_ / "stat").write_text(
            "#!/bin/bash\n"
            '# GNU stat has -c %s natively; BSD (macOS) stat needs -f %z\n'
            'if [ "$1" = "-c" ] && [ "$2" = "%s" ]; then\n'
            '  if /usr/bin/stat --version >/dev/null 2>&1; then exec /usr/bin/stat -c %s "$3";\n'
            '  else exec /usr/bin/stat -f %z "$3"; fi\n'
            "fi\n"
            'exec /usr/bin/stat "$@"\n', encoding="utf-8")
        (bin_ / "curl").write_text(
            "#!/bin/bash\n"
            'log="${CURL_CALL_LOG:?}"; out=""; url=""\n'
            'prev=""; for a in "$@"; do [ "$prev" = "-o" ] && out="$a"; '
            'case "$a" in http*) url="$a";; esac; prev="$a"; done\n'
            'echo "$url" >> "$log"\n'
            'cat "${CURL_PAYLOAD:?}" > "$out"\n'
            "exit 0\n", encoding="utf-8")
        for f in bin_.iterdir():
            f.chmod(0o755)
        return sh, bin_

    def _run_fetch(self, td, out_name, min_bytes, payload):
        sh, bin_ = self._harness(td)
        payload_path = pathlib.Path(td) / "payload"
        payload_path.write_bytes(payload)
        call_log = pathlib.Path(td) / "curl_calls.log"
        out = pathlib.Path(td) / out_name
        env_extra = {
            "PATH": f"{bin_}{pathlib.os.pathsep}{pathlib.os.environ['PATH']}",
            "CURL_PAYLOAD": str(payload_path), "CURL_CALL_LOG": str(call_log),
        }
        r = subprocess.run(
            ["bash", str(sh), "https://example.invalid/x", str(out), str(min_bytes)],
            capture_output=True, text=True, env={**pathlib.os.environ, **env_extra}, timeout=60)
        return r, out, call_log

    @staticmethod
    def _zip_bytes(n):
        import io
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr("payload.bin", b"x" * n)
        return buf.getvalue()

    def test_too_small_download_fails_and_leaves_nothing(self):
        with tempfile.TemporaryDirectory() as td:
            r, out, _ = self._run_fetch(td, "tool.bin", 10000, b"short")
            self.assertNotIn("FETCH_RC=0", r.stdout)
            self.assertIn("download too small", r.stderr)
            self.assertFalse(out.exists())
            self.assertFalse(pathlib.Path(str(out) + ".part").exists())

    def test_non_zip_payload_with_zip_name_fails(self):
        with tempfile.TemporaryDirectory() as td:
            r, out, _ = self._run_fetch(td, "tool.zip", 1000, b"g" * 20000)
            self.assertNotIn("FETCH_RC=0", r.stdout)
            self.assertIn("not a valid zip", r.stderr)
            self.assertFalse(out.exists())
            self.assertFalse(pathlib.Path(str(out) + ".part").exists())

    def test_valid_zip_lands_atomically_and_part_is_consumed(self):
        with tempfile.TemporaryDirectory() as td:
            r, out, _ = self._run_fetch(td, "tool.zip", 1000, self._zip_bytes(20000))
            self.assertIn("FETCH_RC=0", r.stdout)
            self.assertTrue(out.exists())
            self.assertFalse(pathlib.Path(str(out) + ".part").exists())

    def test_undersized_existing_file_triggers_redownload(self):
        with tempfile.TemporaryDirectory() as td:
            payload = self._zip_bytes(20000)
            # 先跑一次拿到合法文件，再截成小于下限的旧缓存重跑：必须重新下载
            r, out, _ = self._run_fetch(td, "tool.zip", 1000, payload)
            self.assertIn("FETCH_RC=0", r.stdout)
            out.write_bytes(b"stale")
            r2, out2, calls = self._run_fetch(td, "tool.zip", 1000, payload)
            self.assertIn("FETCH_RC=0", r2.stdout)
            self.assertEqual(out2.stat().st_size, len(payload))
            self.assertEqual(len(calls.read_text().splitlines()), 2, "旧缓存过小必须触发重下")


if __name__ == "__main__":
    unittest.main(verbosity=2)
