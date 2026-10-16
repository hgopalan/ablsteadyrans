"""Self-tests for check_extern_erf.py, on a throw-away ERF git repository.

    cd tools && python3 -m unittest -v test_check_extern_erf
"""
import pathlib
import subprocess
import tempfile
import unittest

import check_extern_erf as cx

HEADER = """# Code copied from ERF
<!-- manifest:start -->
| Local file | ERF source | ERF SHA | Status | Notes |
|---|---|---|---|---|
{rows}
<!-- manifest:end -->
"""


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          stdout=subprocess.PIPE, text=True).stdout.strip()


class DriftCheck(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self.tmp.name)
        self.erf = base / "erf"
        self.root = base / "abl"
        (self.erf / "Source").mkdir(parents=True)
        (self.root / "extern_erf" / "patches").mkdir(parents=True)
        git(self.erf, "init", "-q")
        git(self.erf, "config", "user.email", "t@t")
        git(self.erf, "config", "user.name", "t")
        (self.erf / "Source" / "A.H").write_text("line 1\nline 2\n")
        (self.erf / "Source" / "B.H").write_text("alpha\nbeta\n")
        git(self.erf, "add", ".")
        git(self.erf, "commit", "-q", "-m", "pin")
        self.sha = git(self.erf, "rev-parse", "--short", "HEAD")

    def tearDown(self):
        self.tmp.cleanup()

    def manifest(self, *rows):
        text = HEADER.format(rows="\n".join(
            f"| `{l}` | `{s}` | {sha} | {st} | |" for l, s, sha, st in rows))
        (self.root / "extern_erf" / "EXTRACTED.md").write_text(text)

    def local(self, name, text):
        (self.root / "extern_erf" / name).write_text(text)

    def run_check(self, report=None):
        lines = []
        errs = cx.check(self.root, self.erf, report, out=lines.append)
        return errs, lines

    def test_empty_manifest_passes_without_erf(self):
        self.manifest()
        self.assertEqual(cx.check(self.root, None), [])

    def test_unchanged_copy_passes(self):
        self.manifest(("extern_erf/A.H", "Source/A.H", self.sha, "unchanged"))
        self.local("A.H", "line 1\nline 2\n")
        self.assertEqual(self.run_check()[0], [])

    def test_unchanged_copy_with_an_edit_fails(self):
        self.manifest(("extern_erf/A.H", "Source/A.H", self.sha, "unchanged"))
        self.local("A.H", "line 1\nline 2 edited\n")
        errs, _ = self.run_check()
        self.assertEqual(len(errs), 1)
        self.assertIn("recorded as unchanged", errs[0])

    def test_modified_copy_with_current_patch_passes(self):
        self.manifest(("extern_erf/B.H", "Source/B.H", self.sha, "modified"))
        self.local("B.H", "alpha\nbeta, local\n")
        patch = cx.make_patch("alpha\nbeta\n", "alpha\nbeta, local\n", "Source/B.H", "extern_erf/B.H")
        (self.root / "extern_erf" / "patches" / "B.H.patch").write_text(patch)
        self.assertEqual(self.run_check()[0], [])

    def test_modified_copy_with_stale_patch_fails(self):
        self.manifest(("extern_erf/B.H", "Source/B.H", self.sha, "modified"))
        self.local("B.H", "alpha\nbeta, local, edited again\n")
        patch = cx.make_patch("alpha\nbeta\n", "alpha\nbeta, local\n", "Source/B.H", "extern_erf/B.H")
        (self.root / "extern_erf" / "patches" / "B.H.patch").write_text(patch)
        errs, _ = self.run_check()
        self.assertEqual(len(errs), 1)
        self.assertIn("does not match the current local changes", errs[0])

    def test_modified_copy_without_patch_fails(self):
        self.manifest(("extern_erf/B.H", "Source/B.H", self.sha, "modified"))
        self.local("B.H", "alpha\nbeta, local\n")
        errs, _ = self.run_check()
        self.assertIn("is missing", errs[0])

    def test_unlisted_file_fails_and_exempt_files_do_not(self):
        self.manifest()
        self.local("Stray.H", "x\n")
        self.local("LICENSE.md", "licence\n")
        (self.root / "extern_erf" / "patches" / "old.patch").write_text("")
        errs = cx.check(self.root, None)
        self.assertEqual(errs, ["extern_erf/Stray.H: in extern_erf/ but has no manifest row"])

    def test_missing_local_file_fails(self):
        self.manifest(("extern_erf/A.H", "Source/A.H", self.sha, "unchanged"))
        errs, _ = self.run_check()
        self.assertIn("does not exist", errs[0])

    def test_bad_status_fails(self):
        self.manifest(("extern_erf/A.H", "Source/A.H", self.sha, "tweaked"))
        self.local("A.H", "line 1\nline 2\n")
        errs, _ = self.run_check()
        self.assertIn("status must be one of", errs[0])

    def test_missing_header_row_fails(self):
        (self.root / "extern_erf" / "EXTRACTED.md").write_text(
            "<!-- manifest:start -->\n"
            f"| `extern_erf/A.H` | `Source/A.H` | {self.sha} | unchanged | |\n"
            f"| `extern_erf/B.H` | `Source/B.H` | {self.sha} | unchanged | |\n"
            "<!-- manifest:end -->\n")
        errs = cx.check(self.root, self.erf)
        self.assertIn("header row", errs[0])

    def test_upstream_report_lists_later_commits_and_never_fails(self):
        self.manifest(("extern_erf/A.H", "Source/A.H", self.sha, "unchanged"))
        self.local("A.H", "line 1\nline 2\n")
        (self.erf / "Source" / "A.H").write_text("line 1\nline 2\nline 3\n")
        git(self.erf, "commit", "-q", "-am", "upstream fix to A")
        errs, lines = self.run_check(report="HEAD")
        self.assertEqual(errs, [])
        self.assertTrue(any("upstream fix to A" in ln for ln in lines))


if __name__ == "__main__":
    unittest.main()
