"""Self-tests for check_versions.py.

    cd tools && python3 -m unittest -v test_check_versions
"""
import pathlib
import subprocess
import tempfile
import unittest

import check_versions as cv

GOOD = """# comment
erf             92410d578f612ab9e844c450bf8e900194becfb5
amrex           347a3a2d26750c0e35f0b949368645b872d5c5e9
fastwindterrain 648f8f1eb0fd833a47b2ffbabcb21e1cf2c3741e
kynema-sgf      8500c568c1053a4c417b7b355e7c383ef4fd0d60
"""


class Versions(unittest.TestCase):
    def test_good_file_parses(self):
        pins, errors = cv.parse(GOOD)
        self.assertEqual(errors, [])
        self.assertEqual(pins["erf"], "92410d578f612ab9e844c450bf8e900194becfb5")

    def test_short_sha_is_rejected(self):
        _, errors = cv.parse(GOOD.replace("92410d578f612ab9e844c450bf8e900194becfb5", "92410d5"))
        self.assertTrue(any("line 2" in e for e in errors))
        self.assertTrue(any("erf is not pinned" in e for e in errors))

    def test_missing_and_repeated_entries(self):
        text = GOOD.replace("kynema-sgf      8500c568c1053a4c417b7b355e7c383ef4fd0d60\n",
                            "amrex 347a3a2d26750c0e35f0b949368645b872d5c5e9\n")
        _, errors = cv.parse(text)
        self.assertTrue(any("amrex is pinned twice" in e for e in errors))
        self.assertTrue(any("kynema-sgf is not pinned" in e for e in errors))


class Submodule(unittest.TestCase):
    """check() compares the amrex pin with the commit checked out."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.tmp.name)
        amrex = self.root / "external" / "amrex"
        amrex.mkdir(parents=True)
        run = lambda *a: subprocess.run(["git", "-C", str(amrex), *a], check=True,
                                        stdout=subprocess.PIPE, text=True).stdout.strip()
        run("init", "-q")
        run("config", "user.email", "t@t")
        run("config", "user.name", "t")
        (amrex / "f").write_text("x")
        run("add", ".")
        run("commit", "-q", "-m", "one")
        self.head = run("rev-parse", "HEAD")

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, amrex_sha):
        text = GOOD.replace("347a3a2d26750c0e35f0b949368645b872d5c5e9", amrex_sha)
        (self.root / "external" / "VERSIONS").write_text(text)

    def test_matching_submodule_passes(self):
        self.write(self.head)
        self.assertEqual(cv.check(self.root), [])

    def test_other_commit_fails(self):
        self.write("0" * 40)
        errors = cv.check(self.root)
        self.assertEqual(len(errors), 1)
        self.assertIn("external/amrex is at " + self.head, errors[0])


if __name__ == "__main__":
    unittest.main()
