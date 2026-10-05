#!/usr/bin/env python3
"""Integration regressions for the repository's Doxygen warning policy."""

import contextlib
import io
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import check_doxygen_warnings as checker


class DoxygenWarningsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = shutil.which(os.environ.get("DOXYGEN", "doxygen"))
        if cls.binary is None:
            raise unittest.SkipTest("Doxygen is required")
        cls.template = Path(checker.REPO, "doc", "Doxygen.in").read_text()
        cls.version = checker.doxygen_version(cls.binary)

    def check_fixture(self, source):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            (repo / "doc").mkdir()
            (repo / "src").mkdir()
            (repo / "examples").mkdir()
            (repo / "src" / "fixture.h").write_text(source)
            (repo / "doc" / "Doxygen.in").write_text(
                self.template + "\nCITE_BIB_FILES =\n")
            baseline = repo / "baseline"
            baseline.write_text(f"# doxygen {self.version}\n")
            with patch.object(checker, "REPO", directory), \
                    patch.object(checker, "BASELINE_PATH", str(baseline)), \
                    patch.object(sys, "argv", ["check", "--doxygen", self.binary]):
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    result = checker.main()
            return result, output.getvalue()

    def test_completely_missing_parameter_documentation_fails(self):
        status, output = self.check_fixture("""/** @file */
/** @brief Sets two values. */
void missing(int first, int second);
""")
        self.assertEqual(status, 1)
        self.assertIn("parameters of member missing are not documented", output)

    def test_partially_documented_parameters_fail(self):
        status, output = self.check_fixture("""/** @file */
/** @brief Sets two values.
 * @param first First value.
 */
void partial(int first, int second);
""")
        self.assertEqual(status, 1)
        self.assertIn("is not documented", output)
        self.assertIn("partial", output)

    def test_missing_return_documentation_fails(self):
        status, output = self.check_fixture("""/** @file */
/** @brief Computes the square.
 * @param value Value to square.
 */
int square(int value);
""")
        self.assertEqual(status, 1)
        self.assertIn("return type of member square is not documented", output)

    def test_complete_documentation_passes(self):
        status, output = self.check_fixture("""/** @file */
/** @brief Computes the square.
 * @param value Value to square.
 * @returns Square of the value.
 */
int square(int value);
""")
        self.assertEqual(status, 0, output)

    def test_failed_doxygen_with_empty_log_cannot_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory, "failed-doxygen")
            binary.write_text(
                f"#!{sys.executable}\n"
                "import pathlib, re, sys\n"
                "if sys.argv[1] == '--version':\n"
                f"    print('{self.version}')\n"
                "    sys.exit(0)\n"
                "config = pathlib.Path(sys.argv[1]).read_text()\n"
                "log = re.findall(r'^WARN_LOGFILE = (.*)$', config, re.M)[-1]\n"
                "pathlib.Path(log).write_text('')\n"
                "sys.exit(7)\n")
            binary.chmod(0o755)
            with patch.object(sys, "argv", ["check", "--doxygen", str(binary)]):
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    status = checker.main()
            self.assertEqual(status, 2)
            self.assertIn("doxygen failed with exit status 7", output.getvalue())


if __name__ == "__main__":
    unittest.main()
