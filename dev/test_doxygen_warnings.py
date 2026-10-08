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
        cls.mcss_template = Path(checker.REPO, "doc", "Doxygen.mcss.in").read_text()
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
            (repo / "doc" / "Doxygen.mcss.in").write_text(self.mcss_template)
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

    def test_mcss_navigation_with_single_line_reference_passes(self):
        status, output = self.check_fixture('''/// @file

/// @brief Does nothing.
void noop();
/**
 * @page guide Guide
 * See @ref target "Target guide" for details.
 * @m_footernavigation
 */
/** @page target Target guide */
''')
        self.assertEqual(status, 0, output)

    def test_mcss_navigation_with_broken_reference_fails(self):
        status, output = self.check_fixture('''/// @file

/// @brief Does nothing.
void noop();
/**
 * @page guide Guide
 * See @ref target "Target
 * guide" for details.
 * @m_footernavigation
 */
/** @page target Target guide */
''')
        self.assertEqual(status, 1, output)
        self.assertIn("unexpected command endxmlonly", output)

    def test_unnamed_parameter_fails_xml_audit(self):
        status, output = self.check_fixture("""/// @file
/** @brief Constant function.
 * @returns One.
 */
int constant(int) { return 1; }
""")
        self.assertEqual(status, 1)
        self.assertIn("is unnamed and lacks parameter documentation", output)

    def test_internal_parameter_and_return_omissions_fail(self):
        status, output = self.check_fixture("""/// @file
/// @internal
int helper(int value) { return value; }
""")
        self.assertEqual(status, 1)
        self.assertIn("lacks parameter documentation", output)
        self.assertIn("lacks return documentation", output)

    def test_private_static_helper_omissions_fail(self):
        status, output = self.check_fixture("""/// @file
/// @brief Cache object.
class Cache {
 private:
  static int helper(int value) { return value; }
};
""")
        self.assertEqual(status, 1)
        self.assertIn("lacks parameter documentation", output)
        self.assertIn("lacks return documentation", output)

    def test_empty_parameter_and_return_descriptions_fail(self):
        status, output = self.check_fixture("""/// @file
/** @brief Identity function.
 * @param value
 * @returns
 */
int identity(int value) { return value; }
""")
        self.assertEqual(status, 1)
        self.assertIn("lacks parameter documentation", output)
        self.assertIn("lacks return documentation", output)

    def test_invalid_private_parameter_documentation_fails(self):
        status, output = self.check_fixture("""/// @file
/// @brief Cache object.
class Cache {
 private:
  /** @brief Identity helper.
   * @param value Input value.
   * @param nonexistent Invalid parameter name.
   * @returns Input value.
   */
  static int helper(int value) { return value; }
};
""")
        self.assertEqual(status, 1)
        self.assertIn("nonexistent", output)

    def test_conversion_operator_requires_return_description(self):
        status, output = self.check_fixture("""/// @file
/// @brief Readiness marker.
struct Ready {
  /// @brief Tests readiness.
  explicit operator bool() const { return true; }
};
""")
        self.assertEqual(status, 1)
        self.assertIn("return type of member operator bool lacks return documentation", output)

    def test_unnamed_deduction_guide_parameter_fails(self):
        status, output = self.check_fixture("""/// @file
/** @brief Value wrapper.
 * @tparam T Value type.
 */
template <class T> struct Box {
  /** @brief Construct a wrapper.
   * @param value Wrapped value.
   */
  Box(const T& value);
};
/// @brief Deduce the wrapped type.
template <class T> Box(const T&) -> Box<T>;
""")
        self.assertEqual(status, 1)
        self.assertIn("is unnamed and lacks parameter documentation", output)

    def test_constructors_void_and_retval_do_not_need_returns(self):
        status, output = self.check_fixture("""/// @file
/** @brief Value wrapper. */
struct Box {
  /// @brief Default constructor.
  constexpr Box() = default;
  /** @brief Deleted copy constructor.
   * @param other Source object.
   */
  Box(const Box& other) = delete;
  /** @brief Deleted assignment.
   * @param other Source object.
   */
  Box& operator=(const Box& other) = delete;
  /// @brief No-op.
  void reset() {}
  /** @brief Gets the status.
   * @retval true Always ready.
   */
  bool ready() const { return true; }
};
""")
        self.assertEqual(status, 0, output)

    def test_successful_doxygen_without_xml_cannot_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory, "no-xml-doxygen")
            binary.write_text(
                f"#!{sys.executable}\n"
                "import pathlib, re, sys\n"
                "if sys.argv[1] == '--version':\n"
                f"    print('{self.version}')\n"
                "else:\n"
                "    cfg = pathlib.Path(sys.argv[1]).read_text()\n"
                "    log = re.findall(r'^WARN_LOGFILE = (.+)$', cfg, re.M)[-1]\n"
                "    pathlib.Path(log).write_text('')\n")
            binary.chmod(0o755)
            with patch.object(sys, "argv", ["check", "--doxygen", str(binary)]):
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    status = checker.main()
            self.assertEqual(status, 2)
            self.assertIn("Doxygen produced no XML index", output.getvalue())

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
