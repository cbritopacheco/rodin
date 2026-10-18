"""Regression checks for the repository's two-space record indentation."""

import os
from pathlib import Path
import contextlib
import io
import subprocess
import tempfile
import unittest
from unittest import mock

from format_cpp import format_source
import check_format


ROOT = Path(__file__).resolve().parent.parent


class FormatConfigTests(unittest.TestCase):
    def check_format(self, expected):
        binary = os.environ.get("CLANG_FORMAT_BINARY", "clang-format")
        # ReflowComments is disabled to preserve mathematical documentation.
        # Keep block-comment layout intact while perturbing C++ indentation.
        source = "".join(
            line if line.lstrip().startswith(("/**", "*"))
            else line.replace("  ", "    ")
            for line in expected.splitlines(keepends=True))
        for text in (source, expected):
            with self.subTest(idempotent=text == expected):
                result = format_source(binary, str(ROOT / "src/FormatFixture.cpp"), text)
                self.assertEqual(result, expected)

    def test_requested_class_and_plain_table_layout(self):
        self.check_format("""class Example
{
  public:
    void run();

  private:
    size_t m_dimension;
};

struct Table
{
  size_t nqp;
  size_t ndof;
};
""")

    def test_plain_struct_and_method_body(self):
        self.check_format("""namespace Example
{
  struct Table
  {
    int value;
    int get() const
    {
      return value;
    }
  };
}
""")

    def test_nested_struct_and_access_labels(self):
        self.check_format("""namespace Example
{
  class Store
  {
    public:
      struct Table
      {
        int value;
      };

    protected:
      int capacity;

    private:
      int next;
  };
}
""")

    def test_documented_members_after_access_labels(self):
        self.check_format("""class Example
{
  public:
    /**
     * @brief Runs the example.
     *
     * Documentation follows the member indentation.
     */
    static void run()
    {
      consume();
    }

  protected:
    int value;

  private:
    int next;
};
""")

    def test_struct_with_explicit_access_labels(self):
        self.check_format("""struct Example
{
  public:
    int value;

  protected:
    int capacity;

  private:
    int next;
};
""")

    def test_nested_records_inside_plain_struct(self):
        self.check_format("""struct Outer
{
  struct Inner
  {
    int value;
  };

  class Store
  {
    public:
      int value;

    private:
      int next;
  };
};
""")

    def test_enum_and_elaborated_return_type(self):
        self.check_format("""enum struct Kind
{
  First,
  Second
};

struct Table;
struct Table* make()
{
  consume();
  return nullptr;
}
""")

    def test_elaborated_parameter_and_variable(self):
        self.check_format("""void consume(struct Table value)
{
  process(value);
}

struct Table value
{
    1, 2,
};
""")

    def test_attributes_inheritance_and_specialization(self):
        self.check_format("""template <class T>
struct Table;

template <>
struct [[nodiscard]] Table<int> : Base<int>
{
  int value;
};
""")

    def test_plain_struct_documentation(self):
        self.check_format("""struct Table
{
  /**
   * @brief Number of quadrature points.
   *
   * Documentation follows the field indentation.
   */
  int nqp;
};
""")

    def test_comment_before_access_section_is_idempotent(self):
        self.check_format("""class Example
{
  /**
   * @brief Public interface.
   *
   * Mathematical indentation after the star is preserved:
   *   value = expression;
   */

  public:
    void run();
};
""")

    def test_inconsistent_comment_star_indentation_is_repaired(self):
        expected = """class Example
{
  public:
    /**
     * @brief Runs the example.
     *
     *   Mathematical body indentation is retained.
     */
    void run();
};
"""
        source = expected.replace("     *", " *")
        binary = os.environ.get("CLANG_FORMAT_BINARY", "clang-format")
        actual = format_source(binary, str(ROOT / "src/FormatFixture.cpp"), source)
        self.assertEqual(actual, expected)
        self.assertEqual(format_source(
            binary, str(ROOT / "src/FormatFixture.cpp"), actual), expected)

    def test_literals_and_comments_do_not_define_access_sections(self):
        # Do not perturb literal contents: whitespace there is application data.
        binary = os.environ.get("CLANG_FORMAT_BINARY", "clang-format")
        expected = '''struct Table
{
  // public: struct Fake { };
  const char* text = R"tag(
    public: struct Fake
    {
        literal data
    };
)tag";
  char brace = '}';
};
'''
        for _ in range(2):
            actual = format_source(binary, str(ROOT / "src/FormatFixture.cpp"), expected)
            self.assertEqual(actual, expected)
            expected = actual

    def test_encoded_raw_literal_contents_are_preserved(self):
        binary = os.environ.get("CLANG_FORMAT_BINARY", "clang-format")
        for prefix in ("u8", "u", "U", "L"):
            with self.subTest(prefix=prefix):
                expected = '''struct Table
{
  const auto text = PREFIXR"tag(
    "quoted" public: struct Fake
    {
        literal data
    };
)tag";
};
'''.replace("PREFIX", prefix)
                self.assertEqual(format_source(
                    binary, str(ROOT / "src/FormatFixture.cpp"), expected), expected)


class ChangedLineFormatTests(unittest.TestCase):
    def run_check(self, content, changed_line):
        binary = os.environ.get("CLANG_FORMAT_BINARY", "clang-format")
        diff = ("--- a/src/Fixture.cpp\n+++ b/src/Fixture.cpp\n"
                f"@@ -{changed_line} +{changed_line} @@\n-old\n+new\n")

        def git(*args, **kwargs):
            self.assertIn(args[0], ("diff", "show"))
            return subprocess.CompletedProcess(
                args, 0, stdout=diff if args[0] == "diff" else content)

        with mock.patch.object(check_format, "git", side_effect=git), \
                mock.patch.object(check_format, "resolve_base", return_value="base"), \
                mock.patch("sys.argv", ["check_format.py", "--head", "HEAD",
                                        "--binary", binary, "--no-color"]), \
                contextlib.redirect_stdout(io.StringIO()):
            return check_format.main()

    def test_plain_struct_four_spaces_rejected(self):
        self.assertEqual(self.run_check("struct Table\n{\n    int nqp;\n};\n", 3), 1)

    def test_plain_struct_two_spaces_accepted(self):
        self.assertEqual(self.run_check("struct Table\n{\n  int nqp;\n};\n", 3), 0)

    def test_explicit_access_layout_accepted(self):
        self.assertEqual(self.run_check(
            "class Example\n{\n  public:\n    void run();\n};\n", 4), 0)

    def test_explicit_access_members_need_four_spaces(self):
        self.assertEqual(self.run_check(
            "class Example\n{\n  public:\n  void run();\n};\n", 4), 1)

    def test_untouched_legacy_layout_is_not_a_failure(self):
        self.assertEqual(self.run_check("struct Table\n{\n    int nqp;\n};\n", 1), 0)

    def test_insertions_include_both_neighbours(self):
        hunks = list(check_format.formatting_hunks("first\nlast\n", "first\n\nlast\n"))
        self.assertEqual(hunks[0][:2], (1, 2))
        first = list(check_format.formatting_hunks("first\n", "\nfirst\n"))
        self.assertEqual(first[0][:2], (0, 1))

    def test_repeated_braces_do_not_hide_an_isolated_edit(self):
        original = "{\n}\n" * 300 + "    int value;\n" + "{\n}\n" * 300
        formatted = original.replace("    int value;", "  int value;")
        hunks = list(check_format.formatting_hunks(original, formatted))
        self.assertEqual(len(hunks), 1)
        self.assertEqual(hunks[0][:2], (601, 601))

    def test_missing_final_newline_is_detected(self):
        hunks = list(check_format.formatting_hunks("int value;", "int value;\n"))
        self.assertEqual(hunks[0][:2], (1, 1))

    def test_fix_leaves_unrelated_legacy_layout_alone(self):
        source = "    first;\n\n    second;\n"
        canonical = "  first;\n\n  second;\n"
        self.assertEqual(check_format.format_changed_source(source, canonical, {1}),
                         "  first;\n\n    second;\n")

    def test_fix_handles_insertions_deletions_and_final_newline(self):
        self.assertEqual(check_format.format_changed_source(
            "first\nlast\n", "first\n\nlast\n", {2}), "first\n\nlast\n")
        self.assertEqual(check_format.format_changed_source(
            "first\n\nlast\n", "first\nlast\n", {2}), "first\nlast\n")
        self.assertEqual(check_format.format_changed_source(
            "int value;", "int value;\n", {1}), "int value;\n")

    def test_fix_is_stable_when_applied_twice(self):
        source = "first\n    int value;\nlast\n"
        canonical = "first\n  int value;\nlast\n"
        fixed = check_format.format_changed_source(source, canonical, {2})
        self.assertEqual(fixed, canonical)
        self.assertEqual(check_format.format_changed_source(fixed, canonical, {2}), fixed)

    def test_repetitive_reindent_does_not_delete_a_statement(self):
        source = "  value;\n    value;\n      value;\n"
        canonical = "    value;\n      value;\n        value;\n"
        self.assertEqual(check_format.format_changed_source(source, canonical, {1}),
                         "    value;\n    value;\n      value;\n")

    def test_fix_rejects_token_and_literal_changes(self):
        for original, canonical in (("int value;\n", "int other;\n"),
                                    ('const char* s = "data";\n',
                                     'const char* s = "other";\n')):
            with self.subTest(original=original), self.assertRaises(ValueError):
                check_format.format_changed_source(original, canonical, {1})

    def test_fix_refuses_historical_head(self):
        with mock.patch("sys.argv", ["check_format.py", "--head", "HEAD", "--fix"]), \
                mock.patch.object(check_format, "git") as git, \
                contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:
                check_format.main()
            self.assertEqual(error.exception.code, 2)
            git.assert_not_called()

    def test_fix_cli_only_writes_changed_regions(self):
        source = ("struct Table\n{\n    int nqp;\n};\n\n"
                  "struct Legacy\n{\n    int ndof;\n};\n")
        diff = "--- a/src/Fixture.cpp\n+++ b/src/Fixture.cpp\n@@ -3 +3 @@\n-old\n+new\n"
        binary = os.environ.get("CLANG_FORMAT_BINARY", "clang-format")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "src").mkdir()
            (root / ".clang-format").write_text((ROOT / ".clang-format").read_text())
            fixture = root / "src/Fixture.cpp"
            fixture.write_text(source)
            with mock.patch.object(check_format, "REPO", directory), \
                    mock.patch.object(check_format, "resolve_base", return_value="base"), \
                    mock.patch.object(check_format, "git", return_value=
                                      subprocess.CompletedProcess([], 0, stdout=diff)), \
                    mock.patch("sys.argv", ["check_format.py", "--fix",
                                            "--binary", binary, "src/Fixture.cpp"]), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(check_format.main(), 0)
            self.assertEqual(fixture.read_text(), source.replace("    int nqp;",
                                                                "  int nqp;"))

    def test_fix_cli_recomputes_ranges_until_stable(self):
        source = ("struct Table\n{\n    int nqp;\n};\n\n"
                  "struct Legacy\n{\n    int ndof;\n};\n")
        diff = "--- a/src/Fixture.cpp\n+++ b/src/Fixture.cpp\n@@ -3 +3 @@\n-old\n+new\n"
        binary = os.environ.get("CLANG_FORMAT_BINARY", "clang-format")
        calls = []

        def git(*args, **kwargs):
            calls.append(args)
            # Model Git exposing another changed line after the first repair.
            changes = diff if len(calls) == 1 else diff + "@@ -8 +8 @@\n-old\n+new\n"
            return subprocess.CompletedProcess(args, 0, stdout=changes)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "src").mkdir()
            (root / ".clang-format").write_text((ROOT / ".clang-format").read_text())
            fixture = root / "src/Fixture.cpp"
            fixture.write_text(source)
            with mock.patch.object(check_format, "REPO", directory), \
                    mock.patch.object(check_format, "resolve_base", return_value="base"), \
                    mock.patch.object(check_format, "git", side_effect=git), \
                    mock.patch("sys.argv", ["check_format.py", "--fix",
                                            "--binary", binary, "src/Fixture.cpp"]), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(check_format.main(), 0)
            self.assertEqual(fixture.read_text(), source.replace("    int", "  int"))
            self.assertEqual(len(calls), 3)


if __name__ == "__main__":
    unittest.main()
