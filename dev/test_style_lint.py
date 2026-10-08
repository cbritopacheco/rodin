"""Regression tests for house-style rules."""

import unittest

from style_lint import check_file, doxygen_style_edits


class ForBracesTests(unittest.TestCase):
    def findings(self, text):
        return [finding for finding in check_file("src/Rodin/Example.cpp",
                                                 text.splitlines())
                if finding.check == "forbraces"]

    def test_one_line_body_and_wrapped_header_pass(self):
        self.assertEqual(self.findings(
            "for (int i = 0;\n     i < size; ++i)\n  consume(i);"), [])

    def test_wrapped_single_statement_requires_braces(self):
        findings = self.findings("for (auto x : xs)\n  consume(\n    x);")
        self.assertEqual(len(findings), 1)
        self.assertEqual(findings[0].line, 1)

    def test_outer_loop_requires_braces_but_inner_body_does_not(self):
        findings = self.findings(
            "for (auto row : rows)\n  for (auto x : row)\n    consume(x);")
        self.assertEqual([f.line for f in findings], [1])

    def test_both_levels_require_braces_when_inner_body_wraps(self):
        findings = self.findings(
            "for (auto row : rows)\n  for (auto x : row)\n    consume(\n      x);")
        self.assertEqual([f.line for f in findings], [1, 2])

    def test_braced_outer_and_one_line_inner_pass(self):
        self.assertEqual(self.findings(
            "for (auto row : rows)\n{\n  for (auto x : row)\n"
            "    consume(x);\n}"), [])

    def test_multiline_lambda_and_initializer_require_braces(self):
        for body in ("consume([] {\n  work();\n});",
                     "consume(Value{\n  1, 2});"):
            with self.subTest(body=body):
                self.assertEqual(len(self.findings("for (auto x : xs)\n" + body)), 1)

    def test_comments_literals_and_directives_do_not_create_loops(self):
        self.assertEqual(self.findings(
            '// for (auto x : xs)\n/* for (auto x : xs)\nconsume(x); */\n'
            'auto s = R"doc(for (auto x : xs)\nconsume(x);)doc";\n'
            '#define LOOP for (auto x : xs) \\\n  consume( \\\n    x);\n'), [])

    def test_multiline_literal_body_requires_braces(self):
        self.assertEqual(len(self.findings(
            'for (auto x : xs)\n  consume(R"doc(first\nsecond)doc");')), 1)

    def test_leading_comment_is_outside_body_extent(self):
        self.assertEqual(self.findings(
            "for (auto x : xs)\n  // Explain consumption.\n  consume(x);"), [])

    def test_nested_if_else_is_measured_as_complete_body(self):
        findings = self.findings(
            "for (auto x : xs)\n  if (ready(x)) consume(x);\n"
            "  else reject(x);")
        self.assertEqual([f.line for f in findings], [1])

    def test_do_while_and_try_catch_are_complete_bodies(self):
        for body in ("do\n  consume(x);\nwhile (ready(x));",
                     "try { consume(x); }\ncatch (...) { reject(x); }"):
            with self.subTest(body=body):
                self.assertEqual(len(self.findings("for (auto x : xs)\n" + body)), 1)

    def test_digit_separators_do_not_hide_following_loops(self):
        findings = self.findings(
            "auto size = 1'000;\nfor (auto x : xs)\n  consume(\n    x);\n"
            "auto other = 2'000;")
        self.assertEqual([f.line for f in findings], [2])

    def test_statement_attributes_preserve_complete_body(self):
        self.assertEqual(self.findings(
            "for (auto x : xs)\n  [[likely]] {\n    consume(x);\n  }"), [])
        self.assertEqual(len(self.findings(
            "for (auto x : xs)\n  [[likely]] if (ready(x)) consume(x);\n"
            "  else reject(x);")), 1)


class DoxygenStyleTests(unittest.TestCase):
    def findings(self, text):
        return [finding for finding in check_file("src/Rodin/Example.cpp",
                                                 text.splitlines())
                if finding.check == "docstyle"]

    def test_multiline_slashes_require_block(self):
        text = "  /// @brief Evaluate.\n  /// @param x Input.\n  void evaluate(int x);"
        self.assertEqual(len(self.findings(text)), 1)
        self.assertEqual(list(doxygen_style_edits(text))[0][2],
                         "  /**\n   * @brief Evaluate.\n   * @param x Input.\n   */")

    def test_single_content_line_requires_slashes(self):
        for text in ("/** @brief Evaluate. */", "/**\n * @brief Evaluate.\n */"):
            with self.subTest(text=text):
                self.assertEqual(len(self.findings(text)), 1)
                self.assertEqual(list(doxygen_style_edits(text))[0][2],
                                 "/// @brief Evaluate.")

    def test_preferred_forms_pass(self):
        self.assertEqual(self.findings(
            "/// @brief Evaluate.\nvoid evaluate();\n"
            "/**\n * @brief Evaluate.\n * @returns Value.\n */"), [])

    def test_strings_and_regular_comments_are_ignored(self):
        self.assertEqual(self.findings(
            'const char* s = R"doc(\n/// First\n/// Second\n/** One */\n)doc";\n'
            'const char* t = "\\\"/** One */";\n'
            '/*\n/// First\n/// Second\n*/\n'
            '// First\n// Second\n//// Divider\n//// Divider'), [])

    def test_trailing_and_inline_docs_are_ignored(self):
        self.assertEqual(self.findings(
            'int x; ///< First\nint y; ///< Second\n'
            '/**< Trailing */\nvoid f(/** Input */ int x);'), [])

    def test_separated_comments_do_not_form_block(self):
        self.assertEqual(self.findings('/// First\n\n/// Second'), [])

    def test_math_and_blank_lines_are_preserved(self):
        text = "/// @f$ x + y @f$\n///\n/// @returns Sum."
        self.assertEqual(list(doxygen_style_edits(text))[0][2],
                         "/**\n * @f$ x + y @f$\n *\n * @returns Sum.\n */")

    def test_code_indentation_is_preserved(self):
        text = "/// Example:\n///     evaluate(x);"
        self.assertEqual(list(doxygen_style_edits(text))[0][2],
                         "/**\n * Example:\n *     evaluate(x);\n */")


if __name__ == "__main__":
    unittest.main()
