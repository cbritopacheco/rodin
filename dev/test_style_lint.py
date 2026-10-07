"""Regression tests for the Doxygen comment-style rule."""

import unittest

from style_lint import check_file, doxygen_style_edits


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
