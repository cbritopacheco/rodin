"""Regression tests for the optional syntactic production-source audit."""

import importlib.util
from pathlib import Path
import tempfile
import unittest

from dev.audit_source import audit


@unittest.skipUnless(importlib.util.find_spec("tree_sitter") and
                     importlib.util.find_spec("tree_sitter_cpp"),
                     "Optional tree-sitter/tree-sitter-cpp dependencies are unavailable")
class SourceAuditTest(unittest.TestCase):
    def inventory(self, source, lines=None):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "src/Rodin/Test").mkdir(parents=True)
            (root / "src/Rodin/A.h").write_text(source)
            (root / "src/Rodin/Test/Helper.h").write_text("void helper() {}")
            return audit(root, lines or {})

    def test_inactive_templates_and_exception_context(self):
        result = self.inventory("""
template<class T> struct A {
  void absent() { throw "bad"; }
  void member() { Alert::Exception(); }
  static void stat() { Alert::Exception(); }
};
#if 0
void inactive() { Alert::Exception(); }
#endif
""")
        self.assertEqual(len(result["files"]), 1)
        self.assertEqual(len(result["definitions"]), 4)
        self.assertEqual(sum(d["template"] for d in result["definitions"]), 3)
        self.assertTrue(any("inactive" in d["signature"] for d in result["definitions"]))
        self.assertEqual([e["text"] for e in result["exceptions"] if e["kind"] == "throw"],
                         ['throw "bad";'])
        self.assertEqual({e["context"] for e in result["exceptions"]
                          if e["kind"] == "Alert::Exception"},
                         {"instance_member", "static", "free"})
        self.assertFalse(result["parse_errors"])

    def test_constructor_initializers_have_separate_attribution(self):
        result = self.inventory("struct A {\n A()\n : value(3)\n {}\n int value;\n};\n",
                                {"src/Rodin/A.h": {"3": True}})
        constructor, = result["definitions"]
        self.assertTrue(constructor["has_initializers"])
        self.assertTrue(constructor["empty_body"])
        self.assertEqual(constructor["measured_definition_lines"], 1)
        self.assertEqual(constructor["hit_definition_lines"], 1)
        self.assertEqual(constructor["measured_body_lines"], 0)

    def test_implicit_candidates_ignore_comments_and_nested_call_prefixes(self):
        result = self.inventory("""
void f() {
  a.at(2); a.value(); std::stoul(s); std::any_cast<int>(v);
  std::get<0>(tuple).getValue();
} // b.at(3)
""")
        self.assertEqual([c["call"] for c in result["implicit_exception_candidates"]],
                         ["a.at", "a.value", "std::stoul", "std::any_cast<int>", "std::get<0>"])
        self.assertFalse(result["parse_errors"])

    def test_parser_recovery_remains_visible(self):
        result = self.inventory("void broken( {\n")
        self.assertTrue(result["files"][0]["parse_error"])
        self.assertTrue(result["parse_errors"])


if __name__ == "__main__":
    unittest.main()
