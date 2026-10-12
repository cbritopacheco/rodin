#!/usr/bin/env python3
"""Regression checks for the production coverage denominator and template audit."""

import tempfile
import unittest
from pathlib import Path

from audit_coverage import audit


class CoverageAuditTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name).resolve()
        self.source = self.root / "src" / "Rodin"
        self.source.mkdir(parents=True)
        (self.source / "Example.h").write_text("inline int answer() { return 42; }\n")

    def report(self, name, records):
        path = self.root / name
        path.write_text(records)
        return path

    def test_union_preserves_zero_hit_lines_and_avoids_duplicate_denominator(self):
        first = self.report("first.info", "SF:src/Rodin/Example.h\nDA:1,3\nDA:2,0\nend_of_record\n")
        second = self.report("second.info", "SF:src/Rodin/Example.h\nDA:1,0\nDA:2,7\nDA:3,0\nend_of_record\n")
        result = audit(self.root, [first, second])
        self.assertEqual(result["lines"], 3)
        self.assertEqual(result["hits"], 2)
        self.assertEqual(result["files"][0]["uncovered_lines"], [3])

    def test_test_benchmark_and_library_helper_records_are_rejected(self):
        records = ""
        for name in ("tests/Test.cpp", "benchmarks/Bench.cpp", "src/Rodin/Test/Random.h"):
            records += f"SF:{name}\nDA:1,1\nend_of_record\n"
        result = audit(self.root, [self.report("foreign.info", records)])
        self.assertEqual(result["lines"], 0)
        self.assertEqual(len(result["nonproduction_records"]), 3)

    def test_absent_template_inventory_preserves_locations_and_ignores_comments(self):
        (self.source / "Template.h").write_text(
            "// template <class Fake>\n/* template <class AlsoFake>\n */\ntemplate <class T>\nstruct Box {};\n")
        helpers = self.source / "Test"
        helpers.mkdir()
        (helpers / "Helper.h").write_text("template <class T> struct Helper {};\n")
        result = audit(self.root, [])
        absent = {file["file"]: file for file in result["absent_sources"]}
        self.assertEqual(absent["src/Rodin/Template.h"]["template_declarations"], [4])
        self.assertNotIn("src/Rodin/Test/Helper.h", absent)
        self.assertIsNone(result["coverage_percent"])

    def test_absolute_and_relative_paths_identify_the_same_source(self):
        report = self.report("paths.info", f"SF:{self.source / 'Example.h'}\nDA:1,1\nend_of_record\nSF:src/Rodin/Example.h\nDA:1,0\nend_of_record\n")
        result = audit(self.root, [report])
        self.assertEqual(result["lines"], 1)
        self.assertEqual(result["hits"], 1)
        self.assertFalse(result["nonproduction_records"])

    def test_shared_line_hits_do_not_hide_zero_count_function_variants(self):
        report = self.report("functions.info",
            "SF:src/Rodin/Example.h\nFNL:0,1,1\n"
            "FNA:0,7,answer<Real>\nFNA:0,0,answer<Complex>\nDA:1,7\nend_of_record\n")
        result = audit(self.root, [report], include_functions=True)
        self.assertEqual(result["coverage_percent"], 100)
        self.assertEqual(result["recorded_functions"], {"observed": 2, "hit": 1})
        self.assertEqual(result["zero_hit_functions"], [
            {"file": "src/Rodin/Example.h", "symbol": "answer<Complex>", "line": 1}])
        self.assertNotIn("zero_hit_functions", audit(self.root, [report]))

    def test_function_union_supports_old_and_new_lcov_and_names_with_commas(self):
        first = self.report("old.info", "SF:src/Rodin/Example.h\n"
            "FN:1,1,answer<Real,2>\nFNDA:0,answer<Real,2>\n"
            "FN:1,answer<Complex,2>\nFNDA:0,answer<Complex,2>\nDA:1,0\nend_of_record\n")
        second = self.report("new.info", f"SF:{self.source / 'Example.h'}\n"
            "FNL:17,1,1\nFNA:17,4,answer<Real,2>\nDA:1,4\nend_of_record\n")
        result = audit(self.root, [first, second], include_functions=True)
        self.assertEqual(result["lines"], 1)
        self.assertEqual(result["hits"], 1)
        self.assertEqual(result["recorded_functions"], {"observed": 2, "hit": 1})
        self.assertEqual(result["zero_hit_functions"][0]["symbol"], "answer<Complex,2>")
        self.assertEqual(result["zero_hit_functions"][0]["line"], 1)
        reverse = audit(self.root, [second, first], include_functions=True)
        self.assertEqual(reverse["recorded_functions"], result["recorded_functions"])
        self.assertEqual(reverse["zero_hit_functions"], result["zero_hit_functions"])

    def test_function_indices_and_names_are_scoped_to_the_production_file(self):
        (self.source / "Other.h").write_text("inline int answer() { return 7; }\n")
        report = self.report("scopes.info", "SF:src/Rodin/Example.h\n"
            "FNL:0,1,1\nFNA:0,2,answer\nDA:1,2\nend_of_record\n"
            "SF:src/Rodin/Other.h\nFNL:0,9,10\nFNA:0,0,answer\nDA:1,0\nend_of_record\n"
            "SF:tests/Test.cpp\nFN:1,answer\nFNDA:9,answer\nDA:1,9\nend_of_record\n")
        result = audit(self.root, [report], include_functions=True)
        self.assertEqual(result["recorded_functions"], {"observed": 2, "hit": 1})
        self.assertEqual(result["zero_hit_functions"], [
            {"file": "src/Rodin/Other.h", "symbol": "answer", "line": 9}])
        self.assertEqual(len(result["nonproduction_records"]), 1)

    def test_function_count_without_a_location_remains_an_audit_candidate(self):
        report = self.report("missing-location.info",
            "SF:src/Rodin/Example.h\nFNDA:0,answer\nDA:1,1\nend_of_record\n")
        result = audit(self.root, [report], include_functions=True)
        self.assertEqual(result["zero_hit_functions"][0]["line"], None)


if __name__ == "__main__":
    unittest.main()
