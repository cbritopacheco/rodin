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


if __name__ == "__main__":
    unittest.main()
