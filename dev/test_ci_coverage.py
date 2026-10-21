"""Regression guards for complete, production-only incremental coverage."""

import importlib.util
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("ci_coverage", Path(__file__).with_name("ci_coverage.py"))
coverage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(coverage)


class ProductionCoverageTest(unittest.TestCase):
    def test_zero_hit_lines_survive_incremental_union(self):
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "initial.info"
            second = Path(directory) / "runtime.info"
            first.write_text("SF:/checkout/src/Rodin/Math/Foo.h\nDA:1,0\nDA:2,0\nDA:3,0\nend_of_record\n")
            second.write_text("SF:src/Rodin/Math/Foo.h\nDA:2,3\nend_of_record\n")
            lines = coverage.production_lines([first, second])
            self.assertEqual(len(lines), 3)
            self.assertEqual(sum(lines.values()), 1)

    def test_test_sources_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "test.info"
            for source in ["tests/unit/foo.cpp", "src/Rodin/Test/Random.h", "third-party/foo.cpp"]:
                report.write_text("SF:" + source + "\nDA:1,1\n")
                with self.assertRaisesRegex(ValueError, "Non-production"):
                    coverage.production_lines([report])

    def test_missing_stale_and_uncompiled_producers_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            definitions = []
            receipts = []
            for profile, build_types in [("unit", ["Debug", "Release"]), ("coverage", ["Release"])]:
                for build_type in build_types:
                    for mt in ["ON", "OFF"]:
                        identity = "-".join([profile, build_type, mt])
                        definition = {"id": identity, "profile": profile, "build_type": build_type,
                                      "multithreaded": mt, "asan": "OFF", "ubsan": "OFF",
                                      "baseline": True, "group": "src/Rodin"}
                        definitions.append(definition)
                        path = folder / identity
                        path.mkdir()
                        receipt = {"definition": definition, "commit": "current", "build": "success",
                                   "tests": "skipped", "coverage_expected": True, "uploaded": True,
                                   "capture": "success", "report_exists": True,
                                   "configured": [{"target": "Library"}, {"target": "TemplateTest"}],
                                   "compiled": [{"target": "Library"}, {"target": "TemplateTest"}]}
                        report = path / "RodinCoverage.info"
                        report.write_text("SF:src/Rodin/Math/Foo.h\nDA:1,0\nend_of_record\n")
                        receipt["report_sha256"] = hashlib.sha256(report.read_bytes()).hexdigest()
                        file = path / "receipt.json"
                        file.write_text(json.dumps(receipt))
                        receipts.append((file, receipt))
            plan = {"baseline": definitions}
            self.assertEqual(len(coverage.verify_receipts(plan, folder, "current")), 6)
            file, receipt = receipts[0]
            report = file.parent / "RodinCoverage.info"
            contents = report.read_bytes()
            report.unlink()
            with self.assertRaisesRegex(ValueError, "upload incomplete"):
                coverage.verify_receipts(plan, folder, "current")
            report.write_bytes(contents + b"# truncated or replaced artifact\n")
            with self.assertRaisesRegex(ValueError, "upload incomplete"):
                coverage.verify_receipts(plan, folder, "current")
            report.write_bytes(contents)
            for field, value, error in [("commit", "old", "Stale"),
                                        ("uploaded", False, "upload incomplete"),
                                        ("compiled", [{"target": "Library"}], "compilation targets")]:
                changed = dict(receipt, **{field: value})
                file.write_text(json.dumps(changed))
                with self.assertRaisesRegex(ValueError, error):
                    coverage.verify_receipts(plan, folder, "current")
            file.unlink()
            with self.assertRaisesRegex(ValueError, "Missing producer"):
                coverage.verify_receipts(plan, folder, "current")


if __name__ == "__main__":
    unittest.main()
