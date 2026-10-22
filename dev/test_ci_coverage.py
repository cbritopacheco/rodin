"""Regression guards for complete, production-only incremental coverage."""

import importlib.util
import hashlib
import json
import shutil
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

    def test_shared_build_lanes_pool_complementary_hits_without_duplicate_denominators(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            initial = folder / "initial.info"
            first = folder / "first.info"
            second = folder / "second.info"
            initial.write_text("SF:/builder/src/Rodin/Math/Foo.h\nDA:1,0\nDA:2,0\nDA:3,0\nDA:4,0\nend_of_record\n")
            first.write_text("SF:/runner/first/src/Rodin/Math/Foo.h\nDA:1,7\nDA:2,0\nDA:3,0\nDA:4,0\nend_of_record\n")
            second.write_text("SF:/runner/second/src/Rodin/Math/Foo.h\nDA:1,0\nDA:2,9\nDA:3,0\nDA:4,0\nend_of_record\n")
            first_lines = coverage.production_lines([first])
            second_lines = coverage.production_lines([second])
            union = coverage.production_lines([initial, first, second])
            self.assertEqual(len(union), 4)
            self.assertEqual(sum(first_lines.values()), 1)
            self.assertEqual(sum(second_lines.values()), 1)
            self.assertEqual(sum(union.values()), 2)
            self.assertFalse(union[("src/Rodin/Math/Foo.h", 3)])
            self.assertFalse(union[("src/Rodin/Math/Foo.h", 4)])
            self.assertEqual(coverage.production_lines([second, initial, first, first]), union)

    def test_test_sources_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "test.info"
            for source in ["tests/unit/foo.cpp", "tests/benchmarks/Rodin/Foo.cpp", "benchmarks/Foo.cpp",
                           "src/Rodin/Test/Random.h", "third-party/foo.cpp"]:
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

    def test_targeted_receipts_require_exact_commit_partition_and_report(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            definition = {"job": "Convergence", "matrix": {"multithreaded": "ON", "suite": "diffusion"}}
            expected = [definition]
            with self.assertRaisesRegex(ValueError, "Missing targeted"):
                coverage.verify_target_receipts(expected, folder, "current")
            report = folder / "RodinCoverage.info"
            report.write_text("SF:src/Rodin/Math/Foo.h\nDA:1,1\nend_of_record\n")
            receipt = {"definition": definition, "commit": "current",
                       "report_sha256": hashlib.sha256(report.read_bytes()).hexdigest()}
            path = folder / "target-receipt.json"
            path.write_text(json.dumps(receipt))
            self.assertEqual(len(coverage.verify_target_receipts(expected, folder, "current")), 1)
            with self.assertRaisesRegex(ValueError, "Stale targeted"):
                coverage.verify_target_receipts(expected, folder, "other")
            with self.assertRaisesRegex(ValueError, "Unexpected"):
                coverage.verify_target_receipts([{**definition, "matrix": {"multithreaded": "OFF"}}], folder, "current")
            report.write_text("changed")
            with self.assertRaisesRegex(ValueError, "changed targeted"):
                coverage.verify_target_receipts(expected, folder, "current")

    def test_shared_executable_keeps_distinct_execution_receipts(self):
        import ci_builds
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            definitions = []
            build_keys = []
            for rank in [1, 2]:
                row = {"multithreaded": "ON", "np": str(rank), "partition": "standard"}
                build_keys.append(ci_builds.identity(ci_builds.definition("MPIH1Restriction", row)))
                definition = {"job": "MPIH1Restriction", "matrix": row}
                definitions.append(definition)
                lane = folder / str(rank)
                lane.mkdir()
                report = lane / "RodinCoverage.info"
                report.write_text("SF:src/Rodin/Math/Foo.h\nDA:1," + str(int(rank == 1))
                                  + "\nDA:2," + str(int(rank == 2)) + "\nend_of_record\n")
                receipt = {"definition": definition, "commit": "current",
                           "report_sha256": hashlib.sha256(report.read_bytes()).hexdigest()}
                (lane / "target-receipt.json").write_text(json.dumps(receipt))
            self.assertEqual(build_keys[0], build_keys[1])
            self.assertNotEqual(coverage.target_identity(definitions[0]), coverage.target_identity(definitions[1]))
            self.assertEqual(len(coverage.verify_target_receipts(definitions, folder, "current")), 2)
            pooled = coverage.production_lines(list(folder.rglob("RodinCoverage.info")))
            self.assertEqual((sum(pooled.values()), len(pooled)), (2, 2))
            shutil.copytree(folder / "2", folder / "duplicate")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                coverage.verify_target_receipts(definitions, folder, "current")
            shutil.rmtree(folder / "duplicate")
            (folder / "2/target-receipt.json").unlink()
            with self.assertRaisesRegex(ValueError, "Missing targeted"):
                coverage.verify_target_receipts(definitions, folder, "current")

    def test_backend_context_combines_execution_lanes_without_combining_build_variants(self):
        definition = {"profile": "unit", "backend": "baseline", "build_type": "Release",
                      "multithreaded": "ON", "asan": "OFF", "ubsan": "OFF"}
        self.assertEqual(coverage.compilation_context(definition),
                         coverage.compilation_context({**definition, "profile": "manufactured"}))
        for field, value in [("backend", "optional"), ("multithreaded", "OFF"), ("build_type", "Debug")]:
            self.assertNotEqual(coverage.compilation_context(definition),
                                coverage.compilation_context({**definition, field: value}))


if __name__ == "__main__":
    unittest.main()
