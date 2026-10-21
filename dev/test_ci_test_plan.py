"""Preserve configuration coverage while changing test scheduling groups."""

import importlib.util
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent))
from ci_test_plan import plan, targeted_plan, native_convergence_groups


class TestPlanTest(unittest.TestCase):
    def test_configuration_matrix_and_job_limits(self):
        matrices = plan()
        rows = [row for entries in matrices.values() for row in entries]
        self.assertEqual(len({row["id"] for row in rows}), len(rows))
        for entries in matrices.values():
            self.assertGreater(len(entries), 0)
            self.assertLessEqual(len(entries), 256)
        required = set(itertools.product(["Debug", "Release"], ["ON", "OFF"], ["ON", "OFF"], ["ON", "OFF"]))
        for profile in ["unit", "manufactured"]:
            selected = [row for row in rows if row["profile"] == profile and not row["baseline"]]
            for group in {row["group"] for row in selected}:
                actual = {(r["build_type"], r["multithreaded"], r["asan"], r["ubsan"])
                          for r in selected if r["group"] == group and r["backend"] != "baseline"}
                self.assertEqual(actual, required)
        baseline = [r for r in rows if r["baseline"] and r["group"] == "src/Rodin"]
        self.assertEqual(len(baseline), 10)
        self.assertNotIn("coverage", matrices)
        self.assertFalse(any(r["profile"] == "coverage" and not r["baseline"] for r in rows))
        self.assertTrue(any(r["group"] == "cmake/ThirdPartyIsolation" for r in rows))

    def test_old_coverage_configuration_moves_to_ordinary_lanes(self):
        rows = [r for entries in plan().values() for r in entries]
        for profile in ["unit", "manufactured"]:
            normal = [r for r in rows if r["profile"] == profile and not r["baseline"]]
            for group in {r["group"] for r in normal if not r["group"].startswith("cmake/")}:
                baseline = [r for r in normal if r["group"] == group and r["backend"] == "baseline"]
                self.assertEqual({(r["build_type"], r["multithreaded"], r["asan"], r["ubsan"])
                                  for r in baseline}, {("Release", mt, "OFF", "OFF") for mt in ["ON", "OFF"]})

    def test_every_targeted_execution_has_independent_collection(self):
        import yaml
        from ci_tests import ROOT
        jobs = yaml.safe_load((ROOT / ".github/workflows/Tests.yml").read_text())["jobs"]
        native = [{"group": group, "multithreaded": mt} for group in native_convergence_groups()
                  for mt in ["ON", "OFF"]]
        definitions = targeted_plan(native)
        self.assertNotIn("Coverage", jobs)
        for name, job in jobs.items():
            steps = job.get("steps", [])
            if not any("cmake -S . -B build" in s.get("run", "") for s in steps):
                continue
            self.assertTrue(any(d["job"] == name for d in definitions))
            for step in steps:
                command = step.get("run", "")
                if "cmake -S . -B build" in command:
                    self.assertIn("-DRODIN_CODE_COVERAGE=ON", command)
                if "apt-get install" in command:
                    self.assertIn("lcov", command)
            publisher = next(s for s in steps if s.get("id") == "coverage_record")
            self.assertEqual(publisher["uses"], "./.github/actions/coverage")
        action = yaml.safe_load((ROOT / ".github/actions/coverage/action.yml").read_text())
        ids = [s.get("id") for s in action["runs"]["steps"]]
        self.assertLess(ids.index("coverage_capture"), ids.index("coverage_upload"))
        self.assertLess(ids.index("coverage_upload"), ids.index("coverage_receipt"))


if __name__ == "__main__":
    unittest.main()
