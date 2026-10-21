"""Preserve configuration coverage while changing test scheduling groups."""

import importlib.util
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent))
from ci_test_plan import plan


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
                          for r in selected if r["group"] == group}
                self.assertEqual(actual, required)
        baseline = [r for r in rows if r["baseline"] and r["group"] == "src/Rodin"]
        self.assertEqual(len(baseline), 6)
        self.assertTrue(any(r["group"] == "cmake/ThirdPartyIsolation" for r in rows))


if __name__ == "__main__":
    unittest.main()
