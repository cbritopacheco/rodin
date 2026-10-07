import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_p1_2d_parameter_campaign import save_case_trace
from summarize_iteration_logs import summarize


class IterationLoggingTest(unittest.TestCase):
    def test_inactive_inner_solves_are_counted_as_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            identity = dict(n=11, lobes=2)
            output = ("wngir geometry: outer=1 phase=accepted geom_sup=0.007 h=0.1 "
                      "min_j=0.2 max_qrel=9 inner_total=0 inner_last=0 "
                      "inner_converged=1 inner_residual=1e-12 seconds=1\n")
            path = save_case_trace(directory, identity, ["example"], output)
            row = summarize(path, {"2": 1})
            self.assertEqual(row["inner_mean"], 0)
            self.assertEqual(row["inner_max"], 0)
            self.assertEqual(row["inner_certified_outer"], 1)

    def trace(self):
        return "\n".join([
            "wngir geometry: outer=0 phase=initial geom_sup=0.1 h=0.1 min_j=1 max_qrel=1 inner_total=0 seconds=0.2",
            "barrier inner=1 outer=0 rel=0.1 alpha=0.5 linear_ok=1 cg_it=7",
            "barrier inner=2 outer=0 rel=0.0001 alpha=1 linear_ok=1 cg_it=5",
            "wngir geometry: outer=1 phase=accepted geom_sup=0.008 h=0.1 min_j=0.2 max_qrel=10 inner_total=2 seconds=1",
            "barrier inner=1 outer=1 rel=0.0001 alpha=1 linear_ok=1 cg_it=5",
            "wngir geometry: outer=2 phase=accepted geom_sup=0.007 h=0.1 min_j=0.2 max_qrel=9 inner_total=3 seconds=2",
            "barrier inner=1 outer=2 linear_ok=0 cg_it=1000",
            "wngir geometry: outer=2 phase=final geom_sup=0.007 h=0.1 min_j=0.2 max_qrel=9 inner_total=3 seconds=3",
        ]) + "\n"

    def test_trace_round_trip_and_first_admissible_hit(self):
        with tempfile.TemporaryDirectory() as directory:
            identity = dict(stage="screen", n=11, lobes=2, kappa_bulk=1e-4, rho=.1, mu_hat=.9)
            output = self.trace()
            path = save_case_trace(directory, identity, ["example", "--trace=1"], output)
            lines = path.read_text().splitlines(keepends=True)
            self.assertEqual(json.loads(lines[0].removeprefix("case=")), identity)
            self.assertEqual("".join(lines[2:]), output)
            row = summarize(path, {"2": 1})
            self.assertEqual(row["first_hit_outer"], 2)
            self.assertEqual(row["first_hit_inner_total"], 3)
            self.assertEqual(row["first_hit_seconds"], 2)
            self.assertEqual(row["inner_attempts"], 4)
            self.assertEqual(row["inner_max"], 2)
            self.assertEqual(row["inner_linear_failures"], 1)
            self.assertEqual(row["inner_damped_steps"], 1)
            self.assertEqual(row["geometry_records"], 3)
            row = summarize(path, {"2": .1})
            self.assertEqual(row["target_hit"], 0)
            self.assertEqual(row["first_hit_outer"], "")

    def test_initial_target_hit_needs_no_newton_steps(self):
        with tempfile.TemporaryDirectory() as directory:
            path = save_case_trace(directory, {"lobes": 0}, [], self.trace())
            self.assertEqual(summarize(path, {"0": 20})["first_hit_outer"], 0)

    def test_canonical_hinge_trace_matches_historical_barrier_trace(self):
        with tempfile.TemporaryDirectory() as directory:
            identity = {"lobes": 2}
            path = save_case_trace(directory, identity, [], self.trace())
            historical = summarize(path, {"2": 1})
            path = save_case_trace(directory, identity, [],
                                   self.trace().replace("barrier inner=", "hinge inner="))
            canonical = summarize(path, {"2": 1})
            self.assertEqual(canonical, historical)

    def test_invalid_reference_coefficient_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = save_case_trace(directory, {"lobes": 0}, [], self.trace())
            with self.assertRaises(ValueError):
                summarize(path, {"0": float("nan")})


if __name__ == "__main__":
    unittest.main()
