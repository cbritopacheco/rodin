import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_p1_2d_parameter_campaign import FIELDS, cases_for, ensure_schema, main, read_done, run_case, parse_responses
from run_p1_3d_screen import command


class CanonicalCampaignTest(unittest.TestCase):
    def test_2d_command_uses_regularization_not_rigid_lift(self):
        args = SimpleNamespace(kappa_f=1, amp=.08, r0=.24, steps=30, barrier_max_iters=15,
                               extra="", log_iterations=False, threads=4,
                               dyld_library_path="", root=Path("/tmp"))
        with patch("run_p1_2d_parameter_campaign.subprocess.run") as run:
            run.return_value = SimpleNamespace(stdout="", returncode=0)
            run_case(args, Path("/tmp/example"), "canonical", 20, 4,
                     1, 1, 90, 1, 1)
        self.check_command(run.call_args.args[0])
        self.assertEqual(run.call_args.kwargs["env"]["OPENBLAS_NUM_THREADS"], "4")

    def test_fitting_grid_and_resume_keys_are_independent(self):
        grid = [1e-4, 1e-3, 1e-2, .1, 1]
        cases = list(cases_for("canonical", [20], [4], grid, grid,
                               [.1, 1, 10, 100, 1000], [1], [1]))
        self.assertEqual(len(cases), 125)
        self.assertEqual(len(set(cases)), 125)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "campaign.csv"
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=FIELDS)
                writer.writeheader()
                for kf in grid:
                    writer.writerow(dict(dataset="canonical", n=20, lobes=4,
                                         kappa_f=kf, kappa_d=.001,
                                         mu_hat=.1, kappa_j=1, kappa_q=1))
            _, done = read_done(path)
            self.assertEqual(len(done), 5)
            self.assertTrue(done.issubset(set(cases)))

    def test_extended_grid_manifest_without_launching_cases(self):
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "example"
            executable.touch()
            argv = ["campaign", "--root", directory, "--exe", str(executable),
                    "--out-dir", directory, "--scratch", directory,
                    "--mu-hat", "0.1,1,10,100,1000",
                    "--deadline", "2000-01-01 00:00:00"]
            with patch.object(sys, "argv", argv), patch(
                    "run_p1_2d_parameter_campaign.run_case") as run:
                main()
                run.assert_not_called()
            manifest = json.loads((Path(directory) / "canonical_p1_2d_manifest.json").read_text())
            self.assertEqual(manifest["expected_cases"], 6875)
            self.assertNotIn("shape_curvature", manifest)
            self.assertEqual(manifest["model"], "F+D-centered-strain-surface12-v10")
            profile = manifest["fixed_profile"]
            self.assertEqual(profile["model"], manifest["model"])
            self.assertEqual(profile["interface_quadrature_order"],
                             "automatic: max(12, 2 * FE order + 2)")
            self.assertEqual(profile["geometric_validation_order"],
                             "automatic: max(14, 2 * FE order + 4)")
            self.assertNotIn("shape_restriction", profile)
            for key in ("kappa_f", "kappa_d"):
                self.assertEqual(manifest[key], [1e-4, 1e-3, 1e-2, .1, 1])

    def test_removed_shape_option_is_rejected(self):
        with patch.object(sys, "argv", ["campaign", "--out-dir", "/tmp",
                                      "--shape-curvature", "full"]):
            with self.assertRaises(SystemExit):
                main()

    def test_primary_responses_preserve_precision_and_inner_counts(self):
        output = ("swift responses: energy=1.2345678901234567e-5 "
                  "geom_sup=0.0012345678901234567 geom_c=1.2345678901234567 geom_sup_target=0.002 "
                  "target_hit=1 quality_ok=1 inner_total=19 inner_max=4 "
                  "inner_last=2 inner_converged=1 inner_residual=1e-9 "
                  "inner_relative_residual=1e-4 inner_residual_tolerance=1e-8 "
                  "min_j=0.0100000000123 max_qrel=9.999999999987")
        fields = parse_responses(output)
        self.assertEqual(fields["geom_sup"], .0012345678901234567)
        self.assertEqual(fields["geom_c"], 1.2345678901234567)
        self.assertEqual(fields["inner_total"], 19)
        self.assertEqual(fields["inner_max"], 4)
        self.assertLessEqual(fields["inner_residual"], fields["inner_residual_tolerance"])
        self.assertLess(fields["max_qrel"], 10)

    def test_3d_command_uses_same_model(self):
        args = SimpleNamespace(exe=Path("/tmp/example"), kappa_f=1, amp=.08, r0=.24,
                               barrier_max_iters=15, cg_rtol=1e-8,
                               log_iterations=False)
        self.check_command(command(args, 20, 4, 1, 90, 30))

    def test_two_metric_coefficients_are_independent(self):
        args = SimpleNamespace(exe=Path("/tmp/example"), kappa_f=2, amp=.08, r0=.24,
                               barrier_max_iters=15, cg_rtol=1e-8,
                               log_iterations=False)
        result = command(args, 20, 4, 5, 90, 30)
        for option in ("--model-fit=2", "--model-distribution-deviatoric=5"):
            self.assertIn(option, result)

    def check_command(self, args):
        self.assertIn("--model-fit=1", args)
        self.assertFalse(any("kappa-s" in arg for arg in args))
        self.assertIn("--model-distribution-deviatoric=1", args)
        self.assertIn("--model-distribution-divergence=1", args)
        self.assertFalse(any("omega-min" in arg or "jls" in arg or "volume-gauge" in arg
                             for arg in args))
        self.assertIn("--linear-solver=mumps", args)
        self.assertIn("--convergence-iterations-inner=15", args)
        self.assertIn("--convergence-iterations-outer=30", args)
        self.assertIn("--convergence-tolerance-step=0", args)
        self.assertIn("--convergence-tolerance-step-over-h=5e-4", args)
        self.assertIn("--convergence-iterations-stagnation=5", args)
        self.assertIn("--globalization-max-step-over-h=0", args)
        self.assertTrue(any(arg.startswith("--convergence-tolerance-geometric=") for arg in args))
        self.assertFalse(any("rms-tol" in arg or "rms-floor" in arg or "descent-fraction" in arg
                             or "positive-shape-curvature" in arg for arg in args))
        self.assertFalse(any("rigid-stabilisation" in arg or "quality-model" in arg
                             or "theta-boundary" in arg or "r-div" in arg or "kappa-bulk" in arg
                             or "kappa-obs" in arg or "kappa-reg" in arg
                             for arg in args))

    def test_old_schema_is_rejected_without_modification(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "old.csv"
            with path.open("w", newline="") as stream:
                csv.writer(stream).writerows([["rho"], ["0.1"]])
            original = path.read_bytes()
            with self.assertRaises(ValueError):
                ensure_schema(path)
            self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
