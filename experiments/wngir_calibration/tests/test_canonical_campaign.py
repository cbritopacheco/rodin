import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_p1_2d_parameter_campaign import FIELDS, cases_for, ensure_schema, main, read_done, run_case
from run_p1_3d_screen import command


class CanonicalCampaignTest(unittest.TestCase):
    def test_2d_command_uses_regularization_not_rigid_lift(self):
        args = SimpleNamespace(kappa_f=1, amp=.08, r0=.24, steps=30, barrier_max_iters=15,
                               extra="", log_iterations=False, threads=4,
                               dyld_library_path="", root=Path("/tmp"))
        with patch("run_p1_2d_parameter_campaign.subprocess.run") as run:
            run.return_value = SimpleNamespace(stdout="", returncode=0)
            run_case(args, Path("/tmp/example"), "canonical", 20, 4,
                     1, 1e-4, 1, 90, 1, 1)
        self.check_command(run.call_args.args[0])
        self.assertEqual(run.call_args.kwargs["env"]["OPENBLAS_NUM_THREADS"], "4")

    def test_fitting_grid_and_resume_keys_are_independent(self):
        grid = [1e-4, 1e-3, 1e-2, .1, 1]
        cases = list(cases_for("canonical", [20], [4], grid, grid, grid,
                               [.1, 1, 10, 100, 1000], [1], [1]))
        self.assertEqual(len(cases), 625)
        self.assertEqual(len(set(cases)), 625)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "campaign.csv"
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=FIELDS)
                writer.writeheader()
                for kf in grid:
                    writer.writerow(dict(dataset="canonical", n=20, lobes=4,
                                         kappa_f=kf, kappa_s=.001, kappa_d=.001,
                                         shape_curvature="full", mu_hat=.1, kappa_j=1, kappa_q=1))
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
            self.assertEqual(manifest["expected_cases"], 68750)
            self.assertEqual(manifest["shape_curvature"], ["full", "psd"])
            for key in ("kappa_f", "kappa_s", "kappa_d"):
                self.assertEqual(manifest[key], [1e-4, 1e-3, 1e-2, .1, 1])

    def test_shape_variants_have_separate_cases_and_commands(self):
        cases = list(cases_for("canonical", [20], [4], [1], [.001], [.001],
                               [.1], [1], [1], ["full", "psd"]))
        self.assertEqual(len(set(cases)), 2)
        args = SimpleNamespace(amp=.08, r0=.24, steps=30, barrier_max_iters=15,
                               extra="", log_iterations=False, threads=4,
                               dyld_library_path="", root=Path("/tmp"))
        for case, expected in zip(cases, (0, 1)):
            with patch("run_p1_2d_parameter_campaign.subprocess.run") as run:
                run.return_value = SimpleNamespace(stdout="", returncode=0)
                row = run_case(args, Path("/tmp/example"), *case)
                self.assertIn(f"--wngir-positive-shape-curvature={expected}", run.call_args.args[0])
                self.assertEqual(row["shape_curvature"], case[-1])

    def test_3d_command_uses_same_model(self):
        args = SimpleNamespace(exe=Path("/tmp/example"), kappa_f=1, amp=.08, r0=.24,
                               barrier_max_iters=15, cg_rtol=1e-8,
                               log_iterations=False)
        self.check_command(command(args, 20, 4, 1e-4, 1, 90, 30))

    def test_three_metric_coefficients_are_independent(self):
        args = SimpleNamespace(exe=Path("/tmp/example"), kappa_f=2, amp=.08, r0=.24,
                               barrier_max_iters=15, cg_rtol=1e-8,
                               log_iterations=False)
        result = command(args, 20, 4, 3, 5, 90, 30)
        for option in ("--wngir-kappa-f=2", "--wngir-kappa-s=3", "--wngir-kappa-d=5"):
            self.assertIn(option, result)

    def check_command(self, args):
        self.assertIn("--wngir-kappa-f=1", args)
        self.assertIn("--wngir-kappa-s=0.0001", args)
        self.assertIn("--wngir-kappa-d=1", args)
        self.assertIn("--wngir-direct-solver=mumps", args)
        self.assertIn("--wngir-primal-barrier-iterations=15", args)
        self.assertIn("--wngir-steps=30", args)
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
