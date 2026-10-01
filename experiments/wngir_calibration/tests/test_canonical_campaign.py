import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_p1_2d_parameter_campaign import ensure_schema, run_case
from run_p1_3d_screen import command


class CanonicalCampaignTest(unittest.TestCase):
    def test_2d_command_uses_regularization_not_rigid_lift(self):
        args = SimpleNamespace(amp=.08, r0=.24, steps=30, barrier_max_iters=15,
                               extra="", log_iterations=False, threads=4,
                               dyld_library_path="", root=Path("/tmp"))
        with patch("run_p1_2d_parameter_campaign.subprocess.run") as run:
            run.return_value = SimpleNamespace(stdout="", returncode=0)
            run_case(args, Path("/tmp/example"), "canonical", 20, 4,
                     1e-4, 1, 90, 1, 1)
        self.check_command(run.call_args.args[0])
        self.assertEqual(run.call_args.kwargs["env"]["OPENBLAS_NUM_THREADS"], "4")

    def test_3d_command_uses_same_model(self):
        args = SimpleNamespace(exe=Path("/tmp/example"), amp=.08, r0=.24,
                               barrier_max_iters=15, cg_rtol=1e-8,
                               log_iterations=False)
        self.check_command(command(args, 20, 4, 1e-4, 1, 90, 30))

    def check_command(self, args):
        self.assertIn("--wngir-kappa-c=1", args)
        self.assertIn("--wngir-kappa-reg=1", args)
        self.assertIn("--wngir-direct-solver=mumps", args)
        self.assertIn("--wngir-primal-barrier-iterations=15", args)
        self.assertIn("--wngir-steps=30", args)
        self.assertFalse(any("rigid-stabilisation" in arg or "quality-model" in arg
                             or "theta-boundary" in arg or "r-div" in arg
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
