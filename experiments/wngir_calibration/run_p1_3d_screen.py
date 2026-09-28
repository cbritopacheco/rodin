#!/usr/bin/env python3
"""Resumable P1-3D classifier preflight and short-cap coefficient screen."""

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import re
import subprocess
import time
from pathlib import Path

from run_p1_2d_parameter_campaign import int_values, real_values


FIELDS = (
    "stage", "n", "elements", "lobes", "kappa_bulk", "rho", "mu_hat",
    "steps", "facets", "inside", "outside", "fit0", "fit", "geom_rms",
    "geom_sup", "normal_rms", "iterations", "min_j", "max_qrel",
    "active_rms_hg", "linear_iterations", "linear_solves", "linear_mean",
    "linear_max", "linear_error", "assembly",
    "solve", "exit", "seconds", "returncode",
)
PAIRS = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")
ELEMENTS = re.compile(r"^\s*elements=(\d+)", re.MULTILINE)
DEFAULT_KAPPA_BULK = (
    "0,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,"
    "1e-2,3e-2,0.1,0.2,0.5,1"
)


def parse_metrics(output, prefix):
    lines = [line.strip() for line in output.splitlines()
             if line.strip().startswith(prefix)]
    return dict(PAIRS.findall(lines[-1])) if lines else {}


def parse_number(values, name, fallback=float("nan")):
    return float(values[name]) if name in values else fallback


def mesh_samples_per_wave(n, lobes, r0, amp):
    if lobes == 0:
        return float("inf")
    h = 1 / (n - 1)
    return 2 * math.pi * (r0 - amp) / (lobes * h)


def command(args, n, lobes, kappa_bulk, rho, mu_hat, steps):
    h = 1 / (n - 1)
    return [
        str(args.exe), f"--n={n}", f"--lobes={lobes}",
        "--cx=0.5", "--cy=0.5", "--cz=0.5", "--phase=0",
        f"--amp={args.amp:.14g}", f"--R0={args.r0:.14g}",
        f"--classifier-eps={1.25*h:.14g}", "--classifier-lambda=0.008",
        f"--wngir-kappa-bulk={kappa_bulk:.14g}",
        f"--wngir-rigid-stabilisation={rho:.14g}",
        f"--wngir-mu-hat={mu_hat:.14g}",
        "--wngir-kappa-j=1", "--wngir-kappa-q=1", "--wngir-r-div=1",
        "--wngir-kappa-obs=1", "--wngir-robust-scale=0",
        "--wngir-jsafe=1e-2", "--wngir-qmax=10",
        f"--wngir-primal-barrier-iterations={args.barrier_max_iters}",
        "--wngir-primal-barrier-relative-tol=1e-3",
        "--wngir-theta-boundary=0.95", "--j-min=1e-8",
        "--wngir-jls=1e-2", "--wngir-armijo=1e-4",
        "--wngir-descent-fraction=1e-4",
        "--wngir-direction-norm-factor=10", "--wngir-alpha-min=1e-4",
        "--wngir-omega-min=0.1", "--wngir-rms-floor=0",
        "--wngir-sup-floor=0", "--wngir-rms-normal-jump-factor=0",
        "--wngir-sup-normal-jump-factor=0", "--wngir-rms-tol=1e-12",
        "--wngir-sup-tol=1e-12", "--wngir-energy-stag-tol=1e-8",
        f"--wngir-step-tol={1e-3*h*h:.14g}",
        f"--wngir-step-h-tol={1e-3*h:.14g}",
        f"--wngir-cg-rtol={args.cg_rtol:.14g}", "--wngir-cg-max-iters=1000",
        f"--wngir-steps={steps}", "--output=0", "--verbose",
    ]


def run_case(args, n, lobes, kappa_bulk, rho, mu_hat, steps):
    env = dict(os.environ)
    env.update({"OMP_NUM_THREADS": str(args.threads),
                "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
                "DYLD_LIBRARY_PATH": args.dyld_library_path})
    start = time.monotonic()
    try:
        proc = subprocess.run(
            command(args, n, lobes, kappa_bulk, rho, mu_hat, steps),
            cwd=args.scratch, env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, timeout=args.timeout)
        output, returncode = proc.stdout, proc.returncode
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or b""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        returncode = 124
    final = parse_metrics(output, "WNGIR it=")
    timing = parse_metrics(output, "wngir timing:")
    geometry = parse_metrics(output, "debug: facets=")
    elements = ELEMENTS.search(output)
    return {
        "stage": args.stage, "n": n,
        "elements": int(elements.group(1)) if elements else 6 * (n - 1) ** 3,
        "lobes": lobes, "kappa_bulk": kappa_bulk, "rho": rho,
        "mu_hat": mu_hat, "steps": steps,
        "facets": int(geometry.get("facets", -1)),
        "inside": int(geometry.get("inside", -1)),
        "outside": int(geometry.get("outside", -1)),
        "fit0": parse_number(geometry, "fit0"),
        "fit": parse_number(final, "fit"),
        "geom_rms": parse_number(final, "geom_rms"),
        "geom_sup": parse_number(final, "geom_sup"),
        "normal_rms": parse_number(final, "normal_rms"),
        "iterations": int(final.get("it", -1)),
        "min_j": parse_number(final, "min_j"),
        "max_qrel": parse_number(final, "max_qrel"),
        "active_rms_hg": parse_number(final, "active_rms_hg"),
        "linear_iterations": int(timing.get("cgIt", -1)),
        "linear_solves": int(timing.get("cgSolves", -1)),
        "linear_mean": parse_number(timing, "cgMean"),
        "linear_max": int(timing.get("cgMax", -1)),
        "linear_error": parse_number(timing, "cgErr"),
        "assembly": parse_number(timing, "assembly"),
        "solve": parse_number(timing, "solve"),
        "exit": "timeout" if returncode == 124 else final.get("exit", "parse-fail"),
        "seconds": time.monotonic() - start, "returncode": returncode,
    }


def read_rows(path):
    if not path.exists():
        return []
    with path.open(newline="") as source:
        reader = csv.DictReader(source)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError(f"unexpected CSV schema: {path}")
        return list(reader)


def key(row):
    return (row["stage"], int(row["n"]), int(row["lobes"]),
            round(float(row["kappa_bulk"]), 14),
            round(float(row["rho"]), 14),
            round(float(row["mu_hat"]), 14))


def append_row(path, row):
    exists = path.exists()
    with path.open("a", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=FIELDS)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("preflight", "screen"))
    parser.add_argument("--n", required=True, help="comma-separated grid-point counts")
    parser.add_argument("--lobes", help="default: all for preflight, 0,5,10 for screen")
    parser.add_argument("--kappa-bulk", default=DEFAULT_KAPPA_BULK)
    parser.add_argument("--rho", default="0.1:1.0:0.1")
    parser.add_argument("--mu-hat", default="0.1:1.0:0.1")
    parser.add_argument("--steps", type=int, default=3)
    parser.add_argument("--barrier-max-iters", type=int, default=10,
                        help="inner Newton cap; use 25 to audit cap failures")
    parser.add_argument("--cg-rtol", type=float, default=1e-6,
                        help="CG relative tolerance; use 1e-9 for finalists")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--timeout", type=int, default=60,
                        help="seconds per case")
    parser.add_argument("--max-cases", type=int, default=0,
                        help="limit new cases for a smoke run; 0 means all")
    parser.add_argument("--min-samples-per-wave", type=float, default=4)
    parser.add_argument("--allow-underresolved", action="store_true")
    parser.add_argument("--amp", type=float, default=0.08)
    parser.add_argument("--r0", type=float, default=0.24)
    parser.add_argument("--exe", type=Path,
                        default=Path("build-p1-3d-clang19/examples/Geometry/"
                                     "LevelSetWNGIRReconstruction3D"))
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--scratch", type=Path,
                        default=Path("/tmp/wngir-p1-3d-screen"))
    parser.add_argument("--dyld-library-path", default="/opt/local/lib")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    args.exe = (root / args.exe).resolve() if not args.exe.is_absolute() else args.exe
    args.out_dir = ((root / args.out_dir).resolve() if not args.out_dir.is_absolute()
                    else args.out_dir)
    ns = int_values(args.n)
    lobes = int_values(args.lobes or
                       ("0:10" if args.stage == "preflight" else "0,5,10"))
    if not ns or any(n < 2 for n in ns) or not lobes or any(l < 0 for l in lobes):
        parser.error("n must be >=2 and lobes must be nonnegative")
    if args.stage == "screen" and not 1 <= args.steps <= 20:
        parser.error("screen steps must be between 1 and 20")
    if args.barrier_max_iters < 1:
        parser.error("barrier-max-iters must be positive")
    if args.threads < 1 or args.timeout < 1 or args.max_cases < 0:
        parser.error("threads and timeout must be positive; max-cases nonnegative")
    if not 0 < args.cg_rtol < 1:
        parser.error("cg-rtol must lie strictly between 0 and 1")
    if not 0 < args.amp < args.r0:
        parser.error("the target requires 0 < amp < R0")
    if not args.exe.is_file():
        parser.error(f"executable not found: {args.exe}")
    controls = (real_values(args.kappa_bulk), real_values(args.rho),
                real_values(args.mu_hat))
    if any(not values for values in controls):
        parser.error("coefficient grids must be nonempty")
    values = [(1e-4, 0.1, 0.9)] if args.stage == "preflight" else list(
        itertools.product(*controls))
    cases = [(n, lobe, *triple) for n in ns for lobe in lobes for triple in values]
    print(f"{args.stage}: {len(cases)} cases; executable={args.exe}", flush=True)
    for n in ns:
        print(f"n={n}: {6*(n-1)**3} tetrahedra; minimum samples/wave="
              f"{min(mesh_samples_per_wave(n, l, args.r0, args.amp) for l in lobes):.2f}",
              flush=True)
    if args.dry_run:
        return
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.scratch.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / f"{args.stage}.csv"
    manifest_path = args.out_dir / f"{args.stage}_manifest.json"
    digest = hashlib.sha256(args.exe.read_bytes()).hexdigest()
    if args.stage == "screen":
        preflight_manifest_path = args.out_dir / "preflight_manifest.json"
        if not preflight_manifest_path.exists():
            parser.error("run preflight in this out-dir before screening")
        preflight_manifest = json.loads(preflight_manifest_path.read_text())
        if (preflight_manifest["sha256"] != digest or
                preflight_manifest["amp"] != args.amp or
                preflight_manifest["r0"] != args.r0):
            parser.error("preflight used a different binary or target")
    manifest = {"stage": args.stage, "n": ns, "lobes": lobes,
                "kappa_bulk": controls[0] if args.stage == "screen" else [1e-4],
                "rho": controls[1] if args.stage == "screen" else [0.1],
                "mu_hat": controls[2] if args.stage == "screen" else [0.9],
                "steps": 0 if args.stage == "preflight" else args.steps,
                "target": "R0 + A/3 sum_i cos(lobes*n_i)",
                "amp": args.amp, "r0": args.r0,
                "classifier": {"epsilon_over_h": 1.25, "lambda_c": 0.008},
                "barrier": {"relative_correction": 1e-3,
                            "max_iterations": args.barrier_max_iters},
                "stopping": {"rms": 1e-12, "sup": 1e-12,
                             "energy_relative": 1e-8,
                             "step_absolute_over_h2": 1e-3,
                             "accepted_step_over_h2": 1e-3,
                             "consecutive_small_steps": 5},
                "linear": {"backend": "CG", "relative_tolerance": args.cg_rtol,
                           "max_iterations_per_solve": 1000},
                "executable": str(args.exe), "sha256": digest}
    if manifest_path.exists():
        if json.loads(manifest_path.read_text()) != manifest:
            parser.error("manifest differs from existing dataset; use a new out-dir")
    else:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    if args.stage == "screen":
        preflight = {(int(row["n"]), int(row["lobes"])): row for row in
                     read_rows(args.out_dir / "preflight.csv")}
        for n in ns:
            for lobe in lobes:
                row = preflight.get((n, lobe))
                if row is None or int(row["facets"]) <= 0 or int(row["returncode"]):
                    parser.error(f"missing or invalid preflight: n={n}, lobes={lobe}")
                samples = mesh_samples_per_wave(n, lobe, args.r0, args.amp)
                if samples < args.min_samples_per_wave and not args.allow_underresolved:
                    parser.error(f"n={n}, lobes={lobe}: {samples:.2f} samples/wave "
                                 "is under-resolved; refine or explicitly override")
    done = {key(row) for row in read_rows(csv_path)}
    new = 0
    for n, lobe, kb, rho, mu in cases:
        case_key = (args.stage, n, lobe, round(kb, 14), round(rho, 14),
                    round(mu, 14))
        if case_key in done:
            continue
        if args.max_cases and new >= args.max_cases:
            break
        row = run_case(args, n, lobe, kb, rho, mu,
                       0 if args.stage == "preflight" else args.steps)
        append_row(csv_path, row)
        new += 1
        print(f"{new} n={n} lobes={lobe} kb={kb:g} rho={rho:g} mu={mu:g} "
              f"facets={row['facets']} D={row['geom_rms']:.4g} "
              f"it={row['iterations']} exit={row['exit']} sec={row['seconds']:.2f}",
              flush=True)


if __name__ == "__main__":
    main()
