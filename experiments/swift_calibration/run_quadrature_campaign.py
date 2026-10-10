#!/usr/bin/env python3
"""Serial quadrature screening with frozen binaries and independent validation."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary-2d", type=Path, required=True)
    parser.add_argument("--binary-3d", type=Path, required=True)
    parser.add_argument("--binary-p2", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--finite-element-only", action="store_true")
    args = parser.parse_args()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / "manifest.json").exists():
        raise SystemExit("Use a fresh directory; existing campaigns are never overwritten.")
    binaries = {}
    for family, source in (("2d", args.binary_2d), ("3d", args.binary_3d),
                           ("p2", args.binary_p2)):
        destination = out / ("solver-" + family)
        shutil.copy2(source.resolve(), destination)
        binaries[family] = destination
    # The 2D experimental binary selects interpolated targets through this
    # environment variable. Classification stays analytic and identical.
    cases = [("2d", n, degree) for n in (8, 16) for degree in (1, 2, 3)]
    if not args.finite_element_only:
        cases += [("2d", 32, 0), ("3d", 8, 0), ("p2", 8, 0)]
    orders = (2, 4, 6, 8, 12, 16, 24, 32)
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       VECLIB_MAXIMUM_THREADS="1")
    manifest = {
        "cases": cases, "orders": orders, "attempts": len(cases) * len(orders),
        "validation_order": 32, "outer_cap": 30, "inner_cap": 15,
        "threads": 1, "timeout_seconds": args.timeout,
        "binary_sha256": {key: hashlib.sha256(value.read_bytes()).hexdigest()
                          for key, value in binaries.items()},
        "protocol": "No warmups/repeats. One process at a time. Final solver "
                    "responses, not the example's subsequently selected best-RMS mesh. "
                    "Order overrides both surface and volume integration. FE target "
                    "degree is independent of P1 displacement degree. Analytic "
                    "classification remains fixed. Validation32 plus vertices is "
                    "sampling, not a certified supremum. Energy rules differ; energy "
                    "is not a common-state integration-error comparison.",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    fields = ["family", "n", "target_degree", "order", "returncode", "seconds",
              "energy", "geom_sup", "geom_c", "target_hit", "quality_ok", "outer",
              "inner_total", "inner_max", "inner_converged", "min_j", "max_qrel", "exit",
              "dense_min_j", "dense_max_qrel", "dense_quality_ok"]
    rows = []
    with (out / "results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for family, n, degree in cases:
            # High-order references come first so a partial campaign is useful.
            for order in (32, 24, *orders[:-2]):
                identity = f"{family}-n{n}-target{degree}-q{order}"
                directory = out / identity
                directory.mkdir()
                command = [str(binaries[family]), f"--n={n}", "--frames=1", "--lobes=4",
                           "--amp=0.08", "--R0=0.24", "--orbitR=0",
                           "--convergence-iterations-outer=30", "--convergence-iterations-inner=15",
                           "--model-fit=1", "--model-distribution-deviatoric=0.0001",
                           "--model-distribution-divergence=0.01", "--model-hinge=10",
                           "--linear-solver=mumps", "--linear-threads=1",
                           "--convergence-tolerance-linear-relative=1e-6", "--trace=1",
                           "--quadrature-validation=32", f"--quadrature-order={order}"]
                env = dict(environment)
                env.pop("SWIFT_TARGET_DEGREE", None)
                if degree:
                    env["SWIFT_TARGET_DEGREE"] = str(degree)
                (directory / "command.json").write_text(json.dumps({
                    "argv": command, "target_degree": degree}, indent=2))
                start = time.monotonic()
                with (directory / "run.log").open("w") as log:
                    try:
                        code = subprocess.run(command, cwd=directory, env=env,
                                              stdout=log, stderr=subprocess.STDOUT,
                                              timeout=args.timeout, check=False).returncode
                    except subprocess.TimeoutExpired:
                        code = "timeout"
                row = {"family": family, "n": n, "target_degree": degree,
                       "order": order, "returncode": code,
                       "seconds": time.monotonic() - start}
                lines = (directory / "run.log").read_text().splitlines()
                if degree and f"quadrature target: degree={degree}" not in lines:
                    row["returncode"] = "target-degree-unverified"
                responses = [line.split("swift responses: ", 1)[1] for line in lines
                             if "swift responses: " in line]
                if responses:
                    response = dict(token.split("=", 1) for token in responses[-1].split()
                                    if "=" in token)
                    row.update({key: response[key] for key in fields if key in response})
                dense = [line.split("quadrature quality: ", 1)[1] for line in lines
                         if "quadrature quality: " in line]
                if dense:
                    row.update(dict(token.split("=", 1) for token in dense[-1].split()
                                    if "=" in token))
                writer.writerow(row)
                stream.flush()
                rows.append(row)
                print(f"{len(rows)}/{manifest['attempts']} {identity} "
                      f"C={row.get('geom_c', 'missing')} exit={row.get('exit', code)}",
                      flush=True)
    comparison = []
    for family, n, degree in cases:
        group = [row for row in rows if (row["family"], row["n"], row["target_degree"])
                 == (family, n, degree)]
        reference = next(row for row in group if row["order"] == 32)
        previous = next(row for row in group if row["order"] == 24)
        def valid(row):
            return (row["returncode"] == 0 and row.get("quality_ok") == "1"
                    and (row["family"] == "3d" or row.get("dense_quality_ok") == "1")
                    and row.get("inner_converged") == "1"
                    and row.get("exit") in ("full-interface-geometric-sup-converged",
                                             "iter-budget", "best-effort-step-stagnation",
                                             "best-effort-energy-stagnation")
                    and math.isfinite(float(row.get("geom_c", "nan"))))
        stable = valid(reference) and valid(previous)
        if stable:
            c = float(reference["geom_c"])
            stable = (abs(float(previous["geom_c"]) - c) <= 0.05 * c + 0.001
                      and previous["target_hit"] == reference["target_hit"])
        accepted = []
        if stable:
            for row in group:
                if (valid(row) and float(row["geom_c"]) <= 1.05 * c + 0.001
                        and (reference["target_hit"] != "1" or row["target_hit"] == "1")
                        and int(row["outer"]) <= int(reference["outer"]) + 2):
                    accepted.append(row["order"])
        comparison.append({"family": family, "n": n, "target_degree": degree,
                           "reference_stable": stable,
                           "reference_C": reference.get("geom_c"),
                           "acceptable_orders": sorted(accepted),
                           "smallest_observed_order": min(accepted) if accepted else None})
    (out / "comparison.json").write_text(json.dumps(comparison, indent=2))
    (out / "completed.json").write_text(json.dumps({"attempts": len(rows)}))


if __name__ == "__main__":
    main()
