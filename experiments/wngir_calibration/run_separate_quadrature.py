"""Serial end-to-end checks with independent surface/volume/validation rules."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--binaries-from", type=Path)
    parser.add_argument("--source-2d", type=Path,
                        help="Experimental 2D driver supporting WNGIR_TARGET_DEGREE; unnecessary with --binaries-from")
    parser.add_argument("--retry-analytic", action="store_true")
    parser.add_argument("--p1-validation", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    build, output = args.build.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / "manifest.json").exists():
        raise RuntimeError("Use a fresh campaign directory")
    entries = json.loads((build / "compile_commands.json").read_text())
    binaries = {}
    for target, family in (("LevelSetWNGIRSweep", "2d"), ("LevelSetWNGIRSweepP2", "p2"), ("LevelSetWNGIRSweep3D", "3d")):
        if args.binaries_from:
            executable = output / ("solver-" + family)
            shutil.copy2(args.binaries_from / executable.name, executable)
            binaries[family] = executable
            continue
        entry = next(e for e in entries if f"CMakeFiles/{target}.dir/" in e["command"])
        command = shlex.split(entry["command"])
        source = Path(command[command.index("-c") + 1])
        if family != "3d":
            if args.source_2d is None:
                parser.error("--source-2d is required when building the FE-target drivers")
            source = args.source_2d.resolve()
            command.insert(1, "-I" + str(root / "examples/Geometry"))
        obj, executable = output / (family + ".o"), output / ("solver-" + family)
        command[command.index("-o") + 1] = str(obj)
        command[command.index("-c") + 1] = str(source)
        link = shlex.split((build / f"examples/Geometry/CMakeFiles/{target}.dir/link.txt").read_text())
        link[link.index("-o") + 1] = str(executable)
        link = [str(obj) if value.endswith(".cpp.o") else value for value in link]
        if args.build_only or not executable.exists():
            (output / (family + ".commands.json")).write_text(json.dumps([command, link], indent=2))
            (output / (family + ".source.sha256")).write_text(hashlib.sha256(source.read_bytes()).hexdigest())
            with (output / (family + ".build.log")).open("w") as log:
                for cmd in (command, link):
                    subprocess.run(cmd, cwd=entry["directory"], stdout=log, stderr=subprocess.STDOUT, check=True)
        binaries[family] = executable
    if args.build_only:
        return
    # Coupled-rule controls preserve the pilot setting. Separate rules change
    # only integration/validation, never model weights or stopping budgets.
    cases = [("2d", 32, 0), ("2d", 16, 1), ("2d", 16, 3), ("3d", 8, 0), ("p2", 8, 0)]
    if args.retry_analytic:
        cases = [("2d", 32, 0), ("p2", 8, 0)]
    variants = {
        "control": ["--quad-order=12"],
        "surface8": ["--surface-quadrature-order=8", "--volume-quadrature-order=2", "--quality-validation-order=16"],
        "surface12": ["--surface-quadrature-order=12", "--volume-quadrature-order=2", "--quality-validation-order=16"],
        "surface24": ["--surface-quadrature-order=24", "--volume-quadrature-order=2", "--quality-validation-order=16"],
    }
    if args.p1_validation:
        cases = [("2d", 32, 0), ("2d", 16, 1), ("2d", 16, 3), ("3d", 8, 0)]
        variants = {key: [value.replace("quality-validation-order=16", "quality-validation-order=2") for value in values]
                    for key, values in variants.items() if key in ("surface8", "surface12")}
    manifest = {"cases": cases, "variants": variants, "attempts": len(cases) * len(variants), "timeout_seconds": 180,
                "binary_sha256": {key: hashlib.sha256(value.read_bytes()).hexdigest() for key, value in binaries.items()},
                "threads": 1, "warmups": 0, "repeats": 0,
                "note": "P2 uses volume8 rather than2; final independent quality32. P1 validation-only subset uses quality2. Geometric validation32. Sampled references, not certificates."}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    rows = []
    for family, n, degree in cases:
        for variant, options in variants.items():
            if family == "p2":
                options = [value.replace("volume-quadrature-order=2", "volume-quadrature-order=8") for value in options]
            directory = output / f"{family}-n{n}-target{degree}-{variant}"
            directory.mkdir()
            command = [str(binaries[family]), f"--n={n}", "--frames=1", "--lobes=4", "--amp=0.08", "--R0=0.24", "--orbitR=0",
                       "--wngir-outer-iterations=30", "--wngir-inner-iterations=15", "--wngir-fit=1",
                       "--wngir-distribution-deviatoric=0.0001", "--wngir-distribution-divergence=0.01", "--wngir-hinge=10",
                       "--wngir-linear-solver=mumps", "--wngir-linear-threads=1", "--wngir-linear-relative-tolerance=1e-6",
                       "--wngir-trace=1", "--geometric-validation-order=32", *options]
            (directory / "command.json").write_text(json.dumps(command, indent=2))
            env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
            env.pop("WNGIR_TARGET_DEGREE", None)
            if degree:
                env["WNGIR_TARGET_DEGREE"] = str(degree)
            start = time.monotonic()
            with (directory / "run.log").open("w") as log:
                try:
                    status = subprocess.run(command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180).returncode
                except subprocess.TimeoutExpired:
                    status = "timeout"
            row = dict(family=family, n=n, target_degree=degree, variant=variant, returncode=status, seconds=time.monotonic() - start)
            for line in (directory / "run.log").read_text().splitlines():
                for marker in ("wngir responses: ", "quadrature quality: "):
                    if marker in line:
                        row.update(dict(token.split("=", 1) for token in line.split(marker)[-1].split() if "=" in token))
            rows.append(row)
            (output / "results.json").write_text(json.dumps(rows, indent=2))
            print(f"{len(rows)}/{manifest['attempts']} {directory.name}: {row.get('exit', status)}, C={row.get('geom_c')}", flush=True)
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with (output / "results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (output / "completed.json").write_text(json.dumps({"attempts": len(rows), "successful": sum(r["returncode"] == 0 for r in rows)}, indent=2))


if __name__ == "__main__":
    main()
