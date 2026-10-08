"""Build and run the serial, fixed-state order calibration without repeats."""

import argparse
import hashlib
import itertools
import json
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
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    build = args.build.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = root / "experiments/swift_calibration/QuadratureCalibration.cpp"
    executable = output / "fixed-quadrature"
    entries = json.loads((build / "compile_commands.json").read_text())
    entry = next(e for e in entries if "CMakeFiles/LevelSetSWIFTSweep.dir/" in e["command"])
    compile_command = shlex.split(entry["command"])
    obj = output / "fixed-quadrature.o"
    compile_command[compile_command.index("-o") + 1] = str(obj)
    compile_command[compile_command.index("-c") + 1] = str(source)
    link = shlex.split((build / "examples/Adaptation/SWIFT/CMakeFiles/LevelSetSWIFTSweep.dir/link.txt").read_text())
    link[link.index("-o") + 1] = str(executable)
    link = [str(obj) if v.endswith(".cpp.o") else v for v in link]
    if (output / "manifest.json").exists():
        raise RuntimeError("Refusing to modify or rerun an existing campaign")
    if args.build_only or not executable.exists():
        shutil.copy2(source, output / source.name)
        shutil.copy2(source.with_name("QuadratureAudit.h"), output / "QuadratureAudit.h")
        (output / "commands.json").write_text(json.dumps([compile_command, link], indent=2))
        with (output / "build.log").open("w") as log:
            for command in (compile_command, link):
                subprocess.run(command, cwd=entry["directory"], stdout=log,
                               stderr=subprocess.STDOUT, check=True)
    if args.build_only:
        return
    if (output / "manifest.json").exists():
        raise RuntimeError("Refusing to overwrite an existing campaign")
    cases = list(itertools.product((2, 3), (1, 2, 3), (0, 1)))
    manifest = {"cases": cases, "binary_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
                "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "audit_sha256": hashlib.sha256(source.with_name("QuadratureAudit.h").read_bytes()).hexdigest(),
                "threads": 1, "repeats": 0, "timeout_seconds": 240}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    import os
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    results = []
    for dimension, degree, stress in cases:
        command = [str(executable), str(dimension), str(degree), str(stress)]
        start = time.monotonic()
        name = f"d{dimension}-p{degree}-s{stress}"
        with (output / (name + ".log")).open("w") as log:
            try:
                process = subprocess.run(command, cwd=output, env=env, stdout=log,
                                         stderr=subprocess.STDOUT, timeout=240)
                status = process.returncode
            except subprocess.TimeoutExpired:
                status = "timeout"
        row = {"case": name, "command": command, "status": status, "seconds": time.monotonic() - start}
        results.append(row)
        (output / "results.json").write_text(json.dumps(results, indent=2))
        print(row, flush=True)
    (output / "completed.json").write_text(json.dumps({"attempts": len(results), "successful": sum(r["status"] == 0 for r in results)}, indent=2))


if __name__ == "__main__":
    main()
