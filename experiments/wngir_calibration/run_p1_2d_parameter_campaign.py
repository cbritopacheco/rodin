#!/usr/bin/env python3
import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path


FINAL_RE = re.compile(
    r"WNGIR it=(?P<it>\d+)\s+fit=(?P<fit>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"alpha=(?P<alpha>(?:[0-9.eE+-]+|inf|nan))\s+step=(?P<step>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"min_j=(?P<min_j>(?:[0-9.eE+-]+|inf|nan))\s+max_j=(?P<max_j>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"max_qrel=(?P<max_qrel>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"active_rms=(?P<active_rms>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"active_sup=(?P<active_sup>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"active_rms_hg=(?P<active_rms_hg>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"act_frac=(?P<act_frac>(?:[0-9.eE+-]+|inf|nan))\s+cR=(?P<cR>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"pb_alpha=(?P<pb_alpha>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"pb_min_alpha=(?P<pb_min_alpha>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"pb_full_steps=(?P<pb_full_steps>\d+)\s+"
    r"rej_j=(?P<rej_j>\d+)\s+rej_q=(?P<rej_q>\d+)\s+rej_e=(?P<rej_e>\d+)\s+"
    r"converged=(?P<converged>\S+)\s+exit=(?P<exit>\S+)\s+"
    r"geom_rms=(?P<geom_rms>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"geom_sup=(?P<geom_sup>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"normal_rms=(?P<normal_rms>(?:[0-9.eE+-]+|inf|nan))")

ELEMENTS_RE = re.compile(r"^\s*elements=(?P<elements>\d+)", re.MULTILINE)

TIMING_RE = re.compile(
    r"wngir timing: it=(?P<timing_it>\d+)\s+"
    r"assembly=(?P<assembly>(?:[0-9.eE+-]+|inf|nan))\s+setup=(?P<setup>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"solve=(?P<solve>(?:[0-9.eE+-]+|inf|nan))\s+cgIt=(?P<cg_it>\d+)\s+"
    r"cgSolves=(?P<cg_solves>\d+)\s+cgMean=(?P<cg_it_mean>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"cgMax=(?P<cg_it_max>\d+)\s+"
    r"cgErr=(?P<cg_err>(?:[0-9.eE+-]+|inf|nan))\s+ls=(?P<ls>(?:[0-9.eE+-]+|inf|nan))\s+"
    r"exit=(?P<timing_exit>\S+)")


RESPONSE_FIELDS = ("energy", "geom_sup_target", "target_hit", "quality_ok",
                   "inner_total", "inner_max", "inner_last", "inner_converged",
                   "inner_residual", "inner_relative_residual", "inner_residual_tolerance")


def parse_responses(output):
    lines = [line for line in output.splitlines() if "wngir responses:" in line]
    if not lines:
        return {key: float("nan") for key in RESPONSE_FIELDS}
    fields = dict(re.findall(r"(\w+)=([^\s]+)", lines[-1]))
    return {key: float(fields.get(key, "nan"))
            for key in (*RESPONSE_FIELDS, "geom_sup", "min_j", "max_qrel")}


FIELDS = [
    "dataset", "n", "elements", "lobes", "kappa_f", "kappa_s", "kappa_d", "shape_curvature", "mu_hat", "kappa_j", "kappa_q",
    "fit", "geom_rms", "geom_sup", "normal_rms", "iterations", "alpha", "step",
    "min_j", "max_j", "max_qrel",
    "active_rms", "active_sup", "active_rms_hg", "act_frac", "cR", "rej_j", "rej_q", "rej_e",
    "pb_alpha", "pb_min_alpha", "pb_full_steps",
    "assembly", "setup", "solve", "cg_it", "cg_solves", "cg_it_mean", "cg_it_max",
    "cg_err", "ls", "converged", "exit",
    "seconds", "returncode", *RESPONSE_FIELDS]


def int_values(text):
    out = []
    for part in [p.strip() for p in text.split(",") if p.strip()]:
        if ":" not in part:
            out.append(int(part))
            continue
        bits = [int(x) for x in part.split(":")]
        if len(bits) == 2:
            a, b = bits
            s = 1
        elif len(bits) == 3:
            a, b, s = bits
        else:
            raise ValueError(f"invalid integer range: {part}")
        out.extend(range(a, b + 1, s))
    return sorted(dict.fromkeys(out))


def real_values(text):
    vals = []
    for part in [p.strip() for p in text.split(",") if p.strip()]:
        if part.startswith("log:"):
            _, a, b, m = part.split(":")
            a = float(a)
            b = float(b)
            m = int(m)
            if m == 1:
                vals.append(a)
            else:
                la = math.log(a)
                lb = math.log(b)
                vals.extend(math.exp(la + (lb - la) * i / (m - 1)) for i in range(m))
        elif ":" in part:
            bits = [float(x) for x in part.split(":")]
            if len(bits) != 3:
                raise ValueError(f"invalid real range: {part}")
            a, b, s = bits
            x = a
            while x <= b + 0.5 * s:
                vals.append(x)
                x += s
        else:
            vals.append(float(part))
    return list(dict.fromkeys(round(v, 14) for v in vals))


def deadline_from(text):
    """A wall-clock stop time, or None to run until the grid is finished."""
    if not text or text.lower() == "none":
        return None
    return dt.datetime.strptime(text, "%Y-%m-%d %H:%M:%S")


def read_done(path):
    done = set()
    rows = []
    if not path.exists():
        return rows, done
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
            done.add((
                row["dataset"],
                int(float(row["n"])),
                int(float(row["lobes"])),
                round(float(row["kappa_f"]), 14),
                round(float(row["kappa_s"]), 14),
                round(float(row["kappa_d"]), 14),
                round(float(row["mu_hat"]), 14),
                round(float(row["kappa_j"]), 14),
                round(float(row["kappa_q"]), 14), row["shape_curvature"]))
    return rows, done


def append_row(path, row):
    exists = path.exists()
    with path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def ensure_schema(path):
    if not path.exists():
        return
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames == FIELDS:
            return
        raise ValueError("CSV schema differs from canonical WNGIR; use a fresh out-dir")


def element_count(n):
    """Element count a run at resolution ``n`` will report.

    The example prints the mesh it actually built and that value is what gets
    logged; this predicts it, so a case can be matched against the already
    finished rows before it is run.  Mesh::UniformGrid(Triangle, {n, n}) takes
    grid *points*, hence 2 (n-1)^2 triangles.
    """
    return 2 * (n - 1) * (n - 1)


def save_case_trace(out_dir, identity, command, output):
    """Retain lossless inner/outer traces separately from terminal CSV records."""
    directory = Path(out_dir) / "iteration_logs"
    directory.mkdir(parents=True, exist_ok=True)
    identity_json = json.dumps(identity, sort_keys=True)
    name = hashlib.sha256(identity_json.encode()).hexdigest()
    path = directory / f"{name}.log"
    path.write_text("case=" + identity_json + "\ncommand=" + json.dumps(command)
                    + "\n" + output)
    return path


def run_case(args, exe, stage, n, lobes, kappa_f, kappa_s, kappa_d, mu_hat, kappa_j, kappa_q,
             shape_curvature="psd"):
    if shape_curvature != "psd":
        raise ValueError("only canonical PSD shape curvature is supported")
    cmd = [
        str(exe),
        f"--n={n}",
        f"--lobes={lobes}",
        "--cx=0.5",
        "--cy=0.5",
        "--phase=0",
        f"--amp={args.amp}",
        f"--R0={args.r0}",
        f"--classifier-eps={1.25 / (n - 1):.14g}",
        "--classifier-lambda=0.008",
        f"--wngir-kappa-s={kappa_s:.14g}",
        f"--wngir-kappa-d={kappa_d:.14g}",
        f"--wngir-mu-hat={mu_hat:.14g}",
        f"--wngir-kappa-j={kappa_j:.14g}",
        f"--wngir-kappa-q={kappa_q:.14g}",
        f"--wngir-kappa-f={kappa_f:.14g}", "--wngir-quality-guard=0.1",
        "--wngir-robust-scale=0",
        "--wngir-direct-solver=mumps", "--wngir-jsafe=1e-2",
        "--wngir-qmax=10",
        f"--wngir-primal-barrier-iterations={args.barrier_max_iters}",
        f"--wngir-primal-barrier-relative-tol={getattr(args, 'inner_rtol', 1e-3):.14g}",

        "--j-min=1e-8",
        "--wngir-jls=1e-2",
        "--wngir-armijo=1e-4",
        "--wngir-alpha-min=1e-4",
        "--wngir-omega-min=0.1",
        f"--wngir-geometric-sup-tol={1 / (n - 1) ** 2:.14g}",
        "--wngir-energy-stag-tol=1e-8",
        # Small steps indicate stagnation, not geometric success.
        f"--wngir-step-tol={1e-3 / (n - 1) ** 2:.14g}",
        f"--wngir-step-h-tol={1e-3 / (n - 1):.14g}",
        f"--wngir-cg-rtol={getattr(args, 'cg_rtol', 1e-9):.14g}",
        "--wngir-cg-max-iters=1000",
        f"--wngir-steps={args.steps}",
        "--output=0",
    ]
    if args.extra:
        cmd.extend(args.extra.split())
    if args.log_iterations:
        cmd.extend(["--trace=1", "--wngir-trace=1"])
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(args.threads)
    env["OPENBLAS_NUM_THREADS"] = str(args.threads)
    env["VECLIB_MAXIMUM_THREADS"] = str(args.threads)
    env["DYLD_LIBRARY_PATH"] = args.dyld_library_path
    t0 = time.time()
    proc = subprocess.run(
        cmd,
        cwd=getattr(args, "work_dir", args.root),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True)
    seconds = time.time() - t0
    if args.log_iterations:
        save_case_trace(args.out_dir,
                        dict(dataset=stage, n=n, lobes=lobes, kappa_f=kappa_f, kappa_s=kappa_s,
                             kappa_d=kappa_d, shape_curvature=shape_curvature,
                             mu_hat=mu_hat, kappa_j=kappa_j, kappa_q=kappa_q),
                        cmd, proc.stdout)
    elements_match = ELEMENTS_RE.search(proc.stdout)
    elements = int(elements_match.group("elements")) if elements_match else element_count(int(n))
    responses = parse_responses(proc.stdout)
    final_matches = list(FINAL_RE.finditer(proc.stdout))
    timing_matches = list(TIMING_RE.finditer(proc.stdout))
    row = {
        "dataset": stage,
        "n": n,
        "elements": elements,
        "lobes": lobes,
        "kappa_f": kappa_f,
        "kappa_s": kappa_s,
        "kappa_d": kappa_d,
        "shape_curvature": shape_curvature,
        "mu_hat": mu_hat,
        "kappa_j": kappa_j,
        "kappa_q": kappa_q,
        "seconds": seconds,
        "returncode": proc.returncode,
    }
    if final_matches:
        g = final_matches[-1].groupdict()
        row.update({
            "fit": float(g["fit"]),
            "geom_rms": float(g["geom_rms"]),
            "geom_sup": float(g["geom_sup"]),
            "normal_rms": float(g["normal_rms"]),
            "iterations": int(g["it"]),
            "alpha": float(g["alpha"]),
            "step": float(g["step"]),
            "min_j": float(g["min_j"]),
            "max_j": float(g["max_j"]),
            "max_qrel": float(g["max_qrel"]),
            "active_rms": float(g["active_rms"]),
            "active_sup": float(g["active_sup"]),
            "active_rms_hg": float(g["active_rms_hg"]),
            "act_frac": float(g["act_frac"]),
            "cR": float(g["cR"]),
            "pb_alpha": float(g["pb_alpha"]),
            "pb_min_alpha": float(g["pb_min_alpha"]),
            "pb_full_steps": int(g["pb_full_steps"]),
            "rej_j": int(g["rej_j"]),
            "rej_q": int(g["rej_q"]),
            "rej_e": int(g["rej_e"]),
            "converged": g["converged"],
            "exit": g["exit"],
        })
    else:
        row.update({
            "fit": float("nan"),
            "geom_rms": float("nan"),
            "geom_sup": float("nan"),
            "normal_rms": float("nan"),
            "iterations": -1,
            "alpha": float("nan"),
            "step": float("nan"),
            "min_j": float("nan"),
            "max_j": float("nan"),
            "max_qrel": float("nan"),
            "active_rms": float("nan"),
            "active_sup": float("nan"),
            "active_rms_hg": float("nan"),
            "act_frac": float("nan"),
            "cR": float("nan"),
            "pb_alpha": float("nan"),
            "pb_min_alpha": float("nan"),
            "pb_full_steps": -1,
            "rej_j": -1,
            "rej_q": -1,
            "rej_e": -1,
            "converged": "parse-fail",
            "exit": "parse-fail",
        })
        sys.stdout.write(f"PARSE_FAIL stage={stage} elements={elements} lobes={lobes} rc={proc.returncode}\n")
        sys.stdout.flush()
    if timing_matches:
        g = timing_matches[-1].groupdict()
        row.update({
            "assembly": float(g["assembly"]),
            "setup": float(g["setup"]),
            "solve": float(g["solve"]),
            "cg_it": int(g["cg_it"]),
            "cg_solves": int(g["cg_solves"]),
            "cg_it_mean": float(g["cg_it_mean"]),
            "cg_it_max": int(g["cg_it_max"]),
            "cg_err": float(g["cg_err"]),
            "ls": float(g["ls"]),
        })
    else:
        row.update({
            "assembly": float("nan"),
            "setup": float("nan"),
            "solve": float("nan"),
            "cg_it": -1,
            "cg_solves": -1,
            "cg_it_mean": float("nan"),
            "cg_it_max": -1,
            "cg_err": float("nan"),
            "ls": float("nan"),
        })
    row.update(responses)
    if not math.isfinite(row["target_hit"]):
        row["converged"] = "parse-fail"
        row["exit"] = "missing-canonical-responses"
    else:
        row["converged"] = "yes" if row["target_hit"] else "best-effort"
    return row


def cases_for(stage, ns, lobes, kappa_f, kappa_s, kappa_d, mu_hat, kappa_j, kappa_q,
              shape_curvatures=("psd",)):
    for n in ns:
        for l in lobes:
            for kf in kappa_f:
                for kb in kappa_s:
                    for r in kappa_d:
                        for mu in mu_hat:
                            for kj in kappa_j:
                                for kq in kappa_q:
                                    for shape in shape_curvatures:
                                        yield stage, n, l, kf, kb, r, mu, kj, kq, shape


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="/Users/carlos/Projects/rodin")
    parser.add_argument("--exe", default="build/examples/Geometry/LevelSetWNGIRReconstruction")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--deadline", default="none",
                        help="wall-clock stop time, or 'none' to run to completion")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--steps", type=int, default=30)
    parser.add_argument("--barrier-max-iters", type=int, default=15)
    parser.add_argument("--log-iterations", action="store_true", default=True,
                        help="save inner Newton and accepted-geometry traces for every case")
    parser.add_argument("--amp", type=float, default=0.08)
    parser.add_argument("--r0", type=float, default=0.24)
    parser.add_argument("--extra", default="")
    parser.add_argument("--scratch", default="/tmp/wngir-p1-2d-canonical")
    parser.add_argument("--dyld-library-path", default="/opt/local/lib/julia:/opt/local/lib")
    parser.add_argument("--dataset", default="canonical_p1_2d",
                        help="dataset name; the campaign reads and appends to <out-dir>/<dataset>.csv")
    parser.add_argument("--n", default="5,10,20,50,100")
    parser.add_argument("--lobes", default="0:10")
    metric_grid = "1e-4,1e-3,1e-2,0.1,1"
    parser.add_argument("--kappa-s", default=metric_grid)
    parser.add_argument("--kappa-d", default=metric_grid)
    parser.add_argument("--kappa-f", default=metric_grid,
                        help="fitting metric coefficient grid")
    parser.add_argument("--mu-hat", default="90")
    parser.add_argument("--shape-curvature", default="psd",
                        help="canonical psd shape curvature")
    parser.add_argument("--kappa-j", default="1")
    parser.add_argument("--kappa-q", default="1")
    args = parser.parse_args()
    shape_curvatures = list(dict.fromkeys(args.shape_curvature.split(",")))
    if not shape_curvatures or any(value != "psd" for value in shape_curvatures):
        parser.error("only canonical psd shape curvature is supported")
    if not 1 <= args.steps <= 200:
        parser.error("--steps must be between 1 and the campaign cap of 200")
    if args.barrier_max_iters < 1:
        parser.error("--barrier-max-iters must be positive")
    if "--initial-mmg" in args.extra:
        sys.stderr.write("element_count assumes the uniform grid; --initial-mmg changes it\n")
        raise SystemExit(2)

    args.root = Path(args.root)
    exe = Path(args.exe)
    if not exe.is_absolute():
        exe = args.root / exe
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = args.root / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir = out_dir
    args.work_dir = Path(args.scratch)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{args.dataset}.csv"
    ensure_schema(csv_path)
    deadline = deadline_from(args.deadline)
    stages = [
        (args.dataset,
         int_values(args.n),
         int_values(args.lobes),
         real_values(args.kappa_f),
         real_values(args.kappa_s),
         real_values(args.kappa_d),
         real_values(args.mu_hat)),
    ]
    if any(not grid or any(not math.isfinite(value) or value <= 0 for value in grid)
           for grid in stages[0][3:6]):
        parser.error("kappa-f, kappa-s and kappa-d grids must be finite and strictly positive")
    kappa_j = real_values(args.kappa_j)
    kappa_q = real_values(args.kappa_q)
    expected_cases = sum(
        len(ns) * len(lobes) * len(kfs) * len(kbs) * len(kappa_ds) * len(mus)
        * len(kappa_j) * len(kappa_q) * len(shape_curvatures)
        for _, ns, lobes, kfs, kbs, kappa_ds, mus in stages)
    manifest = {
        "model": "F+S+D-affine-quadratic-hinges-psd-residual-dinf-v6",
        "executable": str(exe.resolve()),
        "sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
        "dataset": args.dataset,
        "expected_cases": expected_cases,
        "shape_curvature": shape_curvatures,
        "n": stages[0][1],
        "lobes": stages[0][2],
        "kappa_f": stages[0][3],
        "kappa_s": stages[0][4],
        "kappa_d": stages[0][5],
        "mu_hat": stages[0][6],
        "kappa_j": kappa_j,
        "kappa_q": kappa_q,
        "target": {"center": [0.5, 0.5], "R0": args.r0, "amplitude": args.amp, "phase": 0},
        "classifier": {"epsilon_over_h": 1.25, "lambda_c": 0.008},
        "fixed_profile": {
            "model": "F+Splus+D-affine-quadratic-hinges-residual-dinf-v6", "quality_guard": 0.1, "robust_scale": "automatic",
            "j_safe": 1e-2, "q_max": 10,
            "barrier_iterations": args.barrier_max_iters, "barrier_relative_tolerance": 1e-3,
            "j_min": 1e-8, "j_line_search": 1e-2,
            "armijo": 1e-4,
            "alpha_min": 1e-4, "omega_min": 0.1,
            "geometric_sup_target": "h^2 (sampled normalized residual including vertices)",
            "inner_residual_absolute_tolerance": 1e-12,
            "inner_residual_scale": "max(initial residual norm, fitting force norm)",
            "stagnation_iterations": 5, "energy_stagnation": 1e-8,
            "absolute_step": "0.001 h^2", "accepted_step_over_h": "0.001 h",
            "linear_backend": "MUMPS", "cg_relative_tolerance": 1e-9, "cg_max_iterations": 1000,
            "cell_quadrature_order": "automatic: max(2, 2 * FE order)",
            "interface_quadrature_order": "automatic: max(4, 2 * FE order + 2)",
            "geometric_validation_order": "automatic: max(6, 2 * FE order + 4)",
            "max_iterations": args.steps,
        },
    }
    if args.log_iterations:
        manifest["iteration_logging"] = "inner Newton and full-interface geometry trace v1"
    manifest_path = out_dir / f"{args.dataset}_manifest.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        parser.error("manifest differs from existing dataset; use a new out-dir")
    with manifest_path.open("w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    _, done = read_done(csv_path)

    sys.stdout.write(f"out={csv_path}\n")
    sys.stdout.write(f"deadline={deadline or 'none'}\n")
    sys.stdout.write(f"exe={exe}\n")
    sys.stdout.write(f"steps={args.steps} threads={args.threads} amp={args.amp} R0={args.r0}\n")
    sys.stdout.flush()

    completed = 0
    for stage, ns, lobes, kfs, kbs, kappa_ds, mus in stages:
        all_cases = list(cases_for(stage, ns, lobes, kfs, kbs, kappa_ds, mus,
                                  kappa_j, kappa_q, shape_curvatures))
        missing = [c for c in all_cases if (
            c[0], c[1], c[2], round(c[3], 14), round(c[4], 14),
            round(c[5], 14), round(c[6], 14), round(c[7], 14), round(c[8], 14), c[9]) not in done]
        sys.stdout.write(
            f"{stage}: total={len(all_cases)} done={len(all_cases) - len(missing)} "
            f"missing={len(missing)} n={ns} lobes={lobes}\n")
        sys.stdout.flush()
        if deadline is not None and dt.datetime.now() >= deadline:
            sys.stdout.write(f"{stage}: skipped because deadline has been reached\n")
            sys.stdout.flush()
            break
        stage_t0 = time.time()
        for i, case in enumerate(missing, 1):
            if deadline is not None and dt.datetime.now() >= deadline:
                sys.stdout.write(f"{stage}: stopped at deadline after {i - 1} new cases\n")
                sys.stdout.flush()
                break
            row = run_case(args, exe, *case)
            append_row(csv_path, row)
            done.add((
                row["dataset"], int(row["n"]), int(row["lobes"]),
                round(float(row["kappa_f"]), 14),
                round(float(row["kappa_s"]), 14), round(float(row["kappa_d"]), 14),
                round(float(row["mu_hat"]), 14), round(float(row["kappa_j"]), 14),
                round(float(row["kappa_q"]), 14), row["shape_curvature"]))
            completed += 1
            elapsed = time.time() - stage_t0
            eta = elapsed * (len(missing) - i) / i if i else 0
            sys.stdout.write(
                f"{stage} {i}/{len(missing)} elements={row['elements']} lobes={row['lobes']} "
                f"kappa_f={float(row['kappa_f']):.4g} kappa_s={float(row['kappa_s']):.4g} "
                f"kappa_d={float(row['kappa_d']):.3g} "
                f"shape={row['shape_curvature']} "
                f"mu={float(row['mu_hat']):.3g} fit={float(row['fit']):.4g} "
                f"it={row['iterations']} Q={float(row['max_qrel']):.3g} "
                f"sec={float(row['seconds']):.2f} eta={eta/3600:.2f}h\n")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
