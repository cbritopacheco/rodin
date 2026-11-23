"""Run a fixed-calibration Y campaign and compare flow against the full mesh."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np


PERIOD = 0.85
RADIUS_M = np.array([72e-6, 46e-6])
CASES = ("full", "constant", "variable")
COLORS = ("#222222", "#D55E00", "#0072B2")
STYLES = ("-", "--", "-.")
LABELS = ("Fully resolved Y", r"Constant $\mathcal{R}_0$", r"Shear-dependent $\mathcal{R}_\mu$")


def cycle_mean(values):
    # Last cycle sampled at uniform time steps. Its endpoint replaces phase zero.
    return np.mean(values, axis=-1)


def plot_cases(path, waves, steps, dim):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6), sharex=True, layout="constrained")
    phase = np.arange(steps + 1) / steps
    scale = 1e9 if dim == 3 else 1e6
    for col, (condition, cases) in enumerate(waves.items()):
        for case, color, style, label in zip(CASES, COLORS, STYLES, LABELS):
            wide, narrow = cases[case]
            for row, values in enumerate((scale * (wide + narrow), 100 * narrow / (wide + narrow))):
                axes[row, col].plot(phase, np.r_[values[-1], values], c=color, ls=style, lw=2, label=label)
        axes[0, col].set_title(condition.replace("_", " ").capitalize(), pad=7)
        axes[1, col].set_xlabel("Fraction of cycle")
    axes[0, 0].set_ylabel("Total inflow (µL/s)" if dim == 3 else "Inflow per unit depth (mm²/s)")
    axes[1, 0].set_ylabel("Narrow-branch flow share (%)")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside upper center", ncols=3, frameon=False)
    fig.savefig(path, dpi=175)
    plt.close(fig)


def compare(output, conditions, steps, cycles, dim):
    rows = []
    waves = {}
    suffix = f"m{dim}_s"
    flow_fields = (f"wide_{suffix}", f"narrow_{suffix}")
    rate_per_flow = 4 / (np.pi * RADIUS_M**3) if dim == 3 else 3 / (2 * RADIUS_M**2)
    for condition, head in conditions.items():
        folder = output / condition
        waves[condition] = {}
        for case in CASES:
            raw = np.genfromtxt(folder / f"{case}.csv", delimiter=",", names=True, ndmin=1)
            if len(raw) != steps * cycles or not np.all(np.isfinite(raw.view(float))):
                raise ValueError(f"Incomplete or nonfinite waveforms for {condition}/{case}")
            last = raw[-steps:]
            flows = np.array([last[field] for field in flow_fields])
            waves[condition][case] = flows
            total = flows.sum(axis=0)
            if min(total) <= 0:
                raise ValueError(f"The instantaneous split requires positive total inflow for {condition}/{case}")
            mean_flow = cycle_mean(flows)
            reference = waves[condition]["full"]
            reference_mean = cycle_mean(reference)
            mean_share = 100 * mean_flow[1] / mean_flow.sum()
            ref_share = 100 * reference_mean[1] / reference_mean.sum()
            share_difference = 100 * (flows[1] / total - reference[1] / reference.sum(axis=0))
            rms = 100 * np.sqrt(cycle_mean((flows - reference)**2)) / cycle_mean(abs(reference))
            mass_error = 100 * np.max(abs(total - last[f"inlet_{suffix}"])) / np.max(abs(total))
            previous_end = np.array([raw[field][-steps - 1] for field in flow_fields])
            periodicity = np.max(100 * abs(flows[:, -1] - previous_end) / np.max(abs(flows), axis=1))
            if mass_error > 1e-4 or periodicity > 0.01:
                raise ValueError(f"Mass or periodicity check failed for {condition}/{case}: {mass_error=:.4g}, {periodicity=:.4g}")
            shear = rate_per_flow[:, None] * abs(flows)
            rows.append({
                "condition": condition, "mean_head_pa": head, "closure": case,
                "narrow_share_pct": mean_share, "narrow_share_error_pp": mean_share - ref_share,
                "narrow_share_rms_error_pp": np.sqrt(cycle_mean(share_difference**2)),
                f"total_mean_inflow_{suffix}": mean_flow.sum(),
                "total_mean_inflow_error_pct": 100 * (mean_flow.sum() / reference_mean.sum() - 1),
                "wide_flow_rms_error_pct": rms[0], "narrow_flow_rms_error_pct": rms[1],
                "peak_mass_imbalance_pct": mass_error, "periodic_endpoint_change_pct": periodicity,
                "max_picard_iterations": int(last["iterations"].max()),
                "wide_apparent_shear_min_s_1": shear[0].min(), "wide_apparent_shear_max_s_1": shear[0].max(),
                "narrow_apparent_shear_min_s_1": shear[1].min(), "narrow_apparent_shear_max_s_1": shear[1].max(),
            })
    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    plot_cases(output / "flow_split.png", waves, steps, dim)
    for row in rows:
        print(f"{row['condition']:12} {row['closure']:8} narrow share {row['narrow_share_pct']:.3f}% "
              f"({row['narrow_share_error_pp']:+.3f} pp), narrow waveform error {row['narrow_flow_rms_error_pct']:.2f}%")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dim", type=int, choices=(2, 3), default=3)
    parser.add_argument("--head", type=float, default=8.2, help="Baseline mean pressure head, in Pa")
    parser.add_argument("--low-head-ratio", type=float, default=0.20)
    parser.add_argument("--steps", type=int, default=48)
    parser.add_argument("--cycles", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--calibration", type=Path, help="Reuse an external baseline calibration for refinement")
    args = parser.parse_args()
    if not (np.isfinite(args.head) and args.head > 0 and 0 < args.low_head_ratio < 1):
        parser.error("Head must be positive and the low-head ratio must lie between zero and one")
    if args.steps < 4 or args.cycles < 2 or args.jobs < 1:
        parser.error("Use at least four steps, two cycles and one job")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    calibration = output / "calibration.txt"
    conditions = {"baseline": args.head, "reduced_head": args.head * args.low_head_ratio}
    env = os.environ | {"OMP_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    options = ("-ksp_type", "preonly", "-pc_type", "lu", "-pc_factor_mat_solver_type", "mumps")

    def run(mode, head, csv_path, steps):
        mesh = args.mesh / ("full.mesh" if mode == "full" else "cut.mesh")
        command = [args.solver, mesh, mode, csv_path, str(head), str(steps), str(args.cycles), calibration, *options]
        if args.dim == 3 and head == conditions["reduced_head"] and steps and steps % 4 == 0:
            peak_step = (args.cycles - 1) * steps + steps // 4
            command.extend(["-y_snapshot_step", str(peak_step)])
        log_path = csv_path.with_suffix(".log")
        with log_path.open("w") as log:
            run_result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
        if run_result.returncode != 0:
            raise RuntimeError(f"{mode} returned {run_result.returncode}; see {log_path}")

    if args.calibration:
        shutil.copyfile(args.calibration, calibration)
    else:
        run("full", args.head, output / "steady.csv", 0)
    for condition in conditions:
        (output / condition).mkdir(exist_ok=True)
    (output / "protocol.json").write_text(json.dumps({
        "mean_heads_pa": conditions, "period_s": PERIOD, "steps_per_cycle": args.steps,
        "cycles": args.cycles, "mesh": str(args.mesh.resolve()), "dimension": args.dim,
        "full_mesh_sha256": sha256((args.mesh / "full.mesh").read_bytes()).hexdigest(),
        "cut_mesh_sha256": sha256((args.mesh / "cut.mesh").read_bytes()).hexdigest(),
        "calibration_sha256": sha256(calibration.read_bytes()).hexdigest(),
        "calibration_source": str(args.calibration.resolve()) if args.calibration else "steady full Y at baseline head",
    }, indent=2) + "\n")
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run, mode, head, output / condition / f"{mode}.csv", args.steps)
                   for condition, head in conditions.items() for mode in CASES]
        for future in futures:
            future.result()
    compare(output, conditions, args.steps, args.cycles, args.dim)


if __name__ == "__main__":
    main()
