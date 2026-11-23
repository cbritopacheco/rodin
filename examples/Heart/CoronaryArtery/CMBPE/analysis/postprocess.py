"""Compare CMBPEfluid coronary outputs using the paper's third-cycle reporting grid.

The grid is [1.70, 2.55) s at 1 ms by default. Fields and flow are sampled
on their own exported time axes as in the CMBPE reporter. A failed solve or
partial output must never silently substitute the last available samples.
"""
import argparse
from bisect import bisect_left
import csv
import json
import math
import os
from pathlib import Path
import statistics
import tempfile


OUTLET_COLUMNS = [f"outlet{tag}Pressure" for tag in (7, 8, 9, 10, 14, 15)]
FLOW_COLUMNS = ["CoronaryInletFlux", "IntramyoPressure"]
FIELD_COLUMNS = ["dissipation", "meanAbsWSS", "inletPressure", *OUTLET_COLUMNS]


def reporting_grid(period, cycle, dt):
    if min(period, cycle, dt) <= 0:
        raise ValueError("period, cycle and dt must be positive")
    count = round(period / dt)
    if count < 2 or not math.isclose(count * dt, period, abs_tol=1e-10, rel_tol=0):
        raise ValueError("the reporting interval must divide one period")
    return [(cycle - 1) * period + i * dt for i in range(count)]


def sample(path, columns, grid, dt, require_end):
    # CoronaryArtery.csv timestamps are end-of-step. The fields CSV timestamps
    # are the beginning of that step. Both are used as printed by report.py.
    # The flow file must also contain the final accepted step (2.55 by default).
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not {"t", *columns} <= set(reader.fieldnames):
            raise ValueError(f"missing output columns in {path}")
        rows = [(float(row["t"]), [float(row[key]) for key in columns]) for row in reader]
    if not rows:
        raise ValueError(f"empty output file: {path}")
    times = [t for t, _ in rows]
    tol = max(1e-8, dt * 1e-5)
    if not all(math.isfinite(t) for t in times) or any(b < a - tol for a, b in zip(times, times[1:])):
        raise ValueError(f"nonmonotone or invalid time axis: {path}")
    stop = grid[-1] + (dt if require_end else 0)
    if times[0] > grid[0] + tol or times[-1] < stop - tol:
        raise ValueError(f"incomplete cycle in {path}: need {grid[0]:g} through {stop:g} s; last sample {times[-1]:g}")
    selected = {k: [] for k in columns}
    for t in grid:
        idx = bisect_left(times, t - tol)
        if idx < len(times) and abs(times[idx] - t) <= tol:
            values = rows[idx][1]
        elif idx > 0 and idx < len(times) and 0 < times[idx] - times[idx - 1] <= dt + tol:
            # Only interpolate across an accepted step. A missing row on the
            # nominal grid must not be filled by bridging a two-step gap.
            fraction = (t - times[idx - 1]) / (times[idx] - times[idx - 1])
            values = [a + fraction * (b - a) for a, b in zip(rows[idx-1][1], rows[idx][1])]
        else:
            raise ValueError(f"no sample or enclosing accepted step at t={t:g} in {path}")
        for k, v in zip(columns, values):
            if not math.isfinite(v):
                raise ValueError(f"nonfinite {k} at t={t:g} in {path}")
            selected[k].append(v)
    return selected


def flow(run, grid, dt):
    raw = sample(run / "CoronaryArtery.csv", FLOW_COLUMNS, grid, dt, True)
    return [-1e6 * q for q in raw["CoronaryInletFlux"]], raw["IntramyoPressure"]


def absolute_mean(q):
    return statistics.fmean(abs(x) for x in q)


def metrics(run, q, grid, dt, systole):
    fields = sample(run / "CoronaryArtery_fields.csv", FIELD_COLUMNS, grid, dt, False)
    mean_out = statistics.fmean(statistics.fmean(fields[k]) for k in OUTLET_COLUMNS)
    m = {
        "mean_signed_inflow_ml_s": statistics.fmean(q),
        "mean_abs_inflow_ml_s": absolute_mean(q),
        "cycle_dissipation_mJ": 1e3 * dt * sum(fields["dissipation"]),
        "mean_wall_area_mean_WSS_Pa": statistics.fmean(fields["meanAbsWSS"]),
        "mean_inlet_minus_mean_of_six_outlet_pressures_Pa":
            statistics.fmean(fields["inletPressure"]) - mean_out,
    }
    for label, flag in [("systole", True), ("diastole", False)]:
        part = [x for x, state in zip(q, systole) if state == flag]
        m[label + "_signed_volume_ml"] = dt * sum(part)
        m[label + "_abs_volume_ml"] = dt * sum(abs(x) for x in part)
    return m


def atomic_write(path, writer):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, newline="", delete=False) as stream:
        tmp = stream.name
        try:
            writer(stream)
        except BaseException:
            os.unlink(tmp)
            raise
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reference", type=Path, nargs="?")
    ap.add_argument("candidate", type=Path, nargs="?")
    ap.add_argument("--mean", type=Path, help="print the complete-cycle mean absolute inflow in mL/s")
    ap.add_argument("--period", type=float, default=0.85)
    ap.add_argument("--cycle", type=int, default=3)
    ap.add_argument("--dt", type=float, default=0.001)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--waveform", type=Path)
    ap.add_argument("--scale", type=float)
    ap.add_argument("--tol", type=float, default=0.001, help="relative mean-flow matching tolerance")
    args = ap.parse_args()
    grid = reporting_grid(args.period, args.cycle, args.dt)
    if args.mean:
        if args.reference or args.candidate:
            ap.error("specify --mean or a reference/candidate pair")
        print(f"{absolute_mean(flow(args.mean, grid, args.dt)[0]):.12g}")
        return
    if args.reference is None or args.candidate is None:
        ap.error("a reference and candidate are required")
    if args.tol <= 0:
        ap.error("--tol must be positive")
    q0, pim = flow(args.reference, grid, args.dt)
    q1, _ = flow(args.candidate, grid, args.dt)
    if max(pim) <= min(pim):
        raise ValueError("reference intramyocardial pressure is constant")
    threshold = (min(pim) + max(pim)) / 2
    systole = [x >= threshold for x in pim]
    ref = metrics(args.reference, q0, grid, args.dt, systole)
    candidate = metrics(args.candidate, q1, grid, args.dt, systole)
    delta = {k: candidate[k] - v for k, v in ref.items()}
    relative = {k: 100 * delta[k] / abs(v) if v else None for k, v in ref.items()
                if not k.endswith("_volume_ml")}
    scale_q = absolute_mean(q0)
    if scale_q <= 0:
        raise ValueError("mean absolute reference inflow must be positive")
    signed_diffs = [y-x for x, y in zip(q0, q1)]
    absolute_diffs = [abs(y)-abs(x) for x, y in zip(q0, q1)]
    normalized_phase_differences = {
        k: 100 * v / (args.period * scale_q) for k, v in delta.items() if k.endswith("_volume_ml")
    }
    report = {
        "cycle": args.cycle, "period_s": args.period, "dt_s": args.dt, "n": len(grid),
        "reference_dir": str(args.reference.resolve()), "candidate_dir": str(args.candidate.resolve()),
        "scale": args.scale, "systole_fraction": statistics.fmean(systole),
        "reference": ref, "candidate": candidate, "candidate_minus_reference": delta,
        "relative_change_percent": relative,
        "phase_volume_difference_percent_of_reference_cycle_abs_volume": normalized_phase_differences,
        "signed_waveform_rms_difference_pct_mean_abs_reference": 100 * math.sqrt(statistics.fmean(d*d for d in signed_diffs)) / scale_q,
        "abs_waveform_rms_difference_pct_mean_abs_reference": 100 * math.sqrt(statistics.fmean(d*d for d in absolute_diffs)) / scale_q,
        "mean_relative_residual": delta["mean_abs_inflow_ml_s"] / scale_q,
        "matches_mean": abs(delta["mean_abs_inflow_ml_s"]) <= args.tol * scale_q,
    }
    if args.waveform:
        def write_waveform(stream):
            writer = csv.writer(stream)
            writer.writerow(["time_s", "systole", "reference_signed_ml_s", "candidate_signed_ml_s",
                             "candidate_minus_reference_signed_ml_s", "candidate_minus_reference_abs_ml_s"])
            writer.writerows(zip(grid, map(int, systole), q0, q1, signed_diffs, absolute_diffs))
        atomic_write(args.waveform, write_waveform)
    if args.out:
        atomic_write(args.out, lambda stream: json.dump(report, stream, indent=2, allow_nan=False))
    else:
        print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
