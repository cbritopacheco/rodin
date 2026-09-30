#!/usr/bin/env python3
"""Summarize inner Newton work and first D_infty target attainment from saved traces."""

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path

PAIRS = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")


def summarize(path, coefficients, order=2, j_floor=1e-2, q_max=10):
    lines = path.read_text().splitlines()
    identity = json.loads(lines[0].removeprefix("case="))
    coefficient = float(coefficients[str(identity["lobes"])])
    if not math.isfinite(coefficient) or coefficient <= 0:
        raise ValueError("geometric reference coefficients must be finite and positive")
    geometry, inner = [], []
    for line in lines[2:]:
        values = dict(PAIRS.findall(line))
        if line.strip().startswith("wngir geometry:"):
            geometry.append(values)
        elif line.strip().startswith("barrier inner="):
            inner.append(values)
    accepted = [g for g in geometry if g["phase"] in ("initial", "accepted")]
    hits = [g for g in accepted
            if float(g["geom_sup"]) <= coefficient * float(g["h"]) ** order
            and float(g["min_j"]) > j_floor and float(g["max_qrel"]) < q_max]
    hit = min(hits, key=lambda g: int(g["outer"])) if hits else None
    counts = {}
    for row in inner:
        outer = int(row["outer"])
        counts[outer] = counts.get(outer, 0) + 1
    return dict(identity, log=str(path), reference_coefficient=coefficient,
                error_order=order, geometry_records=len(accepted),
                inner_attempts=len(inner),
                inner_mean=statistics.mean(counts.values()) if counts else "",
                inner_median=statistics.median(counts.values()) if counts else "",
                inner_max=max(counts.values(), default=0),
                inner_linear_failures=sum(r.get("linear_ok") == "0" for r in inner),
                inner_damped_steps=sum(0 < float(r.get("alpha", "nan")) < 1 for r in inner),
                target_hit=int(hit is not None),
                first_hit_outer=int(hit["outer"]) if hit else "",
                first_hit_inner_total=int(hit["inner_total"]) if hit else "",
                first_hit_seconds=float(hit["seconds"]) if hit else "",
                final_geom_sup=float(geometry[-1]["geom_sup"]) if geometry else "")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="campaign iteration_logs directory")
    parser.add_argument("--coefficients", required=True, type=Path,
                        help='JSON map of lobe index to independent C_infty reference, e.g. {"0": 0.2}')
    parser.add_argument("--order", type=int, default=2, help="target exponent k+1; P1 uses 2")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.order < 1:
        parser.error("order must be positive")
    coefficients = json.loads(args.coefficients.read_text())
    rows = [summarize(path, coefficients, args.order)
            for path in sorted(args.directory.glob("*.log"))]
    if not rows:
        parser.error("no iteration logs found")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with args.output.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
