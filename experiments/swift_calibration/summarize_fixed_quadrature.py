"""Summarize fixed-state errors against the retained high-order references."""
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    rows = []
    for log in sorted(args.directory.glob("d*-p*-s*.log")):
        kind = "volume"
        records = []
        for line in log.read_text().splitlines():
            if not line.startswith("audit "):
                continue
            category = line.split()[1]
            values = dict(token.split("=", 1) for token in line.split()[2:] if "=" in token)
            if category == "target":
                kind = values["kind"]
            else:
                values.update(case=log.stem, kind=kind, category=category)
                records.append(values)
        quality = [r for r in records if r["category"] == "quality"]
        for row in records:
            if row["category"] == "volume" and quality:
                reference = quality[-1]
                row["j_error"] = abs(float(row["min_j"]) - float(reference["min_j"])) / max(abs(float(reference["min_j"])), 1e-14)
                row["q_error"] = abs(float(row["max_q"]) - float(reference["max_q"])) / max(abs(float(reference["max_q"])), 1e-14)
            if row["category"] == "surface":
                geometry = [r for r in records if r["category"] == "geometry" and r["kind"] == row["kind"]]
                if geometry:
                    row["sup_underestimate"] = max(0, 1 - float(row["geom_sup"]) / float(geometry[-1]["geom_sup"]))
            rows.append(row)
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with (args.directory / "fixed-errors.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for case in sorted({r["case"] for r in rows}):
        print(case)
        for kind in ("analytic", "fe"):
            candidates = [r for r in rows if r["case"] == case and r["kind"] == kind and r["category"] == "surface"]
            stable = next((r for r in candidates if r["order"] == "32"), None)
            stable = stable and all(float(stable[key]) <= limit for key, limit in (("energy_error", 1e-4), ("force_error", 1e-3), ("metric_error", 1e-3)))
            passing = [int(r["order"]) for r in candidates if all(float(r[key]) <= limit for key, limit in (("energy_error", 1e-3), ("force_error", 1e-2), ("metric_error", 1e-2)))]
            print(" surface", kind, "reference_stable", bool(stable), "passing", sorted(passing))
        candidates = [r for r in rows if r["case"] == case and r["category"] == "volume"]
        stable = next((r for r in candidates if r["order"] == "24"), None)
        stable = stable and all(float(stable[key]) <= 1e-3 for key in ("action_error", "hinge_error", "hinge_energy_error"))
        passing = [int(r["order"]) for r in candidates if all(float(r[key]) <= 1e-2 for key in ("action_error", "hinge_error", "hinge_energy_error"))]
        quality = [int(r["order"]) for r in candidates if float(r["j_error"]) <= 1e-2 and float(r["q_error"]) <= 1e-2]
        print(" volume", "reference_stable", bool(stable), "passing", sorted(passing), "quality", sorted(quality))


if __name__ == "__main__":
    main()
