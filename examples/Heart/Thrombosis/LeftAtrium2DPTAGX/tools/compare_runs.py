#!/usr/bin/env python3
"""Compare LeftAtrium2DPTAG / LeftAtrium2DPTAGX runs from their CSV files.

Usage:
  compare_runs.py run_dir_or_csv[:label] ... [--png out.png] [--from-cycle k]

Prints one row per run (Markdown table) and, with --png, plots the domain and
appendage integrals against time.  CSVs of LeftAtrium2DPTAG (no X/integral
columns) are accepted; the missing columns are left blank.
"""
import argparse
import csv
import glob
import math
import os


def load(path):
    if os.path.isdir(path):
        c = sorted(glob.glob(os.path.join(path, "*.csv")))
        if not c:
            raise SystemExit(f"no CSV in {path}")
        path = c[0]
    with open(path) as f:
        rows = list(csv.DictReader(f))
    cols = {k: [float(r[k]) for r in rows] for k in rows[0]} if rows else {}
    return path, cols


def trapz(t, y):
    return sum(0.5 * (y[i] + y[i - 1]) * (t[i] - t[i - 1]) for i in range(1, len(t)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--png")
    ap.add_argument("--from-cycle", type=int, default=0,
                    help="ignore rows before this cycle in the time integrals")
    a = ap.parse_args()

    runs = []
    for spec in a.runs:
        path, label = (spec.split(":", 1) + [None])[:2] if ":" in spec else (spec, None)
        p, c = load(path)
        runs.append((label or os.path.basename(os.path.dirname(os.path.abspath(p))), c))

    head = ["run", "t_end (s)", "max T (nM)", "int T dt", "LAA int T dt", "max X (pM)",
            "int F end", "LAA F end", "LAA F / int F", "max F end (uM)", "max drift"]
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for label, c in runs:
        if not c:
            print(f"| {label} | (empty) |")
            continue
        keep = [i for i, k in enumerate(c["cycle"]) if k >= a.from_cycle]
        t = [c["t"][i] for i in keep]
        g = lambda k: [c[k][i] for i in keep] if k in c else None  # noqa: E731
        fmt = lambda x, s=1.0: "" if x is None else f"{x * s:.3g}"  # noqa: E731
        intT, laaT = g("intT"), g("laaT")
        intF, laaF = (c.get("intF") or [None])[-1], (c.get("laaF") or [None])[-1]
        ratio = laaF / intF if intF and laaF is not None and intF > 0 else None
        row = [label, fmt(c["t"][-1]), fmt(max(c["maxT"]), 1e6),
               fmt(trapz(t, intT) if intT else None), fmt(trapz(t, laaT) if laaT else None),
               fmt(max(c["maxX"]) if "maxX" in c else None, 1e9), fmt(intF), fmt(laaF),
               fmt(ratio), fmt(c["maxF"][-1], 1e3), fmt(max(c["invariantDrift"]))]
        print("| " + " | ".join(row) + " |")
    print("\nUnits: concentrations in mol/m^3 in the CSV; integrals per unit depth "
          "(mol/m in 2D); 'int T dt' in mol s/m.")

    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        panels = [("intT", "∫T dA (mol/m)"), ("laaT", "LAA ∫T dA (mol/m)"),
                  ("intF", "∫F dA (mol/m)"), ("laaF", "LAA ∫F dA (mol/m)"),
                  ("intX", "∫X dA (mol/m)"), ("maxT", "max T (mol/m³)")]
        colors = ["#2a6fdb", "#e0832b", "#2e9e6a", "#c2417a", "#7a5cc7", "#8a8a8a",
                  "#1f9bb0", "#b5892a"]
        fig, ax = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
        for (k, title), axi in zip(panels, ax.flat):
            for j, (label, c) in enumerate(runs):
                if k in c:
                    axi.plot(c["t"], c[k], lw=1.6, color=colors[j % len(colors)], label=label)
            axi.set_title(title, fontsize=10)
            axi.grid(alpha=0.25)
        for axi in ax[-1]:
            axi.set_xlabel("t (s)")
        ax[0][0].legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(a.png, dpi=150)
        print(f"wrote {a.png}")


if __name__ == "__main__":
    main()
