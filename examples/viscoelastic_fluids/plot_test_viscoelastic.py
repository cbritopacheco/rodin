#!/usr/bin/env python3
"""Plot the output of test_viscoelastic_implicit: error curves in h and in dt.

Usage (from the directory where the test was run):
    python3 plot_test_viscoelastic.py [vt_space.csv] [vt_time.csv]
Writes vt_convergence.png. A missing CSV leaves its panel out.

(a) space: Pacheco & Castillo (2023), Sec. 5.1.4, relative errors at t = T
    against h, with dt refined together with h (as their Fig. 5).
(b) time: exact pulsatile Oldroyd-B channel on a fixed mesh, relative errors at
    t = T against dt; solid lines against the dt_min/4 reference (temporal error
    alone), dotted lines against the exact solution (they stall at the spatial
    error of the mesh).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

space = sys.argv[1] if len(sys.argv) > 1 else "vt_space.csv"
time = sys.argv[2] if len(sys.argv) > 2 else "vt_time.csv"

SERIES = [("u_L2", r"$\mathbf{u}$, $L^2$", "o"), ("u_H1", r"$\mathbf{u}$, $H^1$", "s"),
          ("p_L2", r"$p$, $L^2$", "^"), ("sigma_L2", r"$\boldsymbol{\sigma}$, $L^2$", "v")]


def read(path):
    if not os.path.exists(path):
        return None
    return list(csv.DictReader(open(path)))


def slope(ax, x, y0, k, label):
    """Reference line of slope k through (x[0], y0)."""
    ax.loglog(x, [y0 * (xi / x[0]) ** k for xi in x], "k--", lw=0.8, label=label)


panels = [(read(space), "space"), (read(time), "time")]
panels = [(rows, kind) for rows, kind in panels if rows]
if not panels:
    sys.exit("no vt_space.csv or vt_time.csv found")

fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 4.2), squeeze=False)
for ax, (rows, kind) in zip(axes[0], panels):
    if kind == "space":
        x = [float(r["h"]) for r in rows]
        for i, (key, label, marker) in enumerate(SERIES):
            ax.loglog(x, [float(r["e_" + key]) for r in rows], marker + "-", color=f"C{i}",
                      label=label)
        first = min(float(rows[0]["e_" + k]) for k, _, _ in SERIES)
        slope(ax, x, first, 1, r"$O(h)$")
        slope(ax, x, 0.5 * first, 2, r"$O(h^2)$")
        ax.set_xlabel(r"$h$  ($\Delta t \propto h$)")
        ax.set_ylabel(r"relative error at $t = T$")
        ax.set_title("(a) ramped channel, Pacheco & Castillo (2023)", fontsize=9)
    else:
        x = [float(r["dt"]) for r in rows]
        for i, (key, label, marker) in enumerate(SERIES):
            ax.loglog(x, [float(r["ref_" + key]) for r in rows], marker + "-", color=f"C{i}",
                      label=label + " (vs reference)")
            ax.loglog(x, [float(r["ex_" + key]) for r in rows], marker + ":", color=f"C{i}",
                      mfc="none", alpha=0.7)
        first = min(float(rows[0]["ref_" + k]) for k, _, _ in SERIES)
        slope(ax, x, first, 1, r"$O(\Delta t)$")
        ax.set_xlabel(r"$\Delta t$  (fixed mesh)")
        ax.set_ylabel(r"relative error at $t = T$")
        ax.set_title("(b) pulsatile Oldroyd-B channel (exact)\n"
                     "solid: vs reference, dotted: vs exact", fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)

fig.tight_layout()
fig.savefig("vt_convergence.png", dpi=150)
print("wrote vt_convergence.png")
