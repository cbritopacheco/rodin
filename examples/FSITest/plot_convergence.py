#!/usr/bin/env python3
"""Plot the output of FSITest: errors vs h (log-log) and error histories.

Usage (from the directory where FSITest was run):
    python3 plot_convergence.py [fsitest_convergence.csv] [results_fsitest]
Writes fsitest_convergence.png.
"""
import csv
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

conv = sys.argv[1] if len(sys.argv) > 1 else "fsitest_convergence.csv"
hist = sys.argv[2] if len(sys.argv) > 2 else "results_fsitest"

rows = list(csv.DictReader(open(conv)))
h = [1.0 / float(r["N"]) for r in rows]
series = [("e_u", r"$\|u-u_h\|_{L^2}$"), ("e_p", r"$\|p-p_h\|_{L^2}$"),
          ("e_d", r"$\|d-d_h\|_{L^2}$"), ("e_dt", r"$\|\dot d-\dot d_h\|_{L^2}$"),
          ("e_d_sigma", r"$\|d-d_h\|_{L^2(\Sigma)}$")]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
for key, label in series:
    ax1.loglog(h, [float(r[key]) for r in rows], "o-", label=label)
e0 = float(rows[0]["e_u"])
ax1.loglog(h, [e0 * (x / h[0]) for x in h], "k--", lw=1, label=r"$O(h + \Delta t)$")
ax1.set_xlabel(r"$h$  ($\Delta t \propto h$)")
ax1.set_ylabel("error at $t = T$")
ax1.grid(True, which="both", alpha=0.3)
ax1.legend(fontsize=8)

for i, f in enumerate(sorted(glob.glob(os.path.join(hist, "history_L*.csv")))):
    data = list(csv.DictReader(open(f)))
    color = f"C{i}"
    t = [float(r["t"]) for r in data]
    lvl = os.path.basename(f)[len("history_"):-4]
    ax2.semilogy(t, [float(r["e_u"]) for r in data], "-", color=color, label=f"{lvl} fluid")
    ax2.semilogy(t, [float(r["e_d"]) for r in data], ":", color=color, label=f"{lvl} solid")
ax2.set_xlabel("$t$")
ax2.set_ylabel(r"$L^2$ error")
ax2.grid(True, which="both", alpha=0.3)
ax2.legend(fontsize=7, ncol=2)

fig.tight_layout()
fig.savefig("fsitest_convergence.png", dpi=150)
print("wrote fsitest_convergence.png")
