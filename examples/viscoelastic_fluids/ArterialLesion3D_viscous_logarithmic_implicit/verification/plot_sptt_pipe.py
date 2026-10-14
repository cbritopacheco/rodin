#!/usr/bin/env python3
"""Profiles u_x, sigma_xr and sigma_xx across the pipe, discrete (finest mesh of
each case, nodes in the developed region) against the exact solution.

    python3 plot_sptt_pipe.py <OUT dir of run_sptt_pipe.sh>   ->  <OUT>/profiles.pdf
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

out = sys.argv[1]
cases = {}
for f in glob.glob(os.path.join(out, "*_rf*", "profile.csv")):
    c, rf = os.path.basename(os.path.dirname(f)).split("_rf")
    if c not in cases or float(rf) > cases[c][0]:
        cases[c] = (float(rf), f)
if not cases:
    sys.exit("no profile.csv found")

cols = [("u_h", "u_exact", r"$u_x$ (m/s)"), ("sxr_h", "sxr_exact", r"$\sigma_{xr}$ (Pa)"),
        ("sxx_h", "sxx_exact", r"$\sigma_{xx}$ (Pa)")]
fig, ax = plt.subplots(1, 3, figsize=(10, 3.2), constrained_layout=True)
for k, (c, (rf, f)) in enumerate(sorted(cases.items())):
    d = np.genfromtxt(f, delimiter=",", names=True)
    x = np.genfromtxt(f.replace("profile.csv", "profile_exact.csv"), delimiter=",", names=True)
    for j, (h, e, lab) in enumerate(cols):
        ax[j].plot(d["r_over_R"], d[h], ".", ms=2, alpha=0.35, color="C%d" % k)
        ax[j].plot(x["r_over_R"], x[e], "-", lw=1.4, color="C%d" % k,
                   label="%s (exact; dots: h/D = %.3g)" % (c, 0.1 / rf))
        ax[j].set_xlabel(r"$r/R$")
        ax[j].set_ylabel(lab)
ax[0].legend(fontsize=7, frameon=False)
fig.savefig(os.path.join(out, "profiles.pdf"))
fig.savefig(os.path.join(out, "profiles.png"), dpi=150)
