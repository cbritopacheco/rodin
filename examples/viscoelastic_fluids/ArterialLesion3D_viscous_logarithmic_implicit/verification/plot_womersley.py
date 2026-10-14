#!/usr/bin/env python3
"""Convergence plots of run_womersley.sh: <OUT>/summary.csv -> <OUT>/womersley.pdf

Left: time study, d_N = ||x_N - x_2N|| (self-convergence) for BDF1 and BDF2,
with slopes 1 and 2. Right: space study, error against the exact solution
versus h/D, with slope 2.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

out = sys.argv[1]
rows = list(csv.DictReader(open(os.path.join(out, "summary.csv"))))
num = lambda r, k: float(r[k]) if r.get(k) not in (None, "") else np.nan

fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)

# time study: the mesh with the most runs
by_mesh = {}
for r in rows:
    by_mesh.setdefault(r["rf"], []).append(r)
tm = max(by_mesh, key=lambda k: len(by_mesh[k]))
for o, mk in (("1", "o"), ("2", "s")):
    rs = sorted([r for r in by_mesh[tm] if r["order"] == o], key=lambda r: int(r["N"]))
    dt = np.array([num(r, "dt") for r in rs])
    for key, ls, lab in (("dself_u", "-", "u"), ("dself_sigma", "--", r"\sigma")):
        d = np.array([num(r, key) for r in rs])
        ok = np.isfinite(d)
        if ok.sum() >= 1:
            ax[0].loglog(dt[ok], d[ok], ls, marker=mk, label=r"BDF%s, $%s$" % (o, lab))
if rows:
    dts = np.array(sorted({num(r, "dt") for r in by_mesh[tm]}))
    if len(dts) > 1:
        ref = dts / dts.max()
        ax[0].loglog(dts, 0.5 * ref, "k:", lw=0.8)
        ax[0].loglog(dts, 0.5 * ref ** 2, "k-.", lw=0.8)
        ax[0].text(dts.max(), 0.5, "1", fontsize=8)
        ax[0].text(dts.max(), 0.25, "2", fontsize=8)
ax[0].set_xlabel(r"$\Delta t$ (s)")
ax[0].set_ylabel(r"$\|x_{\Delta t} - x_{\Delta t/2}\| / \|x\|$")
ax[0].set_title("time, mesh rf = %s" % tm, fontsize=9)
ax[0].legend(fontsize=7, frameon=False)

# space study: the finest dt of BDF2 on every mesh
n_max = max(int(r["N"]) for r in rows)
sp = sorted([r for r in rows if r["order"] == "2" and int(r["N"]) == n_max], key=lambda r: float(r["rf"]))
if len(sp) >= 1:
    h = np.array([0.1 / float(r["rf"]) for r in sp])
    for key, mk, lab in (("err_u", "o", "u"), ("err_sigma", "s", r"\sigma"),
                         ("err_sxr", "^", r"\sigma_{xr}"), ("err_sxx", "v", r"\sigma_{xx}")):
        ax[1].loglog(h, [num(r, key) for r in sp], "-", marker=mk, label="$%s$" % lab)
    if len(h) > 1:
        e0 = num(sp[0], "err_u")
        ax[1].loglog(h, e0 * (h / h[0]) ** 2, "k-.", lw=0.8, label="slope 2")
ax[1].set_xlabel(r"$h/D$")
ax[1].set_ylabel("relative $L^2$ error (vs exact)")
ax[1].set_title("space, BDF2, %d steps per cycle" % n_max, fontsize=9)
ax[1].legend(fontsize=7, frameon=False)
fig.savefig(os.path.join(out, "womersley.pdf"))
fig.savefig(os.path.join(out, "womersley.png"), dpi=150)
