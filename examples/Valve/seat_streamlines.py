#!/usr/bin/env python3
"""
Streamlines of the time-averaged flow on the plane x = 0, zoomed on the seat
(bottom of the bore, knife at z = 0), for the seat-study runs: one figure
per design with its openings side by side, plus one PNG per case.

Background: time-averaged |u| (log scale); lines: streamlines of the in-plane
components (u_z, u_y). Time average
as in seat_stats.py: the snapshots after t_avg of the run, or the last half
of those saved if the run has not reached it yet (marked "provisional").

Usage:
  python3 seat_streamlines.py runs_seat [-o streamlines] [--designs "cad rampup"]
                              [--openings "5 10 25 100"] [--zoom-z 40] [--zoom-y 25]
                              [--seal-z 12] [--seal-y 4] [--res 0.05] [--dpi 220]
Each figure has the zoom on the seat and, below, a closer detail of the
seal (finer grid, denser lines, and extra lines in pink seeded inside the
pocket). --res sets the grid of the plots in mm (the detail uses res/2).
Writes OUT/streamlines.pdf and OUT/<case>.png. Needs numpy, matplotlib, h5py.
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.tri import Triangulation, LinearTriInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import postprocess_valves as pv  # noqa: E402

BORE = {"FPC8": 101.6e-3, "FPC16": 191.0e-3}  # m, radius
ORDER = ("cad", "fillet", "chamfer", "upramp", "rampup", "rampdn", "lip")


def merged_mesh(meshes):
    """One vertex and tetrahedron list for all ranks (ghost copies removed)."""
    V = np.vstack([m[0] for m in meshes])
    off = np.cumsum([0] + [len(m[0]) for m in meshes[:-1]])
    T = np.vstack([m[1] + o for m, o in zip(meshes, off)])
    key = np.round(V * 1e8).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.reshape(-1)
    T = inv[T]
    _, ut = np.unique(np.sort(T, 1), axis=0, return_index=True)
    return V[first], T[ut], inv, len(first)


def mean_velocity(base, steps, inv, nv):
    U = np.zeros((nv, 3)); n = 0
    for t, ranks in steps:
        Ur = np.vstack([pv.h5(base, r["velocity"]).reshape(-1, 3) for r in ranks])
        Ui = np.zeros((nv, 3)); Ui[inv] = Ur  # ghost copies carry the same values
        U += Ui; n += 1
    return U / n


def plane(V, T, U):
    """Slice x = 0: unique points (z, y), triangles, (u_z, u_y, |u|) at the points."""
    out = []
    for val in (U[:, 2], U[:, 1], np.linalg.norm(U, axis=1)):
        X, Y, F, tri = pv.slice_tets(V, T, val, 0)
        out.append(F)
    ok = np.isfinite(X)
    tri = tri[ok[tri].all(1)]
    key = np.round(np.c_[np.where(ok, X, 0), np.where(ok, Y, 0)] * 1e8).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.reshape(-1)
    tri = inv[tri]
    P = np.c_[X, Y][first]
    a = P[tri[:, 1]] - P[tri[:, 0]]; b = P[tri[:, 2]] - P[tri[:, 0]]
    tri = tri[np.abs(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]) > 1e-16]
    return P, tri, [F[first] for F in out]


def draw(ax, case, zoom, res=0.05e-3, density=1.4, nseed=40, cbar=True):
    P, tri, (uz, uy, um), R, label = case
    (z0, z1), (y0, y1) = zoom
    tr = Triangulation(P[:, 0], P[:, 1], tri)
    s = (P[:, 0] > z0) & (P[:, 0] < z1) & (P[:, 1] > y0) & (P[:, 1] < y1)
    vmax = max(np.nanmax(um[s]), 1e-6) if s.any() else 1.0
    # grid of res [m] (the streamlines and the background are drawn on it)
    zg = np.linspace(z0, z1, int(round((z1 - z0) / res)) + 1)
    yg = np.linspace(y0, y1, int(round((y1 - y0) / res)) + 1)
    Z, Y = np.meshgrid(zg, yg)
    gz = LinearTriInterpolator(tr, uz)(Z, Y); gy = LinearTriInterpolator(tr, uy)(Z, Y)
    gm = LinearTriInterpolator(tr, um)(Z, Y)
    im = ax.pcolormesh(Z, Y, np.ma.masked_invalid(np.clip(gm, vmax * 1e-3, None)), shading="gouraud",
                       cmap="viridis", norm=LogNorm(vmin=vmax * 1e-3, vmax=vmax), rasterized=True)
    GZ, GY = gz.filled(np.nan), gy.filled(np.nan)
    ax.streamplot(Z, Y, GZ, GY, color="w", linewidth=0.6, density=density, arrowsize=0.6)
    # extra lines seeded in the seat pocket (below the bore line), where the
    # flow is slow and the even seeding of streamplot leaves it nearly empty
    R = case[3]
    pk = np.isfinite(GZ) & (Y < -R - 0.3e-3)
    if nseed and pk.sum() > 0:
        idx = np.flatnonzero(pk.ravel())
        pick = idx[np.linspace(0, len(idx) - 1, min(nseed, len(idx))).astype(int)]
        seeds = np.c_[Z.ravel()[pick], Y.ravel()[pick]]
        ax.streamplot(Z, Y, GZ, GY, color="#ffd0d0", linewidth=0.5, start_points=seeds,
                      arrowsize=0.5, integration_direction="both", broken_streamlines=False)
    ax.axhline(-R, color="r", lw=0.6, ls="--")
    ax.axvline(0, color="0.6", lw=0.4, ls=":")
    ax.set_xlim(z0, z1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
    ax.set_facecolor("0.85")
    ax.set_xlabel("z [m]", fontsize=8); ax.set_ylabel("y [m]", fontsize=8); ax.tick_params(labelsize=7)
    ax.set_title(label, fontsize=9)
    if cbar:
        c = plt.colorbar(im, ax=ax, shrink=0.75); c.set_label("|ū| [m/s]", fontsize=7); c.ax.tick_params(labelsize=6)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs")
    ap.add_argument("-o", "--out", default="streamlines")
    ap.add_argument("--designs", default="cad fillet rampup rampdn lip")
    ap.add_argument("--openings", default="5 10 25 100")
    ap.add_argument("--valve", default="FPC8")
    ap.add_argument("--zoom-z", type=float, default=40.0, help="half-width of the zoom along the pipe, mm")
    ap.add_argument("--zoom-y", type=float, default=25.0, help="height of the zoom above the bore bottom, mm")
    ap.add_argument("--seal-z", type=float, default=12.0, help="half-width of the seal detail, mm")
    ap.add_argument("--seal-y", type=float, default=4.0, help="height of the seal detail above the bore bottom, mm")
    ap.add_argument("--res", type=float, default=0.05, help="grid of the plots, mm (seal detail: res/2)")
    ap.add_argument("--dpi", type=int, default=220)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    designs = a.designs.split(); openings = [int(o) for o in a.openings.split()]
    R = BORE[a.valve]
    zoom = ((-a.zoom_z * 1e-3, a.zoom_z * 1e-3), (-R - 0.006, -R + a.zoom_y * 1e-3))
    seal = ((-a.seal_z * 1e-3, a.seal_z * 1e-3), (-R - 0.006, -R + a.seal_y * 1e-3))
    res, sres = a.res * 1e-3, 0.5 * a.res * 1e-3

    cases = {}
    for d_ in designs:
        name = a.valve if d_ == "cad" else a.valve + "-" + d_
        for o in openings:
            tag = "%s_open%03d" % (name, o)
            d = os.path.join(a.runs, tag)
            xd = os.path.join(d, tag + ".xdmf")
            if not os.path.isfile(xd):
                print("  (no %s)" % xd); continue
            env = pv.read_env(d); t0 = env.get("tavg", 0.0)
            base, steps = pv.read_xdmf(xd)
            if not steps:
                print("  (%s: no snapshots)" % tag); continue
            after = [s for s in steps if s[0] > t0 + 1e-12]
            prov = len(after) < 2 or not pv.finished(d)
            steps = after if len(after) >= 2 else steps[len(steps) // 2:]
            V, T, inv, nv = merged_mesh(pv.load_mesh(base, steps[0][1], {}))
            U = mean_velocity(base, steps, inv, nv)
            P, tri, F = plane(V, T, U)
            label = "%s, %d %% (%d snapshots, t = %.3g–%.3g s)%s" % (
                d_, o, len(steps), steps[0][0], steps[-1][0], " — provisorio" if prov else "")
            cases[(d_, o)] = (P, tri, F, R, label)
            fig, axs = plt.subplots(2, 1, figsize=(9, 9))
            draw(axs[0], cases[(d_, o)], zoom, res, nseed=0)
            draw(axs[1], cases[(d_, o)], seal, sres, density=2.5, nseed=60)
            axs[1].set_title("detalle del sello", fontsize=9)
            fig.suptitle("%s — plano x = 0, líneas de corriente del flujo medio en el asiento" % tag, fontsize=10)
            fig.tight_layout(); fig.savefig(os.path.join(a.out, tag + ".png"), dpi=a.dpi); plt.close(fig)
            print("  %s ok%s" % (tag, " (provisional)" if prov else ""))
    if not cases:
        sys.exit("no case found in " + a.runs)

    pdf = PdfPages(os.path.join(a.out, "streamlines.pdf"))
    for d_ in [x for x in ORDER if x in designs] + [x for x in designs if x not in ORDER]:
        ops = [o for o in openings if (d_, o) in cases]
        if not ops:
            continue
        fig, axs = plt.subplots(2, len(ops), figsize=(4.6 * len(ops), 7.6), squeeze=False)
        for j, o in enumerate(ops):
            draw(axs[0, j], cases[(d_, o)], zoom, res, nseed=0)
            draw(axs[1, j], cases[(d_, o)], seal, sres, density=2.5, nseed=60)
            axs[1, j].set_title("detalle del sello", fontsize=8)
        fig.suptitle("Asiento %s — plano x = 0, líneas de corriente del flujo medio (fondo: |ū|, línea roja: diámetro)"
                     % d_, fontsize=11)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    for o in openings:  # the designs side by side at each opening
        ds = [d_ for d_ in designs if (d_, o) in cases]
        if len(ds) < 2:
            continue
        fig, axs = plt.subplots(2, len(ds), figsize=(4.6 * len(ds), 7.6), squeeze=False)
        for j, d_ in enumerate(ds):
            draw(axs[0, j], cases[(d_, o)], zoom, res, nseed=0)
            draw(axs[1, j], cases[(d_, o)], seal, sres, density=2.5, nseed=60)
            axs[1, j].set_title("%s — detalle del sello" % d_, fontsize=8)
        fig.suptitle("%d %% de apertura — diseños del asiento, plano x = 0" % o, fontsize=11)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    pdf.close()
    print("PDF:", os.path.join(a.out, "streamlines.pdf"))


if __name__ == "__main__":
    main()
