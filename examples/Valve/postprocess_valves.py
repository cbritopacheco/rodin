#!/usr/bin/env python3
"""
Post-processing of the KnifeGateValve runs (run_valves.sh): one PDF with

  1. a summary table of every case (time averages after the initial
     transient: F, F_dyn, F_direct, standard deviation and standard error of
     the mean, Q, dp, zeta), also written as CSV;
  2. the mean forces against the opening, both valves;
  3. per case: the force histories after the transient (F_dyn, F, Q, dp),
     and the velocity magnitude on the two longitudinal planes through the
     axis, x = 0 (plane of the stroke) and y = 0, over the whole pipe and
     zoomed on the valve, at several saved times and time-averaged.

Forces come from the case CSV (one row per time step, the same numbers as
the log); if a case has no CSV, the "F_knife = (...)" lines of case.log are
used instead. The transient is cut at t_avg of the run (case.env, i.e.
T_AVG_CONV D/U), or at --skip D/U.

Usage:
  python3 postprocess_valves.py runs/ [-o valves_report.pdf] [--times 4 | --at 0.5,1.0,1.5]
          [--all-times]
          [--skip 8] [--no-mean] [--field mag|uz] [--png]
Needs numpy, matplotlib, h5py (pip install h5py).
"""
import argparse
import csv
import glob
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.tri import Triangulation

BORE = {"FPC8": 0.2032, "FPC16": 0.382}  # m
XDMF_TET = 6


# ----------------------------------------------------------------------------
# Forces
# ----------------------------------------------------------------------------
def read_env(d):
    env = {}
    f = os.path.join(d, "case.env")
    if os.path.isfile(f):
        for line in open(f):
            k, _, v = line.strip().partition("=")
            if v:
                env[k] = float(v)
    return env


def read_history(d, tag):
    """Per-step history: the case CSV, or the log as a fallback."""
    f = os.path.join(d, tag + ".csv")
    if os.path.isfile(f):
        rows = list(csv.DictReader(open(f)))
        if rows:
            return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    log = os.path.join(d, "case.log")
    if not os.path.isfile(log):
        return None
    txt = re.sub(r"\x1b\[[0-9;]*m", "", open(log, errors="ignore").read())
    num = r"([-+0-9.eE]+)"
    vec = r"\(" + num + r", " + num + r", " + num + r"\)"
    t = [float(m) for m in re.findall(r"step \d+\s+t=" + num, txt)]
    dt = [float(m) for m in re.findall(r"\bdt=" + num, txt)]
    q = [float(m) for m in re.findall(r"Q=" + num, txt)]
    F = np.array(re.findall(r"F_knife = " + vec, txt), dtype=float)
    Fd = np.array(re.findall(r"F_dyn \(F - p_out int n\) = " + vec, txt), dtype=float)
    n = min(len(t), len(F))
    if n == 0:
        return None
    h = {"t": np.array(t[:n]), "dt": np.array(dt[:n]), "Q_in": np.array(q[:n])}
    for i, c in enumerate("xyz"):
        h["F" + c] = F[:n, i]
        h["F%s_dyn" % c] = Fd[:n, i] if len(Fd) >= n else F[:n, i]
    return h


def stats(t, w, x, nbatch=6):
    """dt-weighted mean, standard deviation and batch-means standard error."""
    W = w.sum()
    m = (w * x).sum() / W
    s = math.sqrt(max((w * (x - m) ** 2).sum() / W, 0.0))
    sem = float("nan")
    if len(x) >= 2 * nbatch:
        edges = np.linspace(t[0] - w[0], t[-1], nbatch + 1)
        means = []
        for a, b in zip(edges[:-1], edges[1:]):
            k = (t > a) & (t <= b + 1e-12)
            if w[k].sum() > 0:
                means.append((w[k] * x[k]).sum() / w[k].sum())
        if len(means) > 1:
            sem = float(np.std(means, ddof=1) / math.sqrt(len(means)))
    return m, s, sem


# ----------------------------------------------------------------------------
# XDMF / HDF5 (Rodin output: one uniform grid per rank and step)
# ----------------------------------------------------------------------------
def read_xdmf(path):
    """[(time, [ {mesh: (file, topo, geo), field: (file, path)} per rank ])]"""
    base = os.path.dirname(path)
    root = ET.parse(path).getroot()
    steps = []
    for step in root.iter("Grid"):
        if step.get("CollectionType") != "Spatial":
            continue
        t = float(step.find("Time").get("Value"))
        ranks = []
        for g in step.findall("Grid"):
            item = {}
            topo = g.find("Topology/DataItem").text.strip()
            geo = g.find("Geometry/DataItem").text.strip()
            item["topology"] = topo
            item["geometry"] = geo
            for a in g.findall("Attribute"):
                item[a.get("Name")] = a.find("DataItem").text.strip()
            ranks.append(item)
        steps.append((t, ranks))
    return base, steps


def h5(base, ref, cache=None):
    f, _, p = ref.partition(":")
    key = (f, p)
    if cache is not None and key in cache:
        return cache[key]
    with h5py.File(os.path.join(base, f), "r") as fh:
        a = np.array(fh[p])
    if cache is not None:
        cache[key] = a
    return a


def load_mesh(base, ranks, cache):
    """Vertices and tetrahedra of every rank (ghost cells duplicated: harmless)."""
    out = []
    for r in ranks:
        V = h5(base, r["geometry"], cache).reshape(-1, 3)
        topo = h5(base, r["topology"], cache).astype(np.int64)
        if topo.size % 5 == 0 and np.all(topo.reshape(-1, 5)[:, 0] == XDMF_TET):
            T = topo.reshape(-1, 5)[:, 1:]
        else:  # generic mixed list: walk it
            T, i = [], 0
            while i < topo.size:
                assert topo[i] == XDMF_TET, "only tetrahedra are supported"
                T.append(topo[i + 1:i + 5]); i += 5
            T = np.array(T)
        out.append((V, T))
    return out


def load_velocity(base, ranks):
    return [h5(base, r["velocity"]).reshape(-1, 3) for r in ranks]


# ----------------------------------------------------------------------------
# Plane slice of a P1 field on tetrahedra (exact linear interpolation)
# ----------------------------------------------------------------------------
EDGES = np.array([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])


def slice_tets(V, T, val, axis, c=0.0):
    """Triangles of the cut of the tetrahedra by the plane x_axis = c.
    Returns 2D points (z, other), values and triangles for tripcolor."""
    d = V[:, axis] - c
    d = np.where(np.abs(d) < 1e-12, 1e-12, d)  # no vertex exactly on the plane
    dt = d[T]
    cut = (dt.min(1) < 0) & (dt.max(1) > 0)
    T = T[cut]
    if len(T) == 0:
        return None
    dt = dt[cut]
    a, b = T[:, EDGES[:, 0]], T[:, EDGES[:, 1]]          # (n, 6)
    da, db = dt[:, EDGES[:, 0]], dt[:, EDGES[:, 1]]
    cross = (da * db) < 0
    s = np.where(cross, da / np.where(cross, da - db, 1.0), 0.0)
    P = V[a] + s[..., None] * (V[b] - V[a])              # (n, 6, 3)
    F = val[a] + s * (val[b] - val[a])                    # (n, 6)
    other = [k for k in range(3) if k != axis and k != 2][0]
    X, Y = P[..., 2], P[..., other]                       # plot: z horizontal
    X = np.where(cross, X, np.nan); Y = np.where(cross, Y, np.nan)
    # order the 3 or 4 points of each polygon by angle around its centroid
    cx, cy = np.nanmean(X, 1, keepdims=True), np.nanmean(Y, 1, keepdims=True)
    ang = np.where(cross, np.arctan2(Y - cy, X - cx), np.inf)
    o = np.argsort(ang, 1)
    X, Y, F = (np.take_along_axis(A, o, 1) for A in (X, Y, np.where(cross, F, np.nan)))
    npts = cross.sum(1)
    n = len(T)
    base = np.arange(n)[:, None] * 6
    tris = [base + np.array([0, 1, 2])]
    quad = npts == 4
    tris.append((base + np.array([0, 2, 3]))[quad])
    tris = np.vstack(tris)
    return X.ravel(), Y.ravel(), F.ravel(), tris


def plane_field(mesh, vel, axis, field):
    xs, ys, fs, ts, off = [], [], [], [], 0
    for (V, T), U in zip(mesh, vel):
        val = np.linalg.norm(U, axis=1) if field == "mag" else U[:, 2]
        r = slice_tets(V, T, val, axis)
        if r is None:
            continue
        X, Y, F, tri = r
        xs.append(X); ys.append(Y); fs.append(F); ts.append(tri + off); off += len(X)
    X, Y, F, tri = (np.concatenate(a) for a in (xs, ys, fs, ts))
    ok = np.isfinite(X)
    keep = ok[tri].all(1)
    return X, Y, F, tri[keep]


def draw_plane(ax, pf, vmin, vmax, cmap, xlim=None, title=""):
    X, Y, F, tri = pf
    Xs = np.where(np.isfinite(X), X, 0.0); Ys = np.where(np.isfinite(Y), Y, 0.0)
    Fs = np.where(np.isfinite(F), F, 0.0)
    tr = Triangulation(Xs, Ys, tri)
    im = ax.tripcolor(tr, Fs, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    ax.set_aspect("equal")
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=6)
    return im


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs")
    ap.add_argument("-o", "--out", default="valves_report.pdf")
    ap.add_argument("--times", type=int, default=4, help="saved times shown per case (after the transient)")
    ap.add_argument("--skip", type=float, default=None,
                    help="transient cut in D/U (default: t_avg of the run, case.env)")
    ap.add_argument("--at", default="",
                    help="explicit times in s, comma separated (nearest saved step); overrides --times")
    ap.add_argument("--all-times", action="store_true",
                    help="--times spread over the whole run, transient included")
    ap.add_argument("--no-mean", action="store_true", help="skip the time-averaged velocity")
    ap.add_argument("--field", default="mag", choices=["mag", "uz"], help="|u| or axial u_z")
    ap.add_argument("--png", action="store_true", help="also write one PNG per page")
    a = ap.parse_args()

    cases = sorted(d for d in glob.glob(os.path.join(a.runs, "*_open*")) if os.path.isdir(d))
    if not cases:
        sys.exit("no case folder (<valve>_openNNN) in " + a.runs)

    summary = []
    pdf = PdfPages(a.out)
    pngdir = os.path.splitext(a.out)[0] + "_png"
    if a.png:
        os.makedirs(pngdir, exist_ok=True)
    page = [0]

    def save(fig):
        pdf.savefig(fig)
        if a.png:
            page[0] += 1
            fig.savefig(os.path.join(pngdir, "page%03d.png" % page[0]), dpi=130)
        plt.close(fig)

    per_case = []
    for d in cases:
        tag = os.path.basename(d)
        valve, op = tag.split("_open")
        op = int(op)
        env = read_env(d)
        D = BORE.get(valve, env.get("D", 0.0) / 1000.0 or 1.0)
        tc = env.get("tc", 0.0)
        h = read_history(d, tag)
        if h is None:
            print("[%s] no force history, skipped" % tag)
            continue
        if op == 0:
            t0 = h["t"][-1] - 1e-12  # hydrostatic: the last state
        elif a.skip is not None and tc > 0:
            t0 = a.skip * tc
        else:
            t0 = env.get("tavg", 0.5 * h["t"][-1])
        k = h["t"] > t0 if op else np.arange(len(h["t"])) == len(h["t"]) - 1
        if not k.any():
            print("[%s] the run has not passed the transient (t_end = %.4g s < %.4g s): last half used"
                  % (tag, h["t"][-1], t0))
            t0 = 0.5 * h["t"][-1]
            k = h["t"] > t0
        t, w = h["t"][k], h["dt"][k]
        row = {"valve": valve, "opening_pct": op, "t_from_s": t0, "t_end_s": float(h["t"][-1]),
               "window_D_over_U": (w.sum() / tc) if tc > 0 else float("nan"), "steps": int(k.sum())}
        for key in ["Fx", "Fy", "Fz", "Fx_dyn", "Fy_dyn", "Fz_dyn", "Fx_direct", "Fy_direct", "Fz_direct"]:
            if key in h:
                m, s, e = stats(t, w, h[key][k])
                row[key + "_mean"], row[key + "_std"], row[key + "_sem"] = m, s, e
        for key in ["Q_in", "p_in_mean", "p_out_mean", "zeta"]:
            if key in h:
                row[key] = stats(t, w, h[key][k])[0]
        if "p_in_mean" in row:
            row["dp_Pa"] = row["p_in_mean"] - row["p_out_mean"]
        summary.append(row)
        per_case.append((d, tag, valve, op, D, tc, t0, h, k))

    # ---- summary table and CSV ------------------------------------------------
    csvout = os.path.splitext(a.out)[0] + "_summary.csv"
    keys = []
    for r in summary:
        keys += [x for x in r if x not in keys]
    with open(csvout, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        for r in summary:
            wr.writerow({x: ("%.6g" % v if isinstance(v, float) else v) for x, v in r.items()})

    fig = plt.figure(figsize=(11.69, 8.27))
    fig.suptitle("Válvulas de cuchillo — fuerza del fluido sobre la compuerta (promedios tras el transitorio)",
                 fontsize=12)
    cols = ["caso", "ventana\n[D/U]", "F_z [N]", "F_z,dyn [N]", "σ F_z [N]", "F_y,dyn [N]", "σ F_y [N]",
            "e.e. F_y,dyn [N]", "Q [m³/s]", "Δp [kPa]", "ζ"]
    cells = []
    for r in summary:
        g = lambda x, f="%.4g": (f % r[x]) if x in r and np.isfinite(r[x]) else "—"
        cells.append(["%s %d %%" % (r["valve"], r["opening_pct"]), g("window_D_over_U", "%.3g"),
                      g("Fz_mean"), g("Fz_dyn_mean"), g("Fz_std", "%.3g"), g("Fy_dyn_mean"),
                      g("Fy_std", "%.3g"), g("Fy_dyn_sem", "%.2g"), g("Q_in", "%.4g"),
                      ("%.4g" % (r["dp_Pa"] / 1e3)) if "dp_Pa" in r else "—", g("zeta", "%.3g")])
    ax = fig.add_subplot(111); ax.axis("off")
    tb = ax.table(cellText=cells, colLabels=cols, loc="upper center", cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1, 1.6)
    ax.text(0.0, 0.02,
            "F = fuerza del fluido sobre el cuchillo (reacción consistente), ejes: x lateral, y carrera (+ abre), "
            "z flujo.\nF_dyn = F − p_out ∫ n dS: sin la presión de línea uniforme. σ: desviación estándar en el tiempo; "
            "e.e.: error estándar del promedio (medias por bloques).\nTransitorio inicial descartado; 0 %: estado "
            "hidrostático final. Datos: " + os.path.basename(csvout),
            transform=ax.transAxes, fontsize=7, va="bottom")
    save(fig)

    # ---- mean forces vs opening ----------------------------------------------
    fig, axs = plt.subplots(1, 3, figsize=(11.69, 4.6))
    for valve, mk in [("FPC8", "o"), ("FPC16", "s")]:
        rs = sorted([r for r in summary if r["valve"] == valve], key=lambda r: r["opening_pct"])
        if not rs:
            continue
        x = [r["opening_pct"] for r in rs]
        for ax, key, lab in [(axs[0], "Fz", "F_z [kN]"), (axs[1], "Fz_dyn", "F_z,dyn [kN]"),
                             (axs[2], "Fy_dyn", "F_y,dyn [kN]")]:
            y = np.array([r.get(key + "_mean", np.nan) for r in rs]) / 1e3
            e = np.array([r.get(key.replace("_dyn", "") + "_std", np.nan) for r in rs]) / 1e3
            ax.errorbar(x, y, yerr=e, marker=mk, capsize=3, label=valve)
            ax.set_xlabel("apertura [%]"); ax.set_ylabel(lab); ax.grid(alpha=0.3)
    axs[0].legend()
    fig.suptitle("Fuerza media sobre la compuerta vs apertura (barras: ±σ temporal)", fontsize=11)
    fig.tight_layout(); save(fig)

    # ---- per case ---------------------------------------------------------------
    for d, tag, valve, op, D, tc, t0, h, k in per_case:
        # force histories
        fig, axs = plt.subplots(2, 2, figsize=(11.69, 8.27))
        tt = h["t"]
        for ax, keys_, ttl in [(axs[0, 0], ["Fz_dyn", "Fz"], "F_z"), (axs[0, 1], ["Fy_dyn"], "F_y,dyn"),
                               (axs[1, 0], ["Fx_dyn"], "F_x,dyn")]:
            for key in keys_:
                if key in h:
                    ax.plot(tt[k], h[key][k], lw=0.8, label=key)
                    ax.axhline(np.average(h[key][k], weights=h["dt"][k]), ls="--", lw=0.8, color="k")
            ax.set_title(ttl + " [N] (tras el transitorio)", fontsize=9); ax.set_xlabel("t [s]")
            ax.grid(alpha=0.3); ax.legend(fontsize=7)
        ax = axs[1, 1]
        ax.plot(tt, h["Q_in"] * 1e3, lw=0.8, label="Q [L/s]")
        if "p_in_mean" in h:
            ax2 = ax.twinx()
            ax2.plot(tt, (h["p_in_mean"] - h["p_out_mean"]) / 1e3, lw=0.8, color="tab:red", label="Δp [kPa]")
            ax2.set_ylabel("Δp [kPa]", color="tab:red")
        ax.axvline(t0, color="gray", ls=":", lw=1)
        ax.set_title("historia completa (línea punteada: fin del transitorio)", fontsize=9)
        ax.set_xlabel("t [s]"); ax.set_ylabel("Q [L/s]"); ax.grid(alpha=0.3)
        fig.suptitle("%s — %d %% de apertura: fuerzas" % (valve, op), fontsize=12)
        fig.tight_layout(); save(fig)

        # velocity slices
        xd = os.path.join(d, tag + ".xdmf")
        if not os.path.isfile(xd):
            print("[%s] no XDMF, velocity skipped" % tag)
            continue
        base, steps = read_xdmf(xd)
        if not steps:
            continue
        cache = {}
        mesh = load_mesh(base, steps[0][1], cache)
        after = [s for s in steps if s[0] > t0] or steps[-1:]
        if a.at:
            ts = np.array([s[0] for s in steps])
            shown = [steps[int(np.argmin(np.abs(ts - float(x))))] for x in a.at.split(",") if x.strip()]
        else:
            pool = steps if a.all_times else after
            idx = np.unique(np.linspace(0, len(pool) - 1, min(a.times, len(pool))).round().astype(int))
            shown = [pool[i] for i in idx]
        label = "|u| [m/s]" if a.field == "mag" else "u_z [m/s]"
        cmap = "jet" if a.field == "mag" else "RdBu_r"
        frames = []
        for t, ranks in shown:
            vel = load_velocity(base, ranks)
            frames.append(("t = %.4g s" % t, vel))
        if not a.no_mean and len(after) > 1:
            acc = None
            for t, ranks in after:
                vel = load_velocity(base, ranks)
                acc = vel if acc is None else [x + y for x, y in zip(acc, vel)]
            frames.append(("promedio temporal (%d instantes, t > %.4g s)" % (len(after), t0),
                           [x / len(after) for x in acc]))
        planes = []
        for name, vel in frames:
            planes.append((name, plane_field(mesh, vel, 0, a.field), plane_field(mesh, vel, 1, a.field)))
        vals = np.concatenate([np.concatenate([p[2], q[2]]) for _, p, q in planes])
        vals = vals[np.isfinite(vals)]
        if a.field == "mag":
            vmin, vmax = 0.0, float(np.percentile(vals, 99.8))
        else:
            v = float(np.percentile(np.abs(vals), 99.8)); vmin, vmax = -v, v
        zmin = min(float(np.nanmin(p[0])) for _, p, _ in planes)
        zmax = max(float(np.nanmax(p[0])) for _, p, _ in planes)
        for name, px, py in planes:
            fig = plt.figure(figsize=(11.69, 8.27))
            gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 2.2, 2.2], hspace=0.45)
            im = draw_plane(fig.add_subplot(gs[0]), px, vmin, vmax, cmap, (zmin, zmax),
                            "plano x = 0 (plano de la carrera, y vertical): toda la tubería")
            draw_plane(fig.add_subplot(gs[1]), py, vmin, vmax, cmap, (zmin, zmax),
                       "plano y = 0 (x lateral): toda la tubería")
            draw_plane(fig.add_subplot(gs[2]), px, vmin, vmax, cmap, (-1.0 * D, 3.0 * D),
                       "plano x = 0: zoom −1D … 3D")
            draw_plane(fig.add_subplot(gs[3]), py, vmin, vmax, cmap, (-1.0 * D, 3.0 * D),
                       "plano y = 0: zoom −1D … 3D")
            cb = fig.colorbar(im, ax=fig.axes, shrink=0.6, pad=0.02)
            cb.set_label(label)
            fig.suptitle("%s — %d %%: %s   (z [m] a lo largo del flujo, cuchillo en z = 0)"
                         % (valve, op, name), fontsize=11)
            save(fig)
        print("[%s] %d velocity pages" % (tag, len(planes)))

    pdf.close()
    print("PDF:", a.out)
    print("CSV:", csvout)


if __name__ == "__main__":
    main()
