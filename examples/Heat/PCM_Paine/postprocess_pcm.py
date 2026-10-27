#!/usr/bin/env python3
"""Post-processing of the NaturalConvectionPCM thickness sweep (run_pcm_sweep.sh).

For every case OUT/e_XXmm (NaturalConvectionPCM.csv + .xdmf/.h5) it analyses
the last simulated day (periodic regime) and writes to OUT/post:
  summary.csv, summary.md   heat-flux summary per thickness
  flujos_dia.png            q_bottom (PCM -> attic) and q_ceiling (-> dwelling) over the day
  resumen_espesor.png       peaks, daily energies, lag and melting vs thickness
  fraccion_liquida.png      mean liquid fraction over the day
  T_contornos_e_XXmm.png    temperature at several hours (+ melt front Y = 0.5)

Sign convention (W/m^2 of roof): q_top > 0 enters the PCM through the roof
surface, q_bottom > 0 leaves the PCM into the attic, q_ceiling > 0 enters the
dwelling through the ceiling.

    python3 postprocess_pcm.py OUT [--hours 9 13 17 21] [--cases 5 25 50]
"""
import argparse
import csv
import glob
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402

DAY = 86400.0


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def case_info(d):
    info = {}
    p = os.path.join(d, "case.txt")
    if os.path.exists(p):
        for tok in open(p).read().split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                info[k] = float(v)
    if "thickness_mm" not in info:
        info["thickness_mm"] = float(re.findall(r"(\d+)mm", d)[-1])
    return info


def energy(t, q, positive=None):
    """Integral of q dt in Wh/m^2 (trapezoid); positive/negative part if asked."""
    if positive is True:
        q = np.maximum(q, 0.0)
    elif positive is False:
        q = np.minimum(q, 0.0)
    return float(np.sum(0.5 * (q[1:] + q[:-1]) * np.diff(t)) / 3600.0)


def last_day(d):
    t = d["t"]
    nd = max(1, int(math.floor(t[-1] / DAY + 1e-6)))
    t1 = nd * DAY
    t0 = t1 - DAY
    m = (t >= t0 - 1e-6) & (t <= t1 + 1e-6)
    prev = (t >= t0 - DAY - 1e-6) & (t <= t0 + 1e-6)
    return m, prev if nd >= 2 else None, t0, nd


def analyse(d, info):
    m, prev, t0, nd = last_day(d)
    t = d["t"][m]
    h = (t - t0) / 3600.0
    qb, qc, qt = d["qBottom"][m], d["qCeiling"][m], d["qTop"][m]
    ib, ic, it = int(np.argmax(qb)), int(np.argmax(qc)), int(np.argmax(qt))
    r = {
        "e_mm": info["thickness_mm"],
        "malla": "%dx%d" % (info.get("nx", 0), info.get("ny", 0)),
        "dias": nd,
        "dt": info.get("dt", float("nan")),
        "mushy": info.get("mushy", 1.0),
        "qTop_max": qt[it],
        "E_top_neto": energy(t, qt),
        "qBottom_max": qb[ib],
        "h_qBottom_max": h[ib],
        "qBottom_min": float(qb.min()),
        "qBottom_medio": energy(t, qb) / 24.0,
        "E_atico_pos": energy(t, qb, True),
        "E_atico_neto": energy(t, qb),
        "qCeiling_max": qc[ic],
        "h_qCeiling_max": h[ic],
        "E_vivienda_pos": energy(t, qc, True),
        "E_vivienda_neto": energy(t, qc),
        "desfase_h": h[ib] - h[it],
        "Tsup_max": float(d["topT_C"][m].max()),
        "Tinf_max": float(d["bottomT_C"][m].max()),
        "Tatico_max": float(d["Tattic_C"][m].max()),
        "Y_max": float(d["liquidFraction"][m].max()),
        "Y_min": float(d["liquidFraction"][m].min()),
        "horas_fundido": float(np.sum(np.diff(t)[d["liquidFraction"][m][1:] > 0.99]) / 3600.0),
        "u_max_mm_s": 1e3 * float(d["maxU"][m].max()),
        "newton_max": int(d["newtonIterations"][m].max()),
    }
    if prev is not None:
        e_prev = energy(d["t"][prev], d["qBottom"][prev])
        r["periodicidad_pct"] = 100.0 * abs(r["E_atico_neto"] - e_prev) / max(abs(r["E_atico_neto"]), 1e-9)
    else:
        r["periodicidad_pct"] = float("nan")
    return r, h, m


COLUMNS = [
    ("e_mm", "e [mm]", "%.0f"),
    ("malla", "malla", "%s"),
    ("dt", "dt [s]", "%g"),
    ("mushy", "T_l-T_s [K]", "%g"),
    ("qTop_max", "q_sup max [W/m2]", "%.0f"),
    ("E_top_neto", "E_sup neta [Wh/m2]", "%.0f"),
    ("qBottom_max", "q_inf max [W/m2]", "%.1f"),
    ("h_qBottom_max", "hora q_inf max [h]", "%.1f"),
    ("qBottom_min", "q_inf min [W/m2]", "%.1f"),
    ("qBottom_medio", "q_inf medio [W/m2]", "%.1f"),
    ("E_atico_pos", "E al entretecho [Wh/m2]", "%.0f"),
    ("E_atico_neto", "E neta entretecho [Wh/m2]", "%.0f"),
    ("qCeiling_max", "q_cielo max [W/m2]", "%.2f"),
    ("E_vivienda_pos", "E a vivienda [Wh/m2]", "%.1f"),
    ("desfase_h", "desfase [h]", "%.1f"),
    ("Tsup_max", "T_sup max [C]", "%.1f"),
    ("Tinf_max", "T_inf max [C]", "%.1f"),
    ("Tatico_max", "T_entretecho max [C]", "%.1f"),
    ("Y_max", "Y max [-]", "%.2f"),
    ("Y_min", "Y min [-]", "%.2f"),
    ("horas_fundido", "h con Y>0.99", "%.1f"),
    ("u_max_mm_s", "|u| max [mm/s]", "%.3f"),
    ("periodicidad_pct", "cambio dia-1 [%]", "%.1f"),
]


def write_summary(rows, out):
    with open(os.path.join(out, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([c[0] for c in COLUMNS] + ["newton_max"])
        for r in rows:
            w.writerow([r[c[0]] for c in COLUMNS] + [r["newton_max"]])
    ref = rows[0]["qBottom_max"]
    lines = [
        "# Barrido de espesor PCM (RT42), Paine enero, ultimo dia simulado",
        "",
        "q_sup > 0: entra al PCM por el techo; q_inf > 0: sale del PCM al entretecho; "
        "q_cielo > 0: entra a la vivienda por el cielo. Energias por m2 de techo y dia. "
        "Desfase: hora del maximo de q_inf menos hora del maximo de q_sup.",
        "",
        "| " + " | ".join(c[1] for c in COLUMNS) + " | reduccion pico vs %g mm |" % rows[0]["e_mm"],
        "|" + "---|" * (len(COLUMNS) + 1),
    ]
    for r in rows:
        red = 100.0 * (1.0 - r["qBottom_max"] / ref) if ref else float("nan")
        lines.append("| " + " | ".join(c[2] % r[c[0]] for c in COLUMNS) + " | %.0f %% |" % red)
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    # Console version (subset).
    keys = ["e_mm", "qBottom_max", "h_qBottom_max", "E_atico_pos", "E_vivienda_pos",
            "desfase_h", "Tatico_max", "Y_max", "periodicidad_pct"]
    names = {c[0]: c[1] for c in COLUMNS}
    fmts = {c[0]: c[2] for c in COLUMNS}
    width = [max(len(names[k]), 8) for k in keys]
    print("  ".join(names[k].rjust(wd) for k, wd in zip(keys, width)))
    for r in rows:
        print("  ".join((fmts[k] % r[k]).rjust(wd) for k, wd in zip(keys, width)))


def style():
    plt.rcParams.update({
        "figure.dpi": 150, "font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25,
        "legend.frameon": False,
    })


def colors(n):
    cmap = plt.get_cmap("viridis")
    return [cmap(0.05 + 0.85 * i / max(n - 1, 1)) for i in range(n)]


def plot_fluxes(cases, out):
    c = colors(len(cases))
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 8.2), sharex=True)
    for k, (r, d, h, m) in enumerate(cases):
        lab = "%g mm" % r["e_mm"]
        ax[0].plot(h, d["qTop"][m], color=c[k], lw=1.2, label=lab)
        ax[1].plot(h, d["qBottom"][m], color=c[k], lw=1.4, label=lab)
        ax[2].plot(h, d["qCeiling"][m], color=c[k], lw=1.4, label=lab)
    r0, d0, h0, m0 = cases[0]
    axI = ax[0].twinx()
    axI.plot(h0, d0["I"][m0], color="0.6", lw=0.8, ls="--")
    axI.set_ylabel("I techo [W/m$^2$] (- -)", color="0.45")
    axI.spines["top"].set_visible(False)
    axI.grid(False)
    ax[0].set_ylabel("q$_{sup}$ [W/m$^2$]\n(techo $\\rightarrow$ PCM)")
    ax[1].set_ylabel("q$_{inf}$ [W/m$^2$]\n(PCM $\\rightarrow$ entretecho)")
    ax[2].set_ylabel("q$_{cielo}$ [W/m$^2$]\n(entretecho $\\rightarrow$ vivienda)")
    ax[2].set_xlabel("hora solar [h]")
    ax[2].set_xticks(range(0, 25, 3))
    ax[2].set_xlim(0, 24)
    for a in ax:
        a.axhline(0, color="0.3", lw=0.6)
    ax[1].legend(ncol=5, fontsize=7, loc="upper left")
    fig.suptitle("Flujos de calor en la capa de PCM, ultimo dia (Paine, enero)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "flujos_dia.png"))
    plt.close(fig)


def plot_vs_thickness(rows, out):
    e = np.array([r["e_mm"] for r in rows])
    g = lambda k: np.array([r[k] for r in rows])  # noqa: E731
    fig, ax = plt.subplots(2, 2, figsize=(8.5, 6.2))
    a = ax[0, 0]
    a.plot(e, g("qBottom_max"), "o-", color="#3b5ba5", label="maximo")
    a.plot(e, g("qBottom_medio"), "s--", color="#e07b39", label="medio diario")
    a.set_ylabel("q$_{inf}$ [W/m$^2$]")
    a.set_title("Flujo PCM $\\rightarrow$ entretecho", fontsize=9)
    a.legend()
    a = ax[0, 1]
    a.plot(e, g("E_atico_pos"), "o-", color="#3b5ba5", label="al entretecho (q$_{inf}$>0)")
    a.plot(e, g("E_vivienda_pos") * 10, "s--", color="#e07b39", label="a la vivienda x10")
    a.set_ylabel("energia diaria [Wh/m$^2$]")
    a.set_title("Energia diaria transmitida", fontsize=9)
    a.legend()
    a = ax[1, 0]
    a.plot(e, g("desfase_h"), "o-", color="#3b5ba5")
    a.set_ylabel("desfase del pico [h]")
    a.set_xlabel("espesor PCM [mm]")
    a.set_title("Retraso del maximo de q$_{inf}$ respecto de q$_{sup}$", fontsize=9)
    a = ax[1, 1]
    a.plot(e, g("Y_max"), "o-", color="#3b5ba5", label="Y max")
    a.plot(e, g("Y_min"), "s--", color="#e07b39", label="Y min")
    a.set_ylim(-0.05, 1.05)
    a.set_ylabel("fraccion liquida media [-]")
    a.set_xlabel("espesor PCM [mm]")
    a.set_title("Fusion del PCM en el dia", fontsize=9)
    a.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "resumen_espesor.png"))
    plt.close(fig)


def plot_liquid(cases, out):
    c = colors(len(cases))
    fig, ax = plt.subplots(figsize=(7.5, 3.2))
    for k, (r, d, h, m) in enumerate(cases):
        ax.plot(h, d["liquidFraction"][m], color=c[k], lw=1.4, label="%g mm" % r["e_mm"])
    ax.set_xlabel("hora solar [h]")
    ax.set_ylabel("fraccion liquida media [-]")
    ax.set_xticks(range(0, 25, 3))
    ax.set_xlim(0, 24)
    ax.legend(ncol=5, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fraccion_liquida.png"))
    plt.close(fig)


# ---------------------------------------------------------------- XDMF ----
def _dataitem(base, item):
    text = item.text.strip()
    fname, path = text.split(":", 1)
    import h5py
    with h5py.File(os.path.join(base, fname), "r") as f:
        return np.array(f[path])


def read_xdmf(path):
    """[(time, [(xy, tri, fields{name: nodal values}) per rank])]."""
    base = os.path.dirname(path)
    root = ET.parse(path).getroot()
    steps = []
    temporal = root.find("Domain/Grid")
    for step in temporal.findall("Grid"):
        t = float(step.find("Time").get("Value"))
        steps.append((t, step))
    return base, steps


def load_step(base, step, names=("temperature", "liquid_fraction")):
    parts = []
    for g in step.findall("Grid"):
        topo = g.find("Topology")
        conn = _dataitem(base, topo.find("DataItem")).astype(np.int64).ravel()
        if topo.get("TopologyType", "").lower() == "mixed":
            ne = int(topo.get("NumberOfElements"))
            tri = []
            i = 0
            for _ in range(ne):
                typ = conn[i]
                nv = {4: 3, 5: 4, 1: 1, 2: 2}.get(int(typ), 3)
                if typ == 2:  # polyline: type, count, nodes
                    nv = int(conn[i + 1]) + 1
                if typ == 4:
                    tri.append(conn[i + 1:i + 4])
                i += 1 + nv
            tri = np.array(tri)
        else:
            tri = conn.reshape(-1, 3)
        xy = _dataitem(base, g.find("Geometry/DataItem")).reshape(-1, 2)
        fields = {}
        for a in g.findall("Attribute"):
            if a.get("Name") in names:
                fields[a.get("Name")] = _dataitem(base, a.find("DataItem")).ravel()
        parts.append((xy, tri, fields))
    xy = np.vstack([p[0] for p in parts])
    off = np.cumsum([0] + [len(p[0]) for p in parts[:-1]])
    tri = np.vstack([p[1] + o for p, o in zip(parts, off)])
    fields = {n: np.concatenate([p[2][n] for p in parts]) for n in names if all(n in p[2] for p in parts)}
    return xy, tri, fields


def plot_contours(case_dir, r, nd, hours, out):
    xdmf = os.path.join(case_dir, "NaturalConvectionPCM.xdmf")
    if not os.path.exists(xdmf):
        print("  (sin XDMF en %s)" % case_dir)
        return
    base, steps = read_xdmf(xdmf)
    times = np.array([s[0] for s in steps])
    t0 = (nd - 1) * DAY
    sel = []
    for hh in hours:
        k = int(np.argmin(np.abs(times - (t0 + hh * 3600.0))))
        sel.append((hh, times[k], steps[k][1]))
    data = [(hh, tk, load_step(base, st)) for hh, tk, st in sel]
    Tall = np.concatenate([d[2][2]["temperature"] for d in data]) - 273.15
    lo, hi = np.floor(Tall.min()), np.ceil(Tall.max())
    if hi - lo < 1:
        hi = lo + 1
    levels = np.linspace(lo, hi, 41)
    n = len(data)
    fig, ax = plt.subplots(n, 1, figsize=(9, 1.25 * n + 1.1), sharex=True, squeeze=False)
    ax = ax[:, 0]
    cf = None
    for a, (hh, tk, (xy, tri, f)) in zip(ax, data):
        y_mm = xy[:, 1] * 1e3
        tr = mtri.Triangulation(xy[:, 0], y_mm, tri)
        T = f["temperature"] - 273.15
        cf = a.tricontourf(tr, T, levels=levels, cmap="inferno")
        if "liquid_fraction" in f:
            Y = f["liquid_fraction"]
            if Y.min() < 0.5 < Y.max():
                a.tricontour(tr, Y, levels=[0.5], colors="cyan", linewidths=0.9)
        a.set_ylabel("y [mm]")
        a.set_ylim(0, y_mm.max())
        a.grid(False)
        a.set_title("%02d:%02d h solar" % (int((tk - t0) / 3600), int(round((tk - t0) % 3600 / 60))),
                    loc="left", fontsize=8, pad=2)
    ax[-1].set_xlabel("x a lo largo de la pendiente [m] (escala vertical exagerada)")
    cb = fig.colorbar(cf, ax=list(ax), fraction=0.025, pad=0.06)
    cb.set_label("T [$^\\circ$C]")
    fig.suptitle("Temperatura en el PCM, e = %g mm, ultimo dia (cian: frente de fusion Y = 0.5)"
                 % r["e_mm"], fontsize=10)
    fig.savefig(os.path.join(out, "T_contornos_e_%02dmm.png" % int(r["e_mm"])), bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", nargs="?", default="pcm_paine_runs")
    ap.add_argument("--hours", type=float, nargs="+", default=[9, 13, 17, 21],
                    help="solar hours of the last day for the contours (>= 3)")
    ap.add_argument("--cases", type=float, nargs="*", default=None,
                    help="thicknesses [mm] with contours (default: thinnest, middle, thickest; 'all' = every case)")
    ap.add_argument("--all-contours", action="store_true")
    a = ap.parse_args()

    dirs = sorted(d for d in glob.glob(os.path.join(a.out, "e_*mm"))
                  if os.path.exists(os.path.join(d, "NaturalConvectionPCM.csv")))
    if not dirs:
        sys.exit("No cases with NaturalConvectionPCM.csv in %s" % a.out)
    post = os.path.join(a.out, "post")
    os.makedirs(post, exist_ok=True)
    style()

    cases = []
    for d in dirs:
        data = read_csv(os.path.join(d, "NaturalConvectionPCM.csv"))
        if len(data["t"]) < 3:
            continue
        info = case_info(d)
        r, h, m = analyse(data, info)
        r["_dir"] = d
        if not os.path.exists(os.path.join(d, "DONE")):
            print("  aviso: %s no ha terminado (t = %.1f h)" % (d, data["t"][-1] / 3600))
        cases.append((r, data, h, m))
    cases.sort(key=lambda c: c[0]["e_mm"])
    rows = [c[0] for c in cases]

    write_summary(rows, post)
    plot_fluxes(cases, post)
    plot_vs_thickness(rows, post)
    plot_liquid(cases, post)

    es = [r["e_mm"] for r in rows]
    if a.all_contours:
        want = es
    elif a.cases:
        want = a.cases
    else:
        want = sorted(set([es[0], es[len(es) // 2], es[-1]]))
    for r in rows:
        if any(abs(r["e_mm"] - w) < 1e-6 for w in want):
            plot_contours(r["_dir"], r, r["dias"], a.hours, post)
    print("\nResultados en %s:" % post)
    for f in sorted(os.listdir(post)):
        print("  " + f)


if __name__ == "__main__":
    main()
