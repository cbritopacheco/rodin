#!/usr/bin/env python3
"""Vista rapida de Q1: todas las variantes del modelo lado a lado.

    python3 quicklook.py                         # results_pc9, todo lo que haya
    python3 quicklook.py --out results_pc9 --shape CW --rhythm AF --times 5,15,30,57

Escribe en <out>/quicklook/:

  series.png      fibrina en el tiempo, una fila por forma, RS y FA en columnas,
                  una linea por variante. Tres paneles por caso:
                    - fibrina dentro de la orejuela (masa, laaF)
                    - fibrina FUERA de la orejuela (intF - laaF): es la que no
                      deberia estar; en un modelo de gel inmovil, crecer aqui
                      quiere decir que encendio el cuerpo auricular
                    - trombina maxima en la orejuela (laaMaxT), escala log
  campos_<forma>_<ritmo>.png
                  el campo de fibrina en varios instantes, una fila por variante,
                  misma escala de color en toda la figura, con el contorno de la
                  orejuela
  dosis.png       si hay barridos de dosis (dose_<variante>/): tiempo de
                  encendido contra j_0, RS y FA, una linea por variante

Lee run.csv y run.xdmf de cada corrida; no necesita ParaView. Las corridas que
siguen en curso se dibujan hasta donde van.
"""
import argparse
import glob
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

pc = px.pc
HERE = os.path.dirname(os.path.abspath(__file__))
MESH = os.path.join(HERE, "meshes")
SHAPES = ("CW", "WS", "BR", "CA")
RHYTHMS = ("SR", "AF")
RLABEL = {"SR": "RS", "AF": "FA"}
ORDER = ["qureshi", "ptagx_base", "ptagx", "ptagx_notm", "ptagx_mono",
         "ptagx_surface", "ptagx_legacy"]
BURST = 1.0e-4                                   # 100 nM de trombina = encendido


def variants_in(out, seed):
    found = {}
    for d in glob.glob(os.path.join(out, f"*_*_*_{seed}")):
        if not os.path.exists(os.path.join(d, "run.csv")):
            continue
        kin, s, rh, _ = os.path.basename(d).rsplit("_", 3)
        found.setdefault(kin, set()).add((s, rh))
    return sorted(found, key=lambda k: (ORDER.index(k) if k in ORDER else 99, k)), found


def colors(kins):
    # orden fijo por variante, nunca ciclado; qureshi en gris (es la referencia)
    base = ["#2a6fdb", "#d1495b", "#2e933c", "#e09f3e", "#7b4fa0", "#1b998b"]
    c, i = {}, 0
    for k in kins:
        if k == "qureshi":
            c[k] = "#6b6b6b"
        else:
            c[k] = base[i % len(base)]
            i += 1
    return c


def read_run(d, warmup):
    try:
        dd = pc.read_csv(os.path.join(d, "run.csv"))
    except Exception:
        return None
    if len(dd["t"]) == 0:
        return None
    t = dd["t"] - warmup
    m = t >= -1e-9
    return {k: np.asarray(v)[m] for k, v in dd.items()} | {"tau": t[m]}


def series(out, kins, found, cmap, args, plt):
    shapes = [s for s in SHAPES if any(s in {x[0] for x in found[k]} for k in kins)]
    panels = [("laaF", "fibrina en la orejuela (mol/m)", False),
              ("outF", "fibrina fuera de la orejuela (mol/m)", False),
              ("laaMaxT", "trombina max. orejuela (mol/m$^3$)", True)]
    fig, axs = plt.subplots(len(shapes), 2 * len(panels),
                            figsize=(3.0 * 2 * len(panels), 2.3 * len(shapes)),
                            squeeze=False, sharex=True)
    for i, s in enumerate(shapes):
        for j, rh in enumerate(RHYTHMS):
            for p, (key, lab, logy) in enumerate(panels):
                ax = axs[i, 2 * p + j]
                for k in kins:
                    r = read_run(os.path.join(out, f"{k}_{s}_{rh}_{args.seed}"), args.warmup)
                    if r is None:
                        continue
                    y = r["intF"] - r["laaF"] if key == "outF" else r[key]
                    if logy:
                        y = np.maximum(y, 1e-9)
                    ax.plot(r["tau"], y, color=cmap[k], lw=1.6,
                            ls="--" if k == "qureshi" else "-", label=k)
                if logy:
                    ax.set_yscale("log")
                    ax.axhline(BURST, color="#999", lw=0.8, ls=":")
                if i == 0:
                    ax.set_title(f"{lab}\n{RLABEL[rh]}", fontsize=8)
                if 2 * p + j == 0:
                    ax.set_ylabel(s, fontsize=11, fontweight="bold")
                if i == len(shapes) - 1:
                    ax.set_xlabel("t desde la semilla (s)")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=len(kins), bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(args.dest, "series.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def xdmf_steps(path, field):
    """[(t, h5 file, dataset)] of one field, and the mesh h5 of the first step."""
    root = ET.parse(path).getroot()
    steps, meshref = [], None
    for g in root.iter("Grid"):
        tnode = g.find("Time")
        if tnode is None:
            continue
        t = float(tnode.get("Value"))
        for a in g.iter("Attribute"):
            if a.get("Name") == field:
                f, ds = a.find("DataItem").text.strip().split(":")
                steps.append((t, f, ds))
        if meshref is None:
            topo = g.find(".//Topology/DataItem").text.strip().split(":")
            geo = g.find(".//Geometry/DataItem").text.strip().split(":")
            meshref = (topo, geo)
    return steps, meshref


def load_field(d, field, times, warmup):
    import h5py
    x = os.path.join(d, "run.xdmf")
    if not os.path.exists(x):
        return None
    try:
        steps, (topo, geo) = xdmf_steps(x, field)
    except Exception:
        return None
    if not steps:
        return None
    with h5py.File(os.path.join(d, topo[0]), "r") as f:
        tri = np.asarray(f[topo[1]]).reshape(-1, 4)[:, 1:]
    with h5py.File(os.path.join(d, geo[0]), "r") as f:
        xy = np.asarray(f[geo[1]])
    st = np.array([s[0] for s in steps]) - warmup
    out = []
    for tq in times:
        i = int(np.argmin(np.abs(st - tq)))
        if abs(st[i] - tq) > 1.01:              # esa corrida todavia no llega
            out.append((tq, None))
            continue
        t, f, ds = steps[i]
        try:
            with h5py.File(os.path.join(d, f), "r") as h:
                out.append((st[i], np.asarray(h[ds])))
        except Exception:
            out.append((tq, None))
    return xy, tri, out


def fields(out, kins, args, plt, shape, rh):
    import matplotlib.tri as mtri
    times = [float(x) for x in args.times.split(",")]
    rows = []
    for k in kins:
        r = load_field(os.path.join(out, f"{k}_{shape}_{rh}_{args.seed}"), args.field,
                       times, args.warmup)
        if r is not None:
            rows.append((k, r))
    if not rows:
        return None
    vmax = max((np.max(v) for _, (_, _, o) in rows for _, v in o if v is not None),
               default=1.0) or 1.0
    poly = os.path.join(MESH, f"Q_{shape}_laa.txt")
    laa = np.loadtxt(poly) if os.path.exists(poly) else None
    xy0 = rows[0][1][0]
    span = np.ptp(xy0, axis=0)
    w = 2.2 * span[0] / span[1] if span[1] > 0 else 2.2
    fig, axs = plt.subplots(len(rows), len(times), figsize=(w * len(times), 2.2 * len(rows) + 0.6),
                            squeeze=False)
    im = None
    for i, (k, (xy, tri, snaps)) in enumerate(rows):
        T = mtri.Triangulation(xy[:, 0], xy[:, 1], tri)
        for j, (t, v) in enumerate(snaps):
            ax = axs[i, j]
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
            if v is None:
                ax.text(0.5, 0.5, "aun no llega", ha="center", va="center",
                        transform=ax.transAxes, color="#888")
            else:
                im = ax.tripcolor(T, v, shading="gouraud", cmap="viridis", vmin=0, vmax=vmax)
                if laa is not None:
                    ax.plot(laa[:, 0], laa[:, 1], color="white", lw=0.6, alpha=0.8)
            if i == 0:
                ax.set_title(f"t = {t:.0f} s", fontsize=9)
            if j == 0:
                ax.set_ylabel(k, fontsize=9)
                ax.yaxis.set_label_coords(-0.02, 0.5)
    if im is not None:
        cb = fig.colorbar(im, ax=axs, shrink=0.8, pad=0.01)
        cb.set_label(f"{args.field} (mol/m$^3$)")
    fig.suptitle(f"{args.field}: forma {shape}, {RLABEL[rh]} (misma escala en toda la figura)",
                 fontsize=10)
    p = os.path.join(args.dest, f"campos_{shape}_{rh}.png")
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def dose(out, cmap, args, plt):
    dirs = sorted(d for d in glob.glob(os.path.join(out, "dose_*")) if os.path.isdir(d)
                  and not d.endswith("_seeded"))
    if not dirs:
        return None
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
    for d in dirs:
        kin = os.path.basename(d)[len("dose_"):]
        for j, rh in enumerate(RHYTHMS):
            for s in SHAPES:
                pts = []
                for c in glob.glob(os.path.join(d, f"{s}_{rh}_d*", "run.csv")):
                    dj = float(os.path.basename(os.path.dirname(c)).split("_d")[-1])
                    r = read_run(os.path.dirname(c), args.warmup)
                    if r is None:
                        continue
                    ti = pc.first_time(r["tau"], r["maxT"] >= BURST)
                    pts.append((dj, ti if np.isfinite(ti) else np.nan, r["tau"][-1]))
                if not pts:
                    continue
                pts.sort()
                dj, ti, _ = map(np.array, zip(*pts))
                axs[j].plot(dj, ti, "o-", color=cmap.get(kin, "#333"), ms=4, lw=1.2,
                            alpha=0.85, label=f"{kin} {s}")
        axs[j].set_xscale("log")
        axs[j].set_title(f"{RLABEL[rh]}: tiempo de encendido (sin punto = no enciende)", fontsize=9)
        axs[j].set_xlabel("j$_0$ / j$_0$ por defecto")
    axs[0].set_ylabel("t desde el inicio (s)")
    axs[1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    p = os.path.join(args.dest, "dosis.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(HERE, "results_pc9"))
    ap.add_argument("--warmup", type=float, default=3)
    ap.add_argument("--seed", default="P2")
    ap.add_argument("--shape", default="all", help="forma para los campos (o all)")
    ap.add_argument("--rhythm", default="all", help="SR, AF o all")
    ap.add_argument("--times", default="2,10,30,57", help="instantes de los campos (s)")
    ap.add_argument("--field", default="fibrin")
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    args.dest = os.path.join(out, "quicklook")
    os.makedirs(args.dest, exist_ok=True)
    plt = pc.setup_matplotlib()

    kins, found = variants_in(out, args.seed)
    cmap = colors(kins + sorted({os.path.basename(d)[5:] for d in glob.glob(os.path.join(out, "dose_*"))}
                                - set(kins)))
    made = []
    if kins:
        print("variantes:", ", ".join(kins), flush=True)
        made.append(series(out, kins, found, cmap, args, plt))
        shapes = SHAPES if args.shape == "all" else args.shape.split(",")
        rhs = RHYTHMS if args.rhythm == "all" else args.rhythm.split(",")
        for s in shapes:
            for rh in rhs:
                p = fields(out, kins, args, plt, s, rh)
                if p:
                    made.append(p)
    p = dose(out, cmap, args, plt)
    if p:
        made.append(p)
    if not made:
        sys.exit(f"no encontre corridas en {out}")
    print("\n".join(f"  {m}" for m in made))


if __name__ == "__main__":
    main()
