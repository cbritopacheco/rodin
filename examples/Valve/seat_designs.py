#!/usr/bin/env python3
"""
Figure of the seat designs of make_meshes.py --seat (FPC8): the cross-section
of the seat at the bottom of the bore (plane x = 0, flow along z, knife
stroke along y) with the knife seated, and, if meshes are given, the cut of
the fluid meshes by the same plane around the seat at the requested opening.

  python3 seat_designs.py -o seat_designs.pdf [--meshes resources/examples/Valve/seat --opening 25]
"""
import argparse
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_meshes as mm  # noqa: E402

R, Rk, t = 101.6, 105.0, 10.0   # FPC8, mm
d = Rk - R                       # pocket depth, 3.4 mm
Lr = d / math.tan(math.radians(mm.UPRAMP_DEG))
fil = mm.FILLET_MM
cm, cdeg = mm.CHAMFER_MM, mm.CHAMFER_DEG
Lc = cm / math.tan(math.radians(cdeg))
LINER = "#c9b27c"; KNIFE = "#8a8f99"; FLUID = "#dbe9f7"


def knife_patch(ax, seat, y_tip):
    """Knife in the x = 0 plane, t wide, tip at y_tip (local coordinates:
    bore surface at y = 0, pocket below); corners rounded for "fillet"."""
    top = 40.0
    if seat == "fillet":
        pts = [(-0.5 * t, top), (-0.5 * t, y_tip + fil)]
        for th in np.linspace(math.pi, 1.5 * math.pi, 15):
            pts.append((-0.5 * t + fil + fil * math.cos(th), y_tip + fil + fil * math.sin(th)))
        for th in np.linspace(1.5 * math.pi, 2 * math.pi, 15):
            pts.append((0.5 * t - fil + fil * math.cos(th), y_tip + fil + fil * math.sin(th)))
        pts += [(0.5 * t, top)]
        ax.add_patch(Polygon(pts, closed=True, fc=KNIFE, ec="k", lw=0.8, alpha=0.9))
    else:
        ax.add_patch(Rectangle((-0.5 * t, y_tip), t, top - y_tip, fc=KNIFE, ec="k", lw=0.8, alpha=0.9))


def liner_patch(ax, seat, zmax):
    """Liner below the bore surface with the pocket of the design."""
    depth = 14.0
    if seat == "cad":
        pocket = [(-0.5 * t, 0), (-0.5 * t, -d), (0.5 * t, -d), (0.5 * t, 0)]
    elif seat == "fillet":
        pocket = [(-0.5 * t, 0), (-0.5 * t, -d + fil)]
        for th in np.linspace(math.pi, 1.5 * math.pi, 15):
            pocket.append((-0.5 * t + fil + fil * math.cos(th), -d + fil + fil * math.sin(th)))
        for th in np.linspace(1.5 * math.pi, 2 * math.pi, 15):
            pocket.append((0.5 * t - fil + fil * math.cos(th), -d + fil + fil * math.sin(th)))
        pocket += [(0.5 * t, 0)]
    elif seat == "chamfer":
        pocket = [(-0.5 * t - Lc, 0), (-0.5 * t, -cm), (-0.5 * t, -d), (0.5 * t, -d), (0.5 * t, -cm), (0.5 * t + Lc, 0)]
    elif seat == "upramp":
        pocket = [(-0.5 * t - Lr, 0), (-0.5 * t, -d), (0.5 * t, -d), (0.5 * t, 0)]
    pts = [(-zmax, 0)] + pocket + [(zmax, 0), (zmax, -depth), (-zmax, -depth)]
    ax.add_patch(Polygon(pts, closed=True, fc=LINER, ec="k", lw=0.8, hatch="//"))


def schematic(ax, seat, title, note):
    zmax = 36.0
    ax.add_patch(Rectangle((-zmax, 0), 2 * zmax, 40, fc=FLUID, ec="none"))
    liner_patch(ax, seat, zmax)
    knife_patch(ax, seat, -d)
    ax.annotate("", xy=(-14, 13), xytext=(-30, 13), arrowprops=dict(arrowstyle="->", lw=1.2))
    ax.text(-30, 14.5, "flujo", fontsize=8)
    ax.annotate("", xy=(0, 6), xytext=(0, 16), arrowprops=dict(arrowstyle="->", lw=1.2, color="w"))
    ax.text(1.5, 10, "cierre", fontsize=8)
    ax.set_xlim(-zmax, zmax); ax.set_ylim(-7, 18); ax.set_aspect("equal")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.text(0.0, -0.55, note, transform=ax.transAxes, fontsize=8, va="top")
    ax.set_xlabel("z [mm] (a lo largo del flujo)", fontsize=8); ax.set_ylabel("y [mm] (0 = pared del diámetro)", fontsize=8)
    ax.tick_params(labelsize=7)


def mesh_cut(ax, path, title):
    """Edges of the cut of the tetrahedra by x = 0, around the seat."""
    V, T = load_medit(path)
    V = V * 1000.0  # mm
    sel_y = (V[:, 1] < -R + 60) & (np.abs(V[:, 2]) < 60)
    keep = sel_y[T].all(1)
    T = T[keep]
    dx = V[T, 0]
    cut = (dx.min(1) < 0) & (dx.max(1) > 0)
    T = T[cut]
    E = np.array([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
    a, b = T[:, E[:, 0]], T[:, E[:, 1]]
    da, db = V[a, 0], V[b, 0]
    cross = (da * db) < 0
    s = np.where(cross, da / np.where(cross, da - db, 1.0), 0.0)
    P = V[a] + s[..., None] * (V[b] - V[a])
    segs = []
    for i in range(len(T)):
        pts = P[i][cross[i]]
        if len(pts) < 3:
            continue
        c = pts.mean(0)
        o = np.argsort(np.arctan2(pts[:, 1] - c[1], pts[:, 2] - c[2]))
        pts = pts[o]
        for k in range(len(pts)):
            segs.append([(pts[k, 2], pts[k, 1] + R), (pts[(k + 1) % len(pts), 2], pts[(k + 1) % len(pts), 1] + R)])
    from matplotlib.collections import LineCollection
    ax.add_collection(LineCollection(segs, colors="0.25", linewidths=0.25))
    ax.set_xlim(-36, 36); ax.set_ylim(-6, 22); ax.set_aspect("equal")
    ax.axhline(0, color="r", lw=0.6, ls="--")
    ax.set_title(title, fontsize=9); ax.tick_params(labelsize=7)
    ax.set_xlabel("z [mm]", fontsize=8); ax.set_ylabel("y [mm] (0 = pared)", fontsize=8)


def load_medit(fn):
    with open(fn) as fh:
        lines = fh.read().split("\n")
    i = 0; V = T = None
    while i < len(lines):
        s = lines[i].strip()
        if s == "Vertices":
            n = int(lines[i + 1]); V = np.loadtxt(lines[i + 2:i + 2 + n])[:, :3]; i += 2 + n; continue
        if s == "Tetrahedra":
            n = int(lines[i + 1]); T = np.loadtxt(lines[i + 2:i + 2 + n], dtype=int)[:, :4] - 1; i += 2 + n; continue
        i += 1
    return V, T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="seat_designs.pdf")
    ap.add_argument("--meshes", default="", help="folder with FPC8[-seat]_openNNN.mesh (optional)")
    ap.add_argument("--opening", type=int, default=25)
    a = ap.parse_args()
    designs = [
        ("cad", "CAD (referencia)",
         "Ranura rectangular %g × %g mm; el cuchillo\nasienta en el fondo y las paredes verticales\nlo retienen. Bolsa muerta con el cuchillo abierto." % (t, d)),
        ("fillet", "A — Fondo con radios R %g" % fil,
         "Mismo fondo y paredes verticales (%.1f mm);\nesquinas del fondo y de la punta del cuchillo\ncon radio %g mm: sin rincones estancados." % (d - fil, fil)),
        ("chamfer", "B — Chaflanes de entrada %g°" % cdeg,
         "Ranura %g × %g con chaflán de %g mm de hondo\n(%.1f mm de largo) en ambos bordes; paredes\nverticales en los %.1f mm restantes; bidireccional." % (t, d, cm, Lc, d - cm)),
        ("upramp", "C — Rampa aguas arriba %g°" % mm.UPRAMP_DEG,
         "Pared aguas abajo vertical (toma la presión),\nrampa de %.1f mm aguas arriba hasta el fondo:\nel flujo barre el fondo. Válvula unidireccional." % Lr),
    ]
    with_mesh = bool(a.meshes)
    fig, axs = plt.subplots(2 if with_mesh else 1, 4, figsize=(16, 6.2 if with_mesh else 4.4), squeeze=False,
                            gridspec_kw=dict(hspace=0.9))
    for ax, (seat, title, note) in zip(axs[0], designs):
        schematic(ax, seat, title, note)
    if with_mesh:
        for ax, (seat, title, _) in zip(axs[1], designs):
            tag = "FPC8" if seat == "cad" else "FPC8-" + seat
            fn = os.path.join(a.meshes, "%s_open%03d.mesh" % (tag, a.opening))
            if os.path.isfile(fn):
                mesh_cut(ax, fn, "malla, plano x = 0, %d %% de apertura" % a.opening)
            else:
                ax.text(0.5, 0.5, "sin malla\n" + fn, ha="center", va="center", transform=ax.transAxes, fontsize=7)
    fig.suptitle("Diseños del asiento de la FPC8 (corte por el plano x = 0 en el fondo del diámetro; cuchillo cerrado; flujo de izquierda a derecha)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(a.out, dpi=160)
    if a.out.endswith(".pdf"):
        fig.savefig(a.out[:-4] + ".png", dpi=160)
    print(a.out)


if __name__ == "__main__":
    main()
