#!/usr/bin/env python3
"""
Dimensioned drawing of the seat profiles "fillet" and "lipr" of
make_meshes.py for the FPC8 and the FPC16: section by the plane x = 0 at the
bottom of the bore (flow along z, left to right), liner in tan, fluid
pocket in white, knife seated (closed) dashed. Dimensions in mm, taken from
make_meshes.py (the FPC16 ones scaled with its groove depth).

  python3 seat_plan.py -o seat_plan.pdf
"""
import argparse
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_meshes as mm  # noqa: E402

LINER = "#d8c79a"; KNIFE = "#5b6270"; DIM = "#1f4e9a"


def rounded(pts, radii, n=16):
    """Polyline with corners rounded (r > 0), bevelled (r < 0) or sharp."""
    out = [pts[0]]
    for i in range(1, len(pts) - 1):
        P, A, B = (np.asarray(x, float) for x in (pts[i], pts[i - 1], pts[i + 1]))
        r = radii[i]
        if r == 0:
            out.append(tuple(P)); continue
        u = (A - P) / np.linalg.norm(A - P); v = (B - P) / np.linalg.norm(B - P)
        if r < 0:
            out += [tuple(P - r * u), tuple(P - r * v)]; continue
        a_, b_, c_ = mm._fillet(P, A, B, r)
        t0 = math.atan2(a_[1] - c_[1], a_[0] - c_[0]); t1 = math.atan2(b_[1] - c_[1], b_[0] - c_[0])
        dt = (t1 - t0 + math.pi) % (2 * math.pi) - math.pi
        out += [(c_[0] + r * math.cos(t0 + k * dt / n), c_[1] + r * math.sin(t0 + k * dt / n)) for k in range(n + 1)]
    out.append(pts[-1])
    return out


def profile(seat, t, d):
    """Pocket outline (z, y): bore surface y = 0, pocket below (y < 0)."""
    sc = d / mm.SEAT_REF_DEPTH
    h, W = 0.5 * t, 1.5 * t
    if seat == "fillet":
        f = mm.FILLET_MM * sc
        pts = [(-W, 0), (-h, 0), (-h, -d), (h, -d), (h, 0), (W, 0)]
        return rounded(pts, [0, 0, f, f, 0, 0]), dict(f=f)
    land = mm.LIP_LAND * sc
    L = (d - land) / math.tan(math.radians(mm.LIP_DEG))
    pts = [(-W, 0), (-h - L, 0), (-h, -(d - land)), (-h, -d), (h, -d), (h, -(d - land)), (h + L, 0), (W, 0)]
    if seat == "lipb":
        b, l_, fl = mm.LIPB_BORE * sc, mm.LIPB_LAND * sc, mm.LIPB_FLOOR * sc
        return rounded(pts, [0, -b, -l_, -fl, -fl, -l_, -b, 0]), dict(land=land, L=L, b=b, l=l_, fl=fl, ramp=d - land)
    b, l_, fl = mm.LIP_RB * sc, mm.lipr_land_radius(sc), mm.LIPR_FLOOR * sc
    return rounded(pts, [0, b, l_, fl, fl, l_, b, 0]), dict(land=land, L=L, b=b, l=l_, fl=fl, ramp=d - land)


def dim_h(ax, x0, x1, y, txt, off=0.0):
    ax.annotate("", (x0, y), (x1, y), arrowprops=dict(arrowstyle="<->", color=DIM, lw=0.7, shrinkA=0, shrinkB=0))
    ax.text(0.5 * (x0 + x1), y + off, txt, color=DIM, ha="center", va="bottom", fontsize=7)


def dim_v(ax, x, y0, y1, txt, side=1):
    ax.annotate("", (x, y0), (x, y1), arrowprops=dict(arrowstyle="<->", color=DIM, lw=0.7, shrinkA=0, shrinkB=0))
    ax.text(x + 0.15 * side, 0.5 * (y0 + y1), txt, color=DIM, ha="left" if side > 0 else "right",
            va="center", fontsize=7)


def draw(ax, valve, seat):
    v = mm.VALVES[valve]
    t, d = v["t"], v["Rk"] - v["R"]
    out, p = profile(seat, t, d)
    W = 1.5 * t
    z = [q[0] for q in out]; y = [q[1] for q in out]
    ax.add_patch(Polygon(list(zip(z, y)) + [(W, -d - 0.3 * t), (-W, -d - 0.3 * t)], closed=True,
                         fc=LINER, ec="k", lw=0.9))
    h = 0.5 * t
    # knife seated (closed): flat tip on the floor; corners as the floor corners
    kc = p.get("f", p.get("fl", 0.0))
    if seat in ("fillet", "lipr"):
        kn = rounded([(-h, 0.45 * t), (-h, -d), (h, -d), (h, 0.45 * t)], [0, kc, kc, 0])
    else:
        kn = rounded([(-h, 0.45 * t), (-h, -d), (h, -d), (h, 0.45 * t)], [0, -kc, -kc, 0])
    ax.plot([q[0] for q in kn], [q[1] for q in kn], ls="--", color=KNIFE, lw=0.9)
    ax.text(0, 0.25 * t, "cuchillo\n(cerrado)", ha="center", va="center", fontsize=7, color=KNIFE)
    ax.annotate("", (0.95 * W, 0.42 * t), (0.6 * W, 0.42 * t), arrowprops=dict(arrowstyle="->", color="0.4"))
    ax.text(0.78 * W, 0.44 * t, "flujo", fontsize=7, color="0.4", ha="center", va="bottom")
    ax.text(-0.98 * W, 0.03 * t, "diámetro interior (Ø %.1f)" % (2 * v["R"]), fontsize=7, color="0.3", va="bottom")
    # dimensions
    dim_h(ax, -h, h, -d - 0.12 * t, "%.1f (= espesor cuchillo)" % t)
    dim_v(ax, -h - 0.15 * t if seat == "fillet" else -h - p["L"] - 0.12 * t, -d, 0, "%.1f" % d, side=-1)
    if seat == "fillet":
        ax.text(h + 0.12 * t, -d + 0.05 * t, "R %.2g (×2)" % p["f"], fontsize=7, color=DIM, va="center")
        ttl = "fillet: ranura de la CAD, esquinas del fondo R %.2g" % p["f"]
    else:
        dim_v(ax, h + 0.08 * t, -d, -(d - p["land"]), "labio %.2g" % p["land"])
        dim_v(ax, h + p["L"] + 0.1 * t, -(d - p["land"]), 0, "rampa %.2g" % p["ramp"])
        dim_h(ax, -h - p["L"], h + p["L"], 0.06 * t, "%.1f en el diámetro" % (t + 2 * p["L"]))
        ax.text(-h - 0.55 * p["L"], -0.62 * p["ramp"], "%g°" % mm.LIP_DEG, fontsize=7, color=DIM, ha="right")
        if seat == "lipb":
            ax.text(0.98 * W, -d - 0.02 * t,
                    "biseles (catetos):\n  fondo %.2g × %.2g\n  labio-rampa %.2g × %.2g\n  rampa-diámetro %.2g × %.2g"
                    % (p["fl"], p["fl"], p["l"], p["l"], p["b"], p["b"]),
                    fontsize=7, color=DIM, ha="right", va="top")
            ttl = "lipb: labios %.2g + rampas %g° a ambos lados, cambios de dirección biselados" % (p["land"], mm.LIP_DEG)
        else:
            ax.text(0.98 * W, -d - 0.02 * t,
                    "radios:\n  fondo R %.2g\n  labio-rampa R %.2g\n  rampa-diámetro R %.2g"
                    % (p["fl"], p["l"], p["b"]), fontsize=7, color=DIM, ha="right", va="top")
            ttl = "lipr: lip con fondo R %.2g y cambio de dirección R %.2g" % (p["fl"], p["l"])
    ax.set_title("%s — %s" % (valve, ttl), fontsize=9)
    ax.set_xlim(-W, W); ax.set_ylim(-d - 0.32 * t, 0.5 * t)
    ax.set_aspect("equal"); ax.axis("off")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default="seat_plan.pdf")
    a = ap.parse_args()
    fig, axs = plt.subplots(2, 2, figsize=(11.7, 8.3))
    for i, valve in enumerate(("FPC8", "FPC16")):
        for j, seat in enumerate(("fillet", "lipr")):
            draw(axs[j, i], valve, seat)
    fig.suptitle("Asiento: sección por el plano x = 0 en el fondo del diámetro (cotas en mm; FPC16 escalado con la "
                 "profundidad de su ranura, %.1f/%.1f)" % (mm.VALVES["FPC16"]["Rk"] - mm.VALVES["FPC16"]["R"],
                                                            mm.SEAT_REF_DEPTH), fontsize=10)
    fig.tight_layout()
    fig.savefig(a.out)
    fig.savefig(os.path.splitext(a.out)[0] + ".png", dpi=170)
    print(a.out)


if __name__ == "__main__":
    main()
