#!/usr/bin/env python3
"""Meshes of the 2D left atrium of Qureshi et al. (CinC 2021), four appendages.

Sources (geometry/):
  ws_outline.txt, br_outline.txt   boundary of the authors' STL meshes
                                   (github.com/ahmedqureshi1/2D_LAA): exact
  cw_outline.txt                   Fig. 1 of the paper (main panel), digitised
                                   and registered on the common body
                                   (rms 0.42 mm)
  ca_outline.txt                   Fig. 1 (CA thumbnail), digitised; its scale
                                   is not given by the figure and is set so
                                   that the appendage area equals the mean of
                                   the other three
  patches.json                     arclength windows of PV1-4, MV and the LAA
                                   patch on the windsock outline (Fig. 1)
  seeds.json                       injury positions P1-P3 per shape (Fig. 1)

The atrial body is common to all shapes (it is the same in both STL files to
0.15 mm). Labels follow LeftAtrium2DPTAGX: PV1-4 = 6, 8, 10, 12, MV = 14,
LAA wall = 2, rest of the wall = 1.

Output (meshes/): Q_<shape>.mesh (MEDIT), Q_<shape>_laa.txt (appendage
polygon for -ptag_laa_poly), Q_<shape>.json (seeds, areas), and a figure.
"""
import argparse
import json
import os
import sys
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.join(HERE, "geometry")
LABEL = {"PV1": 6, "PV2": 8, "PV3": 10, "PV4": 12, "MV": 14, "LAA": 2, "wall": 1}


def load(name):
    return np.loadtxt(os.path.join(GEO, name))


def ccw(L):
    a = 0.5 * np.sum(L[:, 0] * np.roll(L[:, 1], -1) - np.roll(L[:, 0], -1) * L[:, 1])
    return L if a > 0 else L[::-1]


def densify(L, h, closed=False):
    P = np.vstack([L, L[:1]]) if closed else L
    out = []
    for a, b in zip(P[:-1], P[1:]):
        n = max(1, int(np.ceil(np.linalg.norm(b - a) / h)))
        out.append(a + (b - a) * np.linspace(0, 1, n, endpoint=False)[:, None])
    if not closed:
        out.append(P[-1:])
    return np.vstack(out)


def arclength(L):
    P = np.vstack([L, L[:1]])
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]


def resample(L, h):
    """Uniform resampling of an open polyline at spacing ~h (keeps ends)."""
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(L, axis=0), axis=1))]
    n = max(2, int(np.ceil(s[-1] / h)) + 1)
    si = np.linspace(0, s[-1], n)
    return np.c_[np.interp(si, s, L[:, 0]), np.interp(si, s, L[:, 1])]


def splice(body_ws, appendage, tol=5e-4):
    """Replace the windsock appendage by another one.

    body_ws: closed ccw windsock outline. appendage: points of the other
    outline already in body coordinates. Returns (polygon, mask_appendage).
    """
    ws_body_mask = np.load(os.path.join(GEO, "ws_body_mask.npy"))
    n = len(body_ws)
    # windsock appendage = the run of non-body vertices (contiguous, ccw)
    idx = np.flatnonzero(~ws_body_mask)
    i0, i1 = idx[0], idx[-1]
    assert np.all(np.diff(idx) == 1), "windsock appendage run is not contiguous"
    jR, jL = i0 - 1, (i1 + 1) % n          # body vertices next to the appendage
    # keep the far-from-body part of the other contour, oriented from jR to jL
    d = cKDTree(densify(body_ws[ws_body_mask], 2e-5)).query(appendage)[0]
    far = d > tol
    runs, cur = [], []
    for k, f in enumerate(far):
        if f:
            cur.append(k)
        elif cur:
            runs.append(cur); cur = []
    if cur:
        runs.append(cur)
    # a closed contour may start inside the appendage: join the wrapped run
    closed = np.linalg.norm(appendage[0] - appendage[-1]) < 2e-3
    if closed and len(runs) > 1 and runs[0][0] == 0 and runs[-1][-1] == len(appendage) - 1:
        runs = [runs[-1] + runs[0]] + runs[1:-1]
    run = max(runs, key=len)
    A = appendage[run]
    if np.linalg.norm(A[0] - body_ws[jR]) > np.linalg.norm(A[-1] - body_ws[jR]):
        A = A[::-1]
    # body vertices nearest to the two ends of the appendage
    body_idx = np.flatnonzero(ws_body_mask)
    tb = cKDTree(body_ws[body_idx])
    kR = body_idx[tb.query(A[0])[1]]
    kL = body_idx[tb.query(A[-1])[1]]
    # walk the body ccw from kL to kR, then the appendage from kR to kL
    body_path = [kL]
    k = kL
    while k != kR:
        k = (k + 1) % n
        body_path.append(k)
    P = np.vstack([body_ws[body_path], A])
    is_app = np.r_[np.zeros(len(body_path), bool), np.ones(len(A), bool)]
    return P, is_app


def polygon_area(P):
    return 0.5 * abs(np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]))


def clean_digitised(A, h=0.25e-3, window=0.6e-3):
    """Pixel contour -> smooth polyline: drop duplicates, resample at h/4,
    moving average over `window`, resample at h."""
    keep = np.r_[True, np.linalg.norm(np.diff(A, axis=0), axis=1) > 1e-9]
    A = resample(A[keep], h / 4)
    k = max(1, int(round(window / (h / 4))) | 1)
    pad = k // 2
    Ap = np.vstack([np.repeat(A[:1], pad, 0), A, np.repeat(A[-1:], pad, 0)])
    ker = np.ones(k) / k
    A = np.c_[np.convolve(Ap[:, 0], ker, "valid"), np.convolve(Ap[:, 1], ker, "valid")]
    return resample(A, h)


def ca_points(px, ca_scale=None):
    """CA thumbnail pixels (y up) -> metres, rescaled about the anchor."""
    q = json.load(open(os.path.join(GEO, "ca_transform.json")))
    s = ca_scale if ca_scale is not None else q["scale"]
    (ax, ay), (bx, by) = q["anchor_px"], q["anchor_phys"]
    px = np.atleast_2d(px)
    return np.c_[bx + (px[:, 0] - ax) * s, by + (px[:, 1] - ay) * s]


def build_polygon(shape, ca_scale=None):
    ws = ccw(load("ws_outline.txt"))
    if shape == "WS":
        P = ws.copy()
        return P, ~np.load(os.path.join(GEO, "ws_body_mask.npy"))
    if shape == "BR":
        A = ccw(load("br_outline.txt"))
    elif shape == "CW":
        A = clean_digitised(load("cw_outline.txt"))
    elif shape == "CA":
        A = clean_digitised(ca_points(load("ca_outline_px.txt"), ca_scale))
    else:
        raise ValueError(shape)
    return splice(ws, A)


def label_edges(P, is_app):
    """Label of edge (P[i], P[i+1]) from the windsock arclength windows."""
    patches = json.load(open(os.path.join(GEO, "patches.json")))
    ws = ccw(load("ws_outline.txt"))
    S = arclength(ws)
    per = S[-1]
    tree = cKDTree(ws)
    n = len(P)
    lab = np.full(n, LABEL["wall"])
    mids = 0.5 * (P + np.roll(P, -1, axis=0))
    app_edge = is_app | np.roll(is_app, -1)

    def s_of(x):
        # arclength of the projection on the windsock polyline
        d, i = tree.query(x)
        best = None
        for j in (i - 1, i):
            a, b = ws[j % len(ws)], ws[(j + 1) % len(ws)]
            t = np.clip(np.dot(x - a, b - a) / np.dot(b - a, b - a), 0, 1)
            dd = np.linalg.norm(a + t * (b - a) - x)
            if best is None or dd < best[0]:
                best = (dd, (S[j % len(ws)] + t * np.linalg.norm(b - a)) % per)
        return best[1]

    def inside(s, w):
        a, b = w
        return (a <= s <= b) if a <= b else (s >= a or s <= b)

    for k in range(n):
        if app_edge[k]:
            lab[k] = LABEL["LAA"]
            continue
        s = s_of(mids[k])
        for name in ("PV1", "PV2", "PV3", "PV4", "MV", "LAA"):
            if inside(s, patches[name]):
                lab[k] = LABEL[name]
                break
    return lab


def laa_polygon(P, lab):
    """Appendage region: the LAA-labelled run closed by its chord."""
    on = lab == LABEL["LAA"]
    n = len(P)
    start = next(k for k in range(n) if on[k] and not on[k - 1])
    run = []
    k = start
    while on[k % n]:
        run.append(k % n)
        k += 1
    run.append(k % n)
    return P[run]


def mesh(P, lab, out, h_far, h_laa, laa_poly):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("q")
    geo = gmsh.model.geo
    gp = [geo.addPoint(x, y, 0) for x, y in P]
    lines, bytag = [], {}
    for k in range(len(P)):
        L = geo.addLine(gp[k], gp[(k + 1) % len(P)])
        lines.append(L)
        bytag.setdefault(int(lab[k]), []).append(L)
    surf = geo.addPlaneSurface([geo.addCurveLoop(lines)])
    geo.synchronize()
    for t, ls in bytag.items():
        gmsh.model.addPhysicalGroup(1, ls, t)
    gmsh.model.addPhysicalGroup(2, [surf], 1)
    # size: h_laa inside the appendage and near its walls, h_far elsewhere
    f = gmsh.model.mesh.field
    dist = f.add("Distance")
    f.setNumbers(dist, "CurvesList", bytag[LABEL["LAA"]])
    f.setNumber(dist, "Sampling", 200)
    th = f.add("Threshold")
    f.setNumber(th, "InField", dist)
    f.setNumber(th, "SizeMin", h_laa)
    f.setNumber(th, "SizeMax", h_far)
    f.setNumber(th, "DistMin", 4e-3)
    f.setNumber(th, "DistMax", 20e-3)
    f.setAsBackgroundMesh(th)
    for o in ("MeshSizeFromPoints", "MeshSizeFromCurvature", "MeshSizeExtendFromBoundary"):
        gmsh.option.setNumber("Mesh." + o, 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    ntag, coord, _ = gmsh.model.mesh.getNodes()
    idx = {int(t): i + 1 for i, t in enumerate(ntag)}
    xy = coord.reshape(-1, 3)
    _, tn = gmsh.model.mesh.getElementsByType(2)
    tri = tn.reshape(-1, 3)
    edges = []
    for _, tag in gmsh.model.getPhysicalGroups(1):
        for ent in gmsh.model.getEntitiesForPhysicalGroup(1, tag):
            et, _, en = gmsh.model.mesh.getElements(1, ent)
            for ett, enn in zip(et, en):
                if ett == 1:
                    for e in enn.reshape(-1, 2):
                        edges.append((idx[int(e[0])], idx[int(e[1])], tag))
    gmsh.finalize()
    with open(out, "w") as fo:
        fo.write("MeshVersionFormatted\n2\n\nDimension\n2\n\n")
        fo.write(f"Vertices\n{len(xy)}\n")
        for p in xy:
            fo.write(f"{p[0]:.17g} {p[1]:.17g} 0\n")
        fo.write(f"\nTriangles\n{len(tri)}\n")
        for t in tri:
            A, B, C = (idx[int(v)] for v in t)
            xa, ya = xy[A - 1][:2]; xb, yb = xy[B - 1][:2]; xc, yc = xy[C - 1][:2]
            if (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya) < 0:
                B, C = C, B
            fo.write(f"{A} {B} {C} 1\n")
        fo.write(f"\nEdges\n{len(edges)}\n")
        for A, B, t in edges:
            fo.write(f"{A} {B} {t}\n")
        fo.write("\nEnd\n")
    return len(xy), len(tri), Counter(t for *_, t in edges)


def point_in_poly(x, y, poly):
    inside = False
    n = len(poly)
    for i in range(n):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            inside = not inside
    return inside


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h-far", type=float, default=2.5e-3)
    ap.add_argument("--h-laa", type=float, default=0.5e-3)
    ap.add_argument("--out", default=os.path.join(HERE, "meshes"))
    ap.add_argument("--shapes", default="CW,WS,BR,CA")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from shapely.geometry import Polygon

    seeds = json.load(open(os.path.join(GEO, "seeds.json")))
    # CA scale from the mean appendage area of the other three
    areas = {}
    for sh in ("CW", "WS", "BR"):
        P, app = build_polygon(sh)
        lab = label_edges(P, app)
        areas[sh] = polygon_area(laa_polygon(P, lab))
    target = np.mean(list(areas.values()))
    q = json.load(open(os.path.join(GEO, "ca_transform.json")))
    lo, hi = 0.5 * q["scale"], 1.6 * q["scale"]
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        P, app = build_polygon("CA", mid)
        A = polygon_area(laa_polygon(P, label_edges(P, app)))
        lo, hi = (mid, hi) if A < target else (lo, mid)
    ca_scale = 0.5 * (lo + hi)
    print(f"appendage areas (cm^2): " + ", ".join(f"{k} {v*1e4:.2f}" for k, v in areas.items())
          + f"; CA scale set to {ca_scale/q['scale']:.3f} x its body fit")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    shapes = a.shapes.split(",")
    fig, axs = plt.subplots(1, len(shapes), figsize=(4 * len(shapes), 4.2))
    axs = np.atleast_1d(axs)
    for ax, sh in zip(axs, shapes):
        P, app = build_polygon(sh, ca_scale if sh == "CA" else None)
        # resample: fine along the appendage, keep body vertices
        lab = label_edges(P, app)
        # densify each labelled run for meshing
        Q, L = [], []
        for k in range(len(P)):
            a_, b_ = P[k], P[(k + 1) % len(P)]
            h = a.h_laa if lab[k] == LABEL["LAA"] else a.h_far
            m = max(1, int(np.ceil(np.linalg.norm(b_ - a_) / (0.8 * h))))
            for t in np.linspace(0, 1, m, endpoint=False):
                Q.append(a_ + t * (b_ - a_)); L.append(lab[k])
        Q, L = np.array(Q), np.array(L)
        poly = Polygon(Q)
        assert poly.is_valid, f"{sh}: self-intersecting outline"
        laa = laa_polygon(Q, L)
        area = polygon_area(laa)
        np.savetxt(os.path.join(a.out, f"Q_{sh}_laa.txt"), laa, fmt="%.9g")
        nv, nt, cnt = mesh(Q, L, os.path.join(a.out, f"Q_{sh}.mesh"), a.h_far, a.h_laa, laa)
        sd = dict(seeds[sh])
        if sd.pop("px", False):
            sd = {p: ca_points(np.array(v), ca_scale)[0].tolist() for p, v in sd.items()}
        inside = {p: point_in_poly(*sd[p], Q) for p in sd}
        in_laa = {p: point_in_poly(*sd[p], laa) for p in sd}
        info = dict(shape=sh, vertices=nv, triangles=nt, laa_area_m2=area,
                    edge_labels={int(k): int(v) for k, v in cnt.items()},
                    seeds=sd, seeds_inside_domain=inside, seeds_inside_laa=in_laa,
                    ca_scale_factor=(ca_scale / q["scale"]) if sh == "CA" else None)
        json.dump(info, open(os.path.join(a.out, f"Q_{sh}.json"), "w"), indent=1)
        print(f"{sh}: {nv} vertices, {nt} triangles, LAA area {area*1e4:.2f} cm^2, "
              f"labels {dict(sorted(cnt.items()))}, seeds inside LAA {in_laa}")
        colors = {6: "C0", 8: "C0", 10: "C0", 12: "C0", 14: "C3", 2: "C2", 1: "0.4"}
        for k in range(len(Q)):
            b_ = Q[(k + 1) % len(Q)]
            ax.plot([Q[k, 0] * 1e3, b_[0] * 1e3], [Q[k, 1] * 1e3, b_[1] * 1e3],
                    color=colors[int(L[k])], lw=1.6 if L[k] != 1 else 0.8)
        for p, (x, y) in sd.items():
            ax.plot(x * 1e3, y * 1e3, "o", color="C4", ms=5)
            ax.annotate(p, (x * 1e3, y * 1e3), xytext=(3, 3), textcoords="offset points", fontsize=7)
        ax.set_aspect("equal")
        ax.set_title(f"{sh}: {nt} tri, LAA {area*1e4:.2f} cm$^2$", fontsize=9)
        ax.set_xlabel("x (mm)")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "Q_meshes.png"), dpi=130)
    return 0


if __name__ == "__main__":
    sys.exit(main())
