#!/usr/bin/env python3
"""Structured triangular MEDIT meshes of a stepped axisymmetric pipe.

The pipe is a series of coaxial tubes R_k:L_k (units of the inlet radius),
joined by abrupt contractions, on the meridian half-plane y = r >= 0. Two
meshes are written:

  <prefix>_full.mesh   every tube, a node column at x_T;
  <prefix>_trunc.mesh  the tubes up to the truncation x_T.

Boundary labels (RmuAxi.cpp, Benchmark::Labels):
  1 inlet (x = 0)          2 outlet (last x)       3 wall (other walls, steps)
  4 axis (r = 0)           5 parent cylinder wall  6 stub cylinder wall,
                                                     from the first step to x_T

The radial nodes are shared by every tube (each tube uses the nodes r <= R_k),
so the contractions are conforming. No dependencies beyond the standard
library.

  python3 make_rmu_axi_mesh.py --tubes 1:10,0.6:30,0.35:30 --truncation 13 \\
      --prefix rmu_axi --nr 16 --nmin 8 --cx 0.5
"""
import argparse
import math


def parse_tubes(text):
    tubes = []
    for item in text.split(","):
        if item:
            r, l = item.split(":")
            tubes.append((float(r), float(l)))
    return tubes


def radial_nodes(radii, nr, nmin):
    """Nodes on [0, max R] containing every R_k; spacing at most 1/nr and at
    most R_min/nmin (nmin intervals across the narrowest tube)."""
    breaks = sorted(set([0.0] + radii))
    h = min(1.0 / nr, min(radii) / nmin)
    nodes = [0.0]
    for a, b in zip(breaks[:-1], breaks[1:]):
        n = max(2, int(math.ceil((b - a) / h - 1e-9)))
        nodes += [a + (b - a) * i / n for i in range(1, n + 1)]
    return nodes


def axial_nodes(tubes, truncation, cx, stop_at_truncation):
    """x nodes per tube, spacing cx R_k, with a column at the truncation."""
    columns = []  # (x, radius of the tube to the right of this column's cells)
    x0 = 0.0
    xs = [0.0]
    tube_of_interval = []
    for (R, L) in tubes:
        x1 = x0 + L
        pieces = [(x0, x1)]
        if x0 < truncation < x1:
            pieces = [(x0, truncation), (truncation, x1)]
        for a, b in pieces:
            if stop_at_truncation and a >= truncation - 1e-12:
                break
            n = max(1, int(math.ceil((b - a) / (cx * R) - 1e-9)))
            for i in range(1, n + 1):
                xs.append(a + (b - a) * i / n)
                tube_of_interval.append(R)
        if stop_at_truncation and x1 >= truncation - 1e-12:
            break
        x0 = x1
    return xs, tube_of_interval


def build(tubes, truncation, rnodes, cx, truncated):
    xs, radius_of = axial_nodes(tubes, truncation, cx, truncated)
    first_step = tubes[0][1]
    vid = {}
    verts = []

    def vertex(i, j):
        key = (i, j)
        if key not in vid:
            vid[key] = len(verts) + 1
            verts.append((xs[i], rnodes[j]))
        return vid[key]

    def nr_of(R):
        return max(j for j, r in enumerate(rnodes) if r <= R + 1e-12)

    tris = []
    for i, R in enumerate(radius_of):
        jmax = nr_of(R)
        for j in range(jmax):
            a, b = vertex(i, j), vertex(i + 1, j)
            c, d = vertex(i + 1, j + 1), vertex(i, j + 1)
            # alternate the diagonal to avoid a directional bias
            if (i + j) % 2 == 0:
                tris += [(a, b, c), (a, c, d)]
            else:
                tris += [(a, b, d), (b, c, d)]

    edges = []
    xend = len(radius_of)
    # inlet and outlet
    for j in range(nr_of(radius_of[0])):
        edges.append((vertex(0, j), vertex(0, j + 1), 1))
    for j in range(nr_of(radius_of[-1])):
        edges.append((vertex(xend, j), vertex(xend, j + 1), 2))
    # axis and cylinder walls
    for i, R in enumerate(radius_of):
        edges.append((vertex(i, 0), vertex(i + 1, 0), 4))
        jmax = nr_of(R)
        xm = 0.5 * (xs[i] + xs[i + 1])
        if xm < first_step:
            label = 5
        elif xm < truncation:
            label = 6
        else:
            label = 3
        edges.append((vertex(i, jmax), vertex(i + 1, jmax), label))
    # contraction faces
    for i in range(1, xend):
        jl, jr = nr_of(radius_of[i - 1]), nr_of(radius_of[i])
        lo, hi = sorted((jl, jr))
        for j in range(lo, hi):
            edges.append((vertex(i, j), vertex(i, j + 1), 3))
    return verts, tris, edges


def write_medit(path, verts, tris, edges):
    with open(path, "w") as f:
        f.write("MeshVersionFormatted 2\n\nDimension 2\n\n")
        f.write("Vertices\n%d\n" % len(verts))
        for x, y in verts:
            f.write("%.17g %.17g 0\n" % (x, y))
        f.write("\nTriangles\n%d\n" % len(tris))
        for a, b, c in tris:
            f.write("%d %d %d 1\n" % (a, b, c))
        f.write("\nEdges\n%d\n" % len(edges))
        for a, b, l in edges:
            f.write("%d %d %d\n" % (a, b, l))
        f.write("\nEnd\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tubes", required=True, help="R:L,... (units of the inlet radius)")
    ap.add_argument("--truncation", type=float, required=True, help="x_T, units of R_1")
    ap.add_argument("--prefix", default="rmu_axi")
    ap.add_argument("--nr", type=int, default=16, help="radial intervals per unit radius")
    ap.add_argument("--nmin", type=int, default=8, help="radial intervals in the narrowest tube")
    ap.add_argument("--cx", type=float, default=0.5, help="axial spacing / local radius")
    a = ap.parse_args()

    tubes = parse_tubes(a.tubes)
    if abs(tubes[0][0] - 1.0) > 1e-12:
        raise SystemExit("the first tube must have radius 1 (the inlet radius is the unit)")
    if not (tubes[0][1] < a.truncation < sum(L for _, L in tubes)):
        raise SystemExit("the truncation must lie downstream of the first tube")

    rnodes = radial_nodes([R for R, _ in tubes], a.nr, a.nmin)
    for truncated, name in ((False, "full"), (True, "trunc")):
        verts, tris, edges = build(tubes, a.truncation, rnodes, a.cx, truncated)
        path = "%s_%s.mesh" % (a.prefix, name)
        write_medit(path, verts, tris, edges)
        print("%s: %d vertices, %d triangles, %d boundary edges"
              % (path, len(verts), len(tris), len(edges)))


if __name__ == "__main__":
    main()
