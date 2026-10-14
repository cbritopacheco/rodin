#!/usr/bin/env python3
"""Straight-pipe mesh of the meridian half-plane for ArterialLesionAxi, without gmsh.

    python3 make_pipe_axi_mesh.py <rf> <output.mesh> [--lu 2] [--ld 12]

Structured, quasi-uniform triangles of size h = 0.1 D / rf on
-lu <= x/D <= ld, 0 <= y/D <= 1/2 (y is the radius), in units of D like the
other meshes of the drivers. Each rectangle is split along alternating diagonals.
The boundary tags are those of make_lesion_mesh.py --axi: 1 inlet, 2 outlet,
3 wall, 4 axis, and 5 for the wall with |x| <= 1 (the "lesion" of the healthy
pipe, so that the lesion indices exist). The triangles carry tag 10.
"""
import argparse

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rf", type=float)
    ap.add_argument("out")
    ap.add_argument("--lu", type=float, default=2.0)
    ap.add_argument("--ld", type=float, default=12.0)
    a = ap.parse_args()

    h = 0.1 / a.rf
    nx = max(2, int(round((a.lu + a.ld) / h)))
    ny = max(2, int(round(0.5 / h)))
    x = np.linspace(-a.lu, a.ld, nx + 1)
    y = np.linspace(0.0, 0.5, ny + 1)
    X, Y = np.meshgrid(x, y, indexing="ij")
    V = np.column_stack([X.ravel(), Y.ravel()])
    idx = lambda i, j: i * (ny + 1) + j + 1          # 1-based, MEDIT

    tris = []
    for i in range(nx):
        for j in range(ny):
            a00, a10, a01, a11 = idx(i, j), idx(i + 1, j), idx(i, j + 1), idx(i + 1, j + 1)
            if (i + j) % 2 == 0:
                tris += [(a00, a10, a11), (a00, a11, a01)]
            else:
                tris += [(a00, a10, a01), (a10, a11, a01)]
    tris = np.array(tris)
    # Positive orientation (counter-clockwise in the (x, y) plane).
    P = V[tris - 1]
    area = 0.5 * ((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1])
                  - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
    assert np.all(area > 0)

    edges = []
    for j in range(ny):
        edges.append((idx(0, j), idx(0, j + 1), 1))
        edges.append((idx(nx, j), idx(nx, j + 1), 2))
    for i in range(nx):
        xm = 0.5 * (x[i] + x[i + 1])
        edges.append((idx(i, ny), idx(i + 1, ny), 5 if abs(xm) <= 1.0 else 3))
        edges.append((idx(i, 0), idx(i + 1, 0), 4))

    with open(a.out, "w") as f:
        f.write("MeshVersionFormatted\n2\n\nDimension\n2\n\nVertices\n%d\n" % len(V))
        for p in V:
            f.write("%.16g %.16g 0\n" % (p[0], p[1]))
        f.write("\nTriangles\n%d\n" % len(tris))
        for t in tris:
            f.write("%d %d %d 10\n" % tuple(t))
        f.write("\nEdges\n%d\n" % len(edges))
        for e in edges:
            f.write("%d %d %d\n" % e)
        f.write("\nEnd\n")
    print("axi pipe  rf=%g h/D=%.4g  vertices=%d triangles=%d  x in [%g, %g]"
          % (a.rf, h, len(V), len(tris), -a.lu, a.ld))


if __name__ == "__main__":
    main()
