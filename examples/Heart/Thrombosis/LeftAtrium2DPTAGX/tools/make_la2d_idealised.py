#!/usr/bin/env python3
"""Idealised 2D left atrium + rectangular appendage as a planar MEDIT mesh,
with the boundary labels LeftAtrium2DPTAG expects:
  PV inlets 6,8,10,12 ; MV outlet 14 ; LAA walls 2,3,4 ; body wall 1,5,7,9,11,13.
Boundary built as an explicit tagged polyline; interior triangulated by gmsh.
Usage: python3 make_la2d_idealised.py out.mesh [--h 1.6e-3]"""
import argparse
import math
import sys
from collections import Counter

import gmsh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--h", type=float, default=1.6e-3)
    ap.add_argument("--Rb", type=float, default=0.028)
    ap.add_argument("--Wa", type=float, default=0.015)
    ap.add_argument("--La", type=float, default=0.024)
    a = ap.parse_args()
    Rb, Wa, La, h = a.Rb, a.Wa, a.La, a.h

    # angular windows (deg): 4 PV inlets (top), mitral outlet (bottom), appendage mouth (left)
    PV = [(55, 72), (80, 97), (103, 120), (128, 145)]
    MV = (248, 290)
    APP = (172, 188)                       # appendage mouth on the circle
    pv_tags = [6, 8, 10, 12]
    wall_cycle = [1, 5, 7, 9, 11, 13]

    def tag_of(deg):
        d = deg % 360
        for (lo, hi), t in zip(PV, pv_tags):
            if lo <= d <= hi:
                return t
        if MV[0] <= d <= MV[1]:
            return 14
        return None                        # wall, assigned later

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("la2d")
    geo = gmsh.model.geo

    pts, seg = [], []                      # boundary points (x,y) and segments (i,j,tag)
    step = 1.0

    def addp(x, y):
        pts.append((x, y))
        return len(pts) - 1

    # walk the circle from 0 deg; detour into the appendage box at APP
    first = addp(Rb, 0.0)
    prev, prev_deg = first, 0.0
    d = step
    appended = False
    while d < 360.0 - 1e-9:
        if APP[0] <= d <= APP[1] and not appended:
            # mouth top at APP[0], go out along the box, back at APP[1]
            a0 = math.radians(APP[0]); a1 = math.radians(APP[1])
            m0 = (Rb * math.cos(a0), Rb * math.sin(a0))
            m1 = (Rb * math.cos(a1), Rb * math.sin(a1))
            xt = min(m0[0], m1[0]) - La
            p_m0 = addp(*m0)
            p_c0 = addp(xt, m0[1])                 # corner next to m0 (upper)
            p_c1 = addp(xt, m1[1])                 # corner next to m1 (lower)
            p_m1 = addp(*m1)
            seg.append((prev, p_m0, None))        # wall up to mouth
            seg.append((p_m0, p_c0, 2))           # upper flank
            seg.append((p_c0, p_c1, 3))           # tip
            seg.append((p_c1, p_m1, 4))           # lower flank
            prev = p_m1
            prev_deg = APP[1]
            d = APP[1] + step
            appended = True
            continue
        x, y = Rb * math.cos(math.radians(d)), Rb * math.sin(math.radians(d))
        cur = addp(x, y)
        seg.append((prev, cur, tag_of(0.5 * (prev_deg + d))))
        prev, prev_deg = cur, d
        d += step
    seg.append((prev, first, tag_of(0.5 * (prev_deg + 360.0))))

    # assign wall segments: contiguous None-runs cycle through wall_cycle
    wi = 0
    runprev = "x"
    for k, (i, j, t) in enumerate(seg):
        if t is None:
            if runprev is not None:
                wi += 1
            seg[k] = (i, j, wall_cycle[wi % len(wall_cycle)])
        runprev = t

    gpt = [geo.addPoint(x, y, 0, h) for (x, y) in pts]
    glines, bytag = [], {}
    for (i, j, t) in seg:
        L = geo.addLine(gpt[i], gpt[j])
        glines.append(L)
        bytag.setdefault(t, []).append(L)
    loop = geo.addCurveLoop(glines)
    surf = geo.addPlaneSurface([loop])
    geo.synchronize()
    for t, ls in bytag.items():
        gmsh.model.addPhysicalGroup(1, ls, t)
    gmsh.model.addPhysicalGroup(2, [surf], 1)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.8 * h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.model.mesh.generate(2)

    ntag, coord, _ = gmsh.model.mesh.getNodes()
    idx = {int(t): n + 1 for n, t in enumerate(ntag)}
    xy = coord.reshape(-1, 3)
    _, tri_nodes = gmsh.model.mesh.getElementsByType(2)
    tri = tri_nodes.reshape(-1, 3)
    edges = []
    for dim, tag in gmsh.model.getPhysicalGroups(1):
        for ent in gmsh.model.getEntitiesForPhysicalGroup(1, tag):
            et, _, en = gmsh.model.mesh.getElements(1, ent)
            for ett, enn in zip(et, en):
                if ett == 1:
                    for e in enn.reshape(-1, 2):
                        edges.append((idx[int(e[0])], idx[int(e[1])], tag))
    with open(a.out, "w") as f:
        f.write("MeshVersionFormatted\n2\n\nDimension\n2\n\n")
        f.write(f"Vertices\n{len(xy)}\n")
        for p in xy:
            f.write(f"{p[0]:.17g} {p[1]:.17g} 0\n")
        f.write(f"\nTriangles\n{len(tri)}\n")
        for t in tri:
            A, B, C = idx[int(t[0])], idx[int(t[1])], idx[int(t[2])]
            xa, ya = xy[A - 1][:2]; xb, yb = xy[B - 1][:2]; xc, yc = xy[C - 1][:2]
            if (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya) < 0:
                B, C = C, B
            f.write(f"{A} {B} {C} 1\n")
        f.write(f"\nEdges\n{len(edges)}\n")
        for A, B, t in edges:
            f.write(f"{A} {B} {t}\n")
        f.write("\nEnd\n")
    print("vertices", len(xy), "triangles", len(tri), "edges", len(edges))
    print("edge labels", sorted(Counter(t for *_, t in edges).items()))
    gmsh.finalize()


if __name__ == "__main__":
    sys.exit(main())
