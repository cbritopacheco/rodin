#!/usr/bin/env python3
"""Sanity checks of the axisymmetric MEDIT meshes written by make_lesion_mesh.py --axi."""
import sys, re, numpy as np
from collections import Counter

def read(path):
    tok = open(path).read().split()
    i = 0; sec = {}
    while i < len(tok):
        t = tok[i]
        if t == "Dimension":
            sec["dim"] = int(tok[i+1]); i += 2
        elif t in ("Vertices", "Triangles", "Edges"):
            n = int(tok[i+1]); w = {"Vertices": 3, "Triangles": 4, "Edges": 3}[t]
            arr = np.array(tok[i+2:i+2+n*w], dtype=float).reshape(n, w)
            sec[t] = arr; i += 2 + n*w
        elif t == "End":
            break
        else:
            i += 1
    return sec

def profile(x, case):
    if case.startswith("S"):
        a = np.sqrt(1 - int(case[1:]) / 100) - 1; ell = 1.0
    elif case.startswith("A"):
        a = int(case[1:]) / 100 - 1; ell = 1.5
    else:
        a = 0.0; ell = 1.0
    f = np.where(np.abs(x) <= ell, (1 + np.cos(np.pi * x / ell)) / 2, 0.0)
    return 0.5 * (1 + a * f)

ok_all = True
for path in sys.argv[1:]:
    case = re.search(r"([A-Z]\d*)_axi", path).group(1)
    m = read(path)
    V = m["Vertices"][:, :2]; T = m["Triangles"]; E = m["Edges"]
    tri = T[:, :3].astype(int) - 1; edg = E[:, :2].astype(int) - 1; etag = E[:, 2].astype(int)
    msgs = []
    if m["dim"] != 2: msgs.append("Dimension != 2")
    if V[:, 1].min() < -1e-12: msgs.append("y < 0: %.3e" % V[:, 1].min())
    # orientation and degenerate cells
    a, b, c = V[tri[:, 0]], V[tri[:, 1]], V[tri[:, 2]]
    area = 0.5 * ((b[:, 0]-a[:, 0])*(c[:, 1]-a[:, 1]) - (b[:, 1]-a[:, 1])*(c[:, 0]-a[:, 0]))
    if area.min() <= 0: msgs.append("%d non-positive triangles" % (area <= 0).sum())
    # unused / duplicate vertices
    used = np.zeros(len(V), bool); used[tri.ravel()] = True
    if not used.all(): msgs.append("%d unused vertices" % (~used).sum())
    _, cnt = np.unique(np.round(V, 12), axis=0, return_counts=True)
    if cnt.max() > 1: msgs.append("%d duplicate vertex positions" % (cnt > 1).sum())
    # boundary edges of the triangulation vs tagged edges
    sides = Counter()
    for t in tri:
        for i in range(3):
            sides[tuple(sorted((t[i], t[(i+1) % 3])))] += 1
    bnd = {s for s, k in sides.items() if k == 1}
    tagged = {tuple(sorted(e)) for e in edg}
    if bnd != tagged:
        msgs.append("boundary mismatch: %d untagged boundary sides, %d tagged non-boundary"
                    % (len(bnd - tagged), len(tagged - bnd)))
    tags = set(etag.tolist())
    if tags != {1, 2, 3, 4, 5}: msgs.append("tags %s" % sorted(tags))
    # geometric placement of each tag
    for tag, test, name in [(1, lambda p: np.abs(p[:, 0] + 10) < 1e-9, "inlet x=-10"),
                            (2, lambda p: np.abs(p[:, 0] - 25) < 1e-9, "outlet x=25"),
                            (4, lambda p: np.abs(p[:, 1]) < 1e-9, "axis y=0")]:
        P = V[edg[etag == tag].ravel()]
        if not test(P).all(): msgs.append("tag %d not on %s" % (tag, name))
    for tag in (3, 5):
        P = V[edg[etag == tag].ravel()]
        err = np.abs(P[:, 1] - profile(P[:, 0], case)).max()
        if err > 2e-3: msgs.append("wall tag %d off the profile by %.2e" % (tag, err))
    les = V[edg[etag == 5].ravel()][:, 0]
    ell = 1.5 if case.startswith("A") else 1.0
    if case != "H" and (les.min() < -ell - 1e-6 or les.max() > ell + 1e-6):
        msgs.append("lesion tag outside |x| <= ell")
    # inlet integral of r ds must be R^2/2 = 1/8
    ie = edg[etag == 1]; p, q = V[ie[:, 0]], V[ie[:, 1]]
    irds = (0.5 * (p[:, 1] + q[:, 1]) * np.hypot(*(q - p).T)).sum()
    if abs(irds - 0.125) > 1e-9: msgs.append("inlet int r ds = %.6f != 0.125" % irds)
    # mesh quality
    l = np.stack([np.linalg.norm(b-a, axis=1), np.linalg.norm(c-b, axis=1), np.linalg.norm(a-c, axis=1)], 1)
    q = 4*np.sqrt(3)*area / (l**2).sum(1)   # 1 for equilateral
    hmin = l.min(); hthroat = l[(np.abs(a[:, 0]) < 0.3)].min()
    nin = (etag == 1).sum()
    status = "OK " if not msgs else "BAD"
    ok_all &= not msgs
    print("%s %-5s %-7s V=%6d T=%7d  inlet edges=%2d  h_min=%.4f h_throat=%.4f  q_min=%.2f q<0.3:%d  %s"
          % (status, case, path.split("_axi_")[1].split(".")[0], len(V), len(tri), nin, hmin, hthroat,
             q.min(), (q < 0.3).sum(), "; ".join(msgs)))
sys.exit(0 if ok_all else 1)
