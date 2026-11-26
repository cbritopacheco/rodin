#!/usr/bin/env python3
"""Bounding box of the appendage wall edges of a MEDIT 2D mesh, printed as the
-ptag_laa_box option.  Usage: laa_box.py mesh [labels=2,3,4] [scale=1]"""
import sys

path = sys.argv[1]
labels = {int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "2,3,4").split(",")}
scale = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
tok = open(path).read().split()
i, V, E = 0, [], []
while i < len(tok):
    w = tok[i]
    if w == "Vertices":
        n = int(tok[i + 1]); i += 2
        for _ in range(n):
            V.append((float(tok[i]), float(tok[i + 1]))); i += 3
        continue
    if w == "Edges":
        n = int(tok[i + 1]); i += 2
        for _ in range(n):
            a, b, r = int(tok[i]), int(tok[i + 1]), int(tok[i + 2]); i += 3
            if r in labels:
                E += [a - 1, b - 1]
        continue
    i += 1
if not E:
    sys.exit("no edges with labels %s" % sorted(labels))
xs = [V[k][0] * scale for k in E]; ys = [V[k][1] * scale for k in E]
print("-ptag_laa_box %.6g,%.6g,%.6g,%.6g" % (min(xs), max(xs), min(ys), max(ys)))
