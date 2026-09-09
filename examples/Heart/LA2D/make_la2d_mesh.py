#!/usr/bin/env python3
"""Flatten LA2D_rectLAA.mesh (MEDIT ``Dimension 3``, all z = 0) into a genuine
two-dimensional MEDIT mesh.

Rodin takes the space dimension straight from the ``Dimension`` keyword
(``IO/MEDIT.cpp``, ``readDimension``).  Loading the original file therefore
yields ``getSpaceDimension() == 3`` and ``getDimension() == 2``: a surface
embedded in space, whose velocity space would carry three components and whose
boundary normal would not be the in-plane normal the formulation needs.
Dropping the third coordinate is all that is required; nothing else changes.

Triangles are reoriented to positive area, and the attributes are reported so
the labels used by LeftAtrium2D.h can be checked against the file rather than
assumed:

    PV inlets : 6, 8, 10, 12
    MV outlet : 14
    LAA wall  : 2, 3, 4
    body wall : 1, 5, 7, 9, 11, 13

Usage:
    python3 make_la2d_mesh.py <input .mesh> <output .mesh>
"""
import sys
from collections import Counter


def read_sections(path):
    with open(path) as fh:
        lines = fh.read().split("\n")
    i = 0
    sections = {}
    header = {}
    while i < len(lines):
        word = lines[i].strip().split()
        if not word:
            i += 1
            continue
        key = word[0]
        if key in ("MeshVersionFormatted", "Dimension"):
            if len(word) > 1:
                header[key] = int(word[1])
                i += 1
            else:
                header[key] = int(lines[i + 1].strip())
                i += 2
        elif key in ("Vertices", "Edges", "Triangles", "Quadrilaterals",
                     "Tetrahedra", "Corners", "Ridges", "RequiredVertices"):
            n = int(lines[i + 1].strip())
            sections[key] = [lines[i + 2 + k].split() for k in range(n)]
            i += 2 + n
        elif key == "End":
            break
        else:
            i += 1
    return header, sections


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 1
    src, dst = sys.argv[1], sys.argv[2]

    header, sec = read_sections(src)
    sdim = header.get("Dimension", 3)

    verts = []
    for row in sec["Vertices"]:
        coords = [float(x) for x in row[:sdim]]
        attr = int(float(row[sdim]))
        if sdim == 3 and abs(coords[2]) > 1e-12:
            raise SystemExit("the mesh is not planar: z = %g" % coords[2])
        verts.append((coords[0], coords[1], attr))

    edges = [(int(r[0]), int(r[1]), int(r[2])) for r in sec["Edges"]]
    tris = [(int(r[0]), int(r[1]), int(r[2]), int(r[3])) for r in sec["Triangles"]]

    # Positive orientation, so the Jacobian determinant keeps a single sign.
    flipped = 0
    oriented = []
    for a, b, c, t in tris:
        (xa, ya, _), (xb, yb, _), (xc, yc, _) = verts[a - 1], verts[b - 1], verts[c - 1]
        area2 = (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya)
        if area2 < 0.0:
            a, b = b, a
            flipped += 1
        oriented.append((a, b, c, t))

    with open(dst, "w") as fh:
        fh.write("MeshVersionFormatted\n2\n\nDimension\n2\n\n")
        fh.write("Vertices\n%d\n" % len(verts))
        for x, y, t in verts:
            fh.write("%.17g %.17g %d\n" % (x, y, t))
        fh.write("\nTriangles\n%d\n" % len(oriented))
        for a, b, c, t in oriented:
            fh.write("%d %d %d %d\n" % (a, b, c, t))
        fh.write("\nEdges\n%d\n" % len(edges))
        for a, b, t in edges:
            fh.write("%d %d %d\n" % (a, b, t))
        fh.write("\nEnd\n")

    print("vertices %d  triangles %d (%d reoriented)  edges %d"
          % (len(verts), len(oriented), flipped, len(edges)))
    print("edge attributes : %s" % sorted(Counter(t for _, _, t in edges).items()))
    print("cell attributes : %s" % sorted(Counter(t for *_, t in oriented).items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
