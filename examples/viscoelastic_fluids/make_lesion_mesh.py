#!/usr/bin/env python3
"""Generate the planar lesion meshes of ArterialLesion2D as MEDIT "Dimension 2".

The geometry is lesion_2D.geo (same folder) with planar = 1: a channel of width
D whose walls follow r(z)/R = 1 + a (1 + cos(pi z/ell))/2, a = sqrt(1-S) - 1
(stenosis) or a = Gam - 1 (aneurysm), mirrored about y = 0. Lengths are in units
of D; the driver scales the mesh by Config::diameter.

Rodin reads the space dimension from the MEDIT "Dimension" keyword, so the file
is written with two coordinates per vertex, triangles oriented to positive
area, and every boundary edge carrying its Gmsh physical tag:

    1 inlet   2 outlet   3 wall (parent vessel)   5 lesion wall   (10 fluid)

Usage:
    python3 make_lesion_mesh.py <case> <rf> <output.mesh>
    python3 make_lesion_mesh.py all <rf> <output directory> [suffix]

    case: H, S25, S50, S75, A125, A150, A200    rf: mesh refinement factor
"""
import os
import subprocess
import sys
import tempfile

import gmsh
import numpy as np

CASES = {
    "H":    dict(ltype=0),
    "S25":  dict(ltype=1, S=0.25),
    "S50":  dict(ltype=1, S=0.50),
    "S75":  dict(ltype=1, S=0.75),
    "A125": dict(ltype=2, Gam=1.25),
    "A150": dict(ltype=2, Gam=1.50),
    "A200": dict(ltype=2, Gam=2.00),
}
GEO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lesion_2D.geo")


def build(case, rf):
    # The .geo parameters are DefineConstant's, set the way the command line
    # sets them; the mesh is then read back through the API.
    with tempfile.TemporaryDirectory() as tmp:
        msh = os.path.join(tmp, "lesion.msh")
        cmd = ["gmsh", GEO, "-2", "-o", msh]
        for key, value in dict(CASES[case], planar=1, rf=rf).items():
            cmd += ["-setnumber", key, str(value)]
        subprocess.run(cmd, check=True, capture_output=True)

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(msh)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        xyz = xyz.reshape(-1, 3)
        index = {int(t): i + 1 for i, t in enumerate(tags)}

        triangles, edges = [], []
        for dim, phys in gmsh.model.getPhysicalGroups():
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, phys):
                types, _, nodes = gmsh.model.mesh.getElements(dim, entity)
                for etype, conn in zip(types, nodes):
                    if etype == 1 and dim == 1:
                        edges += [(index[int(a)], index[int(b)], phys)
                                  for a, b in conn.reshape(-1, 2)]
                    elif etype == 2 and dim == 2:
                        triangles += [(index[int(a)], index[int(b)], index[int(c)], phys)
                                      for a, b, c in conn.reshape(-1, 3)]
        gmsh.finalize()

    # Positive orientation.
    out = []
    for a, b, c, ref in triangles:
        pa, pb, pc = xyz[a - 1], xyz[b - 1], xyz[c - 1]
        area = (pb[0] - pa[0]) * (pc[1] - pa[1]) - (pb[1] - pa[1]) * (pc[0] - pa[0])
        out.append((a, c, b, ref) if area < 0 else (a, b, c, ref))
    return xyz, out, edges


def write(path, xyz, triangles, edges):
    with open(path, "w") as fh:
        fh.write("MeshVersionFormatted\n2\n\nDimension\n2\n\n")
        fh.write("Vertices\n%d\n" % len(xyz))
        for x, y, _ in xyz:
            fh.write("%.17g %.17g 0\n" % (x, y))
        fh.write("\nTriangles\n%d\n" % len(triangles))
        for t in triangles:
            fh.write("%d %d %d %d\n" % t)
        fh.write("\nEdges\n%d\n" % len(edges))
        for e in edges:
            fh.write("%d %d %d\n" % e)
        fh.write("\nEnd\n")


def report(case, xyz, triangles, edges):
    refs = {}
    for *_, ref in edges:
        refs[ref] = refs.get(ref, 0) + 1
    print("%-5s vertices=%d triangles=%d edges per tag=%s  y in [%.4f, %.4f]"
          % (case, len(xyz), len(triangles), dict(sorted(refs.items())),
             xyz[:, 1].min(), xyz[:, 1].max()))


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    case, rf, target = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    if case == "all":
        suffix = sys.argv[4] if len(sys.argv) > 4 else "rf%g" % rf
        os.makedirs(target, exist_ok=True)
        for name in CASES:
            xyz, tri, edg = build(name, rf)
            write(os.path.join(target, "%s_planar_%s.mesh" % (name, suffix)), xyz, tri, edg)
            report(name, xyz, tri, edg)
    else:
        xyz, tri, edg = build(case, rf)
        write(target, xyz, tri, edg)
        report(case, xyz, tri, edg)


if __name__ == "__main__":
    main()
