#!/usr/bin/env python3
"""Generate the 3D lesion meshes of ArterialLesion3D as MEDIT "Dimension 3".

The geometry is lesion_3D.geo (same folder): the meridian profile of
lesion_2D.geo, r(x)/R = 1 + a (1 + cos(pi x/ell))/2, a = sqrt(1-S) - 1
(stenosis) or a = Gam - 1 (aneurysm), revolved about the x axis. Lengths are in
units of D; the driver scales the mesh by Config::diameter.

The file is written with three coordinates per vertex, tetrahedra oriented to
positive volume, and every boundary triangle carrying its Gmsh physical tag:

    1 inlet   2 outlet   3 wall (parent vessel)   5 lesion wall   (10 fluid)

Usage:
    python3 make_lesion_mesh_3d.py <case> <rf> <output.mesh> [-setnumber key value ...]
    python3 make_lesion_mesh_3d.py all <rf> <output directory> [suffix]

    case: H, S25, S50, S75, A125, A150, A200    rf: mesh refinement factor
    Extra "-setnumber key value" pairs are passed to gmsh (e.g. Lu, Ld, hb).
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
GEO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lesion_3D.geo")


def build(case, rf, extra=()):
    with tempfile.TemporaryDirectory() as tmp:
        msh = os.path.join(tmp, "lesion.msh")
        cmd = ["gmsh", GEO, "-3", "-o", msh]
        for key, value in dict(CASES[case], rf=rf).items():
            cmd += ["-setnumber", key, str(value)]
        cmd += list(extra)
        subprocess.run(cmd, check=True, capture_output=True)

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(msh)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        xyz = xyz.reshape(-1, 3)
        index = {int(t): i + 1 for i, t in enumerate(tags)}

        tets, faces = [], []
        for dim, phys in gmsh.model.getPhysicalGroups():
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, phys):
                types, _, nodes = gmsh.model.mesh.getElements(dim, entity)
                for etype, conn in zip(types, nodes):
                    if etype == 2 and dim == 2:
                        faces += [(index[int(a)], index[int(b)], index[int(c)], phys)
                                  for a, b, c in conn.reshape(-1, 3)]
                    elif etype == 4 and dim == 3:
                        tets += [(index[int(a)], index[int(b)], index[int(c)], index[int(d)], phys)
                                 for a, b, c, d in conn.reshape(-1, 4)]
        gmsh.finalize()

    if not tets:
        raise RuntimeError("gmsh produced no tetrahedra; check that Lu and Ld exceed the "
                           "lesion half-length ell")

    # Only the referenced vertices, renumbered (gmsh may keep geometry nodes).
    used = sorted({v for t in tets for v in t[:4]})
    renum = {v: i + 1 for i, v in enumerate(used)}
    xyz = xyz[np.array(used) - 1]

    out = []
    for a, b, c, d, ref in tets:
        pa, pb, pc, pd = (xyz[renum[v] - 1] for v in (a, b, c, d))
        vol = np.dot(np.cross(pb - pa, pc - pa), pd - pa)
        t = (renum[a], renum[b], renum[c], renum[d])
        out.append((t[0], t[2], t[1], t[3], ref) if vol < 0 else (*t, ref))
    faces = [(renum[a], renum[b], renum[c], ref) for a, b, c, ref in faces]
    return xyz, out, faces


def write(path, xyz, tets, faces):
    with open(path, "w") as fh:
        fh.write("MeshVersionFormatted\n2\n\nDimension\n3\n\n")
        fh.write("Vertices\n%d\n" % len(xyz))
        for x, y, z in xyz:
            fh.write("%.17g %.17g %.17g 0\n" % (x, y, z))
        fh.write("\nTetrahedra\n%d\n" % len(tets))
        for t in tets:
            fh.write("%d %d %d %d %d\n" % t)
        fh.write("\nTriangles\n%d\n" % len(faces))
        for f in faces:
            fh.write("%d %d %d %d\n" % f)
        fh.write("\nEnd\n")


def report(case, xyz, tets, faces):
    refs = {}
    for *_, ref in faces:
        refs[ref] = refs.get(ref, 0) + 1
    r = np.hypot(xyz[:, 1], xyz[:, 2])
    print("%-5s vertices=%d tetrahedra=%d triangles per tag=%s  x in [%.3f, %.3f]  max r=%.4f"
          % (case, len(xyz), len(tets), dict(sorted(refs.items())),
             xyz[:, 0].min(), xyz[:, 0].max(), r.max()))


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    case, rf, target = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    if case == "all":
        suffix = sys.argv[4] if len(sys.argv) > 4 else "rf%g" % rf
        os.makedirs(target, exist_ok=True)
        for name in CASES:
            xyz, tet, tri = build(name, rf)
            write(os.path.join(target, "%s_pipe_%s.mesh" % (name, suffix)), xyz, tet, tri)
            report(name, xyz, tet, tri)
    else:
        xyz, tet, tri = build(case, rf, sys.argv[4:])
        write(target, xyz, tet, tri)
        report(case, xyz, tet, tri)


if __name__ == "__main__":
    main()
