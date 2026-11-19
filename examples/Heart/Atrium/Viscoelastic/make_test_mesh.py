#!/usr/bin/env python3
"""Small test mesh with the labels of the Atrium mesh, for AtriumViscoelastic.

Not an atrium: a spherical chamber with four venous stubs, a mitral stub and a
blind pocket in the place of the appendage, wrapped in a solid shell. It
exercises what the driver reads from the real mesh: fluid cells of material 1
inside solid cells of material 4, the fluid-solid interface and the stub walls
labelled 4, the venous caps 1, 2, 3, 5 and the mitral cap 6. Lengths in mm.

    python3 make_test_mesh.py out.mesh [h]

h is the target cell size in mm (default 3). The pocket is the ball of centre
(0.034, 0, 0) m and radius 0.012 m: -av_appendage_center 0.034,0,0
-av_appendage_radius 0.012. Requires the gmsh Python module.
"""
import sys

import gmsh
import numpy as np

R = 20.0                      # chamber radius
VEIN_R, VEIN_L = 5.0, 12.0    # venous stubs, from the upper half of the chamber
MITRAL_R, MITRAL_L = 10.0, 8.0
POCKET_R, POCKET_L = 6.0, 22.0
SHELL = 2.0                   # solid thickness around the chamber


def build(path, h):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ

    fluid = [(3, occ.addSphere(0, 0, 0, R))]
    caps = {}
    # Veins: four stubs through the upper hemisphere, tilted by 35 degrees.
    tilt = np.radians(35.0)
    for label, phi in zip((1, 2, 3, 5), np.radians((45.0, 135.0, 225.0, 315.0))):
        d = np.array([np.sin(tilt) * np.cos(phi), np.sin(tilt) * np.sin(phi), np.cos(tilt)])
        start = (R - 3.0) * d
        length = VEIN_L + 3.0
        fluid.append((3, occ.addCylinder(*start, *(length * d), VEIN_R)))
        caps[label] = start + length * d
    # Mitral stub, downwards.
    start = np.array([0.0, 0.0, -(R - 6.0)])
    fluid.append((3, occ.addCylinder(*start, 0, 0, -(MITRAL_L + 6.0), MITRAL_R)))
    caps[6] = start + np.array([0.0, 0.0, -(MITRAL_L + 6.0)])
    # Blind pocket along +x, closed by a hemisphere.
    start = np.array([R - 4.0, 0.0, 0.0])
    fluid.append((3, occ.addCylinder(*start, POCKET_L - POCKET_R, 0, 0, POCKET_R)))
    fluid.append((3, occ.addSphere(start[0] + POCKET_L - POCKET_R, 0, 0, POCKET_R)))
    lumen, _ = occ.fuse(fluid[:1], fluid[1:])
    shell = occ.addSphere(0, 0, 0, R + SHELL)
    occ.synchronize()
    parts, parents = occ.fragment([(3, shell)], lumen)
    occ.synchronize()

    # Volumes that come from the lumen are fluid; the rest of the shell is solid.
    fluid_tags = {tag for dim, tag in parents[1]}
    solid_tags = {tag for dim, tag in parents[0]} - fluid_tags

    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5 * h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)

    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    index = np.zeros(int(tags.max()) + 1, dtype=np.int64)
    index[tags.astype(np.int64)] = np.arange(1, len(tags) + 1)

    tets = []
    for material, volumes in ((1, fluid_tags), (4, solid_tags)):
        for v in volumes:
            types, _, nodes = gmsh.model.mesh.getElements(3, v)
            for t, n in zip(types, nodes):
                if t != 4:
                    raise SystemExit("unexpected volume element type %d" % t)
                c = index[n.astype(np.int64)].reshape(-1, 4)
                # Positive orientation.
                p = xyz[c - 1]
                vol = np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0])
                c[vol < 0] = c[vol < 0][:, [0, 2, 1, 3]]
                tets.append(np.column_stack([c, np.full(len(c), material)]))
    tets = np.vstack(tets)

    # Faces of the fluid: a venous or mitral cap carries its label; any other
    # face with fluid on one side only (exterior or solid) is the wall, 4.
    triangles = []
    for dim, surface in gmsh.model.getEntities(2):
        up, _ = gmsh.model.getAdjacencies(2, surface)
        n_fluid = sum(1 for v in up if v in fluid_tags)
        if n_fluid != 1:
            continue          # inside the fluid, or a face of the solid alone
        com = np.array(gmsh.model.occ.getCenterOfMass(2, surface))
        label = 4
        if not any(v in solid_tags for v in up):
            for cap, centre in caps.items():
                if np.linalg.norm(com - centre) < 1.0e-3 * R:
                    label = cap
        types, _, nodes = gmsh.model.mesh.getElements(2, surface)
        for t, n in zip(types, nodes):
            if t != 2:
                raise SystemExit("unexpected surface element type %d" % t)
            c = index[n.astype(np.int64)].reshape(-1, 3)
            triangles.append(np.column_stack([c, np.full(len(c), label)]))
    triangles = np.vstack(triangles)
    gmsh.finalize()

    with open(path, "w") as f:
        f.write("MeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n%d\n" % len(xyz))
        for p in xyz:
            f.write("%.16g %.16g %.16g 0\n" % tuple(p))
        f.write("\nTetrahedra\n%d\n" % len(tets))
        np.savetxt(f, tets, fmt="%d")
        f.write("\nTriangles\n%d\n" % len(triangles))
        np.savetxt(f, triangles, fmt="%d")
        f.write("\nEnd\n")
    counts = {int(k): int((triangles[:, 3] == k).sum()) for k in np.unique(triangles[:, 3])}
    print("%s  vertices=%d  tetrahedra: fluid=%d solid=%d  triangles per label=%s" % (
        path, len(xyz), int((tets[:, 4] == 1).sum()), int((tets[:, 4] == 4).sum()), counts))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    build(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 3.0)
