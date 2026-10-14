#!/usr/bin/env python3
"""
Fluid meshes of a knife gate valve at several openings.

Fluid = straight pipe (bore of the valve) minus the knife, translated along y
to the requested opening. Flow along +z, knife stroke along +y. CAD in mm,
mesh written in metres (MEDIT .mesh, the format Rodin loads).

Opening f in [0, 1] is the knife stroke measured from the seated position:
    y_c(f) = f * (R + R_k - delta),   delta = 0.0079 D
with y_c the centre of the semicircular knife edge (radius R_k), y_c = 0 the
seated (closed) knife, concentric with the seat groove of the liner, and
delta the residual protrusion of the open knife into the bore (1.6 mm in the
FPC8 CAD, the same fraction of D for FPC16). The gap height on the axis is
h = y_c - R_k + R, so h/D ~ f.

Boundary attributes:  1 inlet (z = -Lup),  2 outlet (z = +Ldown),
                      3 pipe wall,         4 knife (gate).
Volume attribute 1.

Usage:  python3 make_meshes.py FPC8.stp FPC16.stp [--coarse]
"""
import argparse
import math
import os
import sys

import gmsh
import numpy as np

VALVES = {
    # name: (step file key, knife volume tag in the STEP, R bore, R_k knife edge,
    #        knife thickness, knife edge centre y in the CAD)
    "FPC8":  dict(knife_vol=4, R=101.6, Rk=105.0, t=10.0, yc_cad=205.0),
    "FPC16": dict(knife_vol=3, R=191.0, Rk=196.0, t=16.0, yc_cad=-78.0),
}

ATTR_INLET, ATTR_OUTLET, ATTR_WALL, ATTR_KNIFE = 1, 2, 3, 4


def build(step, name, f, lup, ldown, hk_rel, hm_rel, hf_rel, out):
    v = VALVES[name]
    R, Rk, t = v["R"], v["Rk"], v["t"]
    D = 2.0 * R
    delta = 0.0079 * D
    yc = f * (R + Rk - delta)
    gap = yc - Rk + R

    gmsh.clear()
    gmsh.model.add(f"{name}_{int(round(100 * f))}")
    ents = gmsh.model.occ.importShapes(step)
    gmsh.model.occ.synchronize()
    vols = [tg for d, tg in ents if d == 3]
    knife = v["knife_vol"]
    gmsh.model.occ.remove([(3, x) for x in vols if x != knife], recursive=True)
    gmsh.model.occ.translate([(3, knife)], 0.0, yc - v["yc_cad"], 0.0)

    z0, z1 = -lup * D, ldown * D
    pipe = gmsh.model.occ.addCylinder(0, 0, z0, 0, 0, z1 - z0, R)
    # OCC puts the seam of the cylinder at (R, 0): exactly where the bevel of
    # the seated FPC16 knife meets the bore (zero-length edge, sliver tets).
    # Rotate the seam to an angle where nothing happens.
    gmsh.model.occ.rotate([(3, pipe)], 0, 0, 0, 0, 0, 1, 1.1)
    fluid, _ = gmsh.model.occ.cut([(3, pipe)], [(3, knife)])
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    fvols = [tg for d, tg in gmsh.model.getEntities(3)]

    # ---- classify boundary faces --------------------------------------------
    tol = 1e-4 * D
    inlet, outlet, wall, knf = [], [], [], []
    for _, s in gmsh.model.getBoundary([(3, x) for x in fvols], oriented=False):
        s = abs(s)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, s)
        if abs(zmin - z0) < tol and abs(zmax - z0) < tol:
            inlet.append(s); continue
        if abs(zmin - z1) < tol and abs(zmax - z1) < tol:
            outlet.append(s); continue
        (u0, w0), (u1, w1) = gmsh.model.getParametrizationBounds(2, s)
        pts = []
        for a, b in [(0.31, 0.27), (0.5, 0.5), (0.73, 0.81)]:
            pts += [u0 + a * (u1 - u0), w0 + b * (w1 - w0)]
        xyz = np.array(gmsh.model.getValue(2, s, pts)).reshape(-1, 3)
        r = np.hypot(xyz[:, 0], xyz[:, 1])
        if np.all(np.abs(r - R) < 1e-6 * R):
            wall.append(s)
        else:
            if zmin < -0.5 * t - tol or zmax > 0.5 * t + tol:
                raise RuntimeError(f"face {s} unclassified")
            knf.append(s)
    for lst, att, nm in [(inlet, ATTR_INLET, "inlet"), (outlet, ATTR_OUTLET, "outlet"),
                         (wall, ATTR_WALL, "wall"), (knf, ATTR_KNIFE, "knife")]:
        gmsh.model.addPhysicalGroup(2, lst, att, nm)
    gmsh.model.addPhysicalGroup(3, fvols, 1, "fluid")

    # ---- size field --------------------------------------------------------
    hk, hm, hf = hk_rel * D, hm_rel * D, hf_rel * D
    fd = gmsh.model.mesh.field
    dist = fd.add("Distance")
    fd.setNumbers(dist, "SurfacesList", knf)
    fd.setNumber(dist, "Sampling", 60)
    thr = fd.add("Threshold")
    fd.setNumber(thr, "InField", dist)
    fd.setNumber(thr, "SizeMin", hk)
    fd.setNumber(thr, "SizeMax", hf)
    fd.setNumber(thr, "DistMin", 1.5 * hk)
    fd.setNumber(thr, "DistMax", 0.6 * D)
    # jet / separation region downstream of the knife
    box = fd.add("Box")
    fd.setNumber(box, "VIn", hm); fd.setNumber(box, "VOut", hf)
    fd.setNumber(box, "XMin", -R); fd.setNumber(box, "XMax", R)
    fd.setNumber(box, "YMin", -R); fd.setNumber(box, "YMax", R)
    fd.setNumber(box, "ZMin", -0.75 * D); fd.setNumber(box, "ZMax", 3.0 * D)
    fd.setNumber(box, "Thickness", 1.5 * D)
    # knife-wall junction curves (thin wedges at the bevel): finer
    def curves(faces):
        return {abs(c) for _, c in gmsh.model.getBoundary([(2, x) for x in faces],
                                                         combined=False, oriented=False)}
    junction = sorted(curves(knf) & curves(wall))
    fields = [thr, box]
    if junction:
        dj = fd.add("Distance")
        fd.setNumbers(dj, "CurvesList", junction)
        fd.setNumber(dj, "Sampling", 200)
        tj = fd.add("Threshold")
        fd.setNumber(tj, "InField", dj)
        fd.setNumber(tj, "SizeMin", 0.5 * hk)
        fd.setNumber(tj, "SizeMax", hf)
        fd.setNumber(tj, "DistMin", 1.0 * hk)
        fd.setNumber(tj, "DistMax", 6.0 * hk)
        fields.append(tj)
    mn = fd.add("Min")
    fd.setNumbers(mn, "FieldsList", fields)
    fd.setAsBackgroundMesh(mn)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.35 * hk)
    gmsh.option.setNumber("Mesh.MeshSizeMax", hf)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("General.NumThreads", os.cpu_count() or 4)
    gmsh.model.mesh.generate(3)

    # ---- export MEDIT in metres --------------------------------------------
    ntags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3) * 1e-3
    idx = {int(tg): i + 1 for i, tg in enumerate(ntags)}
    tets = []
    for vtag in fvols:
        et, _, en = gmsh.model.mesh.getElements(3, vtag)
        for e, n in zip(et, en):
            assert e == 4, "non-tetrahedral element"
            tets.append(np.array(n, dtype=np.int64).reshape(-1, 4))
    tets = np.vstack(tets)
    tris = []
    for lst, att in [(inlet, ATTR_INLET), (outlet, ATTR_OUTLET),
                     (wall, ATTR_WALL), (knf, ATTR_KNIFE)]:
        for s in lst:
            et, _, en = gmsh.model.mesh.getElements(2, s)
            for e, n in zip(et, en):
                assert e == 2
                a = np.array(n, dtype=np.int64).reshape(-1, 3)
                tris.append(np.hstack([a, np.full((a.shape[0], 1), att)]))
    tris = np.vstack(tris)
    m = np.vectorize(idx.get)
    tets_i = m(tets); tris_i = m(tris[:, :3])

    # orientation: positive volume
    p = xyz
    a, b, c, d = (p[tets_i[:, k] - 1] for k in range(4))
    vol6 = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)
    neg = vol6 < 0
    tets_i[neg, 1], tets_i[neg, 2] = tets_i[neg, 2].copy(), tets_i[neg, 1].copy()
    vol6 = np.abs(vol6)

    with open(out, "w") as fh:
        fh.write("MeshVersionFormatted 2\n\nDimension 3\n\n")
        fh.write(f"Vertices\n{len(p)}\n")
        np.savetxt(fh, np.hstack([p, np.zeros((len(p), 1))]), fmt="%.10g %.10g %.10g %d")
        fh.write(f"\nTriangles\n{len(tris_i)}\n")
        np.savetxt(fh, np.hstack([tris_i, tris[:, 3:4]]), fmt="%d")
        fh.write(f"\nTetrahedra\n{len(tets_i)}\n")
        np.savetxt(fh, np.hstack([tets_i, np.ones((len(tets_i), 1), dtype=np.int64)]), fmt="%d")
        fh.write("\nEnd\n")

    # ---- report ------------------------------------------------------------
    # quality: normalised shape (gamma = 6 sqrt(2) V / l_rms^3)
    q = np.array(gmsh.model.mesh.getElementQualities(
        np.concatenate([gmsh.model.mesh.getElements(3, vt)[1][0] for vt in fvols]), "gamma"))
    knife_area = sum(gmsh.model.occ.getMass(2, s) for s in knf) * 1e-6
    # projected area of the knife on the flow direction = blocked area
    bore = math.pi * R * R * 1e-6
    vol = gmsh.model.occ.getMass(3, fvols[0]) if len(fvols) == 1 else sum(
        gmsh.model.occ.getMass(3, x) for x in fvols)
    info = dict(name=name, opening=f, yc_mm=yc, gap_mm=gap, gap_over_D=gap / D,
                vertices=len(p), tets=len(tets_i), triangles=len(tris_i),
                regions=len(fvols), knife_faces=len(knf), knife_area_m2=knife_area,
                bore_area_m2=bore, volume_m3=vol * 1e-9,
                qmin=float(q.min()), q01=float(np.quantile(q, 0.001)),
                hk_mm=hk, hm_mm=hm, hf_mm=hf)
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fpc8"); ap.add_argument("fpc16")
    ap.add_argument("--openings", default="0,25,50,75,100")
    ap.add_argument("--lup", type=float, default=3.0)
    ap.add_argument("--ldown", type=float, default=8.0)
    ap.add_argument("--hk", type=float, default=0.0125, help="size at the knife / D")
    ap.add_argument("--hm", type=float, default=0.0385, help="size in the jet box / D")
    ap.add_argument("--hf", type=float, default=0.1, help="far size / D")
    ap.add_argument("--only", default="")
    ap.add_argument("--outdir", default="meshes")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    files = {"FPC8": a.fpc8, "FPC16": a.fpc16}
    for name, step in files.items():
        if a.only and name != a.only:
            continue
        for o in [int(x) for x in a.openings.split(",")]:
            out = os.path.join(a.outdir, f"{name}_open{o:03d}.mesh")
            info = build(step, name, o / 100.0, a.lup, a.ldown, a.hk, a.hm, a.hf, out)
            print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                           for k, v in info.items()), flush=True)
    gmsh.finalize()


if __name__ == "__main__":
    main()
