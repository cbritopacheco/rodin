#!/usr/bin/env python3
"""
Fluid meshes of a knife gate valve at several openings.

Fluid = straight pipe (bore of the valve) united with the inner cavity of the
liner (REVESTIMIENTO / CAMISA PU) from the STEP, minus the knife translated
along y to the requested opening. The cavity adds, in the plane of the knife,
the seat groove below the bore (radius R_k > R) and the slot the knife slides
in above it, capped 1 mm below the top of the liner (the gland seals there).
The knife fills the slot above its edge; what stays wet is the groove, the
two side pockets between the bore and the slot walls below the knife edge,
and, for the bevelled FPC16 knife, the wedges between the bevel and the slot
walls. Flow along +z, knife stroke along +y. CAD in mm, mesh written in
metres (MEDIT .mesh, the format Rodin loads).

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
    "FPC8":  dict(knife_vol=4, liner_vol=3, R=101.6, Rk=105.0, t=10.0, yc_cad=205.0),
    "FPC16": dict(knife_vol=3, liner_vol=5, R=191.0, Rk=196.0, t=16.0, yc_cad=-78.0),
}

ATTR_INLET, ATTR_OUTLET, ATTR_WALL, ATTR_KNIFE = 1, 2, 3, 4

# Seat (seal) designs, --seat. The CAD seat is a rectangular groove in the
# lower half of the bore, t wide (the knife thickness) and R_k - R deep
# (3.4 mm in the FPC8), the slot of the same width continuing above: a dead
# pocket where solids settle and get crushed by the knife. The alternatives
# keep what seals and retains the knife, the flat floor at R_k and vertical
# walls that the lips press on and that take the pressure load on the closed
# knife (an inclined wall would give a lifting component), and change the
# pocket only where the flow needs it. Lower half (y < 0) only; the slot and
# the side pockets above the axis stay as in the CAD. Flow along +z.
#   cad      the CAD groove (reference).
#   fillet   floor corners rounded, radius FILLET_MM (the knife tip corners
#            are to be rounded alike): no stagnant corners; vertical walls
#            over the rest of the depth.
#   chamfer  lead-in chamfers of CHAMFER_DEG, CHAMFER_MM deep, on both top
#            edges of the groove (bidirectional): the shear layer is turned
#            into the pocket; vertical walls over the lower depth - CHAMFER_MM.
#   upramp   upstream wall replaced by a ramp of UPRAMP_DEG down to the
#            floor; downstream wall vertical, full depth (it takes the
#            pressure load; unidirectional valve): the flow runs down the
#            ramp over the floor and against the downstream wall, nothing is
#            left behind.
#   rampup   A + C: short ramp of RAMP2_DEG on the upstream side (the flow
#            enters the pocket along it), vertical downstream wall that takes
#            the pressure load, floor corners rounded (FILLET_MM at the wall,
#            RAMP2_FR at the floor-ramp junction) and the ramp-bore edge
#            rounded (RAMP2_FB, tangent to the bore).
#   rampdn   the same mirrored: ramp on the downstream side, vertical wall
#            upstream (flow against the ramp).
#   lip      bidirectional: vertical lands LIP_LAND deep above the floor on
#            both sides (the knife is retained and seals on them and on the
#            floor, whatever the flow direction), then ramps of LIP_DEG to
#            the bore on both sides; floor corners rounded LIP_R, ramp-bore
#            edges rounded LIP_RB (tangent to the bore), land-ramp edges
#            obtuse (150 deg). A chamfer deeper than "chamfer" and without
#            its corners.
SEATS = ("cad", "fillet", "chamfer", "upramp", "rampup", "rampdn", "lip")
LIP_LAND = 1.6
LIP_DEG = 30.0
LIP_R = 0.6
LIP_RB = 3.0
RAMP2_DEG = 25.0
RAMP2_FR = 2.0
RAMP2_FB = 4.0
FILLET_MM = 1.5
CHAMFER_DEG = 30.0
CHAMFER_MM = 1.2
UPRAMP_DEG = 15.0
SEAT_Z = 30.0      # mm, half-length along the pipe of the refined seat zone (--hs)
SEAT_LAYER = 3.0   # mm, layer of the bore above the seat included in it


def _fillet(P, A, B, rho):
    """Fillet of radius rho at the corner P between the segments P-A and P-B
    (2D): tangent points on P-A and P-B and the centre."""
    P, A, B = (np.asarray(x, float) for x in (P, A, B))
    u = (A - P) / np.linalg.norm(A - P)
    v = (B - P) / np.linalg.norm(B - P)
    phi = math.acos(max(-1.0, min(1.0, float(u @ v))))
    s_ = rho / math.tan(0.5 * phi)
    w = (u + v) / np.linalg.norm(u + v)
    return P + s_ * u, P + s_ * v, P + rho / math.sin(0.5 * phi) * w


def ramp_profile_solid(R, Rk, t, sign):
    """Pocket of the "rampup" (sign = -1: ramp upstream, z < 0) or "rampdn"
    (sign = +1) seat: a profile in the (z, r) half-plane revolved about the
    pipe axis. Vertical wall at z = -sign t/2, flat floor at r = R_k, ramp of
    RAMP2_DEG to the bore on the other side; floor corners and the ramp-bore
    edge rounded, the latter tangent to the bore. The profile closes inside
    the bore (r = R - 1) so that it fuses with the pipe."""
    occ = gmsh.model.occ
    h = 0.5 * t
    L = (Rk - R) / math.tan(math.radians(RAMP2_DEG))
    zw = -sign * h  # vertical wall
    # floor-ramp corner: the floor stays flat 0.6 mm past the knife plane
    # (|z| = t/2) before the fillet starts, so that it meets the slot of the
    # upper half (walls at |z| = t/2) at y = 0 without a sliver step.
    s2 = RAMP2_FR / math.tan(0.5 * math.radians(180.0 - RAMP2_DEG))
    zf = sign * (h + 0.6 + s2)
    zb = zf + sign * L
    Q0 = (zw, R - 1.0)
    P1, P2, P3 = (zw, Rk), (zf, Rk), (zb, R)
    Pout = (zb + sign * 20.0, R)
    a1, b1, c1 = _fillet(P1, Q0, P2, FILLET_MM)
    a2, b2, c2 = _fillet(P2, P1, P3, RAMP2_FR)
    a3, b3, c3 = _fillet(P3, P2, Pout, RAMP2_FB)
    Q5 = (b3[0], R - 1.0)

    def pt(p):  # (z, r) -> (x = r, y = 0, z)
        return occ.addPoint(float(p[1]), 0.0, float(p[0]))

    pq0, pa1, pb1, pc1 = pt(Q0), pt(a1), pt(b1), pt(c1)
    pa2, pb2, pc2 = pt(a2), pt(b2), pt(c2)
    pa3, pb3, pc3 = pt(a3), pt(b3), pt(c3)
    pq5 = pt(Q5)
    curves = [occ.addLine(pq0, pa1), occ.addCircleArc(pa1, pc1, pb1), occ.addLine(pb1, pa2),
              occ.addCircleArc(pa2, pc2, pb2), occ.addLine(pb2, pa3), occ.addCircleArc(pa3, pc3, pb3),
              occ.addLine(pb3, pq5), occ.addLine(pq5, pq0)]
    face = occ.addPlaneSurface([occ.addCurveLoop(curves)])
    out = occ.revolve([(2, face)], 0, 0, 0, 0, 0, 1, 2 * math.pi)
    return [e for e in out if e[0] == 3]


def profile_solid(corners, radii, R):
    """Solid of revolution about the pipe axis of a pocket profile in the
    (z, r) half-plane: the corner polyline (it starts and ends on the bore,
    r = R, and runs outside it), each corner rounded with its radius, the
    first and last corners rounded tangent to the bore; closed inside the
    bore at r = R - 1."""
    occ = gmsh.model.occ
    n = len(corners)
    ext = [(corners[0][0] - 20.0 * math.copysign(1.0, corners[1][0] - corners[0][0]), R)] + list(corners) + \
          [(corners[-1][0] + 20.0 * math.copysign(1.0, corners[-1][0] - corners[-2][0]), R)]
    arcs = [_fillet(ext[i + 1], ext[i], ext[i + 2], radii[i]) if radii[i] > 0 else None for i in range(n)]

    def pt(p):  # (z, r) -> (x = r, y = 0, z)
        return occ.addPoint(float(p[1]), 0.0, float(p[0]))

    q_start = pt((arcs[0][0][0], R - 1.0))
    q_end = pt((arcs[-1][1][0], R - 1.0))
    curves, prev = [], q_start
    for i, arc in enumerate(arcs):
        if arc is None:  # sharp corner
            pc_ = pt(corners[i])
            curves.append(occ.addLine(prev, pc_))
            prev = pc_
            continue
        a_, b_, c_ = arc
        pa, pb, pc = pt(a_), pt(b_), pt(c_)
        curves.append(occ.addLine(prev, pa))
        curves.append(occ.addCircleArc(pa, pc, pb))
        prev = pb
    curves += [occ.addLine(prev, q_end), occ.addLine(q_end, q_start)]
    face = occ.addPlaneSurface([occ.addCurveLoop(curves)])
    out = occ.revolve([(2, face)], 0, 0, 0, 0, 0, 1, 2 * math.pi)
    return [e for e in out if e[0] == 3]


def lip_profile_solid(R, Rk, t):
    """Pocket of the "lip" seat (see SEATS)."""
    h = 0.5 * t
    rl = Rk - LIP_LAND  # top of the lands
    L = (rl - R) / math.tan(math.radians(LIP_DEG))
    right = [(h, Rk), (h, rl), (h + L, R)]
    left = [(-z, r) for z, r in right[::-1]]
    corners = left + right
    # The land-ramp corners stay sharp (150 deg, obtuse): rounded, the land
    # would leave the knife face tangentially and a zero-angle fluid wedge
    # would form wherever the knife slides along it (sides of the bore at
    # small openings).
    radii = [LIP_RB, 0.0, LIP_R, LIP_R, 0.0, LIP_RB]
    return profile_solid(corners, radii, R)


def seat_solid(seat, R, Rk, t):
    """The pocket of the modified seat in the lower half of the bore (y < 0),
    overlapping the bore (r >= R - 1) so that it fuses cleanly with the pipe.
    Solids of revolution about the pipe axis: cylinders, cones, tori."""
    occ = gmsh.model.occ
    d = Rk - R
    h = 0.5 * t
    if seat == "fillet":
        f = FILLET_MM
        parts = [occ.addCylinder(0, 0, -h, 0, 0, t, Rk - f),
                 occ.addCylinder(0, 0, -h + f, 0, 0, t - 2 * f, Rk),
                 occ.addTorus(0, 0, h - f, Rk - f, f), occ.addTorus(0, 0, -h + f, Rk - f, f)]
    elif seat == "chamfer":
        c, L = CHAMFER_MM, (CHAMFER_MM + 1.0) / math.tan(math.radians(CHAMFER_DEG))
        parts = [occ.addCylinder(0, 0, -h, 0, 0, t, Rk),
                 occ.addCone(0, 0, h, 0, 0, L, R + c, R - 1.0),
                 occ.addCone(0, 0, -h, 0, 0, -L, R + c, R - 1.0)]
    elif seat == "upramp":
        L = (d + 1.0) / math.tan(math.radians(UPRAMP_DEG))
        parts = [occ.addCylinder(0, 0, -h, 0, 0, t, Rk),
                 occ.addCone(0, 0, -h, 0, 0, -L, Rk, R - 1.0)]
    elif seat in ("rampup", "rampdn"):
        parts = None
        solid = ramp_profile_solid(R, Rk, t, -1.0 if seat == "rampup" else 1.0)
    elif seat == "lip":
        parts = None
        solid = lip_profile_solid(R, Rk, t)
    else:
        raise ValueError(seat)
    if parts is not None:
        solid = occ.fuse([(3, parts[0])], [(3, x) for x in parts[1:]])[0]
    half = occ.addBox(-2 * Rk, -2 * Rk, -2 * Rk, 4 * Rk, 2 * Rk, 4 * Rk)
    solid, _ = occ.intersect(solid, [(3, half)])
    return solid


def build(step, name, f, lup, ldown, hk_rel, hm_rel, hf_rel, out, cavity_mode="full", wake_end=3.0,
          seat="cad", hs_rel=0.0):
    v = VALVES[name]
    R, Rk, t = v["R"], v["Rk"], v["t"]
    D = 2.0 * R
    delta = 0.0079 * D
    yc = f * (R + Rk - delta)
    gap = yc - Rk + R

    gmsh.clear()
    # Auxiliary model: the knife alone, to classify the fluid faces.
    gmsh.model.add("knife")
    ents = gmsh.model.occ.importShapes(step)
    gmsh.model.occ.remove([e for e in ents if e[0] == 3 and e[1] != v["knife_vol"]], recursive=True)
    gmsh.model.occ.translate([(3, v["knife_vol"])], 0.0, yc - v["yc_cad"], 0.0)
    gmsh.model.occ.synchronize()
    knife_aux = v["knife_vol"]

    model = f"{name}_{int(round(100 * f))}"
    gmsh.model.add(model)
    ents = gmsh.model.occ.importShapes(step)
    gmsh.model.occ.synchronize()
    vols = [tg for d, tg in ents if d == 3]
    knife, liner = v["knife_vol"], v["liner_vol"]
    lb = gmsh.model.getBoundingBox(3, liner)
    gmsh.model.occ.remove([(3, x) for x in vols if x not in (knife, liner)], recursive=True)
    gmsh.model.occ.translate([(3, knife)], 0.0, yc - v["yc_cad"], 0.0)

    # Inner cavity of the liner: box minus liner, the component on the axis.
    zc = 30.0  # mm, inside the cylindrical bore of both liners (|z| <= 33)
    # full: groove + slot (the CAD); groove: the slot is left out (as if the
    # liner lips sealed on the knife), only the seat groove below the bore.
    ytop = lb[4] - 1.0 if cavity_mode == "full" else 0.0
    if cavity_mode == "none":
        gmsh.model.occ.remove([(3, liner)], recursive=True)
        cavity = []
    else:
        cbox = gmsh.model.occ.addBox(lb[0] - 5, lb[1] - 5, -zc, lb[3] - lb[0] + 10, ytop - lb[1] + 5, 2 * zc)
        parts, _ = gmsh.model.occ.cut([(3, cbox)], [(3, liner)])
        gmsh.model.occ.synchronize()
        probe = [0.0, -0.5 * R, 0.0]  # in the bore, below the axis
        cavity = [p for p in parts if gmsh.model.isInside(3, p[1], probe)]
        assert len(cavity) == 1, "liner cavity not found"
        gmsh.model.occ.remove([p for p in parts if p not in cavity], recursive=True)

    z0, z1 = -lup * D, ldown * D
    pipe = gmsh.model.occ.addCylinder(0, 0, z0, 0, 0, z1 - z0, R)
    # OCC puts the seam of the cylinder at (R, 0): exactly where the bevel of
    # the seated FPC16 knife meets the bore (zero-length edge, sliver tets).
    # Rotate the seam to an angle where nothing happens.
    gmsh.model.occ.rotate([(3, pipe)], 0, 0, 0, 0, 0, 1, 1.1)
    if seat != "cad" and cavity:
        # the CAD groove goes (lower half of the cavity; the pipe restores
        # the bore there) and the pocket of the design takes its place
        lower = gmsh.model.occ.addBox(-2 * Rk, -2 * Rk, -2 * zc, 4 * Rk, 2 * Rk, 4 * zc)
        cavity, _ = gmsh.model.occ.cut(cavity, [(3, lower)])
    duct = gmsh.model.occ.fuse([(3, pipe)], cavity)[0] if cavity else [(3, pipe)]
    if seat != "cad":
        duct = gmsh.model.occ.fuse(duct, seat_solid(seat, R, Rk, t))[0]
    # The knife is as wide as the slot (2 R_k) and its sides are tangent to its
    # semicircular edge, so each side pocket ends in a cusp at (+-R_k, y_c)
    # (zero angle: sliver tetrahedra). Cut the pocket where it is narrower
    # than g: gap(y) ~ (y_c - y)^2 / (2 R_k) gives a height sqrt(2 R_k g).
    # Always outside the bore, since R < R_k - g.
    g = 0.0075 * D
    hcusp = math.sqrt(2.0 * Rk * g)
    # The cap reaches 1.5 g inwards, into the knife, so that the knife edge
    # crosses its bottom face and not its corner.
    cap_in = 1.5 * g
    caps = [gmsh.model.occ.addBox(sx * Rk - (cap_in if sx > 0 else 0.0), yc - hcusp, -t,
                                  cap_in, hcusp + 2 * t, 2 * t)
            for sx in (1.0, -1.0)]
    duct, _ = gmsh.model.occ.cut(duct, [(3, c) for c in caps])
    fluid, _ = gmsh.model.occ.cut(duct, [(3, knife)])
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    fvols = [tg for d, tg in gmsh.model.getEntities(3)]
    vtot = sum(gmsh.model.occ.getMass(3, x) for x in fvols)
    slivers = [x for x in fvols if gmsh.model.occ.getMass(3, x) < 1e-7 * vtot]
    if slivers:  # coincident knife / groove faces of the seated knife
        gmsh.model.occ.remove([(3, x) for x in slivers], recursive=True)
        gmsh.model.occ.synchronize()
        fvols = [tg for d, tg in gmsh.model.getEntities(3)]

    # ---- classify boundary faces --------------------------------------------
    # Knife faces: on a coarse surface mesh (triangles lie on the trimmed
    # face, which parametric samples do not), the centroid of a triangle
    # pushed a little along its normal falls inside the knife (auxiliary
    # model) on one side.
    tol = 1e-4 * D
    eps = 0.05  # mm
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.06 * D)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.002 * D)
    gmsh.model.mesh.generate(2)
    ntags, ncoord, _ = gmsh.model.mesh.getNodes()
    node = dict(zip(ntags.tolist(), ncoord.reshape(-1, 3)))
    inlet, outlet, wall, knf = [], [], [], []
    for _, s in gmsh.model.getBoundary([(3, x) for x in fvols], oriented=False):
        s = abs(s)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, s)
        if abs(zmin - z0) < tol and abs(zmax - z0) < tol:
            inlet.append(s); continue
        if abs(zmin - z1) < tol and abs(zmax - z1) < tol:
            outlet.append(s); continue
        _, _, en = gmsh.model.mesh.getElements(2, s)
        tri = np.array(en[0], dtype=np.int64).reshape(-1, 3)
        pick = tri[np.linspace(0, len(tri) - 1, min(len(tri), 15)).astype(int)]
        gmsh.model.setCurrent("knife")
        votes = 0
        for a, b, c in pick:
            pa, pb, pc = node[a], node[b], node[c]
            ctr = (pa + pb + pc) / 3.0
            nv = np.cross(pb - pa, pc - pa)
            nv /= np.linalg.norm(nv)
            votes += int(gmsh.model.isInside(3, knife_aux, list(ctr + eps * nv)) > 0
                         or gmsh.model.isInside(3, knife_aux, list(ctr - eps * nv)) > 0)
        gmsh.model.setCurrent(model)
        (knf if 2 * votes > len(pick) else wall).append(s)
    gmsh.model.mesh.clear()
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
    fd.setNumber(box, "ZMin", -0.75 * D); fd.setNumber(box, "ZMax", wake_end * D)
    fd.setNumber(box, "Thickness", 1.5 * D)
    # knife-wall junction curves (thin wedges at the bevel): finer
    def curves(faces):
        return {abs(c) for _, c in gmsh.model.getBoundary([(2, x) for x in faces],
                                                         combined=False, oriented=False)}
    junction = sorted(curves(knf) & curves(wall))
    # liner cavity in the plane of the knife (groove, side pockets, wedges)
    cavity_walls = [x for x in wall
                    if gmsh.model.getBoundingBox(2, x)[2] > -0.5 * t - tol
                    and gmsh.model.getBoundingBox(2, x)[5] < 0.5 * t + tol]
    # seat zone (lower half of the bore, |z| <= SEAT_Z, the pocket and a
    # SEAT_LAYER thick layer of the bore above it): finer, hs (e.g. D/170 =
    # 1.2 mm: 3 cells across a 3.4 mm pocket), growing by 1 mm per mm away
    # from the zone. Geometric, so that the flush seat (no pocket) gets the
    # same near-wall resolution as the others.
    fields = [thr, box]
    if hs_rel > 0.0:
        hs = hs_rel * D
        ms = fd.add("MathEval")
        fd.setString(ms, "F", "Min(%g + 1.0*Sqrt(Max(Abs(z)-%g,0)^2 + Max(%g-Sqrt(x*x+y*y),0)^2 + Max(y,0)^2), %g)"
                     % (hs, SEAT_Z, R - SEAT_LAYER, hf))
        fields.append(ms)
    if cavity_walls:
        dc = fd.add("Distance")
        fd.setNumbers(dc, "SurfacesList", cavity_walls)
        fd.setNumber(dc, "Sampling", 80)
        tc = fd.add("Threshold")
        fd.setNumber(tc, "InField", dc)
        fd.setNumber(tc, "SizeMin", 0.6 * hk)
        fd.setNumber(tc, "SizeMax", hf)
        fd.setNumber(tc, "DistMin", 1.0 * hk)
        fd.setNumber(tc, "DistMax", 8.0 * hk)
        fields.append(tc)
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
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(0.35 * hk, 0.5 * hs_rel * D) if hs_rel > 0 else 0.35 * hk)
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
    info = dict(name=name, seat=seat, opening=f, yc_mm=yc, gap_mm=gap, gap_over_D=gap / D,
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
    ap.add_argument("--wake-end", type=float, default=3.0,
                    help="downstream end of the jet box (size hm) / D, from the knife plane")
    ap.add_argument("--only", default="")
    ap.add_argument("--outdir", default="meshes")
    ap.add_argument("--cavity", default="full", choices=["full", "groove", "none"],
                    help="liner cavity: full (groove + slot, as in the CAD), groove, none (straight pipe)")
    ap.add_argument("--seat", default="cad", help="comma-separated seat designs: " + ",".join(SEATS))
    ap.add_argument("--hs", type=float, default=0.0,
                    help="size in the seat pocket / D (0: none; 0.006 = 1.2 mm in the FPC8)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    files = {"FPC8": a.fpc8, "FPC16": a.fpc16}
    for name, step in files.items():
        if a.only and name != a.only:
            continue
        for seat in a.seat.split(","):
            assert seat in SEATS, seat
            for o in [int(x) for x in a.openings.split(",")]:
                tag = name if seat == "cad" else f"{name}-{seat}"
                out = os.path.join(a.outdir, f"{tag}_open{o:03d}.mesh")
                info = build(step, name, o / 100.0, a.lup, a.ldown, a.hk, a.hm, a.hf, out, a.cavity,
                             a.wake_end, seat, a.hs)
                print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                               for k, v in info.items()), flush=True)
    gmsh.finalize()


if __name__ == "__main__":
    main()
