#!/usr/bin/env python3
"""Publishable results of the 3D lesion campaign, and the quantitative
comparison of the fluids N (Newtonian), CY (Carreau-Yasuda) and R (gsPTT,
viscoelastic) in every (artery, geometry).

Usage:
    python3 postprocess_3d.py <campaign root> [--out <dir>] [--arteries a,b] [--geos S75,A200]
                              [--fluids N,CY,R] [--no-figures] [--cycle-fraction 1.0]

    <campaign root>  $OUT/3D of launch_arteries_pc9.sh, i.e. the folder with
                     <artery>/Re<Re>_Wo<Wo>[_waveform_cca]/<geo>/<fluid>/case.xdmf
                     (a case folder itself, or any folder holding such trees,
                     also works: every case.xdmf below the root is used)
    --out            output folder (default <root>/postprocess)

Per case, over the last cycle of snapshots (the period T from run.log):
  - from case.csv: flow rate, pressure drop dp(t), its cycle mean, peak and
    phase lag against the flow rate; Newton iterations and rejected steps;
  - from the snapshots: the wall shear stress vector tau_w (shearStress),
    TAWSS = (1/T) int |tau_w| dt, OSI = (1 - |int tau_w dt| / int |tau_w| dt)/2,
    their axial profiles (circumferential mean over the wall) and the lesion
    wall areas with TAWSS < 0.4 Pa and OSI > 0.2 (the thresholds of the
    atherosclerosis literature); the recirculation length (axial extent of
    reversed wall shear downstream of the lesion) in time; the shear rate
    gd = sqrt(2 eps:eps) from the P1 velocity (cell gradients, volume-averaged
    to the nodes); for R, the polymer stress |sigma_p|_F and its share of the
    viscous stress 2 eta_CY(gd) |eps|_F;
  - figures (midplane z = 0 cuts at the systolic peak and at the strongest
    deceleration of the inflow): |u|, p, gd, |sigma_p|; axial profiles of
    TAWSS and OSI; dp(t), recirculation length, wall shear at the throat.
Across fluids, in each (artery, geometry), the differences CY - N and R - CY:
  relative L2 norms over the last cycle of u, p and tau_w, and the relative
  change of TAWSS, OSI, dp, recirculation length and low-TAWSS area, in
  comparison.csv and comparison.tex. Needs numpy, h5py, matplotlib.
"""
import argparse
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np

# numpy 2 renamed trapz.
_trapz = getattr(np, "trapezoid", None) or np.trapz

# ---- geometry: the lesion profile r(x)/D of lesion_3D.geo ---------------------------
GEOS = {"H": (0, 0.0), "S25": (1, 0.25), "S50": (1, 0.50), "S75": (1, 0.75),
        "A125": (2, 1.25), "A150": (2, 1.50), "A200": (2, 2.00)}
# Carreau-Yasuda of the campaign (the same steady viscosity for every fluid).
CY = dict(eta0=0.056, etaInf=0.00345, lam=1.902, n=0.22, a=1.25)
# Thresholds of the TAWSS / OSI literature (Malek et al. 1999; Lee et al. 2009).
LOW_TAWSS = 0.4     # Pa
HIGH_OSI = 0.2
# A wall face is a cap when its normal is within 10 degrees of the axis.
CAP_ALIGNMENT = 0.985
# A wall node belongs to the lesion when its radius departs from the parent
# radius by more than this fraction (the lesions change it by 50-100%).
LESION_DEVIATION = 0.05


def wall_radius(x, geo):
    """r_wall(x)/D, x in units of D from the lesion centre."""
    ltype, p = GEOS[geo]
    if ltype == 0:
        return np.full_like(x, 0.5)
    a = np.sqrt(1.0 - p) - 1.0 if ltype == 1 else p - 1.0
    ell = 1.0 if ltype == 1 else 1.5
    return np.where(np.abs(x) <= ell, 0.5 * (1.0 + a * (1.0 + np.cos(np.pi * x / ell)) / 2.0), 0.5)


def lesion_half_length(geo):
    return 1.0 if GEOS[geo][0] != 2 else 1.5


def eta_cy(gd):
    return CY["etaInf"] + (CY["eta0"] - CY["etaInf"]) * (1.0 + (CY["lam"] * gd) ** CY["a"]) ** ((CY["n"] - 1.0) / CY["a"])


# ---- reading a case -------------------------------------------------------------------
class Case:
    """One case folder: mesh (merged over the rank pieces), the snapshots of
    the last cycle, the wall and the axis."""

    def __init__(self, folder, cycle_fraction=1.0, root=None):
        self.folder = folder
        self.name = os.path.relpath(folder, root) if root else os.path.relpath(folder)
        self.read_log()
        self.read_layout()
        self.build_mesh()
        self.select_cycle(cycle_fraction)
        self.classify_boundary()

    # -- run.log: D, Ubar, T, fluid
    def read_log(self):
        self.D = self.Ubar = self.T = None
        self.fluid = "?"
        with open(os.path.join(self.folder, "run.log"), errors="replace") as f:
            for line in f:
                line = re.sub(r"\x1b\[[0-9;]*m", "", line)
                m = re.search(r"\[inflow\] D=([\d.e+-]+) m\s+Ubar=([\d.e+-]+) m/s\s+T=([\d.e+-]+) s", line)
                if m:
                    self.D, self.Ubar, self.T = (float(v) for v in m.groups())
                m = re.search(r"\[groups\].*fluid=(\w+)", line)
                if m:
                    self.fluid = m.group(1)
                if self.D and self.fluid != "?":
                    break
        if self.D is None:
            raise RuntimeError("no [inflow] line with D, Ubar, T in %s/run.log" % self.folder)
        self.geo = next((g for g in sorted(GEOS, key=len, reverse=True) if "/%s/" % g in self.folder + "/"), "S75")
        # The fluid tag of the campaign folder: .../<geo>/<N|CY|R>/
        self.tag = os.path.basename(self.folder.rstrip("/"))

    # -- case.xdmf: times and, per step, the pieces (file:path of each field)
    def read_layout(self):
        tree = ET.parse(os.path.join(self.folder, "case.xdmf"))
        self.times, self.steps = [], []
        for grid in tree.getroot().iter("Grid"):
            if grid.get("CollectionType") != "Spatial":
                continue
            t = grid.find("Time")
            self.times.append(float(t.get("Value")) if t is not None else float(len(self.times)))
            pieces = []
            for g in grid.findall("Grid"):
                d = {"geometry": g.find("Geometry/DataItem").text.strip(),
                     "topology": g.find("Topology/DataItem").text.strip()}
                for a in g.findall("Attribute"):
                    d[a.get("Name")] = a.find("DataItem").text.strip()
                pieces.append(d)
            self.steps.append(pieces)
        self.times = np.array(self.times)
        if len(self.times) < 2:
            raise RuntimeError("fewer than two snapshots in %s" % self.folder)

    def read(self, ref):
        file, path = ref.split(":", 1)
        with h5py.File(os.path.join(self.folder, file), "r") as f:
            return f[path][()]

    # -- the mesh: vertices merged over the pieces, tetrahedra
    def build_mesh(self):
        pieces = self.steps[0]
        coords = [self.read(p["geometry"]).reshape(-1, 3) for p in pieces]
        scale = max(np.abs(c).max() for c in coords)
        q = 1.0e-9 * scale
        lookup, verts, self.maps = {}, [], []
        for c in coords:
            keys = np.round(c / q).astype(np.int64)
            m = np.empty(len(c), dtype=np.int64)
            for i, k in enumerate(map(tuple, keys)):
                j = lookup.get(k)
                if j is None:
                    j = lookup[k] = len(verts)
                    verts.append(c[i])
                m[i] = j
            self.maps.append(m)
        self.X = np.array(verts)
        tets = []
        for p, m in zip(pieces, self.maps):
            topo = self.read(p["topology"]).astype(np.int64)
            if not np.all(topo[::5] == 6):
                raise RuntimeError("non-tetrahedral cells in %s" % self.folder)
            tets.append(m[topo.reshape(-1, 5)[:, 1:]])
        self.tets = np.vstack(tets)
        # Cell volumes and P1 gradient operators (per cell: grad phi_i = -n_i / (3 V)).
        x = self.X[self.tets]                       # (nc, 4, 3)
        e1, e2, e3 = x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]
        self.vol = np.abs(np.einsum("ij,ij->i", e1, np.cross(e2, e3))) / 6.0
        # grad of the barycentric coordinates: rows of the inverse Jacobian.
        J = np.stack([e1, e2, e3], axis=2)         # (nc, 3, 3), columns e_k
        Jinv = np.linalg.inv(J)                     # d(lambda_1..3)/dx
        self.G = np.zeros((len(self.tets), 4, 3))
        self.G[:, 1:, :] = Jinv
        self.G[:, 0, :] = -Jinv.sum(axis=1)
        # Node volumes (lumped) for cell-to-node averages.
        self.nodevol = np.zeros(len(self.X))
        np.add.at(self.nodevol, self.tets.ravel(), np.repeat(self.vol / 4.0, 4))

    def field(self, step, name):
        """Nodal field of a step, merged over the pieces (nodes shared by two
        pieces carry the same value)."""
        pieces = self.steps[step]
        if name not in pieces[0]:
            return None
        v0 = self.read(pieces[0][name])
        out = np.zeros((len(self.X),) + v0.shape[1:])
        for p, m in zip(pieces, self.maps):
            out[m] = self.read(p[name])
        return out

    def gradient(self, u):
        """Nodal gradient of a nodal field (scalar or vector): cell P1
        gradients, volume-averaged to the nodes. Returns (nn, ..., 3)."""
        uc = u[self.tets]                           # (nc, 4, ...)
        gc = np.einsum("cij,ci...->c...j", self.G, uc)   # (nc, ..., 3)
        gn = np.zeros((len(self.X),) + gc.shape[1:])
        w = (self.vol / 4.0)[(slice(None),) + (None,) * (gc.ndim - 1)]
        for k in range(4):
            np.add.at(gn, self.tets[:, k], gc * w)
        return gn / self.nodevol[(slice(None),) + (None,) * (gc.ndim - 1)]

    # -- the last cycle
    def select_cycle(self, fraction):
        t_end = self.times[-1]
        self.sel = np.where(self.times >= t_end - fraction * self.T - 1.0e-9)[0]
        if len(self.sel) < 3:
            self.sel = np.arange(len(self.times))
        self.t = self.times[self.sel]
        self.span = self.t[-1] - self.t[0]

    # -- boundary faces, caps and wall; axial coordinate
    def classify_boundary(self):
        faces = np.vstack([self.tets[:, [1, 2, 3]], self.tets[:, [0, 3, 2]],
                           self.tets[:, [0, 1, 3]], self.tets[:, [0, 2, 1]]])
        key = np.sort(faces, axis=1)
        _, idx, counts = np.unique(key, axis=0, return_index=True, return_counts=True)
        bnd = faces[idx[counts == 1]]
        xf = self.X[bnd]
        n = np.cross(xf[:, 1] - xf[:, 0], xf[:, 2] - xf[:, 0])
        area = 0.5 * np.linalg.norm(n, axis=1)
        n /= (2.0 * area)[:, None]
        # The axis: principal direction of the vertices (the pipe is long).
        c = self.X.mean(axis=0)
        w, v = np.linalg.eigh(np.cov((self.X - c).T))
        self.axis = v[:, -1] * np.sign(v[0, -1] or 1.0)
        cap = np.abs(n @ self.axis) > CAP_ALIGNMENT
        self.wall_faces, self.wall_area = bnd[~cap], area[~cap]
        self.wall_nodes = np.unique(self.wall_faces)
        # Axial coordinate in units of D from the lesion centre (the
        # lesion wall is where r < R or r > R; its centre is the midpoint).
        s = (self.X - c) @ self.axis
        r = np.linalg.norm((self.X - c) - np.outer(s, self.axis), axis=1)
        # The lesion: the wall nodes whose radius departs from the parent
        # radius (the median over the wall) by more than LESION_DEVIATION;
        # the axis through the vertex centroid is tilted by O(1e-4), which
        # moves r by far less than that over the pipe.
        rw = r[self.wall_nodes]
        R = np.median(rw)
        lesion = self.wall_nodes[np.abs(rw - R) > LESION_DEVIATION * R]
        s0 = 0.5 * (s[lesion].min() + s[lesion].max()) if len(lesion) else 0.0
        self.xD = (s - s0) / self.D
        self.rD = r / self.D
        self.inlet_x = self.xD.min()
        self.outlet_x = self.xD.max()
        # Wall area per wall node (a third of each adjacent face).
        self.wall_nodearea = np.zeros(len(self.X))
        np.add.at(self.wall_nodearea, self.wall_faces.ravel(), np.repeat(self.wall_area / 3.0, 3))
        # Lesion wall: |x/D| <= ell.
        self.lesion_nodes = self.wall_nodes[np.abs(self.xD[self.wall_nodes]) <= lesion_half_length(self.geo)]


# ---- per-case analysis ---------------------------------------------------------------
def analyse(case, out, figures=True):
    """Metrics and figures of one case. Returns a dict of scalars and the
    fields used by the comparison."""
    t, sel = case.t, case.sel
    # CSV time series of the last cycle.
    csvrows = []
    with open(os.path.join(case.folder, "case.csv")) as f:
        for row in csv.DictReader(f):
            csvrows.append({k: float(v) for k, v in row.items()})
    tc = np.array([r["t"] for r in csvrows])
    last = tc >= t[0] - 1.0e-9
    q = np.array([r["qIn"] for r in csvrows])[last]
    dp = np.array([r["dp"] for r in csvrows])[last]
    tt = tc[last]
    newton = np.array([r["newtonIts"] for r in csvrows])[last]
    retries = np.array([r["retries"] for r in csvrows]).sum()
    # Phase lag of dp against q: the lag of their fundamental harmonics.
    def phase(y):
        w = 2.0 * np.pi / case.T
        return np.arctan2(_trapz(y * np.sin(w * tt), tt), _trapz(y * np.cos(w * tt), tt))
    lag = (phase(dp) - phase(q)) / (2.0 * np.pi) * case.T
    # Wall shear over the cycle.
    tau = np.array([case.field(s, "shearStress") for s in sel])     # (nt, nn, 3)
    if tau[0] is None:
        raise RuntimeError("no shearStress field in %s" % case.folder)
    mag = np.linalg.norm(tau, axis=2)
    dt = np.gradient(t)
    int_mag = _trapz(mag, t, axis=0)
    int_vec = np.linalg.norm(_trapz(tau, t, axis=0), axis=1)
    tawss = int_mag / case.span
    osi = 0.5 * (1.0 - int_vec / np.maximum(int_mag, 1.0e-30))
    osi[int_mag <= 1.0e-12 * int_mag.max()] = 0.0
    wn, ln = case.wall_nodes, case.lesion_nodes
    A_les = case.wall_nodearea[ln].sum()
    low_area = case.wall_nodearea[ln][tawss[ln] < LOW_TAWSS].sum() / A_les
    high_osi_area = case.wall_nodearea[ln][osi[ln] > HIGH_OSI].sum() / A_les
    # Axial profiles (circumferential area-weighted mean in bins of 0.1 D).
    edges = np.arange(np.floor(case.inlet_x), np.ceil(case.outlet_x) + 0.1, 0.1)
    def axial(f):
        w = case.wall_nodearea[wn]
        num = np.histogram(case.xD[wn], edges, weights=f[wn] * w)[0]
        den = np.histogram(case.xD[wn], edges, weights=w)[0]
        return num / np.maximum(den, 1.0e-30)
    xc = 0.5 * (edges[1:] + edges[:-1])
    prof_tawss, prof_osi = axial(tawss), axial(osi)
    # Recirculation length: downstream of the lesion, the farthest wall node
    # with reversed axial wall shear (tau_w . axis < 0), per snapshot.
    ell = lesion_half_length(case.geo)
    down = wn[case.xD[wn] > 0.0]
    taux = tau[:, down, :] @ case.axis
    rec = np.zeros(len(t))
    for i in range(len(t)):
        rev = case.xD[down][taux[i] < 0.0]
        rec[i] = (rev.max() - ell) if len(rev) and rev.max() > ell else 0.0
    # Throat wall shear in time (|x/D| < 0.1 for a stenosis, the neck for an
    # aneurysm) and the instants of the peak and of the minimum wall shear.
    throat = wn[np.abs(case.xD[wn]) < 0.1]
    tau_throat = mag[:, throat].mean(axis=1)
    # Instants: the systolic peak (largest throat wall shear) and the strongest
    # deceleration of the inflow (most negative dq/dt of the CSV flow rate).
    i_peak = int(np.argmax(tau_throat))
    dq = np.gradient(q, tt)
    i_min = int(np.argmin(np.abs(t - tt[int(np.argmin(dq))])))
    # Fields at the two instants.
    inst = {}
    for label, i in (("peak", i_peak), ("decel", i_min)):
        s = sel[i]
        u = case.field(s, "velocity")
        p = case.field(s, "pressure")
        gu = case.gradient(u)                       # (nn, 3, 3): d u_i / d x_j
        eps = 0.5 * (gu + np.swapaxes(gu, 1, 2))
        gd = np.sqrt(2.0 * np.einsum("nij,nij->n", eps, eps))
        d = {"u": u, "p": p, "gd": gd, "t": t[i]}
        sd, so = case.field(s, "elasticStressDiag"), case.field(s, "elasticStressOff")
        if sd is not None:
            sp = np.sqrt((sd ** 2).sum(1) + 2.0 * (so ** 2).sum(1))
            d["sigma_p"] = sp
            d["share"] = sp / np.maximum(2.0 * eta_cy(gd) * np.sqrt(np.einsum("nij,nij->n", eps, eps)), 1.0e-30)
        inst[label] = d
    # Velocity and pressure over the cycle, for the comparison (coarse in
    # time: every snapshot).
    U = np.array([case.field(s, "velocity") for s in sel])
    P = np.array([case.field(s, "pressure") for s in sel])
    metrics = {
        "case": case.name, "fluid": case.tag, "geo": case.geo, "D": case.D, "Ubar": case.Ubar, "T": case.T,
        "snapshots": len(sel), "dp_mean": float(_trapz(dp, tt) / (tt[-1] - tt[0])), "dp_peak": float(dp.max()),
        "dp_lag_T": float(lag / case.T), "q_mean": float(_trapz(q, tt) / (tt[-1] - tt[0])),
        "tawss_max": float(tawss[wn].max()), "tawss_lesion_mean": float((tawss[ln] * case.wall_nodearea[ln]).sum() / A_les),
        "osi_lesion_max": float(osi[ln].max()), "low_tawss_area": float(low_area), "high_osi_area": float(high_osi_area),
        "recirc_max_D": float(rec.max()), "recirc_mean_D": float(rec.mean()),
        "throat_tau_peak": float(tau_throat.max()), "throat_tau_min": float(tau_throat.min()),
        "t_peak": float(t[i_peak]), "t_decel": float(t[i_min]),
        "umax": float(np.linalg.norm(U, axis=2).max()), "gd_max_peak": float(inst["peak"]["gd"].max()),
        "newton_mean": float(newton.mean()), "retries": int(retries),
    }
    if "sigma_p" in inst["peak"]:
        metrics["sigma_p_max"] = float(max(inst["peak"]["sigma_p"].max(), inst["decel"]["sigma_p"].max()))
        for label in ("peak", "decel"):
            sh = inst[label]["share"]
            metrics["share_p95_" + label] = float(np.percentile(sh[np.isfinite(sh)], 95))
            metrics["share_max_" + label] = float(np.nanmax(sh))
    fields = {"t": t, "U": U, "P": P, "tau": tau, "tawss": tawss, "osi": osi, "rec": rec,
              "tt": tt, "dp": dp, "q": q, "xc": xc, "prof_tawss": prof_tawss, "prof_osi": prof_osi,
              "tau_throat": tau_throat, "inst": inst}
    # Per-case files.
    os.makedirs(out, exist_ok=True)
    np.savetxt(os.path.join(out, "axial_profiles.csv"), np.column_stack([xc, prof_tawss, prof_osi]),
               delimiter=",", header="x/D,TAWSS[Pa],OSI", comments="")
    np.savetxt(os.path.join(out, "timeseries.csv"), np.column_stack([t, rec, tau_throat]),
               delimiter=",", header="t[s],recirculation[D],throat_wall_shear[Pa]", comments="")
    if figures:
        case_figures(case, fields, metrics, out)
    return metrics, fields


# ---- figures --------------------------------------------------------------------------
def plane_basis(case):
    """Unit vectors (e_y, e_z) normal to the axis; the cut plane is z = 0."""
    helper = np.array([0.0, 1.0, 0.0]) if abs(case.axis[1]) < 0.9 else np.array([0.0, 0.0, 1.0])
    ey = helper - helper.dot(case.axis) * case.axis
    ey /= np.linalg.norm(ey)
    return ey, np.cross(case.axis, ey)


def midplane(case, fields):
    """Exact cut of the P1 fields by the midplane z = 0 (x along the axis):
    every tetrahedron crossed by the plane gives a triangle or a quadrilateral
    whose corners are the edge crossings, with the fields interpolated
    linearly along the edges. Returns (x/D, y/D, triangles, [cut fields]),
    for matplotlib's tripcolor."""
    c = case.X.mean(axis=0)
    ey, ez = plane_basis(case)
    z = (case.X - c) @ ez
    z = np.where(z == 0.0, 1.0e-12 * case.D, z)   # no vertex exactly on the plane
    zt = z[case.tets]
    mixed = np.where((zt.min(axis=1) < 0.0) & (zt.max(axis=1) > 0.0))[0]
    pts, vals, tris = [], [[] for _ in fields], []
    xy = np.column_stack([case.xD, (case.X - c) @ ey / case.D])

    def crossing(a, b):
        w = z[a] / (z[a] - z[b])
        return xy[a] + w * (xy[b] - xy[a]), [f[a] + w * (f[b] - f[a]) for f in fields]

    for t in mixed:
        v = case.tets[t]
        pos = [i for i in v if z[i] > 0.0]
        neg = [i for i in v if z[i] < 0.0]
        if len(pos) == 1 or len(neg) == 1:
            # One vertex alone on its side: the three edges from it, a triangle.
            a = pos[0] if len(pos) == 1 else neg[0]
            others = neg if len(pos) == 1 else pos
            corners = [crossing(a, b) for b in others]
            n0 = len(pts)
            tris.append([n0, n0 + 1, n0 + 2])
        else:
            # Two and two: the crossings of the edges ac, ad, bd, bc are a
            # convex quadrilateral in that cyclic order.
            (a, b), (cc, d) = pos, neg
            corners = [crossing(a, cc), crossing(a, d), crossing(b, d), crossing(b, cc)]
            n0 = len(pts)
            tris += [[n0, n0 + 1, n0 + 2], [n0, n0 + 2, n0 + 3]]
        for pt, fv in corners:
            pts.append(pt)
            for k in range(len(fields)):
                vals[k].append(fv[k])
    pts = np.array(pts)
    return pts[:, 0], pts[:, 1], np.array(tris), [np.array(v) for v in vals]


def case_figures(case, F, M, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ttl = "%s %s" % (case.geo, case.tag)
    # 1. Midplane cuts.
    names = [("u", "|u| [m/s]"), ("p", "p [Pa]"), ("gd", "shear rate [1/s]")]
    if "sigma_p" in F["inst"]["peak"]:
        names.append(("sigma_p", "|sigma_p| [Pa]"))
    import matplotlib.tri as mtri
    for label in ("peak", "decel"):
        d = F["inst"][label]
        fields = [np.linalg.norm(d["u"], axis=1) if k == "u" else d[k] for k, _ in names]
        x, y, tris, cut = midplane(case, fields)
        triang = mtri.Triangulation(x, y, tris)
        xlim = (-4.0, 10.0)
        inwin = (x >= xlim[0]) & (x <= xlim[1])
        fig, axes = plt.subplots(len(names), 1, figsize=(12, 1.6 * len(names) + 0.6), constrained_layout=True)
        for ax, (k, lab), Z in zip(np.atleast_1d(axes), names, cut):
            if k in ("gd", "sigma_p"):
                Z = np.log10(np.maximum(Z, 1.0e-3))
                lab = "log10 " + lab
            # Colour limits from the plotted window (the caps carry outliers).
            lo, hi = np.percentile(Z[inwin], [1, 99]) if inwin.any() else (Z.min(), Z.max())
            im = ax.tripcolor(triang, Z, shading="gouraud", cmap="viridis", vmin=lo, vmax=hi, antialiased=False)
            ax.set_xlim(*xlim)
            ax.set_aspect("equal")
            ax.set_ylabel("y/D")
            fig.colorbar(im, ax=ax, label=lab, pad=0.01)
        np.atleast_1d(axes)[-1].set_xlabel("x/D")
        fig.suptitle("%s, %s (t = %.3f s)" % (ttl, "systolic peak" if label == "peak" else "peak deceleration", d["t"]))
        fig.savefig(os.path.join(out, "cut_%s.png" % label), dpi=160)
        plt.close(fig)
    # 2. Axial TAWSS and OSI; time series.
    fig, axes = plt.subplots(2, 2, figsize=(11, 6), constrained_layout=True)
    axes[0, 0].plot(F["xc"], F["prof_tawss"]); axes[0, 0].set_ylabel("TAWSS [Pa]"); axes[0, 0].set_yscale("log")
    axes[0, 0].axhline(LOW_TAWSS, ls="--", c="k", lw=0.8)
    axes[0, 1].plot(F["xc"], F["prof_osi"]); axes[0, 1].set_ylabel("OSI"); axes[0, 1].axhline(HIGH_OSI, ls="--", c="k", lw=0.8)
    for ax in axes[0]:
        ax.set_xlabel("x/D"); ax.set_xlim(-4, 10)
    axes[1, 0].plot(F["tt"], F["dp"], label="dp [Pa]")
    ax2 = axes[1, 0].twinx(); ax2.plot(F["tt"], F["q"] * 1e6, c="C1", label="q [ml/s]"); ax2.set_ylabel("q [ml/s]")
    axes[1, 0].set_xlabel("t [s]"); axes[1, 0].set_ylabel("dp [Pa]")
    axes[1, 1].plot(F["t"], F["rec"], label="recirculation [D]")
    ax3 = axes[1, 1].twinx(); ax3.plot(F["t"], F["tau_throat"], c="C1"); ax3.set_ylabel("throat wall shear [Pa]")
    axes[1, 1].set_xlabel("t [s]"); axes[1, 1].set_ylabel("recirculation length [D]")
    fig.suptitle(ttl)
    fig.savefig(os.path.join(out, "profiles.png"), dpi=160)
    plt.close(fig)


def comparison_figures(cases, out):
    """TAWSS/OSI profiles, dp(t) and recirculation of the fluids of one
    (artery, geometry) on the same axes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 6), constrained_layout=True)
    for tag, (case, M, F) in cases.items():
        axes[0, 0].plot(F["xc"], F["prof_tawss"], label=tag)
        axes[0, 1].plot(F["xc"], F["prof_osi"], label=tag)
        axes[1, 0].plot(F["tt"] - F["tt"][0], F["dp"], label=tag)
        axes[1, 1].plot(F["t"] - F["t"][0], F["rec"], label=tag)
    axes[0, 0].set_yscale("log"); axes[0, 0].set_ylabel("TAWSS [Pa]"); axes[0, 0].axhline(LOW_TAWSS, ls="--", c="k", lw=0.8)
    axes[0, 1].set_ylabel("OSI"); axes[0, 1].axhline(HIGH_OSI, ls="--", c="k", lw=0.8)
    for ax in axes[0]:
        ax.set_xlabel("x/D"); ax.set_xlim(-4, 10)
    axes[1, 0].set_xlabel("t - t0 [s]"); axes[1, 0].set_ylabel("dp [Pa]")
    axes[1, 1].set_xlabel("t - t0 [s]"); axes[1, 1].set_ylabel("recirculation length [D]")
    for ax in axes.ravel():
        ax.legend()
    fig.suptitle(os.path.basename(out))
    fig.savefig(os.path.join(out, "comparison.png"), dpi=160)
    plt.close(fig)


# ---- the comparison of fluids ----------------------------------------------------------
def rel_l2(a, b, w=None):
    """||a - b|| / ||b|| over all entries (time, nodes, components), w the
    node weights (volume or area) broadcast over the node axis."""
    d = a - b
    if w is not None:
        shape = [1] * d.ndim
        shape[1] = -1
        d = d * np.sqrt(w).reshape(shape)
        b = b * np.sqrt(w).reshape(shape)
    return float(np.sqrt((d ** 2).sum() / max((b ** 2).sum(), 1.0e-300)))


def compare(cases):
    """cases: {tag: (case, metrics, fields)} of one (artery, geometry).
    The fields of different fluids are on the same mesh and the same output
    instants (the campaign runs them with the same mesh and steps); the
    comparison uses the common snapshots of the last cycle."""
    rows = []
    pairs = [("CY", "N"), ("R", "CY"), ("R", "N")]
    for a, b in pairs:
        if a not in cases or b not in cases:
            continue
        ca, Ma, Fa = cases[a]
        cb, Mb, Fb = cases[b]
        if len(ca.X) != len(cb.X) or np.abs(ca.X - cb.X).max() > 1.0e-9 * np.abs(ca.X).max():
            print("   %s vs %s: different meshes, only the scalar metrics are compared" % (a, b))
            same = False
        else:
            same = True
        n = min(len(Fa["t"]), len(Fb["t"]))
        row = {"pair": "%s-%s" % (a, b), "case": os.path.dirname(ca.name)}
        if same:
            vol = ca.nodevol
            row["u_relL2"] = rel_l2(Fa["U"][-n:], Fb["U"][-n:], vol)
            pa, pb = Fa["P"][-n:], Fb["P"][-n:]
            pa = pa - pa.mean(axis=1, keepdims=True); pb = pb - pb.mean(axis=1, keepdims=True)
            row["p_relL2"] = rel_l2(pa, pb, vol)
            wa = ca.wall_nodearea
            row["tau_relL2"] = rel_l2(Fa["tau"][-n:], Fb["tau"][-n:], wa)
            row["tawss_relL2"] = rel_l2(Fa["tawss"][None], Fb["tawss"][None], wa)
            ln = ca.lesion_nodes
            row["osi_lesion_absdiff_max"] = float(np.abs(Fa["osi"][ln] - Fb["osi"][ln]).max())
            # Peak-instant fields.
            ua, ub = Fa["inst"]["peak"]["u"], Fb["inst"]["peak"]["u"]
            row["u_peak_relL2"] = rel_l2(ua[None], ub[None], vol)
            row["gd_peak_relL2"] = rel_l2(Fa["inst"]["peak"]["gd"][None], Fb["inst"]["peak"]["gd"][None], vol)
        for k in ("dp_mean", "dp_peak", "tawss_lesion_mean", "tawss_max", "recirc_max_D", "recirc_mean_D",
                  "low_tawss_area", "high_osi_area", "osi_lesion_max", "throat_tau_peak", "umax"):
            den = abs(Mb[k])
            row[k + "_relchange"] = (Ma[k] - Mb[k]) / den if den > 1.0e-30 else float("nan")
        row["dp_lag_diff_T"] = Ma["dp_lag_T"] - Mb["dp_lag_T"]
        rows.append(row)
    return rows


def write_tables(metrics, comparisons, out):
    keys = sorted({k for m in metrics for k in m}, key=lambda k: (k not in ("case", "fluid", "geo"), k))
    with open(os.path.join(out, "summary_cases.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for m in metrics:
            w.writerow(m)
    if comparisons:
        ckeys = sorted({k for r in comparisons for k in r}, key=lambda k: (k not in ("case", "pair"), k))
        with open(os.path.join(out, "comparison.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=ckeys)
            w.writeheader()
            for r in comparisons:
                w.writerow(r)
        # LaTeX: the headline differences per pair and case.
        cols = [("u_relL2", r"$\|\Delta u\|/\|u\|$"), ("tau_relL2", r"$\|\Delta\tau_w\|/\|\tau_w\|$"),
                ("tawss_relL2", r"$\Delta$TAWSS"), ("osi_lesion_absdiff_max", r"max $|\Delta$OSI$|$"),
                ("dp_mean_relchange", r"$\Delta\overline{\Delta p}$"), ("dp_peak_relchange", r"$\Delta\,\Delta p_{\max}$"),
                ("recirc_max_D_relchange", r"$\Delta L_{rec}$"), ("low_tawss_area_relchange", r"$\Delta A_{\rm low}$")]
        with open(os.path.join(out, "comparison.tex"), "w") as f:
            f.write("\\begin{tabular}{ll" + "r" * len(cols) + "}\n\\toprule\ncase & pair & " +
                    " & ".join(c[1] for c in cols) + "\\\\\n\\midrule\n")
            for r in comparisons:
                vals = []
                for k, _ in cols:
                    v = r.get(k)
                    vals.append("--" if v is None or not np.isfinite(v) else "%.1f\\%%" % (100.0 * v))
                f.write("%s & %s & %s\\\\\n" % (r["case"].replace("_", r"\_"), r["pair"], " & ".join(vals)))
            f.write("\\bottomrule\n\\end{tabular}\n")
    # LaTeX of the per-case metrics.
    cols = [("fluid", "fluid"), ("dp_mean", r"$\overline{\Delta p}$ [Pa]"), ("dp_peak", r"$\Delta p_{\max}$ [Pa]"),
            ("dp_lag_T", r"lag/T"), ("tawss_max", r"TAWSS$_{\max}$ [Pa]"), ("osi_lesion_max", r"OSI$_{\max}$"),
            ("low_tawss_area", r"$A_{\rm low}$"), ("recirc_max_D", r"$L_{rec}/D$"), ("share_p95_decel", r"$\sigma_p$ share p95")]
    with open(os.path.join(out, "summary_cases.tex"), "w") as f:
        f.write("\\begin{tabular}{l" + "r" * len(cols) + "}\n\\toprule\ncase & " + " & ".join(c[1] for c in cols) + "\\\\\n\\midrule\n")
        for m in metrics:
            vals = []
            for k, _ in cols:
                v = m.get(k, "--")
                vals.append(v if isinstance(v, str) else ("%.3g" % v))
            f.write("%s & %s\\\\\n" % (os.path.dirname(m["case"]).replace("_", r"\_"), " & ".join(vals)))
        f.write("\\bottomrule\n\\end{tabular}\n")


# ---- main -------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--out", default=None)
    ap.add_argument("--arteries", default="")
    ap.add_argument("--geos", default="")
    ap.add_argument("--fluids", default="")
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--cycle-fraction", type=float, default=1.0)
    args = ap.parse_args()
    root = os.path.abspath(args.root)
    out = args.out or os.path.join(root, "postprocess")
    os.makedirs(out, exist_ok=True)
    # Every case.xdmf below the root, grouped by (artery/Re.../geo).
    groups = {}
    for d, _, files in os.walk(root):
        if "case.xdmf" in files and os.path.abspath(d) != out:
            rel = os.path.relpath(d, root)
            parts = rel.split(os.sep)
            tag = parts[-1]
            group = os.sep.join(parts[:-1]) or "case"
            if args.arteries and not any(a in rel for a in args.arteries.split(",")):
                continue
            if args.geos and not any("/%s/" % g in "/" + rel + "/" for g in args.geos.split(",")):
                continue
            if args.fluids and tag not in args.fluids.split(","):
                continue
            groups.setdefault(group, {})[tag] = d
    if not groups:
        sys.exit("no case.xdmf found below %s" % root)
    all_metrics, all_comp = [], []
    for group, tags in sorted(groups.items()):
        print("== %s: %s" % (group, ", ".join(sorted(tags))))
        cases = {}
        gout = os.path.join(out, group.replace(os.sep, "_"))
        for tag, d in sorted(tags.items()):
            try:
                case = Case(d, args.cycle_fraction, root)
                M, F = analyse(case, os.path.join(gout, tag), figures=not args.no_figures)
            except Exception as e:   # one broken case must not stop the campaign
                print("   %s: FAILED (%s)" % (tag, e))
                continue
            print("   %-3s dp mean %.1f Pa, peak %.1f Pa, lag %.3f T; TAWSS max %.2f Pa, OSI max %.3f, "
                  "low-TAWSS area %.1f%%, recirculation max %.2f D%s" % (
                      tag, M["dp_mean"], M["dp_peak"], M["dp_lag_T"], M["tawss_max"], M["osi_lesion_max"],
                      100 * M["low_tawss_area"], M["recirc_max_D"],
                      (", polymer share p95 %.2f" % M["share_p95_decel"]) if "share_p95_decel" in M else ""))
            cases[tag] = (case, M, F)
            all_metrics.append(M)
        if len(cases) > 1:
            comp = compare(cases)
            for r in comp:
                print("   %s: u %s, tau_w %s, TAWSS %s, dp mean %+.1f%%, recirculation %+.1f%%" % (
                    r["pair"], ("%.2f%%" % (100 * r["u_relL2"])) if "u_relL2" in r else "-",
                    ("%.2f%%" % (100 * r["tau_relL2"])) if "tau_relL2" in r else "-",
                    ("%.2f%%" % (100 * r["tawss_relL2"])) if "tawss_relL2" in r else "-",
                    100 * r["dp_mean_relchange"], 100 * r["recirc_max_D_relchange"]))
            all_comp += comp
            if not args.no_figures:
                comparison_figures(cases, gout)
    write_tables(all_metrics, all_comp, out)
    print("== written to %s: summary_cases.csv/.tex%s" % (out, ", comparison.csv/.tex" if all_comp else ""))


if __name__ == "__main__":
    main()
