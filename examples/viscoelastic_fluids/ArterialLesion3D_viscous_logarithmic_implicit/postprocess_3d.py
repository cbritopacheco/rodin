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
    Only finished cases (a DONE file in the folder) are used unless
    --include-running is given; a group compares the fluids it has.

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
    viscous stress 2 eta_CY(gd) |eps|_F; the first and second normal stress
    differences N1 = s11 - s22, N2 = s22 - s33 of the polymer stress in the
    local (flow, gradient, vorticity) frame -- on the wall (shear direction,
    inward normal) -- (N2 = 0 in steady shear for this sPTT: it flags
    non-viscometric flow), and the recoverable shear N1/(2 tau_w) at the
    lesion wall;
  - the kinematics of the jet region (from the lesion to POST_LESION_LENGTH D
    downstream): the velocity gradient invariants Q = (|W|^2 - |S|^2)/2
    (rotation vs strain) and R = -det(grad u) (stretching vs compression of
    the vortices; sheet vs tube forming strain), their volume-weighted joint
    PDF with the Vieillefosse line, the flow type xi = (|S| - |W|)/(|S| +
    |W|) and the extension rate along the flow; for R: the polymer stretch
    tr c - 3, the Pakdel-McKinley criterion M = sqrt(lambda |u| kappa N1 /
    tau12) of elastic instability on curved streamlines (critical ~5.9) and
    the polymer power sigma_p : eps against the viscous one;
  - the spectral content of the jet fluctuations, from the per-step series
    of case.csv (dp, the m = 1 jet amplitude |(<u_y>, <u_z>)|/Ubar, max|u|):
    the part above the INFLOW_HARMONICS cardiac harmonics, its rms, its
    cycle-to-cycle (irregular) rms, its Welch PSD against St = f d/U_jet
    (d the throat diameter, U_jet the peak jet velocity) and its dominant
    Strouhal number;
  - the cycle statistics of the velocity: the mean u_bar and the rms of the
    fluctuation u' = u - u_bar, the fields on which two fluids can be
    compared in an unsteady jet (the instantaneous fields decorrelate
    between runs, even of the same fluid);
  - figures (midplane z = 0 cuts at the systolic peak and at the strongest
    deceleration of the inflow): |u|, p, gd, |sigma_p|; axial profiles of
    TAWSS and OSI; dp(t), recirculation length, wall shear at the throat.
Across fluids, in each (artery, geometry), the differences CY - N and R - CY:
  relative L2 norms over the last cycle of u, p and tau_w, and the relative
  change of TAWSS, OSI, dp, recirculation length, low-TAWSS area, the jet
  fluctuation rms, its dominant Strouhal number and the L2 distance of the
  jet spectra, in comparison.csv and comparison.tex; each difference is
  set against the noise floor of the pair, the larger cycle-to-cycle
  difference (last cycle against the one before, at the same phase) of its
  two cases, and marked not significant when below it. Needs numpy, scipy,
  h5py, matplotlib.
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
# Spectra: the inflow carries the first INFLOW_HARMONICS harmonics of the
# cardiac frequency (run_3d.sh, HARMONICS=10); content above them is the
# jet's own dynamics (shedding, flapping). Welch segments of SPECTRUM_SEGMENTS
# per window, 50 % overlap, Hann window.
INFLOW_HARMONICS = 10
# Relative rms below which a fluctuation is taken as numerical noise (an
# axisymmetric jet gives A1' ~ 1e-6; flapping jets ~ 1e-3 and above).
FLUCTUATION_FLOOR = 1.0e-4
SPECTRUM_SEGMENTS = 4
# Q-R analysis and the post-lesion region: the jet and its wake, from the
# lesion end to POST_LESION_LENGTH diameters downstream.
POST_LESION_LENGTH = 8.0
# Pakdel-McKinley: elastic instability of curved streamlines when
# M = sqrt(lambda |u| / Rc * N1 / tau12) exceeds about M_CRIT (Pakdel and
# McKinley, PRL 77, 1996; McKinley et al., JNNFM 67, 1996: 5.9 for the
# Taylor-Couette and lid-driven cavity flows of Boger fluids).
M_CRIT = 5.9
# The local flow frame needs a direction: velocities (or vorticities, or wall
# shears) below this fraction of their maximum fall back on the vessel frame.
FRAME_SPEED_FLOOR = 1.0e-6
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
        self.lam = 0.0
        self.fluid = "?"
        with open(os.path.join(self.folder, "run.log"), errors="replace") as f:
            for line in f:
                line = re.sub(r"\x1b\[[0-9;]*m", "", line)
                m = re.search(r"\[inflow\] D=([\d.e+-]+) m\s+Ubar=([\d.e+-]+) m/s\s+T=([\d.e+-]+) s", line)
                if m:
                    self.D, self.Ubar, self.T = (float(v) for v in m.groups())
                m = re.search(r"\[inflow\].*\blambda=([\d.e+-]+) s", line)
                if m:
                    self.lam = float(m.group(1))
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

    def wall_normals(self):
        """Inward unit normal at each wall node (area-weighted mean of the
        normals of its wall faces, oriented towards the axis)."""
        if getattr(self, "_wn", None) is None:
            xf = self.X[self.wall_faces]
            n = np.cross(xf[:, 1] - xf[:, 0], xf[:, 2] - xf[:, 0])   # |n| = 2 area
            acc = np.zeros_like(self.X)
            for k in range(3):
                np.add.at(acc, self.wall_faces[:, k], n)
            c = self.X.mean(axis=0)
            d = self.X - c
            radial = d - np.outer(d @ self.axis, self.axis)
            flip = np.einsum("ij,ij->i", acc, radial) > 0.0           # make it point inwards
            acc[flip] *= -1.0
            nrm = np.linalg.norm(acc, axis=1)
            self._wn = np.where(nrm[:, None] > 0.0, acc / np.maximum(nrm, 1.0e-300)[:, None], 0.0)
        return self._wn

    def flow_frame(self, u, gu, tau):
        """Per node, the local viscometric frame (e1, e2, e3) = (flow,
        gradient, vorticity) directions. Inside: e1 = u/|u|, e3 the vorticity
        orthogonalised to e1, e2 = e3 x e1. On the wall (u = 0): e1 the wall
        shear direction, e2 the inward normal orthogonalised to e1. Nodes with
        neither a defined velocity nor vorticity get the frame of the vessel
        (axis, radial, azimuthal)."""
        n = len(self.X)
        c = self.X.mean(axis=0)
        d = self.X - c
        radial = d - np.outer(d @ self.axis, self.axis)
        r = np.linalg.norm(radial, axis=1)
        helper = np.array([0.0, 1.0, 0.0]) if abs(self.axis[1]) < 0.9 else np.array([0.0, 0.0, 1.0])
        fb = helper - helper.dot(self.axis) * self.axis
        fb /= np.linalg.norm(fb)
        er = np.where((r > 1.0e-6 * self.D)[:, None], radial / np.maximum(r, 1.0e-300)[:, None], fb)
        e1 = np.broadcast_to(self.axis, (n, 3)).copy()
        e2 = er.copy()
        speed = np.linalg.norm(u, axis=1)
        w = np.stack([gu[:, 2, 1] - gu[:, 1, 2], gu[:, 0, 2] - gu[:, 2, 0], gu[:, 1, 0] - gu[:, 0, 1]], axis=1)
        inside = speed > FRAME_SPEED_FLOOR * max(speed.max(), 1.0e-30)
        e1[inside] = u[inside] / speed[inside, None]
        w3 = w - np.einsum("ij,ij->i", w, e1)[:, None] * e1
        wn = np.linalg.norm(w3, axis=1)
        ok = inside & (wn > FRAME_SPEED_FLOOR * max(wn.max(), 1.0e-30))
        e3 = np.cross(e1, e2)
        e3[ok] = w3[ok] / wn[ok, None]
        e2[ok] = np.cross(e3[ok], e1[ok])
        # Wall nodes: shear direction and inward normal.
        wnodes = self.wall_nodes
        tw = tau[wnodes]
        twn = np.linalg.norm(tw, axis=1)
        good = twn > FRAME_SPEED_FLOOR * max(twn.max(), 1.0e-30)
        idx = wnodes[good]
        e1[idx] = tw[good] / twn[good, None]
        nin = self.wall_normals()[idx]
        nin = nin - np.einsum("ij,ij->i", nin, e1[idx])[:, None] * e1[idx]
        e2[idx] = nin / np.maximum(np.linalg.norm(nin, axis=1), 1.0e-300)[:, None]
        e3[idx] = np.cross(e1[idx], e2[idx])
        # Re-orthonormalise the default nodes (e1 = axis, e2 = radial).
        rest = ~(ok)
        rest[idx] = False
        e3[rest] = np.cross(e1[rest], e2[rest])
        return e1, e2, e3

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
        # The cycle before (the noise floor: the same case against itself).
        t0 = self.t[0]
        self.sel_prev = np.where((self.times >= t0 - fraction * self.T - 1.0e-9) & (self.times < t0 - 1.0e-9))[0]
        if len(self.sel_prev) < max(3, len(self.sel) // 2):
            self.sel_prev = np.array([], dtype=int)

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


# ---- spectral content of the jet fluctuations -----------------------------------------
def jet_spectra(case, rows):
    """Spectral content of the jet fluctuations from the per-step series of
    case.csv (the XDMF snapshots, ~20 per cycle, are far too sparse: the
    shedding of the throat jet is at St = f d/U_jet ~ 0.3, i.e. 100-300 Hz).

    Signals: the pressure drop dp, the volume means <u_y>, <u_z> (over Ubar;
    both vanish for an axisymmetric flow, so their modulus A1 = |(<u_y>,
    <u_z>)| is the amplitude of the m = 1 deflected/flapping jet) and max|u|.
    Window: the free cycles (cycle >= 1 in case.csv: cycle 0 is the ramp),
    resampled to a uniform step (the step control makes dt variable).
    Decomposition of each signal y(t):
      y = y_per + y'     y_per its first INFLOW_HARMONICS cardiac harmonics
                         (least squares over the window): the driven part;
                         y' the jet's own dynamics (locked shedding included);
      y' = y_cyc + y''   with >= 2 whole cycles, y_cyc the phase average of y'
                         over the cycles: y'' the cycle-to-cycle (irregular,
                         turbulent) part.
    Returns metrics (rms of y' and y'', the dominant frequency of y' and its
    Strouhal number, the fraction of the variance above the inflow
    harmonics) and the spectra (Welch PSD of y')."""
    from scipy.signal import welch
    t = np.array([r["t"] for r in rows])
    cyc = np.array([r.get("cycle", 0.0) for r in rows])
    free = cyc >= 1
    if free.sum() < 50:
        free = t >= t[-1] - case.T - 1.0e-9   # only the ramp cycle: use the last period
    t = t[free]
    dt = float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt)
    fs = 1.0 / dt
    fc = 1.0 / case.T
    ncyc = int(np.floor((tu[-1] - tu[0]) / case.T + 1.0e-6))
    # The jet scale: throat diameter of a stenosis, D otherwise; the jet
    # velocity at the flow peak.
    ltype, par = GEOS[case.geo]
    d = case.D * (np.sqrt(1.0 - par) if ltype == 1 else 1.0)
    qpk = max(r["qIn"] for r in rows)
    Ujet = qpk / (0.25 * np.pi * d * d)
    # The harmonic basis.
    w = 2.0 * np.pi * fc * (tu - tu[0])
    B = [np.ones_like(tu)]
    for k in range(1, INFLOW_HARMONICS + 1):
        B += [np.cos(k * w), np.sin(k * w)]
    B = np.column_stack(B)
    uy = np.array([r.get("uyMean", 0.0) for r in rows])[free]
    uz = np.array([r.get("uzMean", 0.0) for r in rows])[free]
    signals = {"dp": np.array([r["dp"] for r in rows])[free],
               "A1": np.hypot(uy, uz), "uy": uy, "uz": uz,
               "umax": np.array([r["maxU"] for r in rows])[free]}
    nper = max(16, len(tu) // SPECTRUM_SEGMENTS)
    metrics, spectra = {}, {"St_per_Hz": d / Ujet, "fc": fc, "fs": fs, "window_cycles": ncyc}
    for name, y in signals.items():
        yu = np.interp(tu, t, y)
        coef, *_ = np.linalg.lstsq(B, yu, rcond=None)
        yper = B @ coef
        yp = yu - yper
        var_total = np.var(yu)
        f, P = welch(yp, fs=fs, window="hann", nperseg=nper, noverlap=nper // 2, detrend="constant")
        spectra[name] = (f, P)
        cut = f > (INFLOW_HARMONICS + 0.5) * fc
        fdom = float(f[cut][np.argmax(P[cut])]) if cut.any() and P[cut].max() > 0 else float("nan")
        scale = abs(np.mean(yu)) if name in ("dp", "umax") else 1.0    # A1, uy, uz are already /Ubar
        # A dominant frequency of a fluctuation at the noise level is
        # meaningless (an axisymmetric jet has A1' ~ 1e-6): no St then.
        if np.sqrt(np.mean(yp ** 2)) / max(scale, 1.0e-30) < FLUCTUATION_FLOOR:
            fdom = float("nan")
        m = {"rms": float(np.sqrt(np.mean(yp ** 2))) / max(scale, 1.0e-30),
             "fdom_Hz": fdom, "St": fdom * d / Ujet if np.isfinite(fdom) else float("nan"),
             "hf_fraction": float(_trapz(P[cut], f[cut]) / max(var_total, 1.0e-300)) if cut.any() else 0.0}
        if ncyc >= 2:
            # Phase average of y' over the whole cycles, y'' = y' - <y'>.
            n_in = int(round(case.T / dt))
            stack = np.array([yp[k * n_in:(k + 1) * n_in] for k in range(ncyc) if (k + 1) * n_in <= len(yp)])
            ycyc = np.tile(stack.mean(axis=0), len(stack))
            m["rms_irregular"] = float(np.sqrt(np.mean((yp[:len(ycyc)] - ycyc) ** 2))) / max(scale, 1.0e-30)
        else:
            m["rms_irregular"] = float("nan")
        for k, v in m.items():
            metrics["%s_%s" % (name, k)] = v
    metrics["jet_d"] = d
    metrics["jet_U"] = Ujet
    metrics["spectrum_cycles"] = ncyc
    metrics["nyquist_St"] = 0.5 * fs * d / Ujet
    return metrics, spectra


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
    # Recirculation: the region of reversed axial wall shear attached to the
    # lesion, on the circumferential mean of tau_w . axis in the axial bins.
    # Per snapshot: the first reversed bin downstream of the lesion centre
    # (within the lesion or one diameter past it) starts it; the first
    # forward bin after it is the reattachment x_r; L_rec = x_r - ell (0 when
    # it reattaches within the lesion). The farthest reversed node, used
    # before, measured the near-wall Womersley reversal of the whole parent
    # vessel at high Wo, not the separation.
    ell = lesion_half_length(case.geo)
    w_w = case.wall_nodearea[wn]
    den = np.histogram(case.xD[wn], edges, weights=w_w)[0]
    rec = np.zeros(len(t))
    reach_outlet = np.zeros(len(t), dtype=bool)
    # The sign of forward wall shear: the shearStress field follows the
    # traction convention of the solver, so forward flow may give tau . axis
    # < 0; it is fixed by the cycle mean on the parent wall upstream
    # (x/D < -3, before the lesion), where the flow is forward on average.
    taux_all = tau[:, wn, :] @ case.axis
    upstream = case.xD[wn] < -3.0
    forward = np.sign(_trapz((taux_all[:, upstream] * w_w[upstream]).sum(axis=1), t)) or 1.0
    taux_all = forward * taux_all
    for i in range(len(t)):
        tm = np.histogram(case.xD[wn], edges, weights=taux_all[i] * w_w)[0] / np.maximum(den, 1.0e-30)
        valid = den > 0.0
        cand = np.where(valid & (xc > 0.0) & (xc < ell + 1.0) & (tm < 0.0))[0]
        if len(cand) == 0:
            continue
        k = cand[0]
        while k + 1 < len(xc) and (not valid[k + 1] or tm[k + 1] < 0.0):
            k += 1
        reach_outlet[i] = k + 1 >= len(xc) or xc[k] > case.outlet_x - 1.0
        rec[i] = max(edges[k + 1] - ell, 0.0)
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
        # Velocity gradient invariants (incompressible: P = tr(grad u) = 0):
        #   Q = (|W|^2 - |S|^2)/2  > 0 rotation dominates (vortex cores),
        #                          < 0 strain dominates;
        #   R = -det(grad u)       with Q > 0: R < 0 vortex stretching, R > 0
        #                          vortex compression; with Q < 0: R > 0
        #                          biaxial strain (sheet forming), R < 0
        #                          axial strain (tube forming);
        #   the Vieillefosse line D = 27 R^2/4 + Q^3 = 0 separates the
        #   focal (complex eigenvalues, D > 0) from the nodal topologies.
        # Q is scaled by (Ubar/D)^2, R by (Ubar/D)^3 in the figures.
        W = 0.5 * (gu - np.swapaxes(gu, 1, 2))
        S2 = np.einsum("nij,nij->n", eps, eps)
        W2 = np.einsum("nij,nij->n", W, W)
        d["Q"] = 0.5 * (W2 - S2)
        d["R"] = -np.linalg.det(gu)
        # Flow type xi = (|S| - |W|)/(|S| + |W|): 1 extension, 0 shear, -1 rotation.
        d["xi"] = (np.sqrt(S2) - np.sqrt(W2)) / np.maximum(np.sqrt(S2) + np.sqrt(W2), 1.0e-30)
        # Extension rate along the flow, eps_ss = e1 . eps . e1 with e1 = u/|u|.
        sp_ = np.linalg.norm(u, axis=1)
        e1u = u / np.maximum(sp_, 1.0e-30)[:, None]
        d["ext"] = np.einsum("ni,nij,nj->n", e1u, eps, e1u)
        # Streamline curvature kappa = |u x (u . grad u)| / |u|^3.
        acc = np.einsum("nij,nj->ni", gu, u)
        d["kappa"] = np.linalg.norm(np.cross(u, acc), axis=1) / np.maximum(sp_ ** 3, 1.0e-30)
        sd, so = case.field(s, "elasticStressDiag"), case.field(s, "elasticStressOff")
        # The elastic quantities only for the viscoelastic fluid: N (an sPTT
        # at Wi 0.01, so a viscous polymer) and CY also write the stress
        # fields, but their N1, M or stretch are noise or viscous artefacts.
        if sd is not None and case.fluid == "gsptt":
            sp = np.sqrt((sd ** 2).sum(1) + 2.0 * (so ** 2).sum(1))
            d["sigma_p"] = sp
            d["share"] = sp / np.maximum(2.0 * eta_cy(gd) * np.sqrt(np.einsum("nij,nij->n", eps, eps)), 1.0e-30)
            # The polymer stress in the local viscometric frame (e1, e2, e3)
            # = (flow, gradient, vorticity) directions (Case.flow_frame; on
            # the wall: shear direction, inward normal):
            #   N1 = s11 - s22,  N2 = s22 - s33,  s12 the polymer shear.
            # In a developed pipe flow the frame is (x, r, theta). The sPTT
            # used here (no slip parameter) has N2 = 0 in steady simple
            # shear, so N2 is nonzero only where the flow is not viscometric
            # (the lesion, the recirculation, the unsteady jet); the frame of
            # the vessel instead would mix N1 into N2 wherever the flow is
            # not axial (the inclined lesion wall, a swirling jet).
            S = np.empty((len(sd), 3, 3))
            S[:, 0, 0], S[:, 1, 1], S[:, 2, 2] = sd[:, 0], sd[:, 1], sd[:, 2]
            S[:, 0, 1] = S[:, 1, 0] = so[:, 0]
            S[:, 0, 2] = S[:, 2, 0] = so[:, 1]
            S[:, 1, 2] = S[:, 2, 1] = so[:, 2]
            e1, e2, e3 = case.flow_frame(u, gu, tau[i])
            comp = lambda a, b: np.einsum("ni,nij,nj->n", a, S, b)
            s11, s22, s33, s12 = comp(e1, e1), comp(e2, e2), comp(e3, e3), comp(e1, e2)
            d["N1"], d["N2"], d["s12"] = s11 - s22, s22 - s33, s12
            # Polymer stretch tr(c) - 3, c = exp(psi).
            ld, lo = case.field(s, "logConformationDiag"), case.field(s, "logConformationOff")
            if ld is not None:
                Ps = np.empty((len(ld), 3, 3))
                Ps[:, 0, 0], Ps[:, 1, 1], Ps[:, 2, 2] = ld[:, 0], ld[:, 1], ld[:, 2]
                Ps[:, 0, 1] = Ps[:, 1, 0] = lo[:, 0]
                Ps[:, 0, 2] = Ps[:, 2, 0] = lo[:, 1]
                Ps[:, 1, 2] = Ps[:, 2, 1] = lo[:, 2]
                d["stretch"] = np.exp(np.linalg.eigvalsh(Ps)).sum(axis=1) - 3.0
            # Pakdel-McKinley M = sqrt(lambda |u| kappa N1 / tau12), tau12 the
            # total shear stress in the flow frame (polymer + viscous, the
            # viscous part 2 eta_s(gd) eps12 approximated by eta_CY(gd) gd).
            tau12 = np.abs(s12) + eta_cy(gd) * gd
            M2 = case.lam * sp_ * d["kappa"] * np.maximum(d["N1"], 0.0) / np.maximum(tau12, 1.0e-30)
            d["M"] = np.sqrt(M2)
            # Polymer power sigma_p : eps (> 0 stores or dissipates, < 0 returns
            # energy to the flow) and the viscous one 2 eta_CY eps : eps.
            d["power_p"] = np.einsum("nij,nij->n", S, eps)
            d["power_v"] = 2.0 * eta_cy(gd) * S2
        inst[label] = d
    # Velocity and pressure over the cycle, for the comparison (coarse in
    # time: every snapshot).
    U = np.array([case.field(s, "velocity") for s in sel])
    P = np.array([case.field(s, "pressure") for s in sel])
    # Cycle statistics of the velocity: the mean u_bar = (1/T) int u dt and
    # the rms of the fluctuation u' = u - u_bar, |u'|_rms = sqrt((1/T) int
    # |u'|^2 dt). In an unsteady jet these, not the instantaneous fields, are
    # what two fluids can be compared on.
    def cycle_stats(Us, ts):
        span_ = ts[-1] - ts[0]
        um = _trapz(Us, ts, axis=0) / span_
        ur = np.sqrt(_trapz(((Us - um[None]) ** 2).sum(axis=2), ts, axis=0) / span_)
        return um, ur
    umean, urms = cycle_stats(U, t)
    metrics = {
        "case": case.name, "fluid": case.tag, "geo": case.geo, "D": case.D, "Ubar": case.Ubar, "T": case.T,
        "snapshots": len(sel), "dp_mean": float(_trapz(dp, tt) / (tt[-1] - tt[0])), "dp_peak": float(dp.max()),
        "dp_lag_T": float(lag / case.T), "q_mean": float(_trapz(q, tt) / (tt[-1] - tt[0])),
        "tawss_max": float(tawss[wn].max()), "tawss_lesion_mean": float((tawss[ln] * case.wall_nodearea[ln]).sum() / A_les),
        "osi_lesion_max": float(osi[ln].max()), "low_tawss_area": float(low_area), "high_osi_area": float(high_osi_area),
        "recirc_max_D": float(rec.max()), "recirc_mean_D": float(rec.mean()),
        "recirc_peak_D": float(rec[int(np.argmax(tau_throat))]),
        "recirc_reaches_outlet_fraction": float(reach_outlet.mean()),
        "throat_tau_peak": float(tau_throat.max()), "throat_tau_min": float(tau_throat.min()),
        "t_peak": float(t[i_peak]), "t_decel": float(t[i_min]),
        "umax": float(np.linalg.norm(U, axis=2).max()), "gd_max_peak": float(inst["peak"]["gd"].max()),
        "newton_mean": float(newton.mean()), "retries": int(retries),
    }
    smetrics, spectra = jet_spectra(case, csvrows)
    metrics.update(smetrics)
    # Kinematics of the jet region (lesion end to POST_LESION_LENGTH D downstream):
    # volume fractions of rotation (Q > 0) and of extension (xi > 1/2), the
    # largest Weissenberg number of the extension, lambda eps_ss.
    region = (case.xD > -lesion_half_length(case.geo)) & (case.xD < lesion_half_length(case.geo) + POST_LESION_LENGTH)
    vreg = case.nodevol[region]
    for label in ("peak", "decel"):
        di = inst[label]
        metrics["rotation_fraction_" + label] = float(vreg[di["Q"][region] > 0.0].sum() / vreg.sum())
        metrics["extension_fraction_" + label] = float(vreg[di["xi"][region] > 0.5].sum() / vreg.sum())
        metrics["Wi_ext_max_" + label] = float(case.lam * di["ext"][region].max()) if case.lam > 0.0 else float("nan")
        if "M" in di:
            Mr = di["M"][region]
            metrics["M_max_" + label] = float(Mr.max())
            metrics["M_p99_" + label] = float(np.percentile(Mr, 99))
            metrics["M_supercritical_fraction_" + label] = float(vreg[Mr > M_CRIT].sum() / vreg.sum())
            metrics["stretch_max_" + label] = float(di["stretch"].max()) if "stretch" in di else float("nan")
            # Volume integrals of the powers over the jet region.
            Pp = float((di["power_p"][region] * vreg).sum())
            Pv = float((di["power_v"][region] * vreg).sum())
            metrics["polymer_power_ratio_" + label] = Pp / max(Pv, 1.0e-30)
            metrics["polymer_power_returned_fraction_" + label] = float(
                (np.minimum(di["power_p"][region], 0.0) * vreg).sum() / max(Pv, 1.0e-30))
    if "sigma_p" in inst["peak"]:
        metrics["sigma_p_max"] = float(max(inst["peak"]["sigma_p"].max(), inst["decel"]["sigma_p"].max()))
        for label in ("peak", "decel"):
            sh = inst[label]["share"]
            metrics["share_p95_" + label] = float(np.percentile(sh[np.isfinite(sh)], 95))
            metrics["share_max_" + label] = float(np.nanmax(sh))
            di = inst[label]
            metrics["N1_max_" + label] = float(np.abs(di["N1"]).max())
            metrics["N2_max_" + label] = float(np.abs(di["N2"]).max())
            metrics["N1_lesion_wall_max_" + label] = float(np.abs(di["N1"][ln]).max())
            metrics["N2_lesion_wall_max_" + label] = float(np.abs(di["N2"][ln]).max())
            # Recoverable shear N1/(2 tau): the elastic measure of the
            # literature (the local Weissenberg number of the wall flow), with
            # tau the total wall shear |tau_w| at that instant.
            tw = np.linalg.norm(tau[i_peak if label == "peak" else i_min][ln], axis=1)
            sr = np.abs(di["N1"][ln]) / np.maximum(2.0 * tw, 1.0e-30)
            metrics["recoverable_shear_p95_" + label] = float(np.percentile(sr[tw > 1.0e-3 * tw.max()], 95))
            # |N2|/|N1| at the lesion wall: its median (the bulk of the wall is
            # viscometric, ~0) and its 95th percentile (the non-viscometric spots).
            ratio = np.abs(di["N2"][ln]) / np.maximum(np.abs(di["N1"][ln]), 1.0e-30)
            metrics["N2_over_N1_lesion_wall_median_" + label] = float(np.median(ratio))
            metrics["N2_over_N1_lesion_wall_p95_" + label] = float(np.percentile(ratio, 95))
    # The noise floor: the last cycle against the one before it, at the same
    # phase. In a periodic flow the difference is the convergence to the
    # periodic regime; in an unsteady jet it is the cycle-to-cycle
    # variability. A difference between fluids below it cannot be told from
    # the fluid's own variability.
    noise = {}
    if len(case.sel_prev):
        tp = case.times[case.sel_prev]
        match = np.array([int(np.argmin(np.abs(tp + case.T - ti))) for ti in t])
        Up_all = np.array([case.field(sp, "velocity") for sp in case.sel_prev])
        Up = Up_all[match]
        umean_p, urms_p = cycle_stats(Up_all, tp)
        noise["umean_relL2"] = rel_l2(umean[None], umean_p[None], case.nodevol)
        noise["urms_relL2"] = rel_l2(urms[None], urms_p[None], case.nodevol)
        Pp = np.array([case.field(case.sel_prev[j], "pressure") for j in match])
        taup_all = np.array([case.field(sp, "shearStress") for sp in case.sel_prev])
        magp = np.linalg.norm(taup_all, axis=2)
        span_p = tp[-1] - tp[0]
        tawss_p = _trapz(magp, tp, axis=0) / span_p
        noise["u_relL2"] = rel_l2(U, Up, case.nodevol)
        pc = P - P.mean(axis=1, keepdims=True)
        ppc = Pp - Pp.mean(axis=1, keepdims=True)
        noise["p_relL2"] = rel_l2(pc, ppc, case.nodevol)
        noise["tau_relL2"] = rel_l2(tau, taup_all[match], case.wall_nodearea)
        noise["tawss_relL2"] = rel_l2(tawss[None], tawss_p[None], case.wall_nodearea)
        tcs = np.array([r["t"] for r in csvrows])
        dps = np.array([r["dp"] for r in csvrows])
        prev = (tcs >= tp[0] - 1.0e-9) & (tcs <= t[0] + 1.0e-9)
        if prev.sum() > 2:
            dpm_p = _trapz(dps[prev], tcs[prev]) / (tcs[prev][-1] - tcs[prev][0])
            dpm = _trapz(dp, tt) / (tt[-1] - tt[0])
            noise["dp_mean_relchange"] = (dpm - dpm_p) / max(abs(dpm_p), 1.0e-30)
        for k, v in noise.items():
            metrics["noise_" + k] = v
    metrics["urms_over_umean"] = float(np.sqrt((urms ** 2 * case.nodevol).sum() / max(((umean ** 2).sum(axis=1) * case.nodevol).sum(), 1.0e-300)))
    fields = {"t": t, "U": U, "P": P, "umean": umean, "urms": urms, "tau": tau, "tawss": tawss, "osi": osi, "rec": rec, "noise": noise,
              "tt": tt, "dp": dp, "q": q, "xc": xc, "prof_tawss": prof_tawss, "prof_osi": prof_osi,
              "tau_throat": tau_throat, "inst": inst, "spectra": spectra}
    # Per-case files.
    os.makedirs(out, exist_ok=True)
    np.savetxt(os.path.join(out, "axial_profiles.csv"), np.column_stack([xc, prof_tawss, prof_osi]),
               delimiter=",", header="x/D,TAWSS[Pa],OSI", comments="")
    np.savetxt(os.path.join(out, "timeseries.csv"), np.column_stack([t, rec, tau_throat]),
               delimiter=",", header="t[s],recirculation[D],throat_wall_shear[Pa]", comments="")
    f, P = spectra["A1"]
    cols = [f, f * spectra["St_per_Hz"], P] + [spectra[k][1] for k in ("dp", "uy", "uz", "umax")]
    np.savetxt(os.path.join(out, "spectra.csv"), np.column_stack(cols), delimiter=",",
               header="f[Hz],St,PSD_A1,PSD_dp,PSD_uy,PSD_uz,PSD_umax", comments="")
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
        names += [("sigma_p", "|sigma_p| [Pa]"), ("N1", "N1 = s11 - s22 [Pa]"), ("N2", "N2 = s22 - s33 [Pa]")]
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
            # Colour limits from the plotted window (the caps carry outliers);
            # the signed normal stress differences on a symmetric scale.
            lo, hi = np.percentile(Z[inwin], [1, 99]) if inwin.any() else (Z.min(), Z.max())
            cmap = "viridis"
            if k in ("N1", "N2"):
                hi = max(abs(lo), abs(hi), 1.0e-12)
                lo, cmap = -hi, "RdBu_r"
            im = ax.tripcolor(triang, Z, shading="gouraud", cmap=cmap, vmin=lo, vmax=hi, antialiased=False)
            ax.set_xlim(*xlim)
            ax.set_aspect("equal")
            ax.set_ylabel("y/D")
            fig.colorbar(im, ax=ax, label=lab, pad=0.01)
        np.atleast_1d(axes)[-1].set_xlabel("x/D")
        fig.suptitle("%s, %s (t = %.3f s)" % (ttl, "systolic peak" if label == "peak" else "peak deceleration", d["t"]))
        fig.savefig(os.path.join(out, "cut_%s.png" % label), dpi=160)
        plt.close(fig)
    # 1b. Kinematic and elastic invariants on the midplane, and the Q-R
    #     joint PDF of the jet region.
    kin = [("Q", "Q (D/Ubar)^2"), ("xi", "flow type xi"), ("ext", "lambda eps_ss" if case.lam > 0 else "eps_ss [1/s]")]
    if "M" in F["inst"]["peak"]:
        kin += [("M", "Pakdel-McKinley M"), ("stretch", "tr c - 3"), ("power_p", "sigma_p : eps [W/m^3]")]
    region = (case.xD > -lesion_half_length(case.geo)) & (case.xD < lesion_half_length(case.geo) + POST_LESION_LENGTH)
    sQ, sR = (case.Ubar / case.D) ** 2, (case.Ubar / case.D) ** 3
    for label in ("peak", "decel"):
        d = F["inst"][label]
        fields = []
        for k, _ in kin:
            f = d[k]
            if k == "Q":
                f = f / sQ
            if k == "ext" and case.lam > 0:
                f = case.lam * f
            fields.append(f)
        x, y, tris, cut = midplane(case, fields)
        triang = mtri.Triangulation(x, y, tris)
        xlim = (-4.0, 10.0)
        inwin = (x >= xlim[0]) & (x <= xlim[1])
        fig, axes = plt.subplots(len(kin), 1, figsize=(12, 1.6 * len(kin) + 0.6), constrained_layout=True)
        for ax, (k, lab), Z in zip(np.atleast_1d(axes), kin, cut):
            if k == "stretch":
                Z = np.log10(np.maximum(Z, 1.0e-3)); lab = "log10 " + lab
            lo, hi = np.percentile(Z[inwin], [1, 99]) if inwin.any() else (Z.min(), Z.max())
            cmap = "viridis"
            if k in ("Q", "ext", "power_p"):
                hi = max(abs(lo), abs(hi), 1.0e-12); lo, cmap = -hi, "RdBu_r"
            if k == "xi":
                lo, hi, cmap = -1.0, 1.0, "PuOr_r"
            if k == "M":
                lo, hi, cmap = 0.0, max(hi, M_CRIT), "magma"
            im = ax.tripcolor(triang, Z, shading="gouraud", cmap=cmap, vmin=lo, vmax=hi, antialiased=False)
            if k == "M":
                ax.tricontour(triang, Z, levels=[M_CRIT], colors="c", linewidths=0.8)
            ax.set_xlim(*xlim); ax.set_aspect("equal"); ax.set_ylabel("y/D")
            fig.colorbar(im, ax=ax, label=lab, pad=0.01)
        np.atleast_1d(axes)[-1].set_xlabel("x/D")
        fig.suptitle("%s, %s (t = %.3f s): kinematics%s" % (ttl, "systolic peak" if label == "peak" else "peak deceleration",
                                                          d["t"], " and elasticity" if "M" in d else ""))
        fig.savefig(os.path.join(out, "cut_kinematics_%s.png" % label), dpi=160)
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    for ax, label in zip(axes, ("peak", "decel")):
        d = F["inst"][label]
        Qn, Rn = d["Q"][region] / sQ, d["R"][region] / sR
        w = case.nodevol[region]
        qa, ra = np.percentile(np.abs(Qn), 99), np.percentile(np.abs(Rn), 99)
        H, xe, ye = np.histogram2d(Rn, Qn, bins=120, range=[[-ra, ra], [-qa, qa]], weights=w, density=True)
        im = ax.pcolormesh(xe, ye, np.log10(np.maximum(H.T, H[H > 0].min() if (H > 0).any() else 1.0)), cmap="Greys")
        rr = np.linspace(-ra, ra, 400)
        ax.plot(rr, -np.cbrt(27.0 / 4.0 * rr ** 2), "r-", lw=0.8)       # Vieillefosse line D = 0
        ax.axhline(0.0, c="k", lw=0.5); ax.axvline(0.0, c="k", lw=0.5)
        ax.set_xlabel("R (D/Ubar)^3"); ax.set_ylabel("Q (D/Ubar)^2")
        ax.set_title("%s (t = %.3f s)" % ("systolic peak" if label == "peak" else "peak deceleration", d["t"]))
        fig.colorbar(im, ax=ax, label="log10 joint PDF (volume-weighted)")
    fig.suptitle("%s: Q-R in the jet region (x/D < %.0f past the lesion); red: D = 27R^2/4 + Q^3 = 0" % (ttl, POST_LESION_LENGTH))
    fig.savefig(os.path.join(out, "QR.png"), dpi=160)
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
    # 3. Spectra of the jet fluctuations y' (above the inflow harmonics).
    S = F["spectra"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), constrained_layout=True)
    for ax, keys, lab in ((axes[0], ("A1", "uy", "uz"), "PSD of <u>/Ubar [1/Hz]"), (axes[1], ("dp", "umax"), "PSD [unit^2/Hz]")):
        for k in keys:
            f, P = S[k]
            ax.loglog(f[1:] * S["St_per_Hz"], P[1:], label=k)
        ax.axvline((INFLOW_HARMONICS + 0.5) * S["fc"] * S["St_per_Hz"], ls="--", c="k", lw=0.8)
        ax.set_xlabel("St = f d/U_jet")
        ax.set_ylabel(lab)
        ax.legend()
    fig.suptitle("%s: jet fluctuations (%d free cycles; dashed: last inflow harmonic)" % (ttl, S["window_cycles"]))
    fig.savefig(os.path.join(out, "spectra.png"), dpi=160)
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
    # Mean and rms velocity on the midplane, one row per fluid, common scales.
    import matplotlib.tri as mtri
    tags = list(cases)
    cuts = {}
    for tag in tags:
        case, M, F = cases[tag]
        x, y, tris, cut = midplane(case, [np.linalg.norm(F["umean"], axis=1), F["urms"]])
        cuts[tag] = (mtri.Triangulation(x, y, tris), cut, (x >= -4.0) & (x <= 10.0))
    lims = []
    for k in range(2):
        lims.append(max(np.percentile(cuts[tg][1][k][cuts[tg][2]], 99) for tg in tags))
    fig, axes = plt.subplots(len(tags), 2, figsize=(16, 1.5 * len(tags) + 0.8), constrained_layout=True, squeeze=False)
    for i, tag in enumerate(tags):
        tri, cut, _ = cuts[tag]
        for k, lab in enumerate(("|u mean| [m/s]", "u' rms [m/s]")):
            ax = axes[i, k]
            im = ax.tripcolor(tri, cut[k], shading="gouraud", cmap="viridis" if k == 0 else "magma",
                              vmin=0.0, vmax=lims[k], antialiased=False)
            ax.set_xlim(-4.0, 10.0); ax.set_aspect("equal")
            ax.set_ylabel("%s\ny/D" % tag)
            if i == 0:
                ax.set_title(lab)
            if i == len(tags) - 1:
                ax.set_xlabel("x/D")
                fig.colorbar(im, ax=axes[:, k], label=lab, pad=0.01, shrink=0.8)
    fig.suptitle("%s: cycle-mean velocity and rms of the fluctuation" % os.path.basename(out))
    fig.savefig(os.path.join(out, "comparison_mean_rms.png"), dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), constrained_layout=True)
    for tag, (case, M, F) in cases.items():
        S = F["spectra"]
        for ax, k in zip(axes, ("A1", "dp")):
            f, P = S[k]
            ax.loglog(f[1:] * S["St_per_Hz"], P[1:], label="%s (St* %.2f)" % (tag, M.get(k + "_St", float("nan"))))
    for ax, k in zip(axes, ("m = 1 jet amplitude A1", "pressure drop dp")):
        ax.set_xlabel("St = f d/U_jet"); ax.set_ylabel("PSD"); ax.set_title(k); ax.legend()
    fig.savefig(os.path.join(out, "comparison_spectra.png"), dpi=160)
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
    The fluids are compared at the same Re and Wo, in units of the case's own
    scales: u/Ubar, p and tau/(rho Ubar^2), t/T. N has its own reference
    viscosity (eta_s + eta_p != eta_inf), hence its own Ubar and T; for a
    Newtonian fluid this rescaling is exact, and it is the identity between CY
    and R. Snapshots are paired by phase t/T mod 1, not by index."""
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
        row = {"pair": "%s-%s" % (a, b), "case": os.path.dirname(ca.name)}
        # Scales: velocity Ubar, stress rho Ubar^2 (rho is common).
        sua, sub = 1.0 / ca.Ubar, 1.0 / cb.Ubar
        ssa, ssb = sua ** 2, sub ** 2
        # Phase pairing: for each snapshot of a, the snapshot of b nearest in t/T mod 1.
        pha, phb = np.mod(Fa["t"] / ca.T, 1.0), np.mod(Fb["t"] / cb.T, 1.0)
        dph = np.abs(pha[:, None] - phb[None, :]); dph = np.minimum(dph, 1.0 - dph)
        jb = np.argmin(dph, axis=1)
        dsnap = np.median(np.diff(np.sort(phb))) if len(phb) > 1 else 1.0
        ia = np.where(dph[np.arange(len(pha)), jb] <= 0.5 * dsnap + 1.0e-9)[0]
        jb = jb[ia]
        row["phase_pairs"] = len(ia)
        row["phase_mismatch_max"] = float(dph[ia, jb].max()) if len(ia) else float("nan")
        if same and len(ia):
            vol = ca.nodevol
            row["u_relL2"] = rel_l2(Fa["U"][ia] * sua, Fb["U"][jb] * sub, vol)
            row["umean_relL2"] = rel_l2(Fa["umean"][None] * sua, Fb["umean"][None] * sub, vol)
            row["urms_relL2"] = rel_l2(Fa["urms"][None] * sua, Fb["urms"][None] * sub, vol)
            pa, pb = Fa["P"][ia] * ssa, Fb["P"][jb] * ssb
            pa = pa - pa.mean(axis=1, keepdims=True); pb = pb - pb.mean(axis=1, keepdims=True)
            row["p_relL2"] = rel_l2(pa, pb, vol)
            wa = ca.wall_nodearea
            row["tau_relL2"] = rel_l2(Fa["tau"][ia] * ssa, Fb["tau"][jb] * ssb, wa)
            row["tawss_relL2"] = rel_l2(Fa["tawss"][None] * ssa, Fb["tawss"][None] * ssb, wa)
            ln = ca.lesion_nodes
            row["osi_lesion_absdiff_max"] = float(np.abs(Fa["osi"][ln] - Fb["osi"][ln]).max())
            # Peak instant of a, against b at the same phase.
            ipk = int(np.argmin(np.abs(Fa["t"] - Fa["inst"]["peak"]["t"])))
            k = np.where(ia == ipk)[0]
            if len(k):
                row["u_peak_relL2"] = rel_l2(Fa["U"][ipk][None] * sua, Fb["U"][jb[k[0]]][None] * sub, vol)
            # gd is stored only at each case's own peak: compared only if both peaks share the phase.
            dpk = abs(np.mod(Fa["inst"]["peak"]["t"] / ca.T, 1.0) - np.mod(Fb["inst"]["peak"]["t"] / cb.T, 1.0))
            if min(dpk, 1.0 - dpk) <= 0.5 * dsnap + 1.0e-9:
                row["gd_peak_relL2"] = rel_l2(Fa["inst"]["peak"]["gd"][None] / sua, Fb["inst"]["peak"]["gd"][None] / sub, vol)
        # Scalars in the same units; areas, OSI and lengths in D are already dimensionless
        # (the low-TAWSS area keeps its dimensional 0.4 Pa threshold).
        scale = {"dp_mean": (ssa, ssb), "dp_peak": (ssa, ssb), "tawss_lesion_mean": (ssa, ssb), "tawss_max": (ssa, ssb),
                 "throat_tau_peak": (ssa, ssb), "umax": (sua, sub)}
        for k in ("dp_mean", "dp_peak", "tawss_lesion_mean", "tawss_max", "recirc_max_D", "recirc_mean_D", "recirc_peak_D",
                  "low_tawss_area", "high_osi_area", "osi_lesion_max", "throat_tau_peak", "umax"):
            fa_, fb_ = scale.get(k, (1.0, 1.0))
            va, vb = Ma[k] * fa_, Mb[k] * fb_
            den = abs(vb)
            row[k + "_relchange"] = (va - vb) / den if den > 1.0e-30 else float("nan")
        row["dp_lag_diff_T"] = Ma["dp_lag_T"] - Mb["dp_lag_T"]
        # Jet dynamics: change of the fluctuation rms and of the dominant
        # Strouhal number; the L2 distance of the PSDs of A1 on the St axis.
        for k in ("A1_rms", "dp_rms", "A1_rms_irregular", "dp_rms_irregular"):
            fa_, fb_ = (ssa, ssb) if k.startswith("dp") else (1.0, 1.0)
            den = abs(Mb.get(k, float("nan")) * fb_)
            row[k + "_relchange"] = (Ma[k] * fa_ - Mb[k] * fb_) / den if den > 1.0e-30 else float("nan")
        row["A1_St_diff"] = Ma["A1_St"] - Mb["A1_St"]
        for k in ("rotation_fraction_peak", "rotation_fraction_decel", "extension_fraction_peak", "extension_fraction_decel"):
            row[k + "_diff"] = Ma[k] - Mb[k]
        fa, Pa = Fa["spectra"]["A1"]
        fb, Pb = Fb["spectra"]["A1"]
        Pb_on_a = np.interp(fa * Fa["spectra"]["St_per_Hz"], fb * Fb["spectra"]["St_per_Hz"], Pb)
        row["A1_psd_relL2"] = float(np.linalg.norm(Pa - Pb_on_a) / max(np.linalg.norm(Pb_on_a), 1.0e-300))
        # Noise floor of the pair: the larger cycle-to-cycle difference of
        # the two cases; a difference is significant when it exceeds it.
        for k in ("u_relL2", "umean_relL2", "urms_relL2", "p_relL2", "tau_relL2", "tawss_relL2", "dp_mean_relchange"):
            na, nb = Fa["noise"].get(k), Fb["noise"].get(k)
            if na is None or nb is None or k not in row:
                continue
            floor = max(abs(na), abs(nb))
            row["noise_" + k] = floor
            row[k + "_over_noise"] = abs(row[k]) / max(floor, 1.0e-30)
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
        cols = [("u_relL2", r"$\|\Delta u\|/\|u\|$"), ("umean_relL2", r"$\|\Delta\bar u\|/\|\bar u\|$"),
                ("urms_relL2", r"$\|\Delta u'_{\rm rms}\|/\|u'_{\rm rms}\|$"), ("tau_relL2", r"$\|\Delta\tau_w\|/\|\tau_w\|$"),
                ("tawss_relL2", r"$\Delta$TAWSS"), ("osi_lesion_absdiff_max", r"max $|\Delta$OSI$|$"),
                ("dp_mean_relchange", r"$\Delta\overline{\Delta p}$"), ("dp_peak_relchange", r"$\Delta\,\Delta p_{\max}$"),
                ("recirc_max_D_relchange", r"$\Delta L_{rec}$"), ("low_tawss_area_relchange", r"$\Delta A_{\rm low}$"),
                ("A1_rms_relchange", r"$\Delta A_1'$"), ("A1_psd_relL2", r"$\|\Delta S_{A_1}\|$")]
        with open(os.path.join(out, "comparison.tex"), "w") as f:
            f.write("\\begin{tabular}{ll" + "r" * len(cols) + "}\n\\toprule\ncase & pair & " +
                    " & ".join(c[1] for c in cols) + "\\\\\n\\midrule\n")
            for r in comparisons:
                vals = []
                for k, _ in cols:
                    v = r.get(k)
                    if v is None or not np.isfinite(v):
                        vals.append("--")
                        continue
                    cell = "%.1f\\%%" % (100.0 * v)
                    # Below the cycle-to-cycle noise floor: not significant.
                    if r.get(k + "_over_noise", 2.0) <= 1.0:
                        cell += "$^{\\circ}$"
                    vals.append(cell)
                f.write("%s & %s & %s\\\\\n" % (r["case"].replace("_", r"\_"), r["pair"], " & ".join(vals)))
            f.write("\\bottomrule\n\\end{tabular}\n")
            f.write("% $^\\circ$: below the noise floor (the cycle-to-cycle difference of the cases)\n")
    # LaTeX of the per-case metrics.
    cols = [("fluid", "fluid"), ("dp_mean", r"$\overline{\Delta p}$ [Pa]"), ("dp_peak", r"$\Delta p_{\max}$ [Pa]"),
            ("dp_lag_T", r"lag/T"), ("tawss_max", r"TAWSS$_{\max}$ [Pa]"), ("osi_lesion_max", r"OSI$_{\max}$"),
            ("low_tawss_area", r"$A_{\rm low}$"), ("recirc_max_D", r"$L_{rec}/D$"), ("A1_rms", r"$A_1'$"),
            ("A1_St", r"St$^*$"), ("share_p95_decel", r"$\sigma_p$ share p95"),
            ("N1_lesion_wall_max_peak", r"$N_{1,w}$ [Pa]"), ("N2_lesion_wall_max_peak", r"$N_{2,w}$ [Pa]"),
            ("recoverable_shear_p95_peak", r"$N_1/2\tau_w$"), ("M_max_peak", r"$M_{\max}$"),
            ("polymer_power_ratio_decel", r"$P_p/P_v$")]
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
    ap.add_argument("--include-running", action="store_true",
                    help="also cases without DONE (still running or failed): their last cycle may be incomplete")
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
            if not args.include_running and "DONE" not in files:
                print("   skipping %s (no DONE: running or failed; --include-running to use it)" % rel)
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
                  "low-TAWSS area %.1f%%, recirculation %.2f D at the peak, max %.2f D%s%s" % (
                      tag, M["dp_mean"], M["dp_peak"], M["dp_lag_T"], M["tawss_max"], M["osi_lesion_max"],
                      100 * M["low_tawss_area"], M["recirc_peak_D"], M["recirc_max_D"],
                      (" (reaches the outlet in %.0f%% of the cycle)" % (100 * M["recirc_reaches_outlet_fraction"]))
                      if M["recirc_reaches_outlet_fraction"] > 0 else "",
                      (", polymer share p95 %.2f" % M["share_p95_decel"]) if "share_p95_decel" in M else ""))
            if "N1_lesion_wall_max_peak" in M:
                print("       normal stresses at the lesion wall (peak): N1 max %.3g Pa, N2 max %.3g Pa, |N2/N1| median %.3f p95 %.3f, "
                      "N1/(2 tau_w) p95 %.3g" % (M["N1_lesion_wall_max_peak"], M["N2_lesion_wall_max_peak"],
                                                 M["N2_over_N1_lesion_wall_median_peak"], M["N2_over_N1_lesion_wall_p95_peak"],
                                                 M["recoverable_shear_p95_peak"]))
            if "M_max_peak" in M:
                print("       elasticity (peak/decel): M max %.3g/%.3g (critical %.1f, supercritical volume %.2g/%.2g), "
                      "Wi_ext max %.3g/%.3g, stretch max %.3g/%.3g, polymer/viscous power %.3g/%.3g (returned %.3g/%.3g)" % (
                          M["M_max_peak"], M["M_max_decel"], M_CRIT, M["M_supercritical_fraction_peak"],
                          M["M_supercritical_fraction_decel"], M["Wi_ext_max_peak"], M["Wi_ext_max_decel"],
                          M["stretch_max_peak"], M["stretch_max_decel"], M["polymer_power_ratio_peak"],
                          M["polymer_power_ratio_decel"], M["polymer_power_returned_fraction_peak"],
                          M["polymer_power_returned_fraction_decel"]))
            print("       kinematics (peak/decel): rotation-dominated volume %.2f/%.2f, extensional volume %.2f/%.2f" % (
                M["rotation_fraction_peak"], M["rotation_fraction_decel"], M["extension_fraction_peak"], M["extension_fraction_decel"]))
            print("       jet: A1' rms %.2e (irregular %.2e), St* %s; dp' rms %.2e of mean (irregular %.2e), St* %s; %d cycles, Nyquist St %.1f" % (
                M["A1_rms"], M["A1_rms_irregular"], "%.3f" % M["A1_St"] if np.isfinite(M["A1_St"]) else "- (axisymmetric)",
                M["dp_rms"], M["dp_rms_irregular"], "%.3f" % M["dp_St"] if np.isfinite(M["dp_St"]) else "-",
                M["spectrum_cycles"], M["nyquist_St"]))
            cases[tag] = (case, M, F)
            all_metrics.append(M)
        if len(cases) > 1:
            comp = compare(cases)
            for r in comp:
                def withnoise(k, signed=False):
                    if k not in r:
                        return "-"
                    v = ("%+.2f%%" if signed else "%.2f%%") % (100 * r[k])
                    if "noise_" + k in r:
                        v += " (noise %.2f%%%s)" % (100 * r["noise_" + k], "" if r[k + "_over_noise"] > 1.0 else ", n.s.")
                    return v
                print("   %s: u %s, u mean %s, u rms %s, tau_w %s, TAWSS %s, dp mean %s, recirculation %+.1f%%" % (
                    r["pair"], withnoise("u_relL2"), withnoise("umean_relL2"), withnoise("urms_relL2"),
                    withnoise("tau_relL2"), withnoise("tawss_relL2"),
                    withnoise("dp_mean_relchange", True), 100 * r["recirc_max_D_relchange"]))
            all_comp += comp
            if not args.no_figures:
                comparison_figures(cases, gout)
    write_tables(all_metrics, all_comp, out)
    print("== written to %s: summary_cases.csv/.tex%s" % (out, ", comparison.csv/.tex" if all_comp else ""))


if __name__ == "__main__":
    main()
