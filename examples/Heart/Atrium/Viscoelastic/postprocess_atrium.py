#!/usr/bin/env python3
"""Compares the fluids of an AtriumViscoelastic campaign.

    python3 postprocess_atrium.py <campaign folder> [--out DIR] [--cycles N]
        [--fields] [--appendage x,y,z,r] [--no-figures]

The campaign folder holds one folder per fluid (N, CY, VE, TVE, R; any subset),
each with case.csv, case_probes.csv (optional), case.xdmf and run.log, as
written by run_atrium_rheology.sh. Everything is read from those files.

What it computes

1. Periodicity. Relative L2 difference of K(t), qOut(t) and of the appendage
   mean age rate between the last two cycles, by fluid. A comparison of
   fluids is only meaningful below the difference between cycles: the same
   number is the noise floor quoted with every difference between fluids.

2. Cycle means over the last --cycles cycles (default 1): mitral stroke
   volume, peak flow rate, kinetic energy, solvent and polymer dissipation,
   the dissipation of the convective orthogonal subscale and its ratio to the
   physical one, elastic energy, Weissenberg number (maximum, volume fraction
   above 1, appendage mean), wall indices, and the age and shear rate of the
   appendage. Written as metrics.csv and metrics.tex; the differences between
   pairs of fluids (CY-N: shear thinning, VE-CY: elasticity, TVE-VE: memory
   of the structure) as differences.csv with the noise floor of each.

3. Energy budget of each fluid over the last cycle: dK/dt, the power of the
   openings (from the fluxes and mean pressures), solvent dissipation, polymer
   power and subscale dissipation, and the residual of the balance.

4. Probes. Phase average over the cycles after the first, cycle-to-cycle
   fluctuation u' = u - <u>, its intensity, and the power spectral density of
   u' (Welch) in frequency. Slopes -5/3 and -3 are drawn as guides only: a
   probe in a pulsatile flow without a mean advection has no Taylor
   hypothesis, so a slope in frequency is not a slope in wavenumber.

5. --fields: from the snapshots of the last cycle, by cell: Q criterion,
   shear rate, vorticity, and the normal stress differences in the local
   shear frame (e1 along u, e2 along the part of grad|u| normal to u, e3 =
   e1 x e2): N1 = s11 - s22, N2 = s22 - s33, with s the elastic stress.
   Volume statistics go to fields.csv, the cell fields of the last snapshot
   to <fluid>_fields.vtk, and the shell-averaged kinetic energy spectrum of a
   cube inscribed in the chamber to spectrum_<fluid>.csv, with the largest
   wavenumber the mesh resolves, pi/h. Requires h5py and scipy.

A fluid of an Owens family writes aggregateSize; VE (equilibrium) has
n = n_st(gd) by construction.
"""
import argparse
import csv
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

ORDER = ["N", "CY", "R", "VE", "TVE"]
LABEL = {"N": "Newtonian", "CY": "Carreau-Yasuda", "R": "gsPTT",
         "VE": "Owens, equilibrium structure", "TVE": "Owens"}
PAIRS = [("CY", "N", "shear thinning"), ("VE", "CY", "elasticity"),
         ("TVE", "VE", "structural memory"), ("R", "CY", "elasticity (gsPTT)")]


def read_csv(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    head, data = rows[0], rows[1:]
    data = [r for r in data if len(r) == len(head)]
    a = np.array(data, dtype=float)
    return {name: a[:, i] for i, name in enumerate(head)}


def read_period(folder):
    period = None
    log = os.path.join(folder, "run.log")
    if os.path.exists(log):
        with open(log, errors="replace") as f:
            for line in f:
                m = re.search(r"\[openings\].*\bT=([\d.e+-]+) s", re.sub(r"\x1b\[[0-9;]*m", "", line))
                if m:
                    period = float(m.group(1))
                    break
    return period


def integrate(t, y):
    """Trapezoid rule; the steps are not uniform (step control)."""
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def cycle_slice(d, period, k):
    """Rows of cycle k (0-based) by time, t in (k T, (k + 1) T]."""
    t = d["t"]
    return np.where((t > k * period * (1 + 1e-9)) & (t <= (k + 1) * period * (1 + 1e-9)))[0]


def on_phase(d, period, k, name, phase):
    """Column name of cycle k interpolated at the phases in [0, 1)."""
    i = cycle_slice(d, period, k)
    return np.interp(phase, (d["t"][i] - k * period) / period, d[name][i])


class Fluid:
    def __init__(self, tag, folder, last):
        self.tag, self.folder = tag, folder
        self.d = read_csv(os.path.join(folder, "case.csv"))
        self.T = read_period(folder)
        if self.T is None:
            raise RuntimeError("no [openings] line with the period in %s/run.log" % folder)
        self.cycles = int(round(self.d["t"][-1] / self.T))
        if self.cycles < 1:
            raise RuntimeError("%s: less than one cycle" % folder)
        self.last = min(last, self.cycles)
        self.phase = np.linspace(0.0, 1.0, 400, endpoint=False)

    def series(self, name, k=None):
        k = self.cycles - 1 if k is None else k
        return on_phase(self.d, self.T, k, name, self.phase)

    def mean(self, name):
        """Time mean over the last cycles."""
        out = []
        for k in range(self.cycles - self.last, self.cycles):
            i = cycle_slice(self.d, self.T, k)
            out.append(integrate(self.d["t"][i], self.d[name][i]) / (self.d["t"][i][-1] - self.d["t"][i][0]))
        return float(np.mean(out))

    def peak(self, name):
        i = np.concatenate([cycle_slice(self.d, self.T, k) for k in range(self.cycles - self.last, self.cycles)])
        return float(np.max(self.d[name][i]))

    def periodicity(self, name):
        """Relative L2 difference of the last two cycles; nan with one cycle."""
        if self.cycles < 2:
            return float("nan")
        a, b = self.series(name, self.cycles - 2), self.series(name, self.cycles - 1)
        return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))

    def metrics(self):
        d, m = self.d, {}
        k = self.cycles - 1
        i = cycle_slice(d, self.T, k)
        t = d["t"][i]
        m["cycles"] = self.cycles
        m["stroke_volume_mL"] = 1e6 * integrate(t, d["qOut"][i])
        m["peak_qOut_mL_s"] = 1e6 * self.peak("qOut")
        m["peak_speed_m_s"] = self.peak("maxU")
        m["kinetic_mJ"] = 1e3 * self.mean("kinetic")
        m["solvent_dissipation_mW"] = 1e3 * self.mean("viscousDissipation")
        m["polymer_dissipation_mW"] = 1e3 * self.mean("polymerDissipation")
        m["polymer_power_mW"] = 1e3 * self.mean("polymerPower")
        m["elastic_energy_uJ"] = 1e6 * self.mean("elasticEnergy")
        m["oss_convective_mW"] = 1e3 * self.mean("ossConvective")
        physical = self.mean("viscousDissipation") + self.mean("polymerDissipation")
        m["oss_over_physical"] = self.mean("ossConvective") / max(physical, 1e-300)
        m["subscale_kinetic_fraction"] = self.mean("subscaleKinetic") / max(self.mean("kinetic"), 1e-300)
        m["psi_subscale_ratio"] = self.mean("psiSubscale") / max(self.mean("psiNorm"), 1e-300) \
            if self.mean("psiNorm") > 0 else 0.0
        m["wi_max"] = self.peak("wiMax")
        m["wi_volume_fraction_mean"] = self.mean("wiVolumeFraction")
        m["wi_appendage_mean"] = self.mean("wiAppendage")
        m["shear_appendage_1_s"] = self.mean("shearAppendage")
        m["kinetic_appendage_uJ"] = 1e6 * self.mean("kineticAppendage")
        m["max_tawss_Pa"] = float(d["maxTAWSS"][i][-1])
        m["max_tawss_appendage_Pa"] = float(d["maxTAWSSAppendage"][i][-1])
        m["age_appendage_s"] = float(d["ageAppendage"][i][-1])
        # Growth of the appendage mean age over the last cycle, in periods: 1
        # means that no blood left the appendage, 0 a periodic state.
        m["age_appendage_growth"] = float((d["ageAppendage"][i][-1] - d["ageAppendage"][i][0]) /
                                          max(t[-1] - t[0], 1e-300))
        m["n_min"], m["n_max"] = float(np.min(d["nMin"][i])), float(np.max(d["nMax"][i]))
        m["lambda_max_s"] = float(np.max(d["lambdaMax"][i]))
        m["rejected_steps"] = int(np.sum(d["retries"]))
        m["periodicity_K"] = self.periodicity("kinetic")
        m["periodicity_qOut"] = self.periodicity("qOut")
        return m

    def budget(self):
        """Terms of the kinetic energy balance over the last cycle (W, cycle
        means) and the residual. Opening power: -(p u.n) through the openings
        with the mean pressures, an estimate, so the residual carries it."""
        d = self.d
        i = cycle_slice(d, self.T, self.cycles - 1)
        t = d["t"][i]
        span = t[-1] - t[0]
        dK = (d["kinetic"][i][-1] - d["kinetic"][i][0]) / span
        openings = integrate(t, d["pInletMean"][i] * d["qIn"][i] - d["pOutletMean"][i] * d["qOut"][i]) / span
        solvent = integrate(t, d["viscousDissipation"][i]) / span
        polymer = integrate(t, d["polymerPower"][i]) / span
        oss = integrate(t, d["ossConvective"][i] + d["ossGradDiv"][i]) / span
        return {"dK_dt": dK, "openings": openings, "solvent": solvent, "polymer_power": polymer,
                "subscale": oss, "residual": openings - solvent - polymer - oss - dK}


# ---- Probes --------------------------------------------------------------------
def probe_analysis(fluid, out, figures):
    path = os.path.join(fluid.folder, "case_probes.csv")
    if not os.path.exists(path):
        return []
    p = read_csv(path)
    T, t = fluid.T, p["t"]
    n_probes = sum(1 for k in p if k.startswith("ux"))
    first = 1 if fluid.cycles > 2 else 0          # the ramped cycle is left out
    cycles = list(range(first, fluid.cycles))
    if len(cycles) < 2:
        return []
    steps = max(int(round(T / np.median(np.diff(t)))), 16)
    phase = (np.arange(steps) + 0.5) / steps
    rows = []
    for j in range(n_probes):
        u = []
        for k in cycles:
            i = np.where((t > k * T * (1 + 1e-9)) & (t <= (k + 1) * T * (1 + 1e-9)))[0]
            u.append(np.stack([np.interp(phase, (t[i] - k * T) / T, p["%s%d" % (c, j)][i])
                               for c in ("ux", "uy", "uz")], axis=1))
        u = np.array(u)                             # (cycles, steps, 3)
        mean = u.mean(axis=0)
        fluct = u - mean
        k_mean = 0.5 * np.mean(np.sum(mean ** 2, axis=1))
        k_fluct = 0.5 * np.mean(np.sum(fluct ** 2, axis=2)) * len(cycles) / (len(cycles) - 1)
        # PSD of the fluctuation, concatenated over the cycles, one-sided.
        x = fluct.reshape(-1, 3)
        dt = T / steps
        seg = min(len(x), steps)
        window = np.hanning(seg)
        psd, count = 0.0, 0
        for a in range(0, len(x) - seg + 1, seg // 2):
            block = x[a:a + seg] * window[:, None]
            spec = np.abs(np.fft.rfft(block, axis=0)) ** 2
            psd = psd + spec.sum(axis=1) * 2.0 * dt / np.sum(window ** 2)
            count += 1
        psd = psd / max(count, 1)
        f = np.fft.rfftfreq(seg, dt)
        np.savetxt(os.path.join(out, "probe_psd_%s_%d.csv" % (fluid.tag, j)),
                   np.column_stack([f, psd]), delimiter=",", header="f_Hz,psd_m2_s2_per_Hz", comments="")
        rows.append({"fluid": fluid.tag, "probe": j, "cycles": len(cycles),
                     "mean_kinetic_m2_s2": k_mean, "fluctuation_kinetic_m2_s2": k_fluct,
                     "fluctuation_intensity": np.sqrt(k_fluct / max(k_mean, 1e-300))})
    return rows


# ---- Fields ----------------------------------------------------------------------
class Snapshots:
    """The XDMF series of one fluid: mesh merged over the rank pieces."""

    def __init__(self, folder):
        import h5py
        self.h5py = h5py
        self.folder = folder
        tree = ET.parse(os.path.join(folder, "case.xdmf"))
        self.times, self.steps = [], []
        for grid in tree.getroot().iter("Grid"):
            if grid.get("CollectionType") != "Spatial":
                continue
            time = grid.find("Time")
            self.times.append(float(time.get("Value")))
            pieces = []
            for g in grid.findall("Grid"):
                piece = {"geometry": g.find("Geometry/DataItem").text.strip(),
                         "topology": g.find("Topology/DataItem").text.strip()}
                for a in g.findall("Attribute"):
                    piece[a.get("Name")] = a.find("DataItem").text.strip()
                pieces.append(piece)
            self.steps.append(pieces)
        self.times = np.array(self.times)
        self.build_mesh()

    def read(self, ref):
        name, path = ref.split(":", 1)
        with self.h5py.File(os.path.join(self.folder, name), "r") as f:
            return f[path][()]

    def build_mesh(self):
        pieces = self.steps[0]
        coords = [self.read(p["geometry"]).reshape(-1, 3) for p in pieces]
        q = 1.0e-9 * max(np.abs(c).max() for c in coords)
        allc = np.vstack(coords)
        keys = np.round(allc / q).astype(np.int64)
        _, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
        self.X = allc[first]
        inverse = inverse.reshape(-1)
        self.maps, a = [], 0
        for c in coords:
            self.maps.append(inverse[a:a + len(c)])
            a += len(c)
        tets = []
        for p, m in zip(pieces, self.maps):
            topo = self.read(p["topology"]).astype(np.int64)
            if not np.all(topo[::5] == 6):
                raise RuntimeError("non-tetrahedral cells in %s" % self.folder)
            tets.append(m[topo.reshape(-1, 5)[:, 1:]])
        self.tets = np.vstack(tets)
        x = self.X[self.tets]
        e = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0], x[:, 3] - x[:, 0]], axis=2)
        self.vol = np.abs(np.linalg.det(e)) / 6.0
        inv = np.linalg.inv(e)
        self.G = np.zeros((len(self.tets), 4, 3))
        self.G[:, 1:, :] = inv
        self.G[:, 0, :] = -inv.sum(axis=1)
        self.centroid = x.mean(axis=1)
        self.h = np.cbrt(self.vol)

    def field(self, step, name, components):
        out = np.zeros((len(self.X), components))
        for p, m in zip(self.steps[step], self.maps):
            if name not in p:
                return None
            out[m] = self.read(p[name]).reshape(-1, components)
        return out

    def gradient(self, nodal):
        """Cell gradient of a nodal field: (cells, components, 3)."""
        return np.einsum("cia,cix->cax", nodal[self.tets], self.G)


def field_analysis(fluid, out, ball):
    from scipy.interpolate import griddata
    s = Snapshots(fluid.folder)
    T = fluid.T
    steps = [i for i, t in enumerate(s.times) if t > (fluid.cycles - 1) * T * (1 + 1e-9)]
    if not steps:
        return None
    inside = np.ones(len(s.tets), dtype=bool)
    if ball is not None:
        inside = np.linalg.norm(s.centroid - np.array(ball[:3]), axis=1) < ball[3]
    acc = {k: 0.0 for k in ("q_positive_fraction", "shear_mean", "vorticity_rms", "n1_mean", "n2_mean",
                            "n1_max", "n1_appendage", "shear_appendage")}
    last = None
    for step in steps:
        u = s.field(step, "velocity", 3)
        L = s.gradient(u)                                   # L_ax = d u_a / d x_x
        E = 0.5 * (L + np.swapaxes(L, 1, 2))
        W = 0.5 * (L - np.swapaxes(L, 1, 2))
        gd = np.sqrt(2.0 * np.einsum("cij,cij->c", E, E))
        Q = 0.5 * (np.einsum("cij,cij->c", W, W) - np.einsum("cij,cij->c", E, E))
        vort = np.sqrt(2.0 * np.einsum("cij,cij->c", W, W))
        sd, so = s.field(step, "elasticStressDiag", 3), s.field(step, "elasticStressOff", 3)
        n1 = n2 = np.zeros(len(s.tets))
        if sd is not None and np.abs(sd).max() > 0.0:
            sdc, soc = sd[s.tets].mean(axis=1), so[s.tets].mean(axis=1)
            S = np.zeros((len(s.tets), 3, 3))
            S[:, 0, 0], S[:, 1, 1], S[:, 2, 2] = sdc[:, 0], sdc[:, 1], sdc[:, 2]
            S[:, 0, 1] = S[:, 1, 0] = soc[:, 0]
            S[:, 0, 2] = S[:, 2, 0] = soc[:, 1]
            S[:, 1, 2] = S[:, 2, 1] = soc[:, 2]
            uc = u[s.tets].mean(axis=1)
            speed = np.linalg.norm(uc, axis=1)
            e1 = uc / np.maximum(speed, 1e-300)[:, None]
            # grad |u| = L^T u / |u|; its part normal to u is the gradient direction.
            g = np.einsum("cax,ca->cx", L, e1)
            g = g - np.einsum("cx,cx->c", g, e1)[:, None] * e1
            gn = np.linalg.norm(g, axis=1)
            e2 = g / np.maximum(gn, 1e-300)[:, None]
            e3 = np.cross(e1, e2)
            ok = (speed > 1e-6 * speed.max()) & (gn > 1e-6 * gn.max())
            s11 = np.einsum("ci,cij,cj->c", e1, S, e1)
            s22 = np.einsum("ci,cij,cj->c", e2, S, e2)
            s33 = np.einsum("ci,cij,cj->c", e3, S, e3)
            n1 = np.where(ok, s11 - s22, 0.0)
            n2 = np.where(ok, s22 - s33, 0.0)
        V = s.vol.sum()
        acc["q_positive_fraction"] += s.vol[Q > 0].sum() / V
        acc["shear_mean"] += np.sum(s.vol * gd) / V
        acc["vorticity_rms"] += np.sqrt(np.sum(s.vol * vort ** 2) / V)
        acc["n1_mean"] += np.sum(s.vol * n1) / V
        acc["n2_mean"] += np.sum(s.vol * n2) / V
        acc["n1_max"] = max(acc["n1_max"], float(np.max(n1)))
        if ball is not None and inside.any():
            Va = s.vol[inside].sum()
            acc["n1_appendage"] += np.sum(s.vol[inside] * n1[inside]) / Va
            acc["shear_appendage"] += np.sum(s.vol[inside] * gd[inside]) / Va
        last = (u, Q, gd, vort, n1, n2)
    for k in acc:
        if k != "n1_max":
            acc[k] /= len(steps)
    acc["snapshots"] = len(steps)
    acc["h_mean_mm"] = 1e3 * float(np.mean(s.h))
    acc["cells"] = len(s.tets)

    # Cell fields of the last snapshot, legacy VTK.
    u, Q, gd, vort, n1, n2 = last
    with open(os.path.join(out, "%s_fields.vtk" % fluid.tag), "w") as f:
        f.write("# vtk DataFile Version 3.0\n%s\nASCII\nDATASET UNSTRUCTURED_GRID\n" % fluid.tag)
        f.write("POINTS %d double\n" % len(s.X))
        np.savetxt(f, s.X, fmt="%.9g")
        f.write("CELLS %d %d\n" % (len(s.tets), 5 * len(s.tets)))
        np.savetxt(f, np.column_stack([np.full(len(s.tets), 4), s.tets]), fmt="%d")
        f.write("CELL_TYPES %d\n" % len(s.tets))
        np.savetxt(f, np.full(len(s.tets), 10), fmt="%d")
        f.write("CELL_DATA %d\n" % len(s.tets))
        for name, a in (("Q", Q), ("shearRate", gd), ("vorticity", vort), ("N1", n1), ("N2", n2)):
            f.write("SCALARS %s double 1\nLOOKUP_TABLE default\n" % name)
            np.savetxt(f, a, fmt="%.9g")

    # Kinetic energy spectrum in the largest cube centred at the centroid whose
    # corners are inside the fluid (nearest node closer than 2 h).
    from scipy.spatial import cKDTree
    tree = cKDTree(s.X)
    c = np.average(s.centroid, axis=0, weights=s.vol)
    h = float(np.mean(s.h))
    half = 0.0
    for trial in np.linspace(0.5, 0.02, 49) * np.ptp(s.X, axis=0).max():
        corners = c + trial * np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)])
        if np.all(tree.query(corners)[0] < 2.0 * h):
            half = trial
            break
    if half > 4.0 * h:
        n = max(int(2 * half / h), 8)
        axis = (np.arange(n) + 0.5) / n * 2 * half - half
        gx, gy, gz = np.meshgrid(axis, axis, axis, indexing="ij")
        pts = c + np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
        near = np.all(np.abs(s.X - c) < half + 3 * h, axis=1)
        ug = np.stack([np.nan_to_num(griddata(s.X[near], u[near, k], pts, method="linear")).reshape(n, n, n)
                       for k in range(3)])
        ug -= ug.mean(axis=(1, 2, 3), keepdims=True)
        w1 = np.hanning(n)
        window = w1[:, None, None] * w1[None, :, None] * w1[None, None, :]
        uh = np.fft.fftn(ug * window, axes=(1, 2, 3)) / n ** 3
        e3d = 0.5 * np.sum(np.abs(uh) ** 2, axis=0) / np.mean(window ** 2)
        k1 = 2 * np.pi * np.fft.fftfreq(n, 2 * half / n)
        kk = np.sqrt(k1[:, None, None] ** 2 + k1[None, :, None] ** 2 + k1[None, None, :] ** 2)
        dk = 2 * np.pi / (2 * half)
        shells = np.arange(1, n // 2)
        spectrum = np.array([e3d[(kk >= (i - 0.5) * dk) & (kk < (i + 0.5) * dk)].sum() / dk for i in shells])
        np.savetxt(os.path.join(out, "spectrum_%s.csv" % fluid.tag),
                   np.column_stack([shells * dk, spectrum]), delimiter=",",
                   header="k_rad_per_m,E_m3_s2 (cube of side %.4g m, %d^3 points; the mesh resolves k < %.4g rad/m)"
                          % (2 * half, n, np.pi / h), comments="")
        acc["spectrum_decades"] = float(np.log10((np.pi / h) / dk))
    return acc


# ---- Output -------------------------------------------------------------------------
def write_table(path, rows, columns):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(columns)
        for r in rows:
            w.writerow([r.get(c, "") for c in columns])


def tex_table(path, fluids, metrics):
    keys = [("stroke_volume_mL", "Stroke volume (mL)"), ("peak_speed_m_s", r"Peak speed (m\,s$^{-1}$)"),
            ("kinetic_mJ", "Kinetic energy (mJ)"), ("solvent_dissipation_mW", "Solvent dissipation (mW)"),
            ("polymer_dissipation_mW", "Polymer dissipation (mW)"),
            ("oss_over_physical", "Subscale / physical dissipation"),
            ("wi_max", r"$\max\,\mathit{Wi}$"), ("wi_volume_fraction_mean", r"$V(\mathit{Wi}>1)/V$"),
            ("shear_appendage_1_s", r"Appendage shear rate (s$^{-1}$)"),
            ("age_appendage_s", "Appendage age (s)"), ("max_tawss_Pa", "max TAWSS (Pa)"),
            ("periodicity_K", "Cycle-to-cycle difference of $K$")]
    with open(path, "w") as f:
        f.write("\\begin{tabular}{l" + "r" * len(fluids) + "}\n\\hline\n")
        f.write(" & " + " & ".join(fl.tag for fl in fluids) + " \\\\\n\\hline\n")
        for key, label in keys:
            f.write(label + " & " + " & ".join("%.3g" % metrics[fl.tag][key] for fl in fluids) + " \\\\\n")
        f.write("\\hline\n\\end{tabular}\n")


def figures(fluids, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [("qOut", 1e6, "mitral flow rate (mL/s)"), ("kinetic", 1e3, "kinetic energy (mJ)"),
              ("viscousDissipation", 1e3, "solvent dissipation (mW)"), ("polymerPower", 1e3, "polymer power (mW)"),
              ("wiMax", 1.0, "max Wi"), ("ageAppendage", 1.0, "appendage mean age (s)"),
              ("ossConvective", 1e3, "convective subscale dissipation (mW)"), ("nMax", 1.0, "max aggregate size")]
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.5), constrained_layout=True)
    for ax, (name, scale, label) in zip(axes.ravel(), panels):
        for fl in fluids:
            ax.plot(fl.phase, scale * fl.series(name), label=fl.tag, lw=1.2)
        ax.set_xlabel("t/T (last cycle)")
        ax.set_ylabel(label)
    axes[0, 0].legend(frameon=False)
    fig.savefig(os.path.join(out, "last_cycle.png"), dpi=160)
    plt.close(fig)

    psd = sorted(p for p in os.listdir(out) if p.startswith("probe_psd_"))
    probes = sorted({int(p.rsplit("_", 1)[1][:-4]) for p in psd})
    if probes:
        fig, axes = plt.subplots(1, len(probes), figsize=(4.2 * len(probes), 3.8), constrained_layout=True,
                                 squeeze=False)
        for ax, j in zip(axes[0], probes):
            for fl in fluids:
                path = os.path.join(out, "probe_psd_%s_%d.csv" % (fl.tag, j))
                if os.path.exists(path):
                    a = np.loadtxt(path, delimiter=",", skiprows=1)
                    ax.loglog(a[1:, 0], a[1:, 1], label=fl.tag, lw=1.1)
                    ref = a[1:, 0]
            ymax = ax.get_ylim()[1]
            ax.loglog(ref, 0.1 * ymax * (ref / ref[0]) ** (-5.0 / 3.0), "k--", lw=0.7, label="-5/3")
            ax.loglog(ref, 0.1 * ymax * (ref / ref[0]) ** (-3.0), "k:", lw=0.7, label="-3")
            ax.set_xlabel("f (Hz)")
            ax.set_ylabel("PSD of u' (m$^2$ s$^{-2}$ Hz$^{-1}$)")
            ax.set_title("probe %d" % j)
        axes[0, 0].legend(frameon=False, fontsize=8)
        fig.savefig(os.path.join(out, "probe_spectra.png"), dpi=160)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("campaign")
    ap.add_argument("--out")
    ap.add_argument("--cycles", type=int, default=1)
    ap.add_argument("--fields", action="store_true")
    ap.add_argument("--appendage", help="x,y,z,r in metres, for --fields")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()
    out = args.out or os.path.join(args.campaign, "post")
    os.makedirs(out, exist_ok=True)
    ball = [float(v) for v in args.appendage.split(",")] if args.appendage else None

    fluids = []
    for tag in ORDER + sorted(set(os.listdir(args.campaign)) - set(ORDER)):
        folder = os.path.join(args.campaign, tag)
        if os.path.isfile(os.path.join(folder, "case.csv")):
            try:
                fluids.append(Fluid(tag, folder, args.cycles))
            except Exception as error:          # an unfinished case does not stop the others
                print("[skip] %s: %s" % (tag, error))
    if not fluids:
        raise SystemExit("no case.csv under %s/<fluid>/" % args.campaign)

    metrics = {fl.tag: fl.metrics() for fl in fluids}
    columns = ["fluid"] + list(next(iter(metrics.values())).keys())
    write_table(os.path.join(out, "metrics.csv"), [dict(fluid=t, **m) for t, m in metrics.items()], columns)
    tex_table(os.path.join(out, "metrics.tex"), fluids, metrics)

    rows = []
    for a, b, meaning in PAIRS:
        if a in metrics and b in metrics:
            for key in ("stroke_volume_mL", "kinetic_mJ", "solvent_dissipation_mW", "max_tawss_Pa",
                        "shear_appendage_1_s", "age_appendage_s", "kinetic_appendage_uJ"):
                va, vb = metrics[a][key], metrics[b][key]
                floor = max(metrics[a]["periodicity_K"], metrics[b]["periodicity_K"])
                rel = (va - vb) / max(abs(vb), 1e-300)
                rows.append({"pair": "%s-%s" % (a, b), "effect": meaning, "metric": key, a: va, b: vb,
                             "relative_difference": rel, "cycle_to_cycle_floor": floor,
                             "above_floor": bool(abs(rel) > floor) if np.isfinite(floor) else ""})
    write_table(os.path.join(out, "differences.csv"), rows,
                ["pair", "effect", "metric", "relative_difference", "cycle_to_cycle_floor", "above_floor"])
    write_table(os.path.join(out, "budget.csv"), [dict(fluid=fl.tag, **fl.budget()) for fl in fluids],
                ["fluid", "dK_dt", "openings", "solvent", "polymer_power", "subscale", "residual"])

    probe_rows = []
    for fl in fluids:
        probe_rows += probe_analysis(fl, out, not args.no_figures)
    if probe_rows:
        write_table(os.path.join(out, "probes.csv"), probe_rows, list(probe_rows[0].keys()))

    if args.fields:
        field_rows = []
        for fl in fluids:
            try:
                r = field_analysis(fl, out, ball)
                if r:
                    field_rows.append(dict(fluid=fl.tag, **r))
            except Exception as error:
                print("[fields] %s: %s" % (fl.tag, error))
        if field_rows:
            keys = ["fluid"] + sorted({k for r in field_rows for k in r} - {"fluid"})
            write_table(os.path.join(out, "fields.csv"), field_rows, keys)

    if not args.no_figures:
        try:
            figures(fluids, out)
        except Exception as error:
            print("[figures] %s" % error)

    print("fluid   cycles  SV(mL)  K(mJ)  diss_solvent(mW)  diss_polymer(mW)  OSS/physical  maxWi  age_LAA(s)  dK cycle-to-cycle")
    for fl in fluids:
        m = metrics[fl.tag]
        print("%-6s  %5d  %6.2f  %6.3f  %10.4f  %14.4f  %12.2f  %6.2f  %8.3f  %10.3g" % (
            fl.tag, m["cycles"], m["stroke_volume_mL"], m["kinetic_mJ"], m["solvent_dissipation_mW"],
            m["polymer_dissipation_mW"], m["oss_over_physical"], m["wi_max"], m["age_appendage_s"],
            m["periodicity_K"]))
    print("written to", out)


if __name__ == "__main__":
    main()
