#!/usr/bin/env python3
"""Compare a pulsatile Oldroyd-B pipe run (3D or axisymmetric driver) with the
exact Womersley solution (womersley_ob.py), and optionally with another run.

    python3 compare_womersley.py <case dir> --eta-s .. --eta-p .. --amplitude ..
            [--rho 1060 --diameter 4e-3 --xa 5 --xb 10] [--ref <case dir>]

The case directory holds case.xdmf, the per-rank HDF5 files and run.log; Ubar,
T, lambda and dt are read from the [inflow] line of run.log. The outputs of the
last cycle, (C-1) T < t <= C T, are compared with the exact solution at their
own times, over the elements whose centroid lies in xa <= x/D <= xb (developed
region), with a 4-point (tetrahedra) or 6-point (triangles, weighted by r)
rule on the P1 interpolants:

    err_u     = ||u_h - u||/||u||                 space-time L2 over the cycle
    err_sigma = ||sigma_h - sigma||/||sigma||     full tensor (sigma_tt included)
    err_sxr, err_sxx                               the two non-zero components
    err_u_max = max over the outputs of the relative L2 error of u at one time
    periodicity = ||u(C T) - u((C-1) T)||/||u(C T)||   (transient left)

With --ref B it also prints the relative L2 differences with run B at t = C T,
which both runs output whatever their dt (dself_u, dself_sigma): run B, on the
same mesh with a smaller dt, gives the self-convergence of the time integrator.
One line of key=value pairs.
"""
import argparse
import glob
import math
import os
import re
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from womersley_ob import WomersleyOB  # noqa: E402

TET4 = (np.array([[0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105],
                  [0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105],
                  [0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105],
                  [0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685]]),
        np.full(4, 0.25))
_a, _b = 0.816847572980459, 0.091576213509771
_c, _d = 0.108103018168070, 0.445948490915965
TRI6 = (np.array([[_a, _b, _b], [_b, _a, _b], [_b, _b, _a],
                  [_c, _d, _d], [_d, _c, _d], [_d, _d, _c]]),
        np.array([0.109951743655322] * 3 + [0.223381589678011] * 3))


def run_params(case):
    log = open(os.path.join(case, "run.log"), errors="replace").read()
    log = re.sub(r"\x1b\[[0-9;]*m", "", log)
    line = re.search(r"\[inflow\][^\n]*", log)
    if not line:
        sys.exit("no [inflow] line in %s/run.log" % case)
    get = lambda key: float(re.search(r"\b%s=([-+0-9.eE]+)" % key, line.group(0)).group(1))
    cyc = re.search(r"\[run\] (\d+) cycles of (\d+) steps", log)
    return dict(Ubar=get("Ubar"), T=get("T"), lam=get("lambda"), dt=get("dt"),
                cycles=int(cyc.group(1)), steps=int(cyc.group(2)))


def output_times(case):
    txt = open(os.path.join(case, "case.xdmf")).read()
    steps = re.findall(r'<Grid Name="Step_(\d+)"[^>]*>\s*<Time Value="([^"]+)"', txt)
    return [(int(k), float(t)) for k, t in steps]


def ranks(case):
    return sorted(int(re.search(r"\.r(\d+)\.mesh\.h5$", f).group(1))
                  for f in glob.glob(os.path.join(case, "case.r*.mesh.h5")))


def read(case, name, rank, step):
    with h5py.File(os.path.join(case, "case.%s.r%d.%06d.h5" % (name, rank, step)), "r") as f:
        return f["GridFunction/Values/Data"][()]


class RankMesh:
    def __init__(self, case, rank, D, xa, xb):
        with h5py.File(os.path.join(case, "case.r%d.mesh.h5" % rank), "r") as f:
            X = f["Mesh/Geometry/Vertices"][()]
            topo = f["Mesh/XDMF/Topology"][()].astype(np.int64)
        kind = topo[0]
        npe = {4: 3, 6: 4}.get(int(kind))
        if npe is None or topo.size % (npe + 1) or not np.all(topo[0::npe + 1] == kind):
            sys.exit("unexpected topology (only triangles or tetrahedra)")
        cells = topo.reshape(-1, npe + 1)[:, 1:]
        self.axi = (npe == 3)
        P = X[cells]
        xc = P[:, :, 0].mean(axis=1)
        keep = (xc >= xa * D) & (xc <= xb * D)
        self.cells = cells[keep]
        P = P[keep]
        if self.axi:
            J = P[:, 1:, :2] - P[:, :1, :2]
            meas = 0.5 * np.abs(J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0])
            bary, w = TRI6
        else:
            meas = np.abs(np.linalg.det(P[:, 1:, :] - P[:, :1, :])) / 6.0
            bary, w = TET4
        self.bary = bary
        self.Y = np.einsum("qk,nkd->nqd", bary, P).reshape(-1, X.shape[1])
        self.r = self.Y[:, 1] if self.axi else np.hypot(self.Y[:, 1], self.Y[:, 2])
        dV = (meas[:, None] * w[None, :]).reshape(-1)
        self.dV = dV * (2.0 * math.pi * self.r if self.axi else 1.0)
        self.X = X
        self.keys = None

    def interp(self, nodal):
        return np.einsum("qk,nk...->nq...", self.bary, nodal[self.cells]).reshape(
            (-1,) + nodal.shape[1:])

    def node_keys(self, D):
        if self.keys is None:
            self.keys = [tuple(k) for k in np.round(self.X / (1e-7 * D)).astype(np.int64)]
        return self.keys


def fields(case, rank, step, axi):
    """Nodal u (3 comps, x first) and sigma as (n, 3, 3) in (x, r/y, theta/z)."""
    u = read(case, "velocity", rank, step)
    n = len(u)
    S = np.zeros((n, 3, 3))
    if axi:
        s = read(case, "elasticStress", rank, step)          # Voigt (xx, xr, rr)
        st = read(case, "elasticStressTheta", rank, step).reshape(-1)
        S[:, 0, 0], S[:, 0, 1], S[:, 1, 0], S[:, 1, 1], S[:, 2, 2] = s[:, 0], s[:, 1], s[:, 1], s[:, 2], st
        U = np.zeros((n, 3))
        U[:, :2] = u
    else:
        d = read(case, "elasticStressDiag", rank, step)
        o = read(case, "elasticStressOff", rank, step)
        S[:, 0, 0], S[:, 1, 1], S[:, 2, 2] = d[:, 0], d[:, 1], d[:, 2]
        S[:, 0, 1] = S[:, 1, 0] = o[:, 0]
        S[:, 0, 2] = S[:, 2, 0] = o[:, 1]
        S[:, 1, 2] = S[:, 2, 1] = o[:, 2]
        U = u
    return U, S


def exact(sol, m, t):
    r = np.minimum(m.r, sol.R)
    U = np.zeros((len(r), 3))
    U[:, 0] = sol.u(r, t)
    sxr, sxx = sol.sigma_xr(r, t), sol.sigma_xx(r, t)
    S = np.zeros((len(r), 3, 3))
    S[:, 0, 0] = sxx
    if m.axi:
        S[:, 0, 1] = S[:, 1, 0] = sxr
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            cy = np.where(m.r > 0, m.Y[:, 1] / m.r, 0.0)
            cz = np.where(m.r > 0, m.Y[:, 2] / m.r, 0.0)
        S[:, 0, 1] = S[:, 1, 0] = sxr * cy
        S[:, 0, 2] = S[:, 2, 0] = sxr * cz
    return U, S, sxr, sxx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--eta-s", type=float, required=True)
    ap.add_argument("--eta-p", type=float, required=True)
    ap.add_argument("--amplitude", type=float, required=True)
    ap.add_argument("--rho", type=float, default=1060.0)
    ap.add_argument("--diameter", type=float, default=4e-3)
    ap.add_argument("--xa", type=float, default=5.0)
    ap.add_argument("--xb", type=float, default=10.0)
    ap.add_argument("--ref", default=None)
    a = ap.parse_args()

    D = a.diameter
    prm = run_params(a.case)
    T, C = prm["T"], prm["cycles"]
    sol = WomersleyOB(0.5 * D, a.rho, a.eta_s, a.eta_p, prm["lam"], 2 * math.pi / T,
                      prm["Ubar"], a.amplitude)
    outs = output_times(a.case)
    tol = 0.25 * prm["dt"]
    last = [(k, t) for k, t in outs if t > (C - 1) * T + tol and t <= C * T + tol]
    end = [k for k, t in outs if abs(t - C * T) < tol]
    prev = [k for k, t in outs if abs(t - (C - 1) * T) < tol]
    if not last or not end:
        sys.exit("no outputs in the last cycle of %s" % a.case)

    meshes = [RankMesh(a.case, rk, D, a.xa, a.xb) for rk in ranks(a.case)]
    acc = dict(eu=0.0, nu=0.0, es=0.0, ns=0.0, exr=0.0, nxr=0.0, exx=0.0, nxx=0.0)
    eu_max = 0.0
    for k, t in last:
        e_t = n_t = 0.0
        for rk, m in enumerate(meshes):
            Un, Sn = fields(a.case, ranks(a.case)[rk], k, m.axi)
            Uh, Sh = m.interp(Un), m.interp(Sn)
            Ue, Se, sxr, sxx = exact(sol, m, t)
            du = (m.dV * ((Uh - Ue) ** 2).sum(1)).sum()
            nu = (m.dV * (Ue ** 2).sum(1)).sum()
            acc["eu"] += du
            acc["nu"] += nu
            e_t += du
            n_t += nu
            acc["es"] += (m.dV * ((Sh - Se) ** 2).sum((1, 2))).sum()
            acc["ns"] += (m.dV * (Se ** 2).sum((1, 2))).sum()
            if m.axi:
                sxr_h = Sh[:, 0, 1]
            else:
                with np.errstate(invalid="ignore", divide="ignore"):
                    sxr_h = np.where(m.r > 0, (Sh[:, 0, 1] * m.Y[:, 1] + Sh[:, 0, 2] * m.Y[:, 2]) / m.r, 0.0)
            acc["exr"] += (m.dV * (sxr_h - sxr) ** 2).sum()
            acc["nxr"] += (m.dV * sxr ** 2).sum()
            acc["exx"] += (m.dV * (Sh[:, 0, 0] - sxx) ** 2).sum()
            acc["nxx"] += (m.dV * sxx ** 2).sum()
        eu_max = max(eu_max, math.sqrt(e_t / n_t))

    def norm_diff(caseA, kA, caseB, kB):
        """||A - B||/||B|| for u and sigma at outputs kA of A and kB of B (B mapped on A's nodes)."""
        num_u = den_u = num_s = den_s = 0.0
        B = {}
        for rk in ranks(caseB):
            with h5py.File(os.path.join(caseB, "case.r%d.mesh.h5" % rk), "r") as f:
                XB = f["Mesh/Geometry/Vertices"][()]
            UB, SB = fields(caseB, rk, kB, meshes[0].axi)
            for key, uu, ss in zip(np.round(XB / (1e-7 * D)).astype(np.int64), UB, SB):
                B[tuple(key)] = (uu, ss)
        for rk, m in enumerate(meshes):
            UA, SA = fields(caseA, ranks(caseA)[rk], kA, m.axi)
            keys = m.node_keys(D)
            try:
                UBn = np.array([B[key][0] for key in keys])
                SBn = np.array([B[key][1] for key in keys])
            except KeyError:
                sys.exit("--ref: the two runs do not share the mesh")
            dU, dS = m.interp(UA - UBn), m.interp(SA - SBn)
            num_u += (m.dV * (dU ** 2).sum(1)).sum()
            den_u += (m.dV * (m.interp(UBn) ** 2).sum(1)).sum()
            num_s += (m.dV * (dS ** 2).sum((1, 2))).sum()
            den_s += (m.dV * (m.interp(SBn) ** 2).sum((1, 2))).sum()
        return math.sqrt(num_u / den_u), math.sqrt(num_s / den_s)

    out = dict(cycles=C, steps=prm["steps"], dt=prm["dt"], T=T, lam=prm["lam"], De=prm["lam"] / T,
               outputs=len(last),
               err_u=math.sqrt(acc["eu"] / acc["nu"]), err_sigma=math.sqrt(acc["es"] / acc["ns"]),
               err_sxr=math.sqrt(acc["exr"] / acc["nxr"]), err_sxx=math.sqrt(acc["exx"] / acc["nxx"]),
               err_u_max=eu_max)
    if prev:
        out["periodicity"] = norm_diff(a.case, end[0], a.case, prev[0])[0]
    if a.ref:
        outsB = output_times(a.ref)
        prmB = run_params(a.ref)
        endB = [k for k, t in outsB if abs(t - prmB["cycles"] * prmB["T"]) < 0.25 * prmB["dt"]]
        if not endB or abs(prmB["cycles"] * prmB["T"] - C * T) > 1e-9 * T:
            sys.exit("--ref: the runs do not end at the same time")
        out["dself_u"], out["dself_sigma"] = norm_diff(a.case, end[0], a.ref, endB[0])
    print(" ".join("%s=%.6g" % (k, v) for k, v in out.items()))


if __name__ == "__main__":
    main()
