#!/usr/bin/env python3
"""Compare a steady ArterialLesion3D pipe run with the sPTT Poiseuille solution.

    python3 compare_sptt_pipe.py <case dir> --eta-s .. --eta-p .. --eps .. --re .. --wi ..
                                 [--rho 1060 --diameter 4e-3 --xa 5 --xb 10]

Reads the last step of <case dir>/case.xdmf (every MPI rank), keeps the
tetrahedra whose centroid lies in xa <= x/D <= xb (developed region) and
integrates, with a 4-point rule on each tetrahedron and the P1 interpolant of
the nodal fields, the relative L2 errors

    e_u     = ||u_h - u|| / ||u||                (all three components)
    e_sigma = ||sigma_h - sigma|| / ||sigma||    (full tensor, Frobenius)

The reference solution is evaluated at the flow rate actually carried in that
region, Q_h = (1/L) int u_x dV, so that the comparison isolates the developed
flow from the inlet-flux error, which is reported separately. Also reported:
the pressure gradient G_h (least-squares slope of p over the region) and the
wall shear stress at the wall nodes of the region (TAWSS of the last cycle).
Prints one line of key=value pairs (parsed by run_sptt_pipe.sh). With
--profile out.csv it also writes, for the nodes of the region, r/R and the
discrete and exact u_x, sigma_xx and sigma_xr (for a profile figure).
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
from sptt_pipe import SPTTPipe  # noqa: E402

TET4 = (np.array([[0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105],
                  [0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105],
                  [0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105],
                  [0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685]]),
        np.full(4, 0.25))


def last_step(case):
    files = glob.glob(os.path.join(case, "case.velocity.r0.*.h5"))
    if not files:
        sys.exit("no case.velocity.r0.*.h5 in %s" % case)
    return max(int(re.search(r"\.(\d+)\.h5$", f).group(1)) for f in files)


def load_rank(case, rank, step):
    def field(name):
        with h5py.File(os.path.join(case, "case.%s.r%d.%06d.h5" % (name, rank, step)), "r") as f:
            return f["GridFunction/Values/Data"][()]
    with h5py.File(os.path.join(case, "case.r%d.mesh.h5" % rank), "r") as f:
        X = f["Mesh/Geometry/Vertices"][()]
        topo = f["Mesh/XDMF/Topology"][()].astype(np.int64)
    if topo.size % 5 or not np.all(topo[0::5] == topo[0]):
        sys.exit("unexpected topology layout (expected only tetrahedra)")
    tets = topo.reshape(-1, 5)[:, 1:]
    return dict(X=X, tets=tets, u=field("velocity"), p=field("pressure"),
                sD=field("elasticStressDiag"), sO=field("elasticStressOff"),
                tawss=field("TAWSS"))


def exact_sigma(sol, Y):
    r = np.hypot(Y[:, 1], Y[:, 2])
    sxr = sol.sigma_xr(r)
    sxx = sol.sigma_xx(r)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(r > 0, Y[:, 1] / r, 0.0)
        s = np.where(r > 0, Y[:, 2] / r, 0.0)
    S = np.zeros((len(Y), 3, 3))
    S[:, 0, 0] = sxx
    S[:, 0, 1] = S[:, 1, 0] = sxr * c
    S[:, 0, 2] = S[:, 2, 0] = sxr * s
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--eta-s", type=float, required=True)
    ap.add_argument("--eta-p", type=float, required=True)
    ap.add_argument("--eps", type=float, required=True)
    ap.add_argument("--re", type=float, required=True)
    ap.add_argument("--wi", type=float, required=True)
    ap.add_argument("--rho", type=float, default=1060.0)
    ap.add_argument("--diameter", type=float, default=4e-3)
    ap.add_argument("--xa", type=float, default=5.0)
    ap.add_argument("--xb", type=float, default=10.0)
    ap.add_argument("--profile", default=None, help="write nodal profiles to this CSV")
    a = ap.parse_args()

    D = a.diameter
    eta0 = a.eta_s + a.eta_p
    Ubar = a.re * eta0 / (a.rho * D)
    lam = a.wi * D / Ubar
    xa, xb = a.xa * D, a.xb * D
    step = last_step(a.case)
    ranks = sorted(int(re.search(r"\.r(\d+)\.mesh\.h5$", f).group(1))
                   for f in glob.glob(os.path.join(a.case, "case.r*.mesh.h5")))

    # Pass 1: flow rate carried in the region, and the pressure slope.
    vol = flux = 0.0
    px, pv = [], []
    data = []
    for rk in ranks:
        d = load_rank(a.case, rk, step)
        P = d["X"][d["tets"]]                           # (n, 4, 3)
        xc = P[:, :, 0].mean(axis=1)
        keep = (xc >= xa) & (xc <= xb)
        d["tets"] = d["tets"][keep]
        P = P[keep]
        J = np.einsum("nij->nij", P[:, 1:, :] - P[:, :1, :])
        d["vol"] = np.abs(np.linalg.det(J)) / 6.0
        ux = d["u"][d["tets"], 0].mean(axis=1)          # exact for P1 mean
        vol += d["vol"].sum()
        flux += (d["vol"] * ux).sum()
        nodes = np.unique(d["tets"])
        px.append(d["X"][nodes, 0])
        pv.append(d["p"][nodes])
        data.append(d)
    L = xb - xa
    Qh = flux / L
    px, pv = np.concatenate(px), np.concatenate(pv)
    Gh = -np.polyfit(px, pv, 1)[0]

    sol = SPTTPipe(0.5 * D, a.eta_s, a.eta_p, lam, a.eps)
    sol.set_flow_rate(Qh)

    # Pass 2: L2 errors with a 4-point rule on the P1 interpolants.
    bary, w = TET4
    eu = nu = es = ns = 0.0
    wall_err, wall_val = [], []
    for d in data:
        T = d["tets"]
        P = d["X"][T]
        Y = np.einsum("qk,nkd->nqd", bary, P).reshape(-1, 3)
        dV = (d["vol"][:, None] * w[None, :]).reshape(-1)
        uh = np.einsum("qk,nkd->nqd", bary, d["u"][T]).reshape(-1, 3)
        r = np.hypot(Y[:, 1], Y[:, 2])
        ue = np.zeros_like(uh)
        ue[:, 0] = sol.u(r)
        eu += (dV * ((uh - ue) ** 2).sum(1)).sum()
        nu += (dV * (ue ** 2).sum(1)).sum()
        sDh = np.einsum("qk,nkd->nqd", bary, d["sD"][T]).reshape(-1, 3)
        sOh = np.einsum("qk,nkd->nqd", bary, d["sO"][T]).reshape(-1, 3)
        Sh = np.zeros((len(Y), 3, 3))
        Sh[:, 0, 0], Sh[:, 1, 1], Sh[:, 2, 2] = sDh[:, 0], sDh[:, 1], sDh[:, 2]
        Sh[:, 0, 1] = Sh[:, 1, 0] = sOh[:, 0]
        Sh[:, 0, 2] = Sh[:, 2, 0] = sOh[:, 1]
        Sh[:, 1, 2] = Sh[:, 2, 1] = sOh[:, 2]
        Se = exact_sigma(sol, Y)
        es += (dV * ((Sh - Se) ** 2).sum((1, 2))).sum()
        ns += (dV * (Se ** 2).sum((1, 2))).sum()
        # wall nodes of the region: r within 1e-6 D of R
        nodes = np.unique(T)
        rn = np.hypot(d["X"][nodes, 1], d["X"][nodes, 2])
        wn = nodes[np.abs(rn - 0.5 * D) < 1e-6 * D]
        wall_val.append(d["tawss"][wn])
    tw = np.concatenate(wall_val)
    if a.profile:
        rows = []
        for d in data:
            n = np.unique(d["tets"])
            Xn = d["X"][n]
            rn = np.hypot(Xn[:, 1], Xn[:, 2])
            with np.errstate(invalid="ignore", divide="ignore"):
                sxr = np.where(rn > 0, (d["sO"][n, 0] * Xn[:, 1] + d["sO"][n, 1] * Xn[:, 2]) / rn, 0.0)
            rc = np.minimum(rn, 0.5 * D)
            rows.append(np.column_stack([Xn[:, 0] / D, rn / (0.5 * D), d["u"][n, 0], sol.u(rc),
                                         d["sD"][n, 0], sol.sigma_xx(rc), sxr, sol.sigma_xr(rc)]))
        np.savetxt(a.profile, np.vstack(rows), delimiter=",", comments="",
                   header="x_over_D,r_over_R,u_h,u_exact,sxx_h,sxx_exact,sxr_h,sxr_exact")
        rr = np.linspace(0.0, 0.5 * D, 201)
        np.savetxt(os.path.splitext(a.profile)[0] + "_exact.csv",
                   np.column_stack([rr / (0.5 * D), sol.u(rr), sol.sigma_xx(rr), sol.sigma_xr(rr)]),
                   delimiter=",", comments="", header="r_over_R,u_exact,sxx_exact,sxr_exact")
    # Wall nodes on MPI partition interfaces carry TAWSS = 0 in the 3D driver
    # output (see the README); they are excluded and counted.
    n_zero = int(np.sum(tw == 0.0))
    tw = tw[tw > 0.0]
    tw_ex = sol.wall_shear_stress

    out = dict(
        step=step, region_volume=vol, Q_region=Qh,
        Q_nominal=Ubar * math.pi * D * D / 4.0,
        G_exact=sol.G, G_h=Gh, err_G=Gh / sol.G - 1.0,
        err_u=math.sqrt(eu / nu), err_sigma=math.sqrt(es / ns),
        tauw_exact=tw_ex, err_tauw_mean=tw.mean() / tw_ex - 1.0,
        err_tauw_max=np.max(np.abs(tw / tw_ex - 1.0)), n_wall=len(tw), n_wall_zero=n_zero,
        lam=lam, Ubar=Ubar)
    out["err_Q_region_vs_nominal"] = Qh / out["Q_nominal"] - 1.0
    print(" ".join("%s=%.6g" % (k, v) for k, v in out.items()))


if __name__ == "__main__":
    main()
