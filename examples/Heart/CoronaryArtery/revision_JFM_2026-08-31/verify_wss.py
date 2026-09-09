#!/usr/bin/env python3
"""
Independent check of the wall shear stress written by CoupledLV0DCoronary3D.

Recomputes, straight from the mesh and the velocity field, the element-wise P1
gradient, the shear rate, the Carreau-Yasuda viscosity and tau_w = mu(g) g, and
compares them against the solver's own `viscosity` and `shearStress` outputs.

Before the fix of computeWallShear() this reported tau_w inflated by
mu_0 / mu_inf ~ 65, because the viscosity in the boundary integral was rebuilt
from Jacobian(uSol), which collapses on the no-slip trace and returned mu_0.

Usage:
    python3 verify_wss.py <run-dir> <step> [nranks]
    python3 verify_wss.py hyp2 002500 7
"""
import sys
import glob
import numpy as np
import h5py

MU0, MUINF, LAMBDA, NIDX, AIDX = 0.346, 0.0053, 17.41, 0.22, 0.69


def mu_cy(g):
    return MUINF + (MU0 - MUINF) * (1.0 + (LAMBDA * g) ** AIDX) ** ((NIDX - 1.0) / AIDX)


def load(path, key="GridFunction/Values/Data"):
    with h5py.File(path, "r") as h:
        return h[key][:]


def main(run, step, nranks):
    g_all, tau_all, mu_all, solver_mu, solver_ss = [], [], [], [], []
    for r in range(nranks):
        with h5py.File(f"{run}/CoronaryArtery.r{r}.mesh.h5", "r") as h:
            X = h["Mesh/Geometry/Vertices"][:]
            T = h["Mesh/XDMF/Topology"][:].reshape(-1, 5)[:, 1:]
        U = load(f"{run}/CoronaryArtery.velocity.r{r}.{step}.h5")
        E = X[T][:, 1:, :] - X[T][:, 0:1, :]
        dU = U[T][:, 1:, :] - U[T][:, 0:1, :]
        G = np.einsum("nkj,nki->nij", np.linalg.inv(E), dU)   # G_ij = du_i/dx_j
        D = 0.5 * (G + np.transpose(G, (0, 2, 1)))
        g = np.sqrt(np.maximum(2.0 * np.einsum("nij,nij->n", D, D), 0.0))
        g_all.append(g)
        mu_all.append(mu_cy(g))
        tau_all.append(mu_cy(g) * g)
        solver_mu.append(load(f"{run}/CoronaryArtery.viscosity.r{r}.{step}.h5").ravel())
        solver_ss.append(
            np.linalg.norm(load(f"{run}/CoronaryArtery.shearStress.r{r}.{step}.h5"), axis=1))

    def row(name, x, unit):
        p = [np.percentile(x, q) for q in (50, 90, 99)]
        print(f"  {name:34s} p50={p[0]:10.4g} p90={p[1]:10.4g} p99={p[2]:10.4g} "
              f"max={x.max():10.4g}  [{unit}]")

    print(f"\n{run} @ step {step}\n")
    print("recomputed from mesh + velocity:")
    row("shear rate gamma_dot", np.concatenate(g_all), "1/s")
    row("viscosity mu(gamma)", np.concatenate(mu_all) * 1e3, "mPa.s")
    row("tau = mu(gamma) gamma", np.concatenate(tau_all), "Pa")
    print("\nwritten by the solver:")
    row("viscosity field", np.concatenate(solver_mu) * 1e3, "mPa.s")
    row("shearStress field", np.concatenate(solver_ss), "Pa")

    ratio = np.percentile(np.concatenate(solver_ss), 50) / \
        np.percentile(np.concatenate(tau_all), 50)
    print(f"\n  median(shearStress) / median(mu*gamma) = {ratio:.1f}")
    print(f"  mu_0 / mu_inf                          = {MU0 / MUINF:.1f}")
    print("  -> a ratio near mu_0/mu_inf means the boundary viscosity "
          "collapsed to mu_0 again.\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 7)
