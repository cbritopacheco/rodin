#!/usr/bin/env python3
"""Where does the stress blow up? Reads the XDMF output of an ArterialLesion3D
case (all ranks) and prints, for each output step, the largest |sigma| and
|u|, the extreme eigenvalues of psi, and where they are: x/D from the lesion
centre, r/D, the distance to the wall in units of D, and the region.

Usage:
    python3 locate_blowup.py <case dir> [--geo S75|A200|H|...] [--D 0.004]
                             [--last N] [--top K] [--sigma 1000]

    <case dir>   the folder with case.xdmf (e.g. .../3D/Re300_Wo4_waveform_cca/S75/R)
    --last N     only the last N output steps (default: all)
    --top K      the K nodes of largest |sigma| at each step (default 5)
    --sigma S    also count the nodes with |sigma| > S Pa (default 1000)

|sigma| is the Frobenius norm of the elastic stress; psi = log c. Needs numpy
and h5py.
"""
import argparse
import os
import re
import xml.etree.ElementTree as ET

import h5py
import numpy as np

GEOS = {"H": (0, 0.0), "S25": (1, 0.25), "S50": (1, 0.50), "S75": (1, 0.75),
        "A125": (2, 1.25), "A150": (2, 1.50), "A200": (2, 2.00)}


def wall_radius(x, geo):
    """r_wall(x)/D of lesion_3D.geo (x in units of D)."""
    ltype, p = GEOS[geo]
    if ltype == 0:
        return np.full_like(x, 0.5)
    a = np.sqrt(1.0 - p) - 1.0 if ltype == 1 else p - 1.0
    ell = 1.0 if ltype == 1 else 1.5
    return np.where(np.abs(x) <= ell, 0.5 * (1.0 + a * (1.0 + np.cos(np.pi * x / ell)) / 2.0), 0.5)


def region(x, geo):
    ell = 1.0 if GEOS[geo][0] != 2 else 1.5
    if x < -ell:
        return "upstream"
    if x <= ell:
        return "lesion" if x < 0.0 else "lesion (distal half)"
    if x <= ell + 8.0:
        return "post-lesion (jet / wake)"
    if x >= 20.0:
        return "near outlet"
    return "downstream"


def read(path, root):
    with h5py.File(os.path.join(root, path.split(":")[0]), "r") as f:
        return f[path.split(":")[1]][()]


def steps_of(xdmf):
    """[(time, {rank: {name: h5 path}})] from case.xdmf."""
    tree = ET.parse(xdmf)
    out = []
    for step in tree.getroot().iter("Grid"):
        if step.get("CollectionType") != "Spatial":
            continue
        t = step.find("Time")
        time = float(t.get("Value")) if t is not None else float("nan")
        ranks = {}
        for g in step.findall("Grid"):
            m = re.search(r"_r(\d+)$", g.get("Name", ""))
            r = int(m.group(1)) if m else len(ranks)
            d = {"geometry": g.find("Geometry/DataItem").text.strip()}
            for a in g.findall("Attribute"):
                d[a.get("Name")] = a.find("DataItem").text.strip()
            ranks[r] = d
        out.append((time, ranks))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--geo", default=None)
    ap.add_argument("--D", type=float, default=0.004)
    ap.add_argument("--last", type=int, default=0)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--sigma", type=float, default=1000.0)
    args = ap.parse_args()
    geo = args.geo
    if geo is None:
        for g in sorted(GEOS, key=len, reverse=True):
            if "/%s/" % g in os.path.abspath(args.case) + "/":
                geo = g
                break
    geo = geo or "S75"
    steps = steps_of(os.path.join(args.case, "case.xdmf"))
    if args.last:
        steps = steps[-args.last:]
    print("case %s  geometry %s  D = %g m  %d output steps" % (args.case, geo, args.D, len(steps)))
    for time, ranks in steps:
        X, S, P, U = [], [], [], []
        for r in sorted(ranks):
            d = ranks[r]
            need = ("elasticStressDiag", "elasticStressOff", "logConformationDiag", "logConformationOff", "velocity")
            if not all(k in d for k in need):
                continue
            X.append(read(d["geometry"], args.case))
            S.append(np.hstack([read(d["elasticStressDiag"], args.case), read(d["elasticStressOff"], args.case)]))
            P.append(np.hstack([read(d["logConformationDiag"], args.case), read(d["logConformationOff"], args.case)]))
            U.append(read(d["velocity"], args.case))
        if not X:
            continue
        X = np.vstack(X) / args.D
        S = np.vstack(S)
        P = np.vstack(P)
        U = np.vstack(U)
        sig = np.sqrt((S[:, :3] ** 2).sum(1) + 2.0 * (S[:, 3:] ** 2).sum(1))
        T = np.zeros((len(P), 3, 3))
        T[:, 0, 0], T[:, 1, 1], T[:, 2, 2] = P[:, 0], P[:, 1], P[:, 2]
        T[:, 0, 1] = T[:, 1, 0] = P[:, 3]
        T[:, 0, 2] = T[:, 2, 0] = P[:, 4]
        T[:, 1, 2] = T[:, 2, 1] = P[:, 5]
        ev = np.linalg.eigvalsh(T)
        speed = np.linalg.norm(U, axis=1)
        x = X[:, 0]
        r = np.hypot(X[:, 1], X[:, 2])
        dwall = wall_radius(x, geo) - r
        theta = np.degrees(np.arctan2(X[:, 2], X[:, 1]))

        def where(i):
            return "x/D=%6.2f r/D=%.3f wall dist/D=%.3f theta=%4.0f  %s" % (
                x[i], r[i], dwall[i], theta[i], region(x[i], geo))

        i = int(np.argmax(sig))
        j = int(np.argmax(speed))
        big = sig > args.sigma
        print("\n== t = %.6f s   max|sigma| = %.4g Pa at %s" % (time, sig[i], where(i)))
        print("   max|u| = %.4g m/s at %s" % (speed[j], where(j)))
        print("   psi eigenvalues: max %.3f at %s" % (ev[:, 2].max(), where(int(np.argmax(ev[:, 2])))))
        print("                    min %.3f at %s" % (ev[:, 0].min(), where(int(np.argmin(ev[:, 0])))))
        if big.any():
            print("   %d nodes with |sigma| > %g Pa: x/D in [%.2f, %.2f], wall dist/D in [%.3f, %.3f]"
                  % (big.sum(), args.sigma, x[big].min(), x[big].max(), dwall[big].min(), dwall[big].max()))
        for k in np.argsort(sig)[::-1][:args.top]:
            print("     |sigma| %.4g Pa  |u| %.3g  psi [%.2f %.2f %.2f]  %s"
                  % (sig[k], speed[k], ev[k, 0], ev[k, 1], ev[k, 2], where(int(k))))


if __name__ == "__main__":
    main()
