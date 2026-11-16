#!/usr/bin/env python3
"""Rewrite, in place, the wedge topology of XDMF meshes written by Rodin
before the VTK wedge-order fix, so that ParaView shows them without holes.

Only /Mesh/XDMF/Topology of every *.mesh.h5 found is touched (node indices of
Wedge (8) and Wedge_18 (41) cells are permuted); geometry and fields are
unchanged. Files already in VTK order are detected and left alone.

Usage: fix_xdmf_wedges.py RESULTS_DIR [--dry-run]
"""
import argparse
import glob
import os

import numpy as np

W18 = np.array([0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10, 17, 16, 15])
W6 = np.array([0, 2, 1, 3, 5, 4])
SIZES = {1: 1, 4: 3, 5: 4, 6: 4, 7: 5, 8: 6, 9: 8, 34: 3, 35: 9, 36: 6, 38: 10, 39: 13,
         41: 18, 50: 27}


def fix(path, dry):
    import h5py
    with h5py.File(path, "r+" if not dry else "r") as f:
        if "/Mesh/XDMF/Topology" not in f:
            return "no topology"
        T = np.asarray(f["/Mesh/XDMF/Topology"][()]).astype(np.int64)
        X = np.asarray(f["/Mesh/Geometry/Vertices"][()]).reshape(-1, 3)
        shape = f["/Mesh/XDMF/Topology"].shape
        flat = T.ravel()
        starts, i = [], 0
        while i < len(flat):
            t = int(flat[i])
            if t not in SIZES:
                return f"unsupported cell id {t}"
            n = SIZES[t] + (1 if t == 2 else 0)
            if t in (8, 41):
                starts.append((i + 1, t))
            i += 1 + SIZES[t]
        if not starts:
            return "no wedges"
        s0, t0 = starts[0]
        c = flat[s0:s0 + SIZES[t0]]
        if t0 == 41:
            old = np.allclose(X[c[9]], 0.5 * (X[c[0]] + X[c[3]]), atol=1e-12 * (1 + np.abs(X).max()))
        else:
            n = np.cross(X[c[1]] - X[c[0]], X[c[2]] - X[c[0]])
            old = np.dot(n, X[c[3]] - X[c[0]]) > 0
        if not old:
            return "already VTK order"
        for s, t in starts:
            perm = W18 if t == 41 else W6
            flat[s:s + len(perm)] = flat[s:s + len(perm)][perm]
        if not dry:
            f["/Mesh/XDMF/Topology"][...] = flat.reshape(shape).astype(
                f["/Mesh/XDMF/Topology"].dtype)
        return f"fixed {len(starts)} wedge(s)"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    files = sorted(glob.glob(os.path.join(args.results, "**", "*.mesh.h5"), recursive=True))
    if not files:
        print(f"No *.mesh.h5 under {args.results}")
    for path in files:
        print(f"  {os.path.relpath(path, args.results)}: {fix(path, args.dry_run)}")


if __name__ == "__main__":
    main()
