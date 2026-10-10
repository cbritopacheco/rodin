#!/usr/bin/env python3
"""Check a 3D lesion mesh of ArterialLesion3D (MEDIT, units of D).

    python3 check_3d_mesh.py <case> <file.mesh> [<case> <file.mesh> ...]

case: H, S25, S50, S75, A125, A150, A200. Checks: Dimension 3; every boundary
facet of the tetrahedra tagged and every tagged triangle on the boundary; tags
{1, 2, 3, 5}; positive volumes; volume, inlet and outlet areas and lesion wall
area against the exact revolved profile; wall vertices on r(x); tag 5 exactly
on |x| <= ell; element quality (radius ratio).
"""
import sys
import numpy as np

trapezoid = getattr(np, "trapezoid", None) or np.trapz   # numpy < 2.0

CASES = {"H": (0, 0), "S25": (1, .25), "S50": (1, .50), "S75": (1, .75),
         "A125": (2, 1.25), "A150": (2, 1.50), "A200": (2, 2.00)}


def read(path):
    tok = open(path).read().split()
    i, out = 0, {}
    while i < len(tok):
        k = tok[i]
        if k == "Dimension":
            out["dim"] = int(tok[i + 1]); i += 2
        elif k == "Vertices":
            n = int(tok[i + 1]); a = np.array(tok[i + 2:i + 2 + 4 * n], float).reshape(n, 4)
            out["x"] = a[:, :3]; i += 2 + 4 * n
        elif k == "Tetrahedra":
            n = int(tok[i + 1]); a = np.array(tok[i + 2:i + 2 + 5 * n], int).reshape(n, 5)
            out["tet"] = a; i += 2 + 5 * n
        elif k == "Triangles":
            n = int(tok[i + 1]); a = np.array(tok[i + 2:i + 2 + 4 * n], int).reshape(n, 4)
            out["tri"] = a; i += 2 + 4 * n
        else:
            i += 1
    return out


def profile(case):
    lt, p = CASES[case]
    a = 0.0 if lt == 0 else (np.sqrt(1 - p) - 1 if lt == 1 else p - 1)
    ell = 1.5 if lt == 2 else 1.0
    def r(x):
        x = np.asarray(x, float)
        return np.where(np.abs(x) <= ell, 0.5 * (1 + a * (1 + np.cos(np.pi * x / ell)) / 2), 0.5)
    def drdx(x):
        x = np.asarray(x, float)
        return np.where(np.abs(x) <= ell, -0.5 * a * np.pi / (2 * ell) * np.sin(np.pi * x / ell), 0.0)
    return r, drdx, ell


def check(case, path):
    m = read(path)
    x, tet, tri = m["x"], m["tet"][:, :4] - 1, m["tri"]
    ok = True
    def bad(msg):
        nonlocal ok
        ok = False
        print("  FAIL", msg)
    if m.get("dim") != 3:
        bad("Dimension %s" % m.get("dim"))
    r, drdx, ell = profile(case)
    x0, x1 = x[:, 0].min(), x[:, 0].max()

    # volumes
    P = x[tet]
    vol = np.einsum("ij,ij->i", np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), P[:, 3] - P[:, 0]) / 6
    if (vol <= 0).any():
        bad("%d tetrahedra with non-positive volume" % (vol <= 0).sum())
    xs = np.linspace(x0, x1, 200001)
    vex = trapezoid(np.pi * r(xs) ** 2, xs)

    # boundary facets
    faces = np.sort(np.vstack([tet[:, [0, 1, 2]], tet[:, [0, 1, 3]], tet[:, [0, 2, 3]], tet[:, [1, 2, 3]]]), axis=1)
    u, cnt = np.unique(faces, axis=0, return_counts=True)
    bnd = {tuple(f) for f in u[cnt == 1]}
    if (cnt > 2).any():
        bad("non-manifold faces")
    tagged = {}
    for a, b, c, ref in tri:
        tagged[tuple(sorted((a - 1, b - 1, c - 1)))] = ref
    missing = len(bnd - set(tagged))
    extra = len(set(tagged) - bnd)
    if missing or extra:
        bad("boundary facets untagged=%d, tagged-but-interior=%d" % (missing, extra))
    tags = sorted(set(tri[:, 3]))
    if tags != [1, 2, 3, 5]:
        bad("tags %s" % tags)

    def area(ref):
        T = tri[tri[:, 3] == ref][:, :3] - 1
        Q = x[T]
        return 0.5 * np.linalg.norm(np.cross(Q[:, 1] - Q[:, 0], Q[:, 2] - Q[:, 0]), axis=1), Q
    ain, Qin = area(1)
    aout, Qout = area(2)
    ales, Qles = area(5)
    awall, Qw = area(3)
    xl = np.linspace(-ell, ell, 200001)
    alex = trapezoid(2 * np.pi * r(xl) * np.sqrt(1 + drdx(xl) ** 2), xl)
    awex = 2 * np.pi * 0.5 * ((x1 - x0) - 2 * ell)
    if not np.allclose(Qin[..., 0], x0, atol=1e-9):
        bad("inlet triangles not on x = %g" % x0)
    if not np.allclose(Qout[..., 0], x1, atol=1e-9):
        bad("outlet triangles not on x = %g" % x1)
    xc = Qles[..., 0].mean(axis=1)
    if (np.abs(xc) > ell + 1e-9).any():
        bad("tag 5 outside |x| <= ell")
    xc3 = Qw[..., 0].mean(axis=1)
    if (np.abs(xc3) < ell - 1e-9).any():
        bad("tag 3 inside |x| < ell")

    # wall vertices on r(x)
    wv = np.unique(tri[np.isin(tri[:, 3], [3, 5])][:, :3] - 1)
    rr = np.hypot(x[wv, 1], x[wv, 2])
    dev = np.abs(rr - r(x[wv, 0])).max()
    if dev > 2e-3:
        bad("wall vertices off r(x) by %.2e" % dev)

    # radius ratio quality 3 r_in / r_circ, normalised so that 1 = regular
    e = [P[:, j] - P[:, i] for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))]
    fa = lambda a, b, c: 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    S = fa(P[:, 0], P[:, 1], P[:, 2]) + fa(P[:, 0], P[:, 1], P[:, 3]) + fa(P[:, 0], P[:, 2], P[:, 3]) + fa(P[:, 1], P[:, 2], P[:, 3])
    rin = 3 * np.abs(vol) / S
    L = [np.linalg.norm(v, axis=1) for v in e]
    a_, b_, c_ = L[0] * L[5], L[1] * L[4], L[2] * L[3]
    rc = np.sqrt(np.maximum((a_ + b_ + c_) * (a_ + b_ - c_) * (a_ - b_ + c_) * (-a_ + b_ + c_), 0)) / (24 * np.abs(vol))
    q = 3 * rin / rc
    hmax = np.max(np.vstack(L), axis=0)
    thr = np.abs(P[..., 0].mean(axis=1)) < 0.1
    print("%-5s %-28s V=%d T=%d  x=[%g, %g]" % (case, path.split("/")[-1], len(x), len(tet), x0, x1))
    print("      volume %.5f (exact %.5f, %+.2f%%)   inlet %.5f outlet %.5f (pi/4 = %.5f)"
          % (vol.sum(), vex, 100 * (vol.sum() / vex - 1), ain.sum(), aout.sum(), np.pi / 4))
    print("      lesion wall %.4f (exact %.4f, %+.2f%%)   parent wall %.3f (exact %.3f)   max |r - r(x)| %.1e"
          % (ales.sum(), alex, 100 * (ales.sum() / alex - 1), awall.sum(), awex, dev))
    print("      quality min %.3f, 1%% %.3f, median %.3f;  h at throat/centre %.4f (median), inlet faces %d"
          % (q.min(), np.percentile(q, 1), np.median(q), np.median(hmax[thr]), len(Qin)))
    for name, val, ex in (("inlet", ain.sum(), np.pi / 4), ("outlet", aout.sum(), np.pi / 4)):
        if abs(val / ex - 1) > 0.02:
            bad("%s area off by %.1f%%" % (name, 100 * (val / ex - 1)))
    if abs(vol.sum() / vex - 1) > 0.01:
        bad("volume off by %.2f%%" % (100 * (vol.sum() / vex - 1)))
    if q.min() < 0.05:
        bad("sliver tetrahedra (min quality %.3f)" % q.min())
    print("      " + ("OK" if ok else "PROBLEMS"))
    return ok


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 2 or len(args) % 2:
        sys.exit(__doc__)
    res = [check(args[i], args[i + 1]) for i in range(0, len(args), 2)]
    sys.exit(0 if all(res) else 1)
