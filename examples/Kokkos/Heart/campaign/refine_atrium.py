#!/usr/bin/env python3
"""
refine_atrium.py -- malla adaptada del atrio a partir de Coarse.mesh y de una
corrida previa en ella.

  size   mesh=Coarse.mesh run=c10_TVE target=2200000 out=size.sol [fudge=1]
         Campo de tamaño en los vértices de la malla, para mmg3d:
           h = K * min(h_grad, h_vel, h_pared, h_aberturas, cap), en [hmin, hmax]
         h_grad      = (g / mediana(g))^(-1/2), g = RMS en los dos últimos ciclos
                       de |grad u| por celda: capas de corte y vórtices;
         h_vel       = (U / mediana(U))^(-1/2), U = RMS de |u|: el núcleo de los
                       chorros y donde chocan en el centro, aunque ahí el
                       gradiente sea bajo;
         h_pared     = wall + d_pared / 1 mm (wall = 0.7): capa límite;
         h_aberturas = 0.5 + d_abertura / 3 mm: entrada de los chorros y salida;
         cap (1.2)   = ninguna zona, ni las quietas (orejuela), se engruesa más
                       que 1,2 veces el tamaño típico.
         K se ajusta para que la malla tenga ~target celdas (fudge corrige la
         estimación con el resultado real de mmg3d).
  clean  in=out.mesh out=Fine.mesh   deja Vertices, Tetrahedra y Triangles.
  stats  mesh=Fine.mesh              celdas, volumen, áreas por etiqueta, tamaños.

Unidades de la malla: mm (como Coarse.mesh). Etiquetas: pared 4, venas 1 2 3 5,
mitral 6.
"""
import os
import re
import sys

import numpy as np

WALL = 4
OPENINGS = (1, 2, 3, 5, 6)
args = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
HMIN, HMAX = float(args.get("hmin", 0.1)), float(args.get("hmax", 1.5))   # mm
cmd = sys.argv[1] if len(sys.argv) > 1 else ""


def read_medit(path):
    txt = open(path).read().split("\n")
    out = {}
    i = 0
    while i < len(txt):
        key = txt[i].strip()
        if key in ("Vertices", "Tetrahedra", "Triangles", "Edges", "Ridges", "Corners",
                   "RequiredVertices", "RequiredEdges", "RequiredTriangles"):
            n = int(txt[i + 1])
            out[key] = np.array(" ".join(txt[i + 2:i + 2 + n]).split(), dtype=float)
            out[key] = out[key].reshape(n, -1) if n else out[key].reshape(0, 1)
            i += 2 + n
        else:
            i += 1
    V = out["Vertices"][:, :3]
    T = out["Tetrahedra"][:, :4].astype(np.int64) - 1
    F = out["Triangles"][:, :3].astype(np.int64) - 1
    return dict(V=V, Vref=out["Vertices"][:, 3].astype(int), T=T,
                Tref=out["Tetrahedra"][:, 4].astype(int), F=F,
                Fref=out["Triangles"][:, 3].astype(int))


def geometry(V, T):
    E = V[T[:, 1:]] - V[T[:, [0]]]
    vol = np.abs(np.linalg.det(E)) / 6
    return E, vol


def boundary_faces(T):
    f = np.sort(np.vstack([T[:, [1, 2, 3]], T[:, [0, 2, 3]], T[:, [0, 1, 3]], T[:, [0, 1, 2]]]), axis=1)
    u, c = np.unique(f, axis=0, return_counts=True)
    return u[c == 1]


def stats(m):
    V, T, F, Fr = m["V"], m["T"], m["F"], m["Fref"]
    _, vol = geometry(V, T)
    edges = np.vstack([T[:, [a, b]] for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))])
    L = np.linalg.norm(V[edges[:, 0]] - V[edges[:, 1]], axis=1)
    h = np.cbrt(vol * 6 * np.sqrt(2))
    print(f"  {len(V)} vértices, {len(T)} celdas, volumen {vol.sum() / 1e3:.2f} mL")
    print(f"  arista: min {L.min():.3f}, p5 {np.percentile(L, 5):.3f}, mediana {np.median(L):.3f},"
          f" p95 {np.percentile(L, 95):.3f}, max {L.max():.3f} mm")
    print(f"  celdas con h < 0.1 mm: {(h < 0.1).sum()}, < 0.06 mm: {(h < 0.06).sum()}")
    A = 0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1)
    print("  áreas por etiqueta (mm^2): " + ", ".join(f"{r}: {A[Fr == r].sum():.1f}" for r in np.unique(Fr)))
    bf = boundary_faces(T)
    labelled = {tuple(x) for x in np.sort(F, axis=1)}
    miss = sum(tuple(x) not in labelled for x in bf)
    print(f"  caras de borde sin etiqueta: {miss}")
    return miss


if cmd == "stats":
    stats(read_medit(args["mesh"]))

elif cmd == "clean":
    m = read_medit(args["in"])
    with open(args["out"], "w") as f:
        f.write(f"MeshVersionFormatted 2\n\nDimension 3\n\nVertices\n{len(m['V'])}\n")
        np.savetxt(f, np.column_stack([m["V"], m["Vref"]]), fmt="%.17g %.17g %.17g %d")
        f.write(f"\nTetrahedra\n{len(m['T'])}\n")
        np.savetxt(f, np.column_stack([m["T"] + 1, m["Tref"]]), fmt="%d")
        f.write(f"\nTriangles\n{len(m['F'])}\n")
        np.savetxt(f, np.column_stack([m["F"] + 1, m["Fref"]]), fmt="%d")
        f.write("\nEnd\n")
    if stats(read_medit(args["out"])):
        sys.exit("hay caras de borde sin etiqueta")

elif cmd == "size":
    import h5py
    from scipy.spatial import cKDTree
    m = read_medit(args["mesh"])
    V, T = m["V"], m["T"]
    run = args["run"]
    target = float(args.get("target", 2.2e6))
    fudge = float(args.get("fudge", 1.0))

    _, vol = geometry(V, T)
    if run == "none":
        g = np.ones(len(V))
        U = np.ones(len(V))
        print("  sin corrida: solo pared y aberturas")
    else:
        # Velocidad de la corrida previa en los dos últimos ciclos.
        X = h5py.File(run + ".mesh.h5", "r")["Mesh/Geometry/Vertices"][:]
        TX = h5py.File(run + ".mesh.h5", "r")["Mesh/XDMF/Topology"][:].reshape(-1, 5)[:, 1:].astype(np.int64)
        times = [float(t) for t in re.findall(r'Time Value="([^"]+)"', open(run + ".xdmf").read())]
        cyc = np.loadtxt(run + ".csv", delimiter=",", skiprows=1, usecols=(0, 1))
        last = int(cyc[:, 1].max())
        first = max(1, last - 1) if last > 1 else 1
        t0 = cyc[cyc[:, 1] == first, 0].min() - 1e-9 if last > 1 else 0.0
        ks = [k for k, t in enumerate(times) if t > t0 and os.path.exists(f"{run}.velocity.{k:06d}.h5")]
        if not ks:
            sys.exit("no hay instantes de velocidad en la corrida")
        print(f"  {run}: ciclos {first}-{last}, {len(ks)} instantes")
        EX = X[TX[:, 1:]] - X[TX[:, [0]]]
        Einv = np.linalg.inv(EX)
        G = np.concatenate([-Einv.sum(axis=2, keepdims=True), Einv], axis=2)
        volX = np.abs(np.linalg.det(EX)) / 6
        g2 = np.zeros(len(TX))
        u2 = np.zeros(len(X))
        for k in ks:
            u = h5py.File(f"{run}.velocity.{k:06d}.h5", "r")["GridFunction/Values/Data"][:]
            J = np.einsum("cin,cnk->cki", G, u[TX])
            g2 += np.einsum("cij,cij->c", J, J) / len(ks)
            u2 += np.einsum("vi,vi->v", u, u) / len(ks)
        num = np.zeros(len(X)); den = np.zeros(len(X))
        np.add.at(num, TX.ravel(), np.repeat(volX * g2, 4))
        np.add.at(den, TX.ravel(), np.repeat(volX, 4))
        gX = np.sqrt(num / den)                               # 1/s, en los vértices de la corrida
        # Vértices de la corrida -> vértices de la malla (la corrida está en m).
        s = np.ptp(V, axis=0).max() / np.ptp(X, axis=0).max()
        dist, idx = cKDTree(X * s).query(V)
        if dist.max() > 1e-6 * np.ptp(V):
            sys.exit(f"la corrida {run} no es de esta malla (distancia {dist.max():.3g})")
        g = gX[idx]
        U = np.sqrt(u2)[idx]
    lump = np.zeros(len(V))
    np.add.at(lump, T.ravel(), np.repeat(vol / 4, 4))
    o = np.argsort(g)
    gm = g[o][np.searchsorted(np.cumsum(lump[o]), 0.5 * lump.sum())]
    hg = (np.maximum(g, 1e-3 * gm) / gm) ** -0.5
    o2 = np.argsort(U)
    Um = U[o2][np.searchsorted(np.cumsum(lump[o2]), 0.5 * lump.sum())]
    hu = (np.maximum(U, 1e-3 * Um) / Um) ** -0.5
    F, Fr = m["F"], m["Fref"]
    wall = np.unique(F[Fr == WALL])
    opening = np.unique(F[np.isin(Fr, OPENINGS)])
    dw = cKDTree(V[wall]).query(V)[0]
    do = cKDTree(V[opening]).query(V)[0]
    WALLREL, CAP = float(args.get("wall", 0.7)), float(args.get("cap", 1.2))
    rel = np.minimum.reduce([hg, hu, WALLREL + dw / 1.0, 0.5 + do / 3.0, np.full(len(V), CAP)])

    def count(K):
        h = np.clip(K * rel, HMIN, HMAX)
        hc = h[T].mean(axis=1)
        return fudge * np.sum(vol * 6 * np.sqrt(2) / hc ** 3), h

    lo, hi = 1e-3, 1e2
    for _ in range(60):
        K = np.sqrt(lo * hi)
        n, h = count(K)
        lo, hi = (K, hi) if n > target else (lo, K)
    n, h = count(K)
    print(f"  mediana de |grad u| {gm:.1f} 1/s; K = {K:.4f} mm; celdas estimadas {n:.0f}")
    print(f"  h en la pared: mediana {np.median(h[wall]):.3f} mm; h en el volumen: "
          f"p5 {np.percentile(h, 5):.3f}, mediana {np.median(h):.3f}, p95 {np.percentile(h, 95):.3f} mm")
    for a, b in ((0, 0.5), (0.5, 2), (2, 5), (5, 1e9)):
        k = (dw >= a) & (dw < b)
        if k.any():
            print(f"    a {a:g}-{b:g} mm de la pared: h mediana {np.median(h[k]):.3f} mm "
                  f"({100 * lump[k].sum() / lump.sum():.0f}% del volumen)")
    if run != "none":
        k = U >= np.percentile(U, 90)
        print(f"    10% de vértices más rápidos (chorros): h mediana {np.median(h[k]):.3f} mm")
    with open(args.get("out", "size.sol"), "w") as f:
        f.write(f"MeshVersionFormatted 2\n\nDimension 3\n\nSolAtVertices\n{len(V)}\n1 1\n\n")
        np.savetxt(f, h, fmt="%.6g")
        f.write("\nEnd\n")
    open(args.get("out", "size.sol") + ".K", "w").write(f"{K} {n}\n")

else:
    sys.exit(__doc__)
