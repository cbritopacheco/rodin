#!/usr/bin/env python3
"""
Builds the stenosis-free FSI pair in the frame of the 80% stenosis meshes.

The stenotic pair (Stenosis80 fluid, Prismatic80 solid) was obtained from the
healthy FSIFluid/FSISolid pair by a similarity transform (scale 25/66) and a
local inward deformation of the lumen.  This script

  1. fits that similarity on the FSI vertices far from the lesion (residual
     ~1e-8 mm), and maps the healthy fluid mesh with it;
  2. takes the stenotic solid (same wedge connectivity, same labels) and
     replaces only its lesion vertices: inner surface = mapped healthy lumen,
     outer surface = inner + the stenotic wall thickness along the stenotic
     extrusion direction.

Away from the lesion both pairs coincide; both fluid/solid pairs are
conforming.  Usage:

  make_nostenosis_mesh.py FSIFluid.mesh FSISolid.mesh Prismatic80.mesh OUTDIR
"""
import sys

import numpy as np
from scipy.spatial import cKDTree

SECTIONS = ('Vertices', 'Tetrahedra', 'Triangles', 'Wedges', 'Edges', 'Quadrilaterals')
LESION_CENTER = np.array([-33.0, 17.1, -153.4])  # mm, stenotic frame
LESION_RADIUS = 10.0  # mm, excluded from the fit


def parse(path):
    lines = open(path).read().split('\n')
    sections, order, i = {}, [], 0
    while i < len(lines):
        key = lines[i].strip()
        if key in SECTIONS:
            n = int(lines[i + 1].split()[0])
            sections[key] = np.array([l.split() for l in lines[i + 2:i + 2 + n]], dtype=float)
            order.append(key)
            i += 2 + n
        else:
            i += 1
    return sections, order


def write(path, sections, order):
    with open(path, 'w') as f:
        f.write('MeshVersionFormatted 2\n\nDimension\n3\n\n')
        for key in order:
            rows = sections[key]
            f.write(f'{key}\n{len(rows)}\n')
            for r in rows:
                if key == 'Vertices':
                    f.write(f'{r[0]:.17g} {r[1]:.17g} {r[2]:.17g} {int(r[3])}\n')
                else:
                    f.write(' '.join(str(int(x)) for x in r) + '\n')
            f.write('\n')
        f.write('End\n')


def similarity(A, B):
    """Umeyama: B ~ s R A + t."""
    ca, cb = A.mean(0), B.mean(0)
    Ac, Bc = A - ca, B - cb
    U, S, Vt = np.linalg.svd(Ac.T @ Bc / len(A))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    s = np.trace(np.diag(S) @ D) / (Ac ** 2).sum(1).mean()
    return s, R, cb - s * R @ ca


def wedge_volume(X, W):
    a, b, c, d = (X[W[:, k]] for k in range(4))
    return np.einsum('ij,ij->i', np.cross(b - a, c - a), d - a)


def tet_volume(X, T):
    a, b, c, d = (X[T[:, k]] for k in range(4))
    return np.einsum('ij,ij->i', np.cross(b - a, c - a), d - a) / 6


def main(fluid_path, solid_path, stenotic_solid_path, outdir):
    F, forder = parse(fluid_path)
    S, _ = parse(solid_path)
    P, porder = parse(stenotic_solid_path)
    XS, XP = S['Vertices'][:, :3], P['Vertices'][:, :3]
    W = P['Wedges'][:, :6].astype(int) - 1
    assert np.array_equal(W, S['Wedges'][:, :6].astype(int) - 1), 'solid connectivity differs'

    tri = P['Triangles']
    inner = np.unique(tri[tri[:, 3] == 1][:, :3].astype(int) - 1)
    bottom, top = W[:, :3].ravel(), W[:, 3:].ravel()
    _, first = np.unique(bottom, return_index=True)
    ia, ob = bottom[first], top[first]
    assert set(ia) == set(inner), 'wedge bottoms are not the FSI surface'

    far = inner[np.linalg.norm(XP[inner] - LESION_CENTER, axis=1) > LESION_RADIUS]
    s, R, t = similarity(XS[far], XP[far])
    T = lambda X: s * X @ R.T + t
    res = np.linalg.norm(T(XS[far]) - XP[far], axis=1).max()
    print(f'similarity scale {s:.12f}, max residual away from the lesion {res:.2e} mm')
    assert res < 1e-5

    thickness = np.linalg.norm(XP[ob] - XP[ia], axis=1)
    direction = (XP[ob] - XP[ia]) / thickness[:, None]
    XH = XP.copy()
    XH[ia] = T(XS[ia])
    XH[ob] = XH[ia] + thickness[:, None] * direction
    vol = wedge_volume(XH, W)
    assert (vol > 0).all(), 'inverted wedge'
    moved = (np.linalg.norm(XH - XP, axis=1) > 1e-6).sum()
    print(f'solid: {moved} vertices moved, min wedge volume ratio '
          f'{(vol / wedge_volume(XP, W)).min():.3f}')

    XF = T(F['Vertices'][:, :3])
    tets = F['Tetrahedra'][:, :4].astype(int) - 1
    assert (tet_volume(XF, tets) > 0).all(), 'inverted tetrahedron'
    ftri = F['Triangles']
    fsi = np.unique(ftri[ftri[:, 3] == 2][:, :3].astype(int) - 1)
    gap = cKDTree(XH[inner]).query(XF[fsi])[0].max()
    print(f'fluid: lumen volume {tet_volume(XF, tets).sum():.1f} mm^3, '
          f'interface gap {gap:.1e} mm')
    assert gap < 1e-9

    F['Vertices'][:, :3] = XF
    P['Vertices'][:, :3] = XH
    write(f'{outdir}/NoStenosisFluid.mesh', F, [k for k in forder if k != 'Edges'])
    write(f'{outdir}/NoStenosisSolid.mesh', P, porder)


if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
