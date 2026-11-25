#!/usr/bin/env python3
"""Volume rendering of the Carreau--Yasuda straight-tube validation (RmuAxi).

Reads the last snapshot of the axisymmetric RmuAxi fields (XDMF + HDF5) of one
or more steady runs, revolves the (x, r) solution about the axis and renders a
cut-away segment of the developed tube: a ray-cast volume of u/U (U the mean
velocity), a glass wall, and on the cut plane the computed velocity profile
(tube), the exact Carreau--Yasuda profile at the same pressure gradient (beads)
and the Newtonian (Poiseuille) profile with the same flow (dashed).

Fields (pc9, next to rmu_axi_validation/, about 15 min with the three runs in parallel):

  EXE=~/rodin-cmbpe2/build/examples/Heart/CoronaryArtery/RmuAxi/RmuAxi
  CY="-rmu_cy_mu0 0.0712616 -rmu_cy_muinf 0.0041 -rmu_cy_lambda 9.2780 -rmu_cy_n 0.3572 -rmu_cy_a 1.8773"
  C="-rmu_full_tubes 1:10,1:30 -rmu_reduced_tubes 1:20 -rmu_truncation 20 -rmu_scale 1e-4 -rmu_period 0.85 $CY"
  for g in 30 100 400; do d=rmu_axi_fields/g$g; mkdir -p $d
    $EXE -rmu_mesh rmu_axi_validation/mesh/tube32_full.mesh -rmu_output_prefix $d -rmu_closure none \
      -rmu_gamma_ref $g -rmu_amplitude 0 -rmu_steps_per_cycle 50 -rmu_cycles 2 -rmu_output_every 100 \
      $C > $d/log.txt 2>&1 &
  done; wait

Rendering (Mac):

  pip install pyvista scipy h5py
  python3 render_cy_tube.py rmu_axi_fields/g30 rmu_axi_fields/g100 rmu_axi_fields/g400 --interactive
  python3 render_cy_tube.py rmu_axi_fields/g* --scale 3 --out cy_tube.png
  python3 render_cy_tube.py rmu_axi_fields/g* --field viscosity       # mu(gamma_dot) instead of u/U

Each folder holds the RmuAxi.xdmf written with -rmu_output_every; the label is
taken from the folder name (g30 -> gamma_0 = 30 1/s).
"""
import argparse
import glob
import os
import re
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pyvista as pv

# Healthy blood (Table 1)
CY = dict(mu0=0.0712616, muinf=0.0041, lam=9.2780, n=0.3572, a=1.8773)

# XDMF mixed-topology codes -> number of nodes (None: count follows the code)
XDMF_NODES = {1: None, 2: None, 4: 3, 5: 4, 6: 4, 34: 3, 35: 9, 36: 6, 37: 8, 38: 10}


def mu_cy(g, p=CY):
    return p["muinf"] + (p["mu0"] - p["muinf"]) * (1 + (p["lam"] * g) ** p["a"]) ** ((p["n"] - 1) / p["a"])


# ----------------------------------------------------------------------------- reading
def _dataset(item, folder):
    fname, path = item.text.strip().split(":", 1)
    with h5py.File(os.path.join(folder, fname), "r") as f:
        return np.asarray(f[path])


def _uniform_grids(node):
    if node.get("GridType", "Uniform") == "Uniform":
        return [node]
    out = []
    for g in node.findall("Grid"):
        out += _uniform_grids(g)
    return out


def read_last(folder):
    """Last snapshot of RmuAxi.xdmf as a 2-D mesh in the (x, r, 0) plane."""
    xdmf = sorted(glob.glob(os.path.join(folder, "*.xdmf")))
    if not xdmf:
        raise SystemExit(f"no .xdmf in {folder}: rerun RmuAxi with -rmu_output_every")
    root = ET.parse(xdmf[0]).getroot()
    temporal = [g for g in root.iter("Grid") if g.get("CollectionType") == "Temporal"][0]
    last = temporal.findall("Grid")[-1]
    parts = []
    for g in _uniform_grids(last):
        pts = _dataset(g.find("Geometry/DataItem"), folder).reshape(-1, 2)
        top = _dataset(g.find("Topology/DataItem"), folder)
        ttype = g.find("Topology").get("TopologyType", "Mixed")
        if top.ndim == 2:                                   # uniform
            tri = top[:, :3].astype(np.int64)
        else:                                               # mixed stream
            tri, i, top = [], 0, top.astype(np.int64)
            while i < len(top):
                code = int(top[i]); i += 1
                k = XDMF_NODES.get(code)
                if k is None:
                    k = int(top[i]); i += 1
                if code in (4, 36):
                    tri.append(top[i:i + 3])
                i += k
            tri = np.array(tri)
        del ttype
        cells = np.hstack([np.full((len(tri), 1), 3), tri]).ravel()
        m = pv.UnstructuredGrid(cells, np.full(len(tri), pv.CellType.TRIANGLE),
                                np.c_[pts, np.zeros(len(pts))])
        for a in g.findall("Attribute"):
            v = _dataset(a.find("DataItem"), folder)
            if a.get("Center", "Node") == "Node":
                m.point_data[a.get("Name")] = v
            else:
                m.cell_data[a.get("Name")] = v
        parts.append(m)
    m = parts[0] if len(parts) == 1 else pv.merge(parts)
    return m.clean(tolerance=1e-12 * m.length)


# ----------------------------------------------------------------------------- physics
def developed(m, R, L):
    """Mean velocity, pressure gradient and computed profile at mid-length."""
    xc = 0.5 * (m.bounds[0] + m.bounds[1])
    r = np.linspace(0, R, 201)
    line = pv.PolyData(np.c_[np.full_like(r, xc), r, np.zeros_like(r)]).sample(m)
    u = line.point_data["velocity"][:, 0]
    Q = np.trapezoid(2 * np.pi * r * u, r)
    U = Q / (np.pi * R ** 2)
    dx = 0.1 * L
    pp = pv.PolyData(np.array([[xc - dx, 0.3 * R, 0], [xc + dx, 0.3 * R, 0]])).sample(m)
    G = -(pp.point_data["pressure"][1] - pp.point_data["pressure"][0]) / (2 * dx)
    return xc, r, u, U, G


def exact_profile(G, R, r):
    """Exact developed CY profile at pressure gradient G: tau = G r / 2."""
    from scipy.optimize import brentq
    tw = 0.5 * G * R
    s = np.linspace(0, R, 2001)
    gd = np.array([0.0 if t == 0 else brentq(lambda g: mu_cy(g) * g - t, 0, t / CY["muinf"] * 1.01)
                   for t in 0.5 * G * s])
    # u(r) = int_r^R gd ds
    cum = np.concatenate([[0], np.cumsum(0.5 * (gd[1:] + gd[:-1]) * np.diff(s))])
    u = cum[-1] - cum
    del tw
    return np.interp(r, s, u)


# ----------------------------------------------------------------------------- volume
def revolve(m, xc, half, R, spacing, cut, field):
    """Revolve the (x, r) field into a voxel volume of the segment |x-xc| < half (units of R)."""
    nx = int(np.ceil(2 * half / spacing)) + 1
    ny = int(np.ceil(2.0 / spacing)) + 1
    img = pv.ImageData(dimensions=(nx, ny, ny), spacing=(spacing,) * 3,
                       origin=(-half, -1.0, -1.0))
    P = img.points
    rr = np.sqrt(P[:, 1] ** 2 + P[:, 2] ** 2)
    mn = m.copy()
    mn.points = (mn.points - [xc, 0, 0]) / R      # units of R: VTK locators misbehave at 1e-4 m
    q = pv.PolyData(np.c_[P[:, 0], np.minimum(rr, 1.0 - 1e-9), np.zeros(len(P))])
    s = q.sample(mn, tolerance=1e-6)
    val = s.point_data[field].copy()
    val[rr > 1.0] = 0.0
    if cut == "quarter":
        val[(P[:, 1] > 0) & (P[:, 2] > 0)] = 0.0
    elif cut == "half":
        val[P[:, 2] > 0] = 0.0
    img.point_data["v"] = val
    return img


def wall_surface(half, cut):
    w = pv.Cylinder(center=(0, 0, 0), direction=(1, 0, 0), radius=1.0, height=2 * half,
                    resolution=180, capping=False).triangulate().subdivide(2)
    P = w.points
    if cut == "quarter":
        w = w.extract_points((P[:, 1] <= 1e-9) | (P[:, 2] <= 1e-9), adjacent_cells=False)
    elif cut == "half":
        w = w.extract_points(P[:, 2] <= 1e-9, adjacent_cells=False)
    return w.extract_surface(algorithm="dataset_surface")


def fix_sampling(actor, step):
    """Fixed ray step: VTK's automatic step adaptation shows banding off-screen."""
    mp = actor.mapper
    for name, arg in (("SetAutoAdjustSampleDistances", 0), ("SetSampleDistance", step),
                      ("SetInteractiveSampleDistance", step)):
        try:
            getattr(mp, name)(arg)
        except (AttributeError, TypeError):
            pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folders", nargs="+")
    ap.add_argument("--out", default="cy_tube_rendering.png")
    ap.add_argument("--field", choices=["speed", "viscosity"], default="speed")
    ap.add_argument("--length", type=float, default=8.0, help="segment length (radii)")
    ap.add_argument("--spacing", type=float, default=0.02, help="voxel size (radii)")
    ap.add_argument("--cut", choices=["half", "quarter", "none"], default="half")
    ap.add_argument("--cmap", default=None)
    ap.add_argument("--opacity", type=float, default=0.3)
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--mapper", default="smart", help="smart (GPU) or fixed_point (CPU)")
    ap.add_argument("--dark", action="store_true", help="dark background (default white)")
    ap.add_argument("--scale", type=int, default=2)
    ap.add_argument("--interactive", action="store_true")
    a = ap.parse_args()

    bg, fg = ("#0b0f14", "#e8edf2") if a.dark else ("white", "#111111")
    wallc = "#cfe3f0" if a.dark else "#7d93a3"
    comp_c = "#7fd3ff" if a.dark else "#1f6fb2"     # computed profile
    bead_c = "#ffd166" if a.dark else "#e07a00"     # exact profile
    half = 0.5 * a.length
    runs = []
    for f in a.folders:
        m = read_last(f)
        R = m.bounds[3]
        L = m.bounds[1] - m.bounds[0]
        xc, r, u, U, G = developed(m, R, L)
        m.point_data["speed"] = np.linalg.norm(m.point_data["velocity"], axis=1) / U
        if a.field == "viscosity":
            d = m.compute_derivative(scalars="velocity", gradient="grad")
            g = d.point_data["grad"]          # du/dx du/dy dv/dx dv/dy ...
            m.point_data["viscosity"] = 1e3 * mu_cy(np.abs(g[:, 1]))      # mPa s
        ue = exact_profile(G, R, r)
        gm = re.search(r"g(\d+)", os.path.basename(os.path.normpath(f)))
        label = (f"Reference shear rate {gm.group(1)} 1/s" if gm else os.path.basename(f))
        eQ = 100 * (np.trapezoid(r * u, r) / np.trapezoid(r * ue, r) - 1)
        print(f"{f}: R = {R:.3g} m, U = {U:.4g} m/s, dp/dx = {G:.4g} Pa/m, "
              f"u_max/U = {u.max() / U:.3f} (Poiseuille 2), flow vs exact {eQ:+.2f} %")
        runs.append((m, R, xc, r / R, u / U, ue / U, label))

    if a.field == "speed":
        clim, cmap, title = (0.0, 2.0), a.cmap or ("inferno" if a.dark else "YlOrRd"), "u / U"
    else:
        vals = np.concatenate([x[0].point_data["viscosity"] for x in runs])
        clim = (float(np.percentile(vals, 1)), float(np.percentile(vals, 99)))
        cmap, title = a.cmap or "viridis", "μ (mPa s)"
    opacity = np.linspace(0, 1, 64) ** a.gamma * a.opacity
    if a.field == "viscosity":
        opacity = opacity[::-1] * 0.6 + 0.08 * a.opacity   # viscous core denser
        opacity[0] = 0.0                                    # outside the lumen (value 0)

    pv.global_theme.font.family = "times"
    n = len(runs)
    pl = pv.Plotter(shape=(n, 1), window_size=(1800, 520 * n), off_screen=not a.interactive,
                    border=False)
    for k, (m, R, xc, rh, uh, ueh, label) in enumerate(runs):
        pl.subplot(k, 0)
        pl.set_background(bg)
        vol = revolve(m, xc, half, R, a.spacing, a.cut, a.field)
        act = pl.add_volume(vol, scalars="v", clim=clim, cmap=cmap, opacity=opacity, shade=False,
                      mapper=a.mapper, show_scalar_bar=(k == n - 1),
                      scalar_bar_args=dict(title=title, vertical=True, position_x=0.90,
                                           position_y=0.15, width=0.035, height=0.7,
                                           color=fg, title_font_size=20, label_font_size=16,
                                           n_labels=5, fmt="%.2g"))
        fix_sampling(act, a.spacing / 3)
        pl.add_mesh(wall_surface(half, a.cut), color=wallc, opacity=0.10, smooth_shading=True,
                    specular=0.6, specular_power=40)
        for x0 in (-half, half):
            pl.add_mesh(pv.Circle(radius=1.0, resolution=180).rotate_y(90).translate((x0, 0, 0))
                        .extract_feature_edges(), color=wallc, line_width=2)
        # profiles on the cut plane z = 0, drawn as an axial displacement u/U * amp
        x0, amp = -0.75 * half, 0.55 * half
        ys = np.r_[-rh[::-1], rh[1:]]
        uc = np.r_[uh[::-1], uh[1:]]
        ux = np.r_[ueh[::-1], ueh[1:]]
        z = 0.02 if a.cut == "half" else 0.0
        base = pv.Line((x0, ys.min(), z), (x0, ys.max(), z))
        pl.add_mesh(base, color=fg, line_width=2)
        comp = pv.Spline(np.c_[x0 + amp * uc, ys, np.full_like(ys, z)], 400)
        pl.add_mesh(comp.tube(radius=0.018), color=comp_c, smooth_shading=True)
        pois = 2 * (1 - ys ** 2)
        dash = np.c_[x0 + amp * pois, ys, np.full_like(ys, z)]
        for i in range(0, len(dash) - 6, 12):
            pl.add_mesh(pv.Spline(dash[i:i + 6], 10), color=fg, line_width=3)
        beads = pv.PolyData(np.c_[x0 + amp * ux, ys, np.full_like(ys, z)][::14])
        pl.add_mesh(beads.glyph(geom=pv.Sphere(radius=0.035), scale=False, orient=False), color=bead_c,
                    smooth_shading=True)
        pl.add_text(label, position="upper_left", font_size=16, color=fg)
        if k == 0:
            pl.add_text("tube: computed   beads: exact Carreau–Yasuda   dashes: Poiseuille, same flow",
                        position="upper_right", font_size=11, color=fg)
        pl.camera_position = [(0.0, -1.1 * half, 1.9 * half), (0.0, 0.0, 0.0), (0, 1, 0)]
        pl.camera.zoom(1.25)
        pl.renderer.ResetCameraClippingRange()
    pl.link_views()
    if a.interactive:
        pl.show(auto_close=False)
    pl.enable_anti_aliasing("fxaa")   # ssaa drops ray-cast volumes on some drivers
    pl.screenshot(a.out, scale=a.scale)
    print("written", a.out)


if __name__ == "__main__":
    main()
