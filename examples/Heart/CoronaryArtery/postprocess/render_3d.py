#!/usr/bin/env python3
"""3D results of the coronary FSI: WSS, TAWSS, OSI, streamlines and wall motion.

Reads the Rodin XDMF/HDF5 output of every case
<results>/<geometry>/<scenario>/CoronaryArtery_Explicit_PoroElastic{fluid,solid}.xdmf
(sequential or MPI, one HDF5 piece per rank), merges the rank pieces, and over
the LAST cardiac cycle computes on the wall

  TAWSS = 1/T int_T |tau_w| dt,   OSI = 1/2 (1 - |int_T tau_w dt| / int_T |tau_w| dt).

Snapshots are taken at peak diastolic inflow and at peak LV pressure (from
the CSV). Panels are rendered off-screen with PyVista (no scalar bars) and
assembled with matplotlib into CMAME-width figures with vector colour bars:

  fig3d_wss_diastole   |tau_w| at peak diastolic inflow       (geometry x scenario)
  fig3d_tawss          TAWSS (log scale)
  fig3d_osi            OSI
  fig3d_streamlines    instantaneous streamlines at peak diastolic inflow
  fig3d_wall_motion    cyclic wall displacement at peak systole (deformed, exaggerated)
  fig3d_lesion         close-up of the lesion (80% stenosis): TAWSS, OSI, streamlines
  metrics_3d.csv       wall-averaged indices per case

Usage: render_3d.py RESULTS_DIR FIGURE_DIR [--period 0.85] [--width 2400]
Requires: numpy, h5py, pyvista (pip install h5py pyvista), matplotlib.
"""
import argparse
import glob
import os
import sys
import warnings
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmame_style as st  # noqa: E402

PREFIX = "CoronaryArtery_Explicit_PoroElastic"
LOW_TAWSS = 0.4   # Pa (4 dyn/cm^2)
HIGH_OSI = 0.2

METRICS3D = {
    "TAWSS_mean": (r"$\langle \mathrm{TAWSS}\rangle$", "Pa"),
    "TAWSS_p99": (r"$\mathrm{TAWSS}_{99}$", "Pa"),
    "low_TAWSS_area": (r"$A(\mathrm{TAWSS}<0.4\,\mathrm{Pa})$", "%"),
    "OSI_mean": (r"$\langle \mathrm{OSI}\rangle$", "-"),
    "high_OSI_area": (r"$A(\mathrm{OSI}>0.2)$", "%"),
    "WSS_dia_p99": (r"$|\tau_w|_{99}^{dia}$", "Pa"),
    "umax_dia": (r"$\max|\mathbf{u}|^{dia}$", "m/s"),
    "disp_amp_max": (r"$\max|\mathbf{d}-\mathbf{d}_0|$", "um"),
}

# XDMF mixed-topology ids -> (VTK cell type, nodes)
XDMF_CELLS = {4: (5, 3), 6: (10, 4), 8: (13, 6), 9: (12, 8), 38: (24, 10),
              41: (32, 18)}
# Rodin's former Wedge_18 node order -> VTK order (files written before the fix).
WEDGE18_OLD_TO_VTK = [0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10, 17, 16, 15]
WEDGE6_OLD_TO_VTK = [0, 2, 1, 3, 5, 4]


def surface(mesh, **kw):
    """extract_surface with the classic dataset-surface filter on any PyVista."""
    try:
        return mesh.extract_surface(algorithm="dataset_surface", **kw)
    except TypeError:
        return mesh.extract_surface(**kw)


# ---------------------------------------------------------------------------
# XDMF / HDF5 reading
# ---------------------------------------------------------------------------
class XdmfSeries:
    """Temporal collection of spatial collections (one uniform grid per rank)."""

    def __init__(self, path):
        import h5py
        self.h5py = h5py
        self.dir = os.path.dirname(os.path.abspath(path))
        root = ET.parse(path).getroot()
        temporal = root.find("Domain/Grid")
        self.steps = []
        for step in temporal.findall("Grid"):
            t = float(step.find("Time").get("Value"))
            pieces = step.findall("Grid") or [step]
            self.steps.append((t, pieces))
        if not self.steps:
            raise ValueError(f"{path}: no time steps")
        self.times = np.array([t for t, _ in self.steps])
        self._cache = {}
        self._build_mesh()

    def _read(self, item):
        ref = item.text.strip()
        fname, dset = ref.split(":", 1)
        fpath = os.path.join(self.dir, fname)
        key = (fpath, dset)
        if key in self._cache:
            return self._cache[key]
        with self.h5py.File(fpath, "r") as f:
            data = np.asarray(f[dset][()])
        if "Mesh" in dset:
            self._cache[key] = data
        return data

    @staticmethod
    def _parse_mixed(stream):
        stream = np.asarray(stream, np.int64)
        t0 = int(stream[0])
        if t0 in XDMF_CELLS:
            n = XDMF_CELLS[t0][1]
            if len(stream) % (n + 1) == 0 and np.all(stream[:: n + 1] == t0):
                return {t0: stream.reshape(-1, n + 1)[:, 1:]}
        out, i = {}, 0
        while i < len(stream):
            t = int(stream[i])
            if t not in XDMF_CELLS:
                raise ValueError(f"unsupported XDMF cell id {t}")
            n = XDMF_CELLS[t][1]
            out.setdefault(t, []).append(stream[i + 1:i + 1 + n])
            i += n + 1
        return {t: np.array(v) for t, v in out.items()}

    def _fix_wedges(self, cells, X):
        if 41 in cells:
            c = cells[41][:1]
            mid03 = 0.5 * (X[c[:, 0]] + X[c[:, 3]])
            if np.allclose(X[c[:, 9]], mid03, atol=1e-12 * (1 + np.abs(X).max())):
                cells[41] = cells[41][:, WEDGE18_OLD_TO_VTK]
                self.fixed_wedges = True
        if 8 in cells:
            c = cells[8][:64]
            n = np.cross(X[c[:, 1]] - X[c[:, 0]], X[c[:, 2]] - X[c[:, 0]])
            if np.mean(np.einsum("ij,ij->i", n, X[c[:, 3]] - X[c[:, 0]]) > 0) > 0.5:
                cells[8] = cells[8][:, WEDGE6_OLD_TO_VTK]
                self.fixed_wedges = True

    def _build_mesh(self):
        import pyvista as pv
        self.fixed_wedges = False
        _, pieces = self.steps[0]
        Xs, cellsets, offset = [], {}, 0
        self.piece_sizes = []
        for g in pieces:
            X = self._read(g.find("Geometry/DataItem")).reshape(-1, 3)
            cells = self._parse_mixed(self._read(g.find("Topology/DataItem")))
            self._fix_wedges(cells, X)
            for t, c in cells.items():
                cellsets.setdefault(t, []).append(c + offset)
            Xs.append(X)
            self.piece_sizes.append(len(X))
            offset += len(X)
        X = np.vstack(Xs)
        # Merge the rank pieces (and the per-cell nodes of curved output).
        scale = np.ptp(X, axis=0).max()
        key = np.round(X / (1e-9 * scale)).astype(np.int64)
        _, first, inverse = np.unique(key, axis=0, return_index=True, return_inverse=True)
        self.inverse = inverse.ravel()
        self.points = X[first]
        # Linear cells on the corner nodes (P1/Q1 geometry has straight edges;
        # the surface filters return exact point ids only for linear cells).
        linear = {41: (13, 6), 40: (13, 6), 38: (10, 4)}
        cells = {}
        for t, c in cellsets.items():
            c = self.inverse[np.vstack(c)]
            if t in linear:
                vtk_type, n = linear[t]
                cells[vtk_type] = c[:, :n]
            else:
                cells[XDMF_CELLS[t][0]] = c
        self.grid = pv.UnstructuredGrid(cells, self.points.copy())

    def field(self, step, name):
        """Point field of a step on the merged points (or None)."""
        _, pieces = self.steps[step]
        vals = []
        for g in pieces:
            item = None
            for a in g.findall("Attribute"):
                if a.get("Name") == name and a.get("Center", "Node") == "Node":
                    item = a.find("DataItem")
            if item is None:
                return None
            vals.append(np.asarray(self._read(item), float))
        v = np.concatenate([x.reshape(len(x), -1) for x in vals])
        if len(v) != len(self.inverse):
            raise ValueError(f"{name}: {len(v)} values for {len(self.inverse)} nodes")
        out = np.zeros((len(self.points), v.shape[1]))
        out[self.inverse] = v
        return out if out.shape[1] > 1 else out[:, 0]

    def nearest_step(self, t):
        return int(np.argmin(np.abs(self.times - t)))


# ---------------------------------------------------------------------------
# Case analysis
# ---------------------------------------------------------------------------
def read_csv(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return {n: np.asarray(d[n], float) for n in d.dtype.names}


def cycle_times(csv, period):
    t = csv["t"]
    n = int(np.floor(t[-1] / period + 1e-9))
    t0 = (n - 1) * period if n >= 1 else t[0]
    m = (t >= t0 - 1e-12) & (t <= t0 + period + 1e-12)
    tc, pv, qin = t[m], csv["lv_pv"][m], -csv["q_in"][m]
    sys_m = pv > pv.min() + 0.1 * (pv.max() - pv.min())
    t_sys = tc[np.argmax(pv)]
    dia = ~sys_m if (~sys_m).any() else np.ones_like(sys_m)
    t_dia = tc[dia][np.argmax(qin[dia])]
    return t0, t0 + period, t_sys, t_dia


class Case:
    def __init__(self, results, geom, scen, period):
        self.geom, self.scen = geom, scen
        self.dir = os.path.join(results, geom, scen)
        self.csv = read_csv(os.path.join(self.dir, PREFIX + ".csv"))
        self.t0, self.t1, self.t_sys, self.t_dia = cycle_times(self.csv, period)
        self.fluid = XdmfSeries(os.path.join(self.dir, PREFIX + "fluid.xdmf"))
        spath = os.path.join(self.dir, PREFIX + "solid.xdmf")
        self.solid = XdmfSeries(spath) if os.path.exists(spath) else None
        self.metrics = {}
        self.inlet_center = None
        tt = self.fluid.times
        self.cycle_steps = np.flatnonzero((tt >= self.t0 - 1e-9) & (tt <= self.t1 + 1e-9))
        if len(self.cycle_steps) < 2:
            self.cycle_steps = np.arange(len(tt))
            print(f"  [{geom}/{scen}] warning: fewer than 2 snapshots in the last cycle; "
                  f"using all {len(tt)} snapshot(s) (t = {tt[0]:.3f}..{tt[-1]:.3f} s)")
        self.s_dia = self.fluid.nearest_step(self.t_dia)
        self.s_sys = self.fluid.nearest_step(self.t_sys)

    def label(self):
        return f"{self.geom}/{self.scen}"

    # --- fluid ------------------------------------------------------------
    def analyse_fluid(self):
        import pyvista as pv
        fl = self.fluid
        surf = surface(fl.grid, pass_pointid=True, pass_cellid=False)
        vid = np.asarray(surf.point_data["vtkOriginalPointIds"])
        surf.point_data["vol_id"] = vid
        steps = self.cycle_steps
        taus = np.stack([fl.field(s, "shearStress")[vid] for s in steps])  # (n, p, 3)
        times = fl.times[steps]
        mags = np.linalg.norm(taus, axis=2)
        wall = mags.max(axis=0) > 0.0
        if len(steps) > 1:
            T = times[-1] - times[0]
            integ = getattr(np, "trapezoid", None) or np.trapz
            tawss = integ(mags, times, axis=0) / T
            mean_vec = integ(taus, times, axis=0) / T
        else:
            tawss, mean_vec = mags[0], taus[0]
        osi = np.where(tawss > 1e-14,
                       0.5 * (1 - np.linalg.norm(mean_vec, axis=1) / np.maximum(tawss, 1e-30)),
                       0.0)
        osi = np.clip(osi, 0.0, 0.5)
        tau_dia = np.linalg.norm(fl.field(self.s_dia, "shearStress")[vid], axis=1)
        surf.point_data["TAWSS"] = tawss
        surf.point_data["OSI"] = osi
        surf.point_data["WSS_dia"] = tau_dia
        surf.point_data["wall"] = wall.astype(float)
        # Wall only: drop the inlet/outlet caps (no shear data there).
        cell_wall = surf.point_data_to_cell_data(pass_point_data=False)["wall"] > 0.999
        self.wall = surface(surf.extract_cells(np.flatnonzero(cell_wall)))
        self.caps = surface(surf.extract_cells(np.flatnonzero(~cell_wall)))
        area = self.wall.compute_cell_sizes(length=False, volume=False)["Area"]
        cd = self.wall.point_data_to_cell_data()
        A = area.sum()
        self.metrics.update({
            "TAWSS_mean": float((cd["TAWSS"] * area).sum() / A),
            "TAWSS_p99": float(np.percentile(self.wall.point_data["TAWSS"], 99)),
            "low_TAWSS_area": float(100 * area[cd["TAWSS"] < LOW_TAWSS].sum() / A),
            "OSI_mean": float((cd["OSI"] * area).sum() / A),
            "high_OSI_area": float(100 * area[cd["OSI"] > HIGH_OSI].sum() / A),
            "WSS_dia_p99": float(np.percentile(self.wall.point_data["WSS_dia"], 99)),
        })
        # Current configuration and velocity at peak diastolic inflow.
        g = fl.grid.copy()
        disp = fl.field(self.s_dia, "ALEDisp")
        if disp is not None:
            g.points = fl.points + disp
            self.wall_dia = self.wall.copy()
            # wall point coordinates follow the volume ones
            self.wall_dia.points = g.points[np.asarray(self.wall.point_data["vol_id"])]
        else:
            self.wall_dia = self.wall
        u = fl.field(self.s_dia, "FluidVelocity")
        g.point_data["U"] = u
        g.point_data["Umag"] = np.linalg.norm(u, axis=1)
        self.grid_dia = g
        self.metrics["umax_dia"] = float(g.point_data["Umag"].max())
        self.streamlines = self._streamlines(g)

    def _inlet_seeds(self, g, n_seeds=600):
        import pyvista as pv
        caps = self.caps.connectivity(extraction_mode="all")
        best, best_area = None, -1.0
        for rid in np.unique(caps.cell_data["RegionId"]):
            part = caps.extract_cells(np.flatnonzero(caps.cell_data["RegionId"] == rid))
            part = surface(part)
            a = part.compute_cell_sizes(length=False, volume=False)["Area"].sum()
            if a > best_area:
                best, best_area = part, a
        if best is None:
            return None
        c = best.points.mean(axis=0)
        self.inlet_center = c
        normal = best.compute_normals(cell_normals=True, point_normals=False)["Normals"].mean(0)
        normal /= np.linalg.norm(normal) + 1e-30
        r = np.sqrt(best_area / np.pi)
        rng = np.random.default_rng(1)
        rho = 0.92 * r * np.sqrt(rng.random(4 * n_seeds))
        th = 2 * np.pi * rng.random(4 * n_seeds)
        e1 = np.cross(normal, [1.0, 0.0, 0.0])
        if np.linalg.norm(e1) < 1e-6:
            e1 = np.cross(normal, [0.0, 1.0, 0.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(normal, e1)
        for sign in (-1.0, 1.0):
            P = (c + sign * 0.05 * r * normal + rho[:, None] * np.cos(th)[:, None] * e1
                 + rho[:, None] * np.sin(th)[:, None] * e2)
            inside = g.find_containing_cell(P) >= 0
            if inside.mean() > 0.5:
                return pv.PolyData(P[inside][:n_seeds])
        return None

    def _streamlines(self, g):
        seeds = self._inlet_seeds(g)
        if seeds is None:
            print(f"  [{self.label()}] warning: no inlet seeds, streamlines skipped")
            return None
        diag = np.linalg.norm(np.ptp(g.points, axis=0))
        lines = g.streamlines_from_source(
            seeds, vectors="U", integration_direction="forward", initial_step_length=0.2,
            step_unit="cl", min_step_length=0.01, max_step_length=0.5, max_steps=20000,
            max_length=4 * diag, terminal_speed=1e-8, compute_vorticity=False)
        return lines if lines.n_points else None

    # --- solid ------------------------------------------------------------
    def analyse_solid(self):
        if self.solid is None:
            self.solid_surf = None
            return
        so = self.solid
        steps = np.flatnonzero((so.times >= self.t0 - 1e-9) & (so.times <= self.t1 + 1e-9))
        if len(steps) == 0:
            steps = np.arange(len(so.times))
        d0 = so.field(steps[0], "Displacement")
        amp = np.zeros(len(so.points))
        for s in steps:
            amp = np.maximum(amp, np.linalg.norm(so.field(s, "Displacement") - d0, axis=1))
        s_sys = so.nearest_step(self.t_sys)
        dsys = so.field(s_sys, "Displacement")
        surf = surface(so.grid, pass_pointid=True)
        pid = np.asarray(surf.point_data["vtkOriginalPointIds"])
        surf.point_data["dcyc"] = 1e6 * amp[pid]  # cyclic amplitude, um
        surf.point_data["d0"] = d0[pid]
        surf.point_data["dc"] = (dsys - d0)[pid]
        self.solid_surf = surf
        self.metrics["disp_amp_max"] = float(1e6 * np.percentile(amp, 99.5))


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
class View:
    """Common camera: PCA of the fluid geometry, longest axis horizontal."""

    def __init__(self, points, inlet_hint=None):
        c = points.mean(axis=0)
        _, _, vt = np.linalg.svd(points - c, full_matrices=False)
        e_long, e_mid, e_normal = vt
        if inlet_hint is not None and np.dot(inlet_hint - c, e_long) > 0:
            e_long = -e_long  # inlet on the left
        if np.dot(np.cross(e_long, e_mid), e_normal) < 0:
            e_normal = -e_normal
        self.center = c
        self.up = e_mid
        self.normal = e_normal
        proj = (points - c) @ np.stack([e_long, e_mid]).T
        self.half = np.abs(proj).max(axis=0)

    def apply(self, plotter, focus=None, scale=None, aspect=16 / 9):
        cam = plotter.camera
        f = self.center if focus is None else focus
        cam.focal_point = f
        cam.position = f + self.normal * 10 * max(self.half)
        cam.up = self.up
        plotter.enable_parallel_projection()
        if scale is None:
            scale = 1.06 * max(self.half[1], self.half[0] / aspect)
        cam.parallel_scale = scale
        plotter.reset_camera_clipping_range()


def make_plotter(width, height):
    import pyvista as pv
    p = pv.Plotter(off_screen=True, window_size=(width, height), lighting="light kit")
    p.set_background("white")
    try:
        p.enable_anti_aliasing("ssaa")
    except Exception:
        pass
    return p


SURF_KW = dict(smooth_shading=True, ambient=0.25, diffuse=0.75, specular=0.15,
               specular_power=20, show_scalar_bar=False)


def render_scalar(path, mesh, scalars, cmap, clim, view, width, height, log=False,
                  context=None, focus=None, scale=None):
    p = make_plotter(width, height)
    if context is not None:
        p.add_mesh(context, color="#d9d9d9", opacity=0.25, smooth_shading=True)
    p.add_mesh(mesh, scalars=scalars, cmap=cmap, clim=clim, log_scale=log, **SURF_KW)
    view.apply(p, focus, scale, width / height)
    p.screenshot(path, transparent_background=False)
    p.close()


def render_streamlines(path, lines, wall, clim, view, width, height, focus=None,
                       scale=None):
    p = make_plotter(width, height)
    p.add_mesh(wall, color="#cfcfcf", opacity=0.18, smooth_shading=True)
    if lines is not None:
        lw = max(1.5, width / 900)
        p.add_mesh(lines, scalars="Umag", cmap="turbo", clim=clim, line_width=lw,
                   render_lines_as_tubes=True, show_scalar_bar=False)
    view.apply(p, focus, scale, width / height)
    p.screenshot(path, transparent_background=False)
    p.close()


def render_wall_motion(path, surf, clim, view, width, height, exaggeration):
    p = make_plotter(width, height)
    s = surf.copy()
    s.points = s.points + s.point_data["d0"] + exaggeration * s.point_data["dc"]
    p.add_mesh(s, scalars="dcyc", cmap="viridis", clim=clim, **SURF_KW)
    view.apply(p, None, None, width / height)
    p.screenshot(path, transparent_background=False)
    p.close()


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------
def crop_box(paths, pad=12):
    from PIL import Image
    box, size = None, None
    for path in paths:
        im = np.asarray(Image.open(path).convert("L"))
        size = im.shape
        ys, xs = np.nonzero(im < 250)
        if len(xs) == 0:
            continue
        b = [xs.min(), ys.min(), xs.max(), ys.max()]
        box = b if box is None else [min(box[0], b[0]), min(box[1], b[1]),
                                     max(box[2], b[2]), max(box[3], b[3])]
    if box is None:
        return None
    return (max(box[0] - pad, 0), max(box[1] - pad, 0), min(box[2] + pad, size[1]),
            min(box[3] + pad, size[0]))


def assemble(grid_paths, row_labels, col_labels, cbar, out, width=st.DOUBLE,
             note=None):
    """grid_paths[r][c] -> PNG (or None); cbar = (cmap, norm, label)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    flat = [p for row in grid_paths for p in row if p]
    if not flat:
        return
    box = crop_box(flat)
    ims = [[np.asarray(Image.open(p).convert("RGB").crop(box)) if p else None for p in row]
           for row in grid_paths]
    h, w = next(im for row in ims for im in row if im is not None).shape[:2]
    nr, nc = len(ims), len(ims[0])
    cb_w = 0.09
    lab_w = 0.05 if any(row_labels) else 0.0
    panel_w = width * (1 - cb_w - lab_w) / nc
    fig_h = nr * panel_w * h / w * 1.1 + 10 * st.MM
    fig = plt.figure(figsize=(width, fig_h))
    top = 1 - 8 * st.MM / fig_h
    gs = fig.add_gridspec(nr, nc, left=lab_w, right=1 - cb_w, top=top, bottom=0.01,
                          wspace=0.03, hspace=0.10)
    letters = "abcdefghijklmnopqrstuvwxyz"
    k = 0
    for r in range(nr):
        for c in range(nc):
            ax = fig.add_subplot(gs[r, c])
            ax.set_axis_off()
            ax.imshow(ims[r][c] if ims[r][c] is not None else np.full((h, w, 3), 255, np.uint8))
            if ims[r][c] is not None:
                ax.text(0.0, 1.0, f"({letters[k]})", transform=ax.transAxes,
                        ha="left", va="bottom", fontsize=8)
                k += 1
            if r == 0 and col_labels:
                ax.set_title(col_labels[c], fontsize=8, pad=9)
            if c == 0 and row_labels and row_labels[r]:
                ax.text(-0.03, 0.5, row_labels[r], transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=8)
    cmap, norm, label = cbar
    cax = fig.add_axes([1 - cb_w + 0.018, 0.12, 0.012, top - 0.2])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(label)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2.5, which="major")
    cb.ax.tick_params(width=0.4, length=1.5, which="minor")
    if note:
        fig.text(1 - 0.01, 0.005, note, ha="right", va="bottom", fontsize=6.5, color="0.3")
    st.save(fig, out)
    plt.close(fig)
    print(f"  wrote {out}.pdf")


def robust_range(arrays, lo=1.0, hi=99.5, positive=False):
    v = np.concatenate([np.ravel(a) for a in arrays if a is not None and np.size(a)])
    if positive:
        v = v[v > 0]
    if v.size == 0:
        return (1e-3, 1.0) if positive else (0.0, 1.0)
    a, b = np.percentile(v, [lo, hi])
    if b <= a:
        b = a + (abs(a) + 1e-12)
    return float(a), float(b)


def nice_log_range(a, b):
    """Widen [a, b] to half decades."""
    lo = 10 ** (np.floor(2 * np.log10(max(a, 1e-6))) / 2)
    hi = 10 ** (np.ceil(2 * np.log10(max(b, 1e-6))) / 2)
    return float(lo), float(max(hi, 3 * lo))


# ---------------------------------------------------------------------------
def discover(results):
    out = []
    for path in glob.glob(os.path.join(results, "*", "*", PREFIX + "fluid.xdmf")):
        d = os.path.dirname(path)
        if os.path.exists(os.path.join(d, PREFIX + ".csv")):
            out.append((os.path.basename(os.path.dirname(d)), os.path.basename(d)))

    def key(c):
        g = st.GEOMETRY_ORDER.index(c[0]) if c[0] in st.GEOMETRY_ORDER else 99
        s = st.SCENARIO_ORDER.index(c[1]) if c[1] in st.SCENARIO_ORDER else 99
        return (g, s) + c

    return sorted(out, key=key)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results")
    ap.add_argument("figures")
    ap.add_argument("--period", type=float, default=0.85)
    ap.add_argument("--width", type=int, default=2400, help="panel width in pixels")
    ap.add_argument("--exaggeration", type=float, default=0.0,
                    help="wall-motion warp factor (0: automatic)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import colors
    import pyvista as pv
    pv.OFF_SCREEN = True
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", message=".*algorithm.*")
    st.apply()

    found = discover(args.results)
    if not found:
        sys.exit(f"No {PREFIX}fluid.xdmf + CSV under {args.results}/<geometry>/<scenario>/")
    panels = os.path.join(args.figures, "panels3d")
    os.makedirs(panels, exist_ok=True)

    cases = []
    for geom, scen in found:
        print(f"  [{geom}/{scen}] reading ...", flush=True)
        c = Case(args.results, geom, scen, args.period)
        if c.fluid.fixed_wedges or (c.solid is not None and c.solid.fixed_wedges):
            print(f"  [{geom}/{scen}] wedge node order converted to VTK (pre-fix output)")
        c.analyse_fluid()
        c.analyse_solid()
        print(f"  [{geom}/{scen}] snapshots in last cycle: {len(c.cycle_steps)}, "
              f"t_dia={c.fluid.times[c.s_dia]:.3f} s, t_sys={c.fluid.times[c.s_sys]:.3f} s",
              flush=True)
        cases.append(c)

    # metrics
    with open(os.path.join(args.figures, "metrics_3d.csv"), "w") as f:
        names = list(METRICS3D)
        f.write("geometry,scenario," + ",".join(names) + "\n")
        for c in cases:
            f.write(f"{c.geom},{c.scen}," +
                    ",".join(f"{c.metrics.get(n, np.nan):.6g}" for n in names) + "\n")

    geoms = list(dict.fromkeys(c.geom for c in cases))
    scens = list(dict.fromkeys(c.scen for c in cases))
    by = {(c.geom, c.scen): c for c in cases}
    ref = by.get((geoms[0], scens[0]))
    # inlet (largest cap) of the reference case on the left
    view = View(ref.fluid.points, ref.inlet_center)
    W = args.width
    H = int(W * 9 / 16)

    def grid_of(fn, tag):
        out = []
        for g in geoms:
            row = []
            for s in scens:
                c = by.get((g, s))
                if c is None:
                    row.append(None)
                    continue
                path = os.path.join(panels, f"{tag}_{g}_{s}.png")
                fn(c, path)
                row.append(path)
            out.append(row)
        return out

    rows = [st.geometry_label(g) for g in geoms] if len(geoms) > 1 else [""]
    cols = [st.scenario_label(s) for s in scens]

    # WSS at peak diastole, TAWSS: shared log scale
    wl = robust_range([c.wall.point_data["WSS_dia"] for c in cases]
                      + [c.wall.point_data["TAWSS"] for c in cases], 2, 99.5, positive=True)
    wl = nice_log_range(*wl)
    norm_w = colors.LogNorm(*wl)
    g = grid_of(lambda c, p: render_scalar(p, c.wall_dia, "WSS_dia", "magma", wl, view, W, H,
                                           log=True), "wss_dia")
    assemble(g, rows, cols, ("magma", norm_w, r"$|\boldsymbol{\tau}_w|$ [Pa]"),
             os.path.join(args.figures, "fig3d_wss_diastole"),
             note="Peak diastolic inflow, last cycle")
    g = grid_of(lambda c, p: render_scalar(p, c.wall, "TAWSS", "magma", wl, view, W, H,
                                           log=True), "tawss")
    assemble(g, rows, cols, ("magma", norm_w, "TAWSS [Pa]"),
             os.path.join(args.figures, "fig3d_tawss"))

    ol = (0.0, max(0.05, min(0.5, robust_range([c.wall.point_data["OSI"] for c in cases],
                                               0, 99.5)[1])))
    g = grid_of(lambda c, p: render_scalar(p, c.wall, "OSI", "cividis", ol, view, W, H),
                "osi")
    assemble(g, rows, cols, ("cividis", colors.Normalize(*ol), "OSI [-]"),
             os.path.join(args.figures, "fig3d_osi"))

    ul = (0.0, robust_range([c.grid_dia.point_data["Umag"] for c in cases], 0, 99.9)[1])
    g = grid_of(lambda c, p: render_streamlines(p, c.streamlines, c.wall_dia, ul, view, W, H),
                "streamlines")
    assemble(g, rows, cols, ("turbo", colors.Normalize(*ul), r"$|\mathbf{u}|$ [m/s]"),
             os.path.join(args.figures, "fig3d_streamlines"),
             note="Instantaneous streamlines from the inlet, peak diastolic inflow")

    solids = [c for c in cases if c.solid_surf is not None]
    if solids:
        dl = (0.0, robust_range([c.solid_surf.point_data["dcyc"] for c in solids], 0, 99.5)[1])
        diag = np.linalg.norm(np.ptp(ref.fluid.points, axis=0))
        dmax = max(1e-12, dl[1] * 1e-6)
        if args.exaggeration:
            ex = args.exaggeration
        elif dl[1] < 1e-3:  # below a nanometre: nothing to exaggerate
            ex = 1.0
        else:
            ex = max(1.0, float(10 ** np.round(np.log10(0.02 * diag / dmax))))
        g = grid_of(lambda c, p: render_wall_motion(p, c.solid_surf, dl, view, W, H, ex)
                    if c.solid_surf is not None else None, "wall_motion")
        assemble(g, rows, cols, ("viridis", colors.Normalize(*dl),
                                 r"$\max_t|\mathbf{d}(t)-\mathbf{d}(t_0)|$ [$\mu$m]"),
                 os.path.join(args.figures, "fig3d_wall_motion"),
                 note=("Cyclic displacement amplitude over the last cycle; wall shown at "
                       f"peak LV pressure, motion exaggerated {ex:g}$\\times$"))

    # Lesion close-up (80% stenosis): focus on the TAWSS maximum.
    lesion = [c for c in cases if c.geom == "stenosis80"]
    if lesion:
        h = lesion[0]
        healthy_geom = [c for c in cases if c.geom == "nostenosis"]
        if healthy_geom:
            # where the stenotic wall departs from the stenosis-free one
            dist = np.abs(h.wall.compute_implicit_distance(healthy_geom[0].wall)
                          .point_data["implicit_distance"])
            top = h.wall.points[dist >= 0.5 * dist.max()]
        else:
            tw = h.wall.point_data["TAWSS"]
            top = h.wall.points[tw >= np.percentile(tw, 99.7)]
        focus = top.mean(axis=0)
        extent = np.linalg.norm(np.ptp(top, axis=0)) if len(top) > 1 else 0.0
        scale = min(max(0.12 * max(view.half), 1.6 * extent), 0.35 * max(view.half))
        Wc, Hc = W, int(W * 9 / 16)
        rows_c, grid_c = [], []
        for tag, fn in (
            ("TAWSS", lambda c, p: render_scalar(p, c.wall, "TAWSS", "magma", wl, view, Wc,
                                                 Hc, log=True, focus=focus, scale=scale)),
            ("OSI", lambda c, p: render_scalar(p, c.wall, "OSI", "cividis", ol, view, Wc, Hc,
                                               focus=focus, scale=scale)),
            ("Streamlines", lambda c, p: render_streamlines(p, c.streamlines, c.wall_dia, ul,
                                                            view, Wc, Hc, focus=focus,
                                                            scale=scale))):
            row = []
            for s in scens:
                c = by.get(("stenosis80", s))
                if c is None:
                    row.append(None)
                    continue
                path = os.path.join(panels, f"lesion_{tag}_{s}.png")
                fn(c, path)
                row.append(path)
            grid_c.append(row)
            rows_c.append(tag)
        for (tag, cmap, norm, label), row in zip(
                (("TAWSS", "magma", norm_w, "TAWSS [Pa]"),
                 ("OSI", "cividis", colors.Normalize(*ol), "OSI [-]"),
                 ("Streamlines", "turbo", colors.Normalize(*ul), r"$|\mathbf{u}|$ [m/s]")),
                grid_c):
            assemble([row], [""], cols, (cmap, norm, label),
                     os.path.join(args.figures, f"fig3d_lesion_{tag.lower()}"),
                     note="80% stenosis, close-up of the lesion")
    print(f"  wrote {os.path.join(args.figures, 'metrics_3d.csv')}")


if __name__ == "__main__":
    main()
