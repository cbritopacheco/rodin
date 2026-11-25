#!/usr/bin/env python3
"""Volume rendering of the speed in the three-dimensional Y at peak inlet pressure.

Reads the three peak-pressure fields of the lower-load run (fully resolved Y,
shear-dependent outlets, constant outlets), clips the resolved Y to the region
common to all three, and renders the speed with the same colour and opacity
scales in every panel: a ray-cast volume of |u|, a glass vessel wall and,
optionally, streamline tubes seeded at the inlet.

  pip install pyvista scipy
  python3 render_y_volume.py data/ --streamlines --interactive   # rotate, press q to save
  python3 render_y_volume.py data/ --streamlines --scale 3       # 7200x2700 PNG, dark
  python3 render_y_volume.py data/ --streamlines --dark          # dark background
  python3 render_y_volume.py data/ --mapper fixed_point          # CPU mapper if the GPU one fails

Inputs (in the folder given): truncation_y_low_head_peak.vtu,
truncation_y_low_head_peak_variable.vtu, truncation_y_low_head_peak_constant.vtu,
each with the point field velocity_m_s on quadratic tetrahedra.
"""
import argparse
import os

import numpy as np
import pyvista as pv

CASES = [("truncation_y_low_head_peak.vtu", "Fully resolved"),
         ("truncation_y_low_head_peak_variable.vtu", "Shear-dependent outlets"),
         ("truncation_y_low_head_peak_constant.vtu", "Constant outlets")]


def load(path, common=None):
    """Linear copy of the field with the speed in mm/s.

    With `common` (a truncated mesh), keep only the cells of the resolved Y that
    the truncated mesh also contains: it is a subset of the full mesh, so its
    cell centres coincide with those of the retained cells.
    """
    m = pv.read(path).linear_copy()
    m.points = m.points * 1e6  # m -> um: VTK ray casting and depth peeling fail at sub-mm scales
    if common is not None:
        from scipy.spatial import cKDTree
        ref = 1e6 * common.linear_copy().cell_centers().points
        d, _ = cKDTree(ref).query(m.cell_centers().points)
        h = 1e-3 * (m.bounds[5] - m.bounds[4])
        m = m.extract_cells(np.where(d < h)[0])
    m.point_data["speed"] = 1e3 * np.linalg.norm(m.point_data["velocity_m_s"], axis=1)  # mm/s
    return m


def resample(mesh, spacing):
    """Speed on a uniform grid; points outside the lumen get speed 0."""
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    dims = np.ceil(np.array([xmax - xmin, ymax - ymin, zmax - zmin]) / spacing).astype(int) + 1
    img = pv.ImageData(dimensions=dims, spacing=(spacing,) * 3, origin=(xmin, ymin, zmin))
    s = img.sample(mesh, tolerance=spacing * 0.5)
    mask = s.point_data["vtkValidPointMask"].astype(bool)
    speed = np.where(mask, s.point_data["speed"], 0.0)
    img.point_data["speed"] = speed
    return img


def smooth(img, std):
    """Gaussian-smoothed speed: removes the facets of the piecewise-linear field."""
    if std <= 0:
        return img
    out = img.gaussian_smooth(radius_factor=3.0, std_dev=std, scalars="speed")
    out.point_data["speed"] = np.maximum(out.point_data["speed"], 0.0)
    return out


def streamlines(mesh, n=60):
    """Streamline tubes seeded on a disc near the inlet (x = x_min)."""
    m = mesh.copy()
    m.point_data.set_vectors(m.point_data["velocity_m_s"], "velocity_m_s")
    x0 = m.bounds[0] + 0.02 * (m.bounds[1] - m.bounds[0])
    r = 0.9 * 0.5 * (m.bounds[5] - m.bounds[4])
    seed = pv.Disc(center=(x0, 0.0, 0.0), inner=0.0, outer=r, normal=(1, 0, 0),
                   r_res=max(2, int(np.sqrt(n / 6))), c_res=12)
    sl = m.streamlines_from_source(seed, vectors="velocity_m_s", integration_direction="forward",
                                   max_steps=4000, initial_step_length=0.2,
                                   max_step_length=0.5, step_unit="cl")
    return sl.tube(radius=0.006 * (m.bounds[5] - m.bounds[4]), n_sides=12) if sl.n_points else None


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
    ap.add_argument("folder", help="folder with the three .vtu files")
    ap.add_argument("--out", default="y_volume_rendering.png")
    ap.add_argument("--spacing", type=float, default=1.5, help="voxel size (um)")
    ap.add_argument("--vmax", type=float, default=None, help="upper end of the speed scale (mm/s)")
    ap.add_argument("--cmap", default=None)
    ap.add_argument("--streamlines", action="store_true")
    ap.add_argument("--scale", type=int, default=3, help="screenshot magnification")
    ap.add_argument("--interactive", action="store_true")
    ap.add_argument("--zoom", type=float, default=1.15)
    ap.add_argument("--mapper", default="smart",
                    help="volume mapper: smart (GPU), fixed_point (CPU, works without GPU)")
    ap.add_argument("--opacity", type=float, default=0.35, help="peak volume opacity")
    ap.add_argument("--gamma", type=float, default=1.6,
                    help="opacity ramp exponent: >1 lets the fast core glow, <1 shows slow flow")
    ap.add_argument("--shade", action="store_true", help="gradient shading of the volume")
    ap.add_argument("--smooth", type=float, default=1.0,
                    help="Gaussian smoothing of the voxel field (std in voxels, 0 = off)")
    ap.add_argument("--no-peel", action="store_true", help="disable depth peeling")
    ap.add_argument("--dark", action="store_true", help="dark background (default white)")
    ap.add_argument("--seeds", type=int, default=60, help="streamlines per panel")
    a = ap.parse_args()

    # common region: the bounds of the truncated meshes
    trunc = pv.read(os.path.join(a.folder, CASES[1][0]))
    meshes = [load(os.path.join(a.folder, CASES[0][0]), trunc)] + \
             [load(os.path.join(a.folder, f)) for f, _ in CASES[1:]]
    print("cells:", [m.n_cells for m in meshes])
    vmax = a.vmax or max(float(m.point_data["speed"].max()) for m in meshes)
    clim = (0.0, vmax)
    # opacity: transparent at rest, rising smoothly with speed (64 points: pyvista
    # takes a 256-entry array as raw uint8 and interpolates shorter float ones)
    opacity = np.clip(np.linspace(0.0, 1.0, 64) ** a.gamma, 0.0, 1.0) * a.opacity

    pv.global_theme.font.family = "times"
    bg, fg = ("#0b0f14", "#e8edf2") if a.dark else ("white", "black")
    a.light = not a.dark
    cmap = a.cmap or ("inferno" if a.dark else "YlOrRd")
    pl = pv.Plotter(shape=(1, 3), window_size=(2400, 900), off_screen=not a.interactive,
                    border=False)
    for k, (mesh, (_, title)) in enumerate(zip(meshes, CASES)):
        pl.subplot(0, k)
        pl.set_background(bg)
        vol = smooth(resample(mesh, a.spacing), a.smooth)
        act = pl.add_volume(vol, scalars="speed", clim=clim, cmap=cmap, opacity=opacity,
                      shade=a.shade, ambient=0.35, diffuse=0.75, specular=0.25,
                      show_scalar_bar=(k == len(meshes) - 1), mapper=a.mapper,
                      scalar_bar_args=dict(title="Speed (mm/s)", vertical=False,
                                           position_x=0.2, position_y=0.05, height=0.05,
                                           width=0.6, color=fg, title_font_size=20,
                                           label_font_size=16, n_labels=5, fmt="%.1f"))
        wall = mesh.extract_surface(algorithm="dataset_surface").smooth(n_iter=30)
        fix_sampling(act, a.spacing / 3)
        pl.add_mesh(wall, color="#9fb3c2" if a.light else "#cfe3f0", opacity=0.08, smooth_shading=True,
                    specular=0.6, specular_power=40)
        rims = wall.extract_feature_edges(feature_angle=60, boundary_edges=True,
                                          feature_edges=True, manifold_edges=False,
                                          non_manifold_edges=False)
        if rims.n_points:
            pl.add_mesh(rims, color="#4a5a66" if a.light else "#8aa0b0", line_width=2)
        if a.streamlines:
            tubes = streamlines(mesh, a.seeds)
            if tubes is not None:
                pl.add_mesh(tubes, scalars="speed", clim=clim, cmap=cmap, smooth_shading=True,
                            specular=0.4, show_scalar_bar=False)
        pl.add_text(title, position="upper_edge", font_size=14, color=fg)
        pl.view_xy()
        pl.camera.elevation = 25
        pl.camera.azimuth = -15
        pl.reset_camera()
        pl.camera.zoom(a.zoom)
        if not a.no_peel:
            pl.enable_depth_peeling(number_of_peels=8)
    pl.link_views()
    if a.interactive:
        pl.show(auto_close=False)
    pl.enable_anti_aliasing("fxaa")   # ssaa drops ray-cast volumes on some drivers
    pl.screenshot(a.out, scale=a.scale, transparent_background=False)
    print("written", a.out, "speed scale 0 -", round(vmax, 3), "mm/s")


if __name__ == "__main__":
    main()
