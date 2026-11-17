#!/usr/bin/env python3
"""Velocity profiles on cross-sections along the coronary branches.

For every case <results>/<geometry>/<scenario>/ (Rodin XDMF/HDF5 output):

  1. centerlines (once per geometry, cached in <figures>/centerlines_<geometry>.npz):
     Dijkstra paths on the fluid mesh from the inlet to every outlet, with
     edge weights |e| / r^2 (r the distance to the wall) so that the paths
     follow the lumen centre; smoothed and resampled by arc length. The
     longest path is the main vessel; each other outlet contributes the part
     of its path after it departs from the vessels already kept.
  2. cross-sections normal to the centerline at s/L = 0.1, 0.37, 0.63, 0.9
     of the main vessel and of the N longest branches, coloured by the
     velocity magnitude at peak diastolic inflow (or peak LV pressure), all
     at the same scale so the calibres compare;
  3. figures (CMAME widths): per case a locator + sections grid, and per
     geometry the main-vessel sections of all scenarios side by side.

Usage: sections_3d.py RESULTS FIGURES [--branches 2] [--stations 4]
                      [--time diastole|systole] [--period 0.85]
"""
import argparse
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmame_style as st  # noqa: E402
import render_3d as r3  # noqa: E402

STATION_COLORS = ["#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#56B4E9"]
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII"]


# ---------------------------------------------------------------------------
# Geometry: wall, caps, centerlines
# ---------------------------------------------------------------------------
def wall_and_caps(series, step):
    surf = r3.surface(series.grid, pass_pointid=True, pass_cellid=False)
    vid = np.asarray(surf.point_data["vtkOriginalPointIds"])
    tau = series.field(step, "shearStress")
    wallp = (np.linalg.norm(tau, axis=1)[vid] > 0.0) if tau is not None else None
    if wallp is None or not wallp.any():
        raise RuntimeError("shearStress is zero: cannot tell the wall from the caps")
    surf.point_data["wall"] = wallp.astype(float)
    cw = surf.point_data_to_cell_data(pass_point_data=False)["wall"] > 0.999
    wall = r3.surface(surf.extract_cells(np.flatnonzero(cw)))
    caps = surf.extract_cells(np.flatnonzero(~cw)).connectivity(extraction_mode="all")
    parts = []
    for rid in np.unique(caps.cell_data["RegionId"]):
        part = r3.surface(caps.extract_cells(np.flatnonzero(caps.cell_data["RegionId"] == rid)))
        area = part.compute_cell_sizes(length=False, volume=False)["Area"].sum()
        parts.append((area, part.points.mean(axis=0)))
    parts.sort(key=lambda a: -a[0])
    return wall, parts  # parts[0] = inlet (largest cap)


def smooth_resample(P, h, passes=3, window=7):
    P = np.asarray(P, float)
    for _ in range(passes):
        Q = P.copy()
        k = window // 2
        for i in range(1, len(P) - 1):
            a, b = max(0, i - k), min(len(P), i + k + 1)
            Q[i] = P[a:b].mean(axis=0)
        P = Q
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    n = max(int(s[-1] / h) + 1, 2)
    t = np.linspace(0.0, s[-1], n)
    return np.stack([np.interp(t, s, P[:, c]) for c in range(3)], axis=1)


def centerlines(series, step, cache=None):
    if cache and os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return list(z["vessels"]), float(z["h"])
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    from scipy.spatial import cKDTree

    g = series.grid
    wall, caps = wall_and_caps(series, step)
    X = np.asarray(g.points)
    dist = np.abs(np.asarray(g.compute_implicit_distance(wall)["implicit_distance"]))
    tets = g.cells_dict.get(10)
    if tets is None:
        raise RuntimeError("centerlines need a tetrahedral fluid mesh")
    pairs = np.vstack([tets[:, [a, b]] for a, b in
                       ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))])
    pairs = np.unique(np.sort(pairs, axis=1), axis=0)
    L = np.linalg.norm(X[pairs[:, 0]] - X[pairs[:, 1]], axis=1)
    rr = np.maximum(0.5 * (dist[pairs[:, 0]] + dist[pairs[:, 1]]), 0.05 * dist.max())
    w = L / rr ** 2
    n = len(X)
    A = coo_matrix((np.r_[w, w], (np.r_[pairs[:, 0], pairs[:, 1]],
                                  np.r_[pairs[:, 1], pairs[:, 0]])), shape=(n, n)).tocsr()
    tree = cKDTree(X)
    src = int(tree.query(caps[0][1])[1])
    _, pred = dijkstra(A, indices=src, return_predecessors=True)
    h = float(np.median(L))

    paths = []
    for _, c in caps[1:]:
        j = int(tree.query(c)[1])
        idx = []
        while j >= 0 and j != src:
            idx.append(j)
            j = pred[j]
        idx.append(src)
        P = smooth_resample(X[idx[::-1]], h)
        rad = dist[tree.query(P)[1]]
        paths.append((P, rad))
    paths.sort(key=lambda pr: -np.sum(np.linalg.norm(np.diff(pr[0], axis=0), axis=1)))

    vessels = [{"name": "Main vessel", "P": paths[0][0], "r": paths[0][1]}]
    for P, rad in paths[1:]:
        kept = np.vstack([v["P"] for v in vessels])
        krad = np.concatenate([v["r"] for v in vessels])
        kd, ki = cKDTree(kept).query(P)
        outside = kd > 1.2 * krad[ki]
        if not outside.any():
            continue
        # departure: first index after which the path stays outside
        last_in = np.flatnonzero(~outside)
        j0 = int(last_in[-1] + 1) if len(last_in) else 0
        if len(P) - j0 < 8:
            continue
        vessels.append({"name": None, "P": P[j0:], "r": rad[j0:]})
    main, rest = vessels[0], vessels[1:]
    rest.sort(key=lambda v: -np.sum(np.linalg.norm(np.diff(v["P"], axis=0), axis=1)))
    for k, v in enumerate(rest):
        v["name"] = f"Branch {ROMAN[k] if k < len(ROMAN) else k + 1}"
    vessels = [main] + rest
    if cache:
        np.savez(cache, vessels=np.array(vessels, dtype=object), h=h)
    return vessels, h


def stations(vessel, fractions):
    P = vessel["P"]
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    out = []
    for f in fractions:
        i = int(np.clip(np.searchsorted(s, f * s[-1]), 1, len(P) - 2))
        t = P[i + 1] - P[i - 1]
        t /= np.linalg.norm(t)
        out.append((f, P[i], t, float(vessel["r"][i])))
    return out


# ---------------------------------------------------------------------------
# Sections and rendering
# ---------------------------------------------------------------------------
def section(grid, c, t, radius, owner=None, vessel_id=None):
    sl = grid.slice(normal=t, origin=c)
    if sl.n_points == 0:
        return None
    # A plane near a bifurcation also cuts the neighbouring lumen, which may
    # touch the section of interest: keep the cells whose nearest centerline
    # point belongs to this vessel (and not too far from the station).
    cc = sl.cell_centers().points
    keep = np.linalg.norm(cc - c, axis=1) <= 4.0 * radius
    if owner is not None:
        tree, ids = owner
        keep &= ids[tree.query(cc)[1]] == vessel_id
    if not keep.any():
        return None
    # fraction of the near-plane cells given to other vessels (bifurcation)
    near = np.linalg.norm(cc - c, axis=1) <= 2.0 * radius
    foreign = 1.0 - (keep & near).sum() / max(near.sum(), 1)
    sl = sl.extract_cells(np.flatnonzero(keep)).extract_surface()
    sl.field_data["foreign"] = [foreign]
    try:
        sl = sl.connectivity(extraction_mode="closest", closest_point=c)
    except TypeError:
        sl = sl.connectivity(largest=False).extract_largest()
    return sl if sl.n_points else None


def section_up(view, t):
    up = view.up - np.dot(view.up, t) * t
    if np.linalg.norm(up) < 1e-6:
        up = view.normal - np.dot(view.normal, t) * t
    return up / np.linalg.norm(up)


def render_section(path, sl, c, t, up, scale, clim, size=700):
    p = r3.make_plotter(size, size)
    p.add_mesh(sl, scalars="Umag", cmap="viridis", clim=clim, show_scalar_bar=False,
               smooth_shading=True, ambient=0.6, diffuse=0.4, specular=0.0)
    cam = p.camera
    cam.focal_point = c
    cam.position = c + t * 10 * scale
    cam.up = up
    p.enable_parallel_projection()
    cam.parallel_scale = scale
    p.reset_camera_clipping_range()
    p.screenshot(path, transparent_background=False)
    p.close()


def render_locator(path, wall, view, vessel, sts, width=900, height=620):
    import pyvista as pv
    p = r3.make_plotter(width, height)
    p.add_mesh(wall, color="#cfcfcf", opacity=0.45, smooth_shading=True, ambient=0.35,
               diffuse=0.65, specular=0.0)
    line = pv.lines_from_points(vessel["P"])
    diag = np.linalg.norm(np.ptp(wall.points, axis=0))
    p.add_mesh(line.tube(radius=0.004 * diag), color="#303030")
    for k, (_, c, _, _) in enumerate(sts):
        p.add_mesh(pv.Sphere(radius=0.021 * diag, center=c),
                   color=STATION_COLORS[k % len(STATION_COLORS)], smooth_shading=True)
    view.apply(p, None, None, width / height)
    p.screenshot(path, transparent_background=False)
    p.close()


def crop(paths, pad=10):
    return r3.crop_box([q for q in paths if q], pad)


def assemble_case(rows, out, clim, note):
    """rows: [(name, locator_png, [(png, label, colour), ...]), ...]"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from PIL import Image

    sec_box = crop([q for _, _, secs in rows for q, _, _ in secs])
    loc_box = crop([loc for _, loc, _ in rows])
    nst = max(len(secs) for _, _, secs in rows)
    nr = len(rows)
    cb_w = 0.08
    width = st.DOUBLE
    ratios = [1.7] + [1.0] * nst
    unit = width * (1 - cb_w) / sum(ratios)
    row_h = unit * 1.05 + 4.5 * st.MM
    fig_h = nr * row_h + 6 * st.MM
    fig = plt.figure(figsize=(width, fig_h))
    gs = fig.add_gridspec(nr, nst + 1, width_ratios=ratios, left=0.0, right=1 - cb_w,
                          top=1 - 6 * st.MM / fig_h, bottom=0.03,
                          wspace=0.04, hspace=0.32)
    for r, (name, loc, secs) in enumerate(rows):
        ax = fig.add_subplot(gs[r, 0])
        ax.set_axis_off()
        ax.imshow(np.asarray(Image.open(loc).convert("RGB").crop(loc_box)))
        ax.text(0.0, 1.02, name, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=8, fontweight="bold")
        for k in range(nst):
            ax = fig.add_subplot(gs[r, k + 1])
            ax.set_axis_off()
            if k < len(secs):
                png, label, col = secs[k]
                ax.imshow(np.asarray(Image.open(png).convert("RGB").crop(sec_box)))
                ax.text(0.5, -0.04, label, transform=ax.transAxes, ha="center", va="top",
                        fontsize=7, color=col)
    cax = fig.add_axes([1 - cb_w + 0.02, 0.2, 0.012, 0.6])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=colors.Normalize(*clim),
                                                   cmap="viridis"), cax=cax)
    cb.set_label(r"$|\mathbf{u}|$ [m/s]")
    cb.outline.set_linewidth(0.5)
    if note:
        fig.text(1 - 0.01, 1.0, note, ha="right", va="top", fontsize=6.5, color="0.3")
    st.save(fig, out)
    plt.close(fig)
    print(f"  wrote {out}.pdf")


def assemble_compare(grid, row_labels, col_labels, out, clim, note):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from PIL import Image

    box = crop([q for row in grid for q in row])
    nr, nc = len(grid), len(grid[0])
    cb_w, lab_w = 0.08, 0.06
    width = st.DOUBLE * 0.85
    unit = width * (1 - cb_w - lab_w) / nc
    fig_h = nr * unit * 1.02 + 6 * st.MM
    fig = plt.figure(figsize=(width, fig_h))
    top = 1 - 5 * st.MM / fig_h
    gs = fig.add_gridspec(nr, nc, left=lab_w, right=1 - cb_w, top=top, bottom=0.01,
                          wspace=0.04, hspace=0.06)
    for r in range(nr):
        for c in range(nc):
            ax = fig.add_subplot(gs[r, c])
            ax.set_axis_off()
            if grid[r][c]:
                ax.imshow(np.asarray(Image.open(grid[r][c]).convert("RGB").crop(box)))
            if r == 0:
                ax.set_title(col_labels[c][0], fontsize=7.5, color=col_labels[c][1], pad=3)
            if c == 0:
                ax.text(-0.06, 0.5, row_labels[r], transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=8)
    cax = fig.add_axes([1 - cb_w + 0.02, 0.15, 0.014, top - 0.25])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=colors.Normalize(*clim),
                                                   cmap="viridis"), cax=cax)
    cb.set_label(r"$|\mathbf{u}|$ [m/s]")
    cb.outline.set_linewidth(0.5)
    if note:
        fig.text(1 - 0.01, 0.0, note, ha="right", va="bottom", fontsize=6.5, color="0.3")
    st.save(fig, out)
    plt.close(fig)
    print(f"  wrote {out}.pdf")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results")
    ap.add_argument("figures")
    ap.add_argument("--period", type=float, default=0.85)
    ap.add_argument("--branches", type=int, default=2)
    ap.add_argument("--stations", type=int, default=4)
    ap.add_argument("--time", choices=["diastole", "systole"], default="diastole")
    args = ap.parse_args()

    import pyvista as pv
    pv.OFF_SCREEN = True
    warnings.filterwarnings("ignore")
    st.apply()

    found = r3.discover(args.results)
    if not found:
        sys.exit(f"No {r3.PREFIX}fluid.xdmf + CSV under {args.results}/<geometry>/<scenario>/")
    panels = os.path.join(args.figures, "panels_sections")
    os.makedirs(panels, exist_ok=True)
    fr = np.linspace(0.1, 0.9, args.stations)

    data = {}   # (geom, scen) -> list of rows (name, stations, sections)
    walls, views, vessels_of, owners = {}, {}, {}, {}
    umax = []
    for geom, scen in found:
        d = os.path.join(args.results, geom, scen)
        csv = r3.read_csv(os.path.join(d, r3.PREFIX + ".csv"))
        t0, t1, t_sys, t_dia = r3.cycle_times(csv, args.period)
        fl = r3.XdmfSeries(os.path.join(d, r3.PREFIX + "fluid.xdmf"))
        s = fl.nearest_step(t_dia if args.time == "diastole" else t_sys)
        print(f"  [{geom}/{scen}] t = {fl.times[s]:.3f} s ({args.time})", flush=True)
        if geom not in vessels_of:
            vs, _ = centerlines(fl, len(fl.times) - 1,
                                os.path.join(args.figures, f"centerlines_{geom}.npz"))
            vessels_of[geom] = vs[:1 + args.branches]
            from scipy.spatial import cKDTree
            owners[geom] = (cKDTree(np.vstack([v["P"] for v in vs])),
                            np.concatenate([np.full(len(v["P"]), k) for k, v in enumerate(vs)]))
            wall, caps = wall_and_caps(fl, len(fl.times) - 1)
            walls[geom] = wall
            views[geom] = r3.View(fl.points, caps[0][1])
            print(f"  [{geom}] centerlines: " + ", ".join(
                f"{v['name']} ({1e3 * np.sum(np.linalg.norm(np.diff(v['P'], axis=0), axis=1)):.1f} mm)"
                for v in vessels_of[geom]), flush=True)
        g = fl.grid.copy()
        u = fl.field(s, "FluidVelocity")
        g.point_data["U"] = u
        g.point_data["Umag"] = np.linalg.norm(u, axis=1)
        rows = []
        for vid, v in enumerate(vessels_of[geom]):
            sts = stations(v, fr)
            secs = []
            for (f, c, t, rad) in sts:
                # slide off bifurcations: first clean cut within +-9% of L
                best = None
                for df in (0.0, 0.03, -0.03, 0.06, -0.06, 0.09, -0.09):
                    ff = float(np.clip(f + df, 0.02, 0.98))
                    (_, cc_, tt_, rr_), = stations(v, [ff])
                    sl = section(g, cc_, tt_, rr_, owners[geom], vid)
                    if sl is None:
                        continue
                    fo = float(sl.field_data["foreign"][0])
                    if best is None or fo < best[0]:
                        best = (fo, (ff, cc_, tt_, rr_, sl))
                    if fo < 0.02:
                        break
                ff, c, t, rad, sl = best[1] if best else (f, c, t, rad, None)
                secs.append((ff, c, t, rad, sl))
                if sl is not None:
                    umax.append(np.percentile(sl.point_data["Umag"], 99.5))
            rows.append((v, sts, secs))
        data[(geom, scen)] = (rows, fl.times[s])

    clim = (0.0, float(max(umax)) if umax else 1.0)
    geoms = list(dict.fromkeys(gm for gm, _ in found))
    for geom in geoms:
        # common scale of all the sections of the geometry (calibres compare)
        rad = [sec[3] for (gm, _), (rows, _) in data.items() if gm == geom
               for _, _, secs in rows for sec in secs]
        scale = 1.25 * max(rad)
        view = views[geom]
        locs = {}
        first = next(rows for (gm, _), (rows, _) in data.items() if gm == geom)
        for v, _, secs in first:
            path = os.path.join(panels, f"locator_{geom}_{v['name'].replace(' ', '_')}.png")
            render_locator(path, walls[geom], view, v,
                           [(f, c, t, r) for (f, c, t, r, _) in secs])
            locs[v["name"]] = path
        for (gm, scen), (rows, tsnap) in data.items():
            if gm != geom:
                continue
            fig_rows = []
            for v, sts, secs in rows:
                cells = []
                for k, (f, c, t, _, sl) in enumerate(secs):
                    png = os.path.join(panels, f"sec_{geom}_{scen}_{v['name'].replace(' ', '_')}_{k}.png")
                    if sl is not None:
                        render_section(png, sl, c, t, section_up(view, t), scale, clim)
                    else:
                        png = None
                    cells.append((png, f"{k + 1}   $s/L = {f:.2f}$",
                                  STATION_COLORS[k % len(STATION_COLORS)]))
                fig_rows.append((v["name"], locs[v["name"]], [q for q in cells if q[0]]))
            assemble_case(fig_rows, os.path.join(args.figures, f"fig3d_sections_{geom}_{scen}"),
                          clim, f"{st.geometry_label(geom)}, {st.scenario_label(scen)}: "
                                f"velocity magnitude at peak {args.time} ($t = {tsnap:.3f}$ s)")
        # scenarios side by side, main vessel
        scens = [sc for gm, sc in found if gm == geom]
        if len(scens) > 1:
            grid = []
            for sc in scens:
                rows, _ = data[(geom, sc)]
                grid.append([os.path.join(panels, f"sec_{geom}_{sc}_Main_vessel_{k}.png")
                             if rows[0][2][k][4] is not None else None
                             for k in range(len(fr))])
            secs0 = data[(geom, scens[0])][0][0][2]
            cols = [(f"{k + 1}   $s/L = {sec[0]:.2f}$", STATION_COLORS[k % len(STATION_COLORS)])
                    for k, sec in enumerate(secs0)]
            assemble_compare(grid, [st.scenario_label(sc) for sc in scens], cols,
                             os.path.join(args.figures, f"fig3d_sections_main_{geom}"), clim,
                             f"{st.geometry_label(geom)}, main vessel, peak {args.time}")


if __name__ == "__main__":
    main()
