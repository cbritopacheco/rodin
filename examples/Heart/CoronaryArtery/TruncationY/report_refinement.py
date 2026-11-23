"""Report branch-flow redistribution in the circular and planar Y benchmarks."""

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection, TriMesh
from matplotlib.colors import Normalize, to_rgb
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from matplotlib.tri import Triangulation
import numpy as np

from run_benchmark import CASES, COLORS, LABELS, STYLES


def mesh_block(lines, name):
    start = lines.index(name) + 1
    return lines[start + 1:start + 1 + int(lines[start])]


def collect(runs, destination):
    rows = []
    calibrations = set()
    for directory in runs:
        protocol = json.loads((directory / "protocol.json").read_text())
        if protocol["dimension"] != 3:
            raise ValueError(f"Expected a 3D campaign in {directory}")
        calibrations.add(tuple(np.loadtxt(directory / "calibration.txt")))
        mesh_dir = directory / "mesh"
        if not mesh_dir.is_dir():
            mesh_dir = Path(protocol["mesh"])
        mesh = (mesh_dir / "full.mesh").read_bytes()
        cells = len(mesh_block(mesh.decode().splitlines(), "Tetrahedra"))
        mesh_digest = sha256(mesh).hexdigest()
        if protocol.get("full_mesh_sha256", mesh_digest) != mesh_digest:
            raise ValueError(f"Full mesh differs from its recorded solver input in {directory}")
        steps = protocol["steps_per_cycle"]
        for condition, head in protocol["mean_heads_pa"].items():
            for closure in CASES:
                with (directory / condition / f"{closure}.csv").open() as stream:
                    trace = list(csv.DictReader(stream))
                if len(trace) != steps * protocol["cycles"]:
                    raise ValueError(f"Incomplete {condition}/{closure} in {directory}")
                for i, sample in enumerate(trace, 1):
                    if int(sample["step"]) != i or not np.isclose(float(sample["time_s"]), i * protocol["period_s"] / steps):
                        raise ValueError(f"Unexpected time grid in {directory}/{condition}/{closure}")
                    rows.append(dict(full_tets=cells, full_mesh_sha256=mesh_digest,
                                     steps_per_cycle=steps, mean_head_pa=head,
                                     closure=closure, **sample))
    if len(calibrations) != 1:
        raise ValueError("The outlet calibration differs between runs")
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def project_volume(volume, ticks=None, *, opacity_gain=4.5):
    """Composite orthographic samples of the interpolated finite-element velocity."""
    speed = volume["speed_mm_s"]
    lo, hi = volume["lower"], volume["upper"]
    if ticks is None:
        ticks = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]).tick_values(0, np.nanmax(speed))
    scale = plt.cm.ScalarMappable(norm=Normalize(0, ticks[-1]), cmap="viridis")
    image = np.ones((*speed.shape[1:], 3), dtype=np.float32)
    dz = (hi[2] - lo[2]) / (len(speed) - 1)
    for layer in speed:  # Back to front; zero speed has zero opacity.
        inside = np.isfinite(layer)
        fraction = scale.norm(layer[inside])
        alpha = -np.expm1(-opacity_gain * fraction**1.05 * dz)[:, None]
        image[inside] = (1 - alpha) * image[inside] + alpha * scale.cmap(fraction)[:, :3]
    dx, dy = (hi[:2] - lo[:2]) / (np.array(speed.shape[:0:-1]) - 1)
    extent = (lo[0] - dx/2, hi[0] + dx/2, lo[1] - dy/2, hi[1] + dy/2)
    return image, extent, scale, ticks


def draw_geometry(axis, lines, volume_image=None, *, schematic=False):
    nodes = np.loadtxt(mesh_block(lines, "Vertices"))[:, :3]
    tets = np.loadtxt(mesh_block(lines, "Tetrahedra"), dtype=int)[:, :4] - 1
    triangles = np.loadtxt(mesh_block(lines, "Triangles"), dtype=int)
    # In the full mesh tags 3 and 4 are interior. They close the truncated mesh.
    full = np.any(triangles[:, 3] == 5)
    exterior = triangles[~np.isin(triangles[:, 3], [3, 4])] if full else triangles
    indices = exterior[:, :3] - 1
    faces = nodes[indices]
    eye = np.array([.6, 0., .8]) if schematic else np.array([0.36, -0.04, 0.93])
    eye /= np.linalg.norm(eye)
    horizontal = np.array([1, 0, 0]) - eye[0] * eye
    horizontal /= np.linalg.norm(horizontal)
    view = np.column_stack((horizontal, np.cross(eye, horizontal)))
    centres = faces.mean(axis=1)
    directions = np.deg2rad([30, -38])
    downstream = (centres[:, :2] @ np.array([np.cos(directions), np.sin(directions)])).max(axis=1) > 3
    colors = np.where(downstream[:, None], to_rgb("#ebbf91"), to_rgb("#82adbb"))
    if schematic:
        colors = np.where((full & downstream)[:, None], to_rgb("#e9b180"), to_rgb("#81afbf"))
        colors[exterior[:, 3] != 2] *= .70
    # Orient boundary normals out of their owning tetrahedron. Average the
    # wall normals for lighting only; retain all mesh vertices and triangles.
    opposite = {}
    for i in range(4):
        for face, node in zip(np.sort(np.delete(tets, i, axis=1), axis=1), tets[:, i]):
            opposite[tuple(face)] = node
    inside = nodes[[opposite[tuple(face)] for face in np.sort(indices, axis=1)]]
    areas = np.cross(faces[:, 1] - faces[:, 0], faces[:, 2] - faces[:, 0])
    areas[np.einsum("ij,ij->i", areas, inside - centres) > 0] *= -1
    normals = areas / np.linalg.norm(areas, axis=1)[:, None]
    wall = exterior[:, 3] == 2
    smooth = np.zeros_like(nodes)
    for corner in range(3):
        np.add.at(smooth, indices[wall, corner], areas[wall])
    smooth /= np.maximum(np.linalg.norm(smooth, axis=1)[:, None], 1e-12)
    corner_normals = np.where(wall[:, None, None], smooth[indices], normals[:, None, :])
    light = 0.72 * eye + 0.67 * view[:, 1] - 0.16 * horizontal
    illumination = np.maximum(0, corner_normals @ light)[:, :, None]
    if schematic:
        shine = np.maximum(0, corner_normals @ ((light + eye) / np.linalg.norm(light + eye))) ** 11
        paints = np.clip(colors[:, None, :] * (.52 + .48 * illumination) + .17 * shine[:, :, None], 0, 1)
    else:
        paints = colors[:, None, :] * (.71 + .29 * illumination)
    front = np.flatnonzero(normals @ eye > 0)
    order = front[np.argsort(centres[front] @ eye)]
    projected = faces[order] @ view
    points = projected.reshape(-1, 2)
    triangles_2d = Triangulation(*points.T, triangles=np.arange(len(points)).reshape(-1, 3))
    if volume_image is None:
        axis.add_collection(TriMesh(triangles_2d, facecolors=paints[order].reshape(-1, 3)))
    if schematic:
        for tag in (3, 4, 5, 6):
            cap = triangles[triangles[:, 3] == tag, :3] - 1
            if not len(cap):
                continue
            edges = np.sort(cap[:, [(0, 1), (1, 2), (2, 0)]].reshape(-1, 2), axis=1)
            edges, counts = np.unique(edges, axis=0, return_counts=True)
            rim = nodes[np.unique(edges[counts == 1])]
            center = rim.mean(axis=0)
            # Sort the projected rim into one closed path so dashes continue across edges.
            projected_rim = (rim - center) @ view
            angles = np.arctan2(projected_rim[:, 1], projected_rim[:, 0])
            loop = np.argsort(angles)
            ring = projected_rim[np.r_[loop, loop[0]]]
            if tag in (5, 6) or not full:
                axis.plot(*(ring + center @ view).T, color="#45555b", lw=.48)
            if tag in (3, 4):
                axis.plot(*(ring * 1.28 + center @ view).T, color="#995027",
                          lw=.83, dashes=(2.4, 1.4))
    elif volume_image is None:
        axis.text(3.1, 3.5, r"$72\,\mu$m", ha="center")
        axis.text(2.8, -5.05, r"$46\,\mu$m", ha="center")
        axis.annotate("", xy=(-4.05, 0), xytext=(-5.9, 0),
                      arrowprops=dict(arrowstyle="-|>", lw=1.1, color="#232323"))
    else:
        image, extent = volume_image
        axis.imshow(image, extent=extent, origin="lower", interpolation="bilinear", zorder=0)
    if not schematic:
        axis.add_collection(PolyCollection(projected, facecolors="none",
                                           edgecolors=(0.14, 0.22, 0.25, 0.20), linewidths=0.18))
    edges = indices[front][:, [(0, 1), (1, 2), (2, 0)]].reshape(-1, 2)
    unique, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    # In a schematic, draw the outline underneath the filled surface so the
    # far-side creases at the bifurcation do not show through the near wall.
    axis.add_collection(LineCollection(nodes[unique[counts == 1]] @ view, colors="#485057",
                                       linewidths=.70 if schematic else .48, zorder=0 if schematic else 2))
    axis.set(xlim=(-6, 7.4), ylim=(-5.8, 4.6), aspect="equal")
    axis.set_axis_off()


def analyse(data):
    last_cycles = {}
    metrics = []
    planar = "full_triangles" in data.dtype.names
    cells_column = "full_triangles" if planar else "full_tets"
    units = "m2_s" if planar else "m3_s"
    for cells, steps in sorted(set(zip(data[cells_column], data["steps_per_cycle"]))):
        grid = data[(data[cells_column] == cells) & (data["steps_per_cycle"] == steps)]
        for head in sorted(set(grid["mean_head_pa"]), reverse=True):
            reference = None
            for closure in CASES:
                run = grid[(grid["mean_head_pa"] == head) & (grid["closure"] == closure)]
                if len(run) != 2 * steps or not np.array_equal(run["step"], np.arange(1, 2 * steps + 1)):
                    raise ValueError(f"Missing or repeated steps for {cells}/{steps}/{head}/{closure}")
                flows = np.column_stack((run[f"wide_{units}"], run[f"narrow_{units}"]))
                previous, last = np.split(flows, 2)
                if not np.isfinite(flows).all() or min(last.sum(axis=1)) <= 0:
                    raise ValueError("Invalid branch flows")
                if closure == "full":
                    reference = last
                share = 100 * last[:, 1].sum() / last.sum()
                error = share - 100 * reference[:, 1].sum() / reference.sum()
                wave_error = 100 * np.sqrt(np.mean((last[:, 1] - reference[:, 1])**2)) / np.mean(abs(reference[:, 1]))
                repeat = 100 * np.linalg.norm(last - previous, axis=0) / np.linalg.norm(last, axis=0)
                balance = 100 * np.max(abs(last.sum(axis=1) - run[f"inlet_{units}"][-steps:])) / np.max(last.sum(axis=1))
                metrics.append((cells, steps, head, closure, share, error, wave_error, repeat.max(), balance))
                last_cycles[cells, steps, head, closure] = last
    return last_cycles, metrics


def draw_volume_detail(lines, raster, extent, scale, ticks, volume, detail):
    fig = plt.figure(figsize=(6, 5.3))
    fig.suptitle("Flow inside the fully resolved Y", y=0.985, fontsize=14)
    fig.text(0.5, 0.922, f"Peak inlet pressure, second cycle   ·   "
             f"{float(volume['inlet_pressure_pa']):.3f} Pa   ·   {float(volume['time_s']):.4f} s",
             ha="center", fontsize=9)
    view = fig.add_axes([0.02, 0.00, 0.96, 0.75])
    draw_geometry(view, lines, (raster, extent))
    color_ax = fig.add_axes([0.29, 0.80, 0.42, 0.019])
    fig.colorbar(scale, cax=color_ax, orientation="horizontal", ticks=ticks)
    color_ax.set_title(r"Speed (mm s$^{-1}$)", pad=5)
    fig.savefig(detail, dpi=250, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def set_figure_fonts():
    plt.rcParams.update({"font.size": 10, "font.family": "serif", "text.usetex": True,
                         "font.serif": ["Computer Modern Roman"],
                         "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}",
                         "axes.spines.top": False, "axes.spines.right": False})


def draw_matched_fields(fig, volume, cut_mesh, data):
    cases = ("full", "variable", "constant")
    cases_text = ("(a) Fully resolved", "(b) Shear-dependent", "(c) Constant")
    mesh_bytes = cut_mesh.read_bytes()
    if list(volume["cases"]) != list(cases) or volume["cut_mesh_sha256"].item() != sha256(mesh_bytes).hexdigest():
        raise ValueError("Matched volume cases or retained mesh do not match the figure")
    finest = data[(data["full_tets"] == max(data["full_tets"])) &
                  (data["steps_per_cycle"] == max(data["steps_per_cycle"]))]
    if set(finest["full_mesh_sha256"]) != {volume["full_mesh_sha256"].item()}:
        raise ValueError("The velocity snapshots and cycle-volume shares have different reference meshes")
    steps = int(max(finest["steps_per_cycle"]))
    head = min(finest["mean_head_pa"])
    peak = finest[(finest["mean_head_pa"] == head) & (finest["step"] == steps + steps//4)]
    if (len(peak) != len(cases) or volume["source_step"] != steps + steps//4 or
            not np.isclose(volume["mean_head_pa"], head) or
            not np.allclose(peak["time_s"], volume["time_s"]) or
            not np.allclose(peak["inlet_pressure_pa"], volume["inlet_pressure_pa"])):
        raise ValueError("Expected three simultaneous snapshots at the second-cycle low-load pressure peak")
    lines = mesh_bytes.decode().splitlines()
    valid = volume["valid"]
    lower, upper = volume["lower"], volume["upper"]
    fields = [volume[case] for case in cases]
    ticks = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]).tick_values(0, max(field.max() for field in fields))
    for i, (case, label, field) in enumerate(zip(cases, cases_text, fields)):
        flux = peak["narrow_m3_s"][peak["closure"] == case].item()
        if not np.isclose(flux, volume["narrow_m3_s"][i], atol=0., rtol=1e-7):
            raise ValueError(f"The snapshot and plotted narrow-branch flux differ for {case}")
        speed = np.full(valid.shape, np.nan, dtype=np.float32)
        speed[valid] = field
        raster, extent, scale, _ = project_volume(dict(speed_mm_s=speed, lower=lower, upper=upper),
                                                   ticks, opacity_gain=8)
        x = .016 + .329*i
        ax = fig.add_axes([x, .594, .31, .225])
        draw_geometry(ax, lines, (raster, extent))
        ax.set(xlim=(lower[0]-.06, upper[0]+.06), ylim=(lower[1]-.07, upper[1]+.07))
        ax.set_title(label, fontsize=10, pad=5)
        fig.text(x+.155, .558, rf"$Q_{{\rm n}}={1e12*flux:.3f}\,\mathrm{{nL\,s^{{-1}}}}$",
                 ha="center", va="center", fontsize=9)
    color_ax = fig.add_axes([.43, .944, .33, .018])
    fig.colorbar(scale, cax=color_ax, orientation="horizontal", ticks=ticks)
    fig.text(.40, .951, r"Speed (mm s$^{-1}$)", ha="right", va="center")
    color_ax.tick_params(labelsize=9, length=2, pad=2)


def draw_redistribution(circular, planar, output, comparison, cut_mesh, data):
    """Compare the change in narrow-branch volume share on a common signed scale."""
    rows = ("full", "variable", "constant")
    names = ("Fully resolved", "Shear-dependent", "Constant")
    colors = dict(zip(CASES, COLORS))
    changes = []
    pressures = None
    for metrics in (circular, planar):
        grid = max((row[0], row[1]) for row in metrics)
        shares = {(head, case): share for cells, steps, head, case, share, *_ in metrics
                  if (cells, steps) == grid}
        heads = sorted({head for head, _ in shares}, reverse=True)
        if len(heads) != 2 or (pressures is not None and heads != pressures):
            raise ValueError("The two geometries must share the same two mean pressures")
        pressures = heads
        high, low = pressures
        changes.append([shares[low, case] - shares[high, case] for case in rows])

    set_figure_fonts()
    fig = plt.figure(figsize=(7.0, 4.86))
    draw_matched_fields(fig, comparison, cut_mesh, data)
    fig.text(.50, .477, f"Mean inlet pressure lowered from {high:g} to {low:g} Pa",
             ha="center", va="center")
    lower = min(0., np.min(changes))
    upper = max(0., np.max(changes))
    padding = .27 * (upper-lower) if upper > lower else .27
    limits = (lower-padding, upper+padding)
    ticks = MaxNLocator(nbins=4, steps=[1, 2, 5, 10]).tick_values(*limits)
    ticks = [tick for tick in ticks if limits[0] <= tick <= limits[1]]
    tick_labels = [f"${tick:g}$" for tick in ticks]
    titles = ("(d) Circular tubes (3D)", "(e) Planar channels (2D)")
    for panel, values in enumerate(changes):
        ax = fig.add_axes((.17 + .43*panel, .112, .35, .228))
        ax.set_title(titles[panel], pad=11, fontsize=10)
        ax.axvline(0, color="#67727a", linewidth=.9, zorder=1)
        for y, case, change in zip((2, 1, 0), rows, values):
            ax.hlines(y, 0, change, color=colors[case], linewidth=2.4)
            ax.plot(change, y, "o", color=colors[case], markersize=5)
            offset = -5 if change < 0 else 5
            ax.annotate(f"${change:+.2f}$", (change, y), xytext=(offset, 0),
                        textcoords="offset points", ha="right" if change < 0 else "left",
                        va="center", fontsize=9)
        ax.set_xticks(ticks, tick_labels)
        ax.set(xlim=limits, ylim=(-.40, 2.45), yticks=(2, 1, 0))
        ax.set_yticklabels(names if panel == 0 else [])
        ax.tick_params(axis="y", length=0, pad=4)
        ax.tick_params(axis="x", length=3, pad=4)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_color("#829095")
        ax.grid(axis="x", color="#e2e6e6", linewidth=.55)
        ax.set_axisbelow(True)
    fig.text(.55, .025, "Change in narrow-branch share of cycle volume (percentage points)",
             ha="center", va="center")
    fig.savefig(output, dpi=310, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)


def draw_diagnostics(waves, data, mesh, volume, output, detail=None):
    grid = max((cells, steps) for cells, steps, _, _ in waves)
    mesh_bytes = mesh.read_bytes()
    finest = data[(data["full_tets"] == grid[0]) & (data["steps_per_cycle"] == grid[1])]
    digest = sha256(mesh_bytes).hexdigest()
    if set(finest["full_mesh_sha256"]) != {digest} or volume["mesh_sha256"].item() != digest:
        raise ValueError("The mesh, volume and finest waveforms must use the same full mesh")
    heads = sorted({head for cells, steps, head, _ in waves if (cells, steps) == grid}, reverse=True)
    if len(heads) != 2:
        raise ValueError("The figure compares two pressure heads")
    reference = finest[(finest["closure"] == "full") & (finest["mean_head_pa"] == heads[-1])]
    peak = reference[reference["step"] == grid[1] + grid[1]//4][0]
    if (volume["source_step"] != peak["step"] or not np.isclose(volume["time_s"], peak["time_s"])
            or not np.isclose(volume["inlet_pressure_pa"], peak["inlet_pressure_pa"])):
        raise ValueError("The volume must come from peak inlet pressure in the second low-load cycle")
    set_figure_fonts()
    fig = plt.figure(figsize=(7.0, 5.65))
    upper = [fig.add_axes([x, 0.557, 0.41, 0.31]) for x in [0.06, 0.55]]
    axes = [fig.add_axes([x, 0.09, 0.37, 0.245]) for x in [0.105, 0.60]]
    lines = mesh_bytes.decode().splitlines()
    raster, extent, scale, ticks = project_volume(volume)
    draw_geometry(upper[0], lines)
    draw_geometry(upper[1], lines, (raster, extent))
    upper[0].set_title("(a) Refined Y mesh", fontsize=10, pad=7)
    upper[1].set_title("(b) Resolved fluid velocity", fontsize=10, pad=7)
    fig.legend(handles=[Patch(facecolor="#82adbb", label="Retained vessel"),
                        Patch(facecolor="#ebbf91", label="Distal extensions")],
               loc="center", bbox_to_anchor=(0.265, 0.957), ncols=2,
               handlelength=1, columnspacing=1, frameon=False, fontsize=10)
    color_ax = fig.add_axes([0.60, 0.951, 0.31, 0.015])
    fig.colorbar(scale, cax=color_ax, orientation="horizontal", ticks=ticks)
    color_ax.set_title(r"Speed (mm s$^{-1}$)", fontsize=9, pad=3)
    color_ax.tick_params(labelsize=9, length=2, pad=1)
    phases = np.linspace(0, 1, grid[1] + 1)
    for ix, (ax, head) in enumerate(zip(axes, heads)):
        for closure, color, style, label in zip(CASES, COLORS, STYLES, LABELS):
            flow = waves[grid[0], grid[1], head, closure]
            split = 100 * flow[:, 1] / flow.sum(axis=1)
            ax.plot(phases, np.r_[split[-1], split], color=color, ls=style, lw=1.6, label=label)
        ax.set_title(f"({'cd'[ix]}) Mean inlet pressure {head:g} Pa", fontsize=10, pad=8)
        ax.set(xlabel="Fraction of cycle", xlim=(0, 1), xticks=(0, 0.5, 1))
        ax.grid(axis="y", color="#dedede", lw=0.6)
    limits = [ax.get_ylim() for ax in axes]
    for ax in axes:
        ax.set_ylim(min(x[0] for x in limits), max(x[1] for x in limits))
    axes[0].set_ylabel(r"Narrow-branch share (\%)")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="center", bbox_to_anchor=(0.52, 0.459),
               ncols=3, frameon=False, columnspacing=2)
    if output:
        fig.savefig(output, dpi=310, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)
    if detail:
        draw_volume_detail(lines, raster, extent, scale, ticks, volume, detail)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path, help="Combined circular-tube time-series CSV")
    parser.add_argument("--planar-data", type=Path, default=Path(__file__).resolve().parent / "data/truncation_y_2d.csv",
                        help="Combined planar-channel time-series CSV")
    parser.add_argument("--runs", nargs="+", type=Path, help="Rebuild the CSV from complete campaign directories")
    parser.add_argument("--mesh", type=Path, help="Finest full Y mesh for the field and waveform diagnostics")
    parser.add_argument("--volume", type=Path, help="Sampled quadratic field from sample_velocity.py")
    parser.add_argument("--comparison", type=Path, help="Three aligned velocity fields from sample_comparison.py")
    parser.add_argument("--cut-mesh", type=Path, default=Path(__file__).resolve().parent / "data/truncation_y_refined/cut.mesh",
                        help="Retained common mesh for the main figure")
    parser.add_argument("--figure", type=Path, help="Peak speed fields and pressure-induced changes in flow split")
    parser.add_argument("--diagnostics", type=Path, help="Full 3D field, mesh and branch-flow waveforms")
    parser.add_argument("--detail", type=Path, help="Also save the velocity rendering at larger scale")
    args = parser.parse_args()
    if args.figure and not args.comparison:
        parser.error("The main figure requires the three sampled fields (--comparison)")
    if (args.diagnostics or args.detail) and not (args.mesh and args.volume):
        parser.error("The velocity and waveform diagnostics require --mesh and --volume")
    if args.runs:
        collect(args.runs, args.data)
    data = np.genfromtxt(args.data, names=True, delimiter=",", dtype=None, encoding="utf-8")
    waveform, metrics = analyse(data)
    print("full_cells,steps_per_cycle,mean_head_pa,closure,share_pct,error_pp,narrow_rms_error_pct,max_branch_cycle_rms_change_pct,mass_imbalance_pct")
    for result in metrics:
        print(",".join(str(item) for item in result))
    if args.figure:
        planar = np.genfromtxt(args.planar_data, names=True, delimiter=",", dtype=None, encoding="utf-8")
        _, planar_metrics = analyse(planar)
        with np.load(args.comparison) as volume:
            draw_redistribution(metrics, planar_metrics, args.figure, volume, args.cut_mesh, data)
    if args.diagnostics or args.detail:
        with np.load(args.volume) as volume:
            draw_diagnostics(waveform, data, args.mesh, volume, args.diagnostics, args.detail)
