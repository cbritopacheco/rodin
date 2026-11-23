"""Map matched-minus-variable cycle-mean wall shear stress magnitude from Rodin XDMFs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

import numpy as np


def mesh_data(path):
    # CMBPE's Medit medium mesh is in mm; the solver scales it to m.
    lines = iter(line.strip() for line in path.read_text().splitlines() if line.strip())
    blocks = {}
    for line in lines:
        if line in ("Vertices", "Triangles", "Tetrahedra"):
            count = int(next(lines))
            blocks[line] = np.array([next(lines).split() for _ in range(count)])
    coordinates = blocks["Vertices"][:, :3].astype(float) * 1e-3
    triangles = blocks["Triangles"].astype(np.int64)
    wall = triangles[triangles[:, 3] == 2, :3] - 1
    tetrahedra = blocks["Tetrahedra"][:, :4].astype(np.int64) - 1
    if len(wall) == 0 or len(tetrahedra) == 0:
        raise ValueError("expected tag-2 wall triangles and tetrahedra")
    return coordinates, wall, tetrahedra


def read_item(item, root):
    if item is None or item.get("Format") != "HDF":
        raise ValueError("expected an HDF-backed Rodin XDMF DataItem")
    filename, dataset = item.text.strip().split(":", 1)
    filename = root / filename
    dims = tuple(map(int, item.attrib["Dimensions"].split()))
    try:
        import h5py
    except ImportError:
        dtype = {"Float": "f", "UInt": "u", "Int": "i"}[item.get("NumberType")]
        dtype = np.dtype("<" + dtype + item.attrib["Precision"])
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / "data.bin"
            subprocess.run(["h5dump", "-d", dataset, "-b", "LE", "-o", str(raw), str(filename)],
                           check=True, stdout=subprocess.DEVNULL)
            values = np.fromfile(raw, dtype=dtype)
    else:
        with h5py.File(filename) as h5:
            values = h5[dataset][...]
    return values.reshape(dims)


def xdmf_path(path):
    if path.is_file():
        return path
    candidates = list(path.glob("*.xdmf"))
    if len(candidates) != 1:
        raise ValueError(f"expected one XDMF in {path}, found {candidates}")
    return candidates[0]


def sample_grids(path, times):
    document = ET.parse(path)
    timeline = document.findall(".//Grid[@CollectionType='Temporal']")
    if len(timeline) != 1:
        raise ValueError(f"expected one temporal grid in {path}")
    steps = timeline[0].findall("Grid")
    recorded = np.array([float(step.find("Time").attrib["Value"]) for step in steps])
    sampled = []
    for t in times:
        hits = np.flatnonzero(np.isclose(recorded, t, atol=1e-7, rtol=0))
        if len(hits) != 1:
            raise ValueError(f"{path}: expected one saved field at {t:.6f} s, found {len(hits)}")
        snapshot = steps[hits[0]]
        pieces = snapshot.findall("Grid") if snapshot.get("CollectionType") == "Spatial" else [snapshot]
        if not pieces:
            raise ValueError(f"{path}: empty snapshot at {t}")
        sampled.append(pieces)
    return sampled


def coord_keys(coords):
    return [tuple(row) for row in np.rint(coords * 1e11).astype(np.int64)]


def cells_for_piece(grid, root, ids):
    topology = grid.find("Topology")
    cells = read_item(topology.find("DataItem"), root)
    if topology.get("TopologyType") == "Mixed":
        cells = cells.reshape(-1, 5)
        if not np.all(cells[:, 0] == 6):
            raise ValueError("expected a linear tetrahedral Rodin volume mesh")
        cells = cells[:, 1:]
    elif topology.get("TopologyType") == "Tetrahedron":
        cells = cells.reshape(-1, 4)
    else:
        raise ValueError("expected tetrahedral XDMF topology")
    return ids[cells]


def canonical_cells(rows):
    normalized = np.sort(rows, axis=1)
    return normalized[np.lexsort(normalized.T[::-1])]


def static_key(grid):
    return tuple(grid.find(s).text.strip() for s in ("Geometry/DataItem", "Topology/DataItem"))


def sample_wall(path, coordinates, wall_ids, tetrahedra, times):
    frames = sample_grids(path, times)
    global_by_coordinate = dict(zip(coord_keys(coordinates), range(len(coordinates))))
    if len(global_by_coordinate) != len(coordinates):
        raise ValueError("reference mesh contains spatially ambiguous vertices")
    wall_index = np.full(len(coordinates), -1, dtype=int)
    wall_index[wall_ids] = np.arange(len(wall_ids))
    projection, cells = {}, []
    for grid in frames[0]:
        key = static_key(grid)
        if key in projection:
            raise ValueError("duplicate MPI mesh piece")
        coords = read_item(grid.find("Geometry/DataItem"), path.parent)
        ids = np.array([global_by_coordinate[x] for x in coord_keys(coords)])
        if not np.allclose(coords, coordinates[ids], rtol=0, atol=2e-11):
            raise ValueError("mesh vertices differ from the campaign mesh")
        cells.append(cells_for_piece(grid, path.parent, ids))
        wall_local = np.flatnonzero(wall_index[ids] >= 0)
        projection[key] = (len(coords), wall_local, wall_index[ids[wall_local]])
    if not np.array_equal(canonical_cells(np.concatenate(cells)), canonical_cells(tetrahedra)):
        raise ValueError(f"{path}: tetrahedral connectivity differs from the tagged campaign mesh")
    series_mean = np.zeros(len(wall_ids))
    for frame in frames:
        if {static_key(g) for g in frame} != set(projection) or len(frame) != len(projection):
            raise ValueError("XDMF mesh changed within the cycle")
        sample = np.full(len(wall_ids), np.nan)
        for grid in frame:
            attr = grid.find("Attribute[@Name='shearStress']")
            if attr is None or attr.get("Center") != "Node" or attr.get("AttributeType") != "Vector":
                raise ValueError("missing nodal shearStress vector")
            n, local, global_wall = projection[static_key(grid)]
            vec = read_item(attr.find("DataItem"), path.parent)
            if vec.shape != (n, 3):
                raise ValueError("invalid nodal vector dimensions")
            magnitudes = np.linalg.norm(vec[local], axis=1)
            if not np.all(np.isfinite(magnitudes)):
                raise ValueError("nonfinite wall stress in exported vector")
            duplicated = np.isfinite(sample[global_wall])
            if not np.allclose(sample[global_wall[duplicated]], magnitudes[duplicated],
                               rtol=1e-7, atol=1e-9):
                raise ValueError("MPI copies of a wall vertex give inconsistent stress magnitudes")
            sample[global_wall] = magnitudes
        if not np.all(np.isfinite(sample)):
            raise ValueError("wall has missing vertices in this snapshot")
        series_mean += sample
    return series_mean / len(frames)


def save_surface(path, coords, triangles, fields):
    with path.open("w") as f:
        f.write("# vtk DataFile Version 3.0\nCoronary tag-2 wall, coordinates m, stresses Pa\n"
                "ASCII\nDATASET POLYDATA\n")
        f.write(f"POINTS {len(coords)} double\n")
        np.savetxt(f, coords, fmt="%.14g")
        f.write(f"POLYGONS {len(triangles)} {len(triangles)*4}\n")
        np.savetxt(f, np.column_stack((np.full(len(triangles), 3), triangles)), fmt="%d")
        f.write(f"POINT_DATA {len(coords)}\n")
        for name, values in fields.items():
            f.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
            np.savetxt(f, values, fmt="%.14g")


def render(path, coords, triangles, delta, limit):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    # Equal axes and opposite views; the color on each rendered face is its
    # mean nodal difference. The companion VTK stores the full nodal field.
    center = (coords.max(axis=0) + coords.min(axis=0)) / 2
    r = 0.52 * np.ptp(coords, axis=0).max()
    norm, cmap = Normalize(-limit, limit), plt.get_cmap("RdBu_r")
    fig = plt.figure(figsize=(9, 5))
    for index, azimuth in enumerate([-60, 120]):
        ax = fig.add_axes([.01 + index*.49, .02, .49, .76], projection="3d")
        surf = Poly3DCollection(coords[triangles], linewidths=0, antialiaseds=False)
        surf.set_facecolor(cmap(norm(delta[triangles].mean(axis=1))))
        ax.add_collection3d(surf)
        ax.set(xlim=(center[0]-r,center[0]+r), ylim=(center[1]-r,center[1]+r),
               zlim=(center[2]-r,center[2]+r))
        ax.set_box_aspect((1,1,1)); ax.set_proj_type("ortho")
        ax.view_init(elev=16 if index==0 else -16, azim=azimuth); ax.set_axis_off()
    fig.suptitle("Cycle-mean wall shear stress magnitude", y=.98)
    cax = fig.add_axes([.28,.81,.44,.025])
    bar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap), cax=cax,
                       orientation="horizontal")
    bar.set_label("Matched constant minus shear-dependent outlet (Pa)")
    fig.savefig(path,dpi=200); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["reference","matched","mesh","out"]:
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--cycle", type=int, default=3)
    parser.add_argument("--period", type=float, default=.85)
    parser.add_argument("--sample-dt", type=float, default=.025)
    args = parser.parse_args()
    if args.cycle < 1 or min(args.period,args.sample_dt) <= 0:
        parser.error("cycle, period and sample interval must be positive")
    n = round(args.period/args.sample_dt)
    if n < 2 or not np.isclose(n*args.sample_dt,args.period,rtol=0,atol=1e-10):
        parser.error("the snapshot interval must divide the cycle")
    ref, matched = map(xdmf_path, [args.reference,args.matched])
    if ref.resolve() == matched.resolve():
        parser.error("two distinct simulation directories are required")
    times = (args.cycle-1)*args.period+np.arange(n)*args.sample_dt
    xyz, wall, tets = mesh_data(args.mesh)
    wall_ids, tris = np.unique(wall,return_inverse=True)
    tris = tris.reshape(-1,3)
    a,b = [sample_wall(p,xyz,wall_ids,tets,times) for p in (ref,matched)]
    delta = b-a
    coords = xyz[wall_ids]
    ep1 = coords[tris[:,1]] - coords[tris[:,0]]
    ep2 = coords[tris[:,2]] - coords[tris[:,0]]
    areas = np.linalg.norm(np.cross(ep1,ep2),axis=1)/2
    meta = dict(reference=str(ref.resolve()), matched=str(matched.resolve()),
                mesh_sha256=hashlib.sha256(args.mesh.read_bytes()).hexdigest(), wall_tag=2,
                coordinate_unit="m", stress_unit="Pa",
                statistic="Mean of saved nodal shearStress vector magnitudes on the wall; matched constant minus nonlinear.",
                cycle=args.cycle, period_s=args.period, snapshot_interval_s=args.sample_dt,
                first_snapshot_time_s=float(times[0]), last_snapshot_time_s=float(times[-1]),
                snapshots=len(times), wall_vertices=len(wall_ids), wall_triangles=len(tris),
                min_difference_pa=float(delta.min()), max_difference_pa=float(delta.max()),
                mean_difference_area_weighted_pa=float(np.average(delta[tris].mean(axis=1),weights=areas)),
                color_limit_pa=float(np.max(np.abs(delta))))
    fields = dict(variable_mean_abs_WSS_Pa=a,matched_mean_abs_WSS_Pa=b,matched_minus_variable_Pa=delta)
    args.out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.out) as stage:
        stage = Path(stage)
        save_surface(stage/"wall-wss-difference.vtk",coords,tris,fields)
        render(stage/"wall-wss-difference.png",coords,tris,delta,meta["color_limit_pa"] or 1e-12)
        (stage/"wall-wss-summary.json").write_text(json.dumps(meta,indent=2)+"\n")
        for output in stage.iterdir():
            os.replace(output,args.out/output.name)
    print(json.dumps(meta,indent=2))


if __name__ == "__main__":
    main()
