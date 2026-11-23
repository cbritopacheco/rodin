"""Validate the quadratic Y velocity field and sample it for volume rendering."""

import argparse
from hashlib import sha256
from pathlib import Path

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.util.vtkConstants import VTK_QUADRATIC_TETRA
from vtkmodules.vtkCommonCore import vtkSMPTools
from vtkmodules.vtkFiltersCore import vtkResampleToImage, vtkStaticCleanUnstructuredGrid
from vtkmodules.vtkIOLegacy import vtkUnstructuredGridReader
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader, vtkXMLUnstructuredGridWriter


EDGE_NODES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]


def read_block(lines, name, dtype):
    start = lines.index(name) + 1
    return np.loadtxt(lines[start + 1:start + 1 + int(lines[start])], dtype=dtype)


def validate(field, mesh, series):
    seconds = field.GetFieldData().GetArray("time_s").GetValue(0)
    pressure = field.GetFieldData().GetArray("inlet_pressure_pa").GetValue(0)
    matched = series[np.isclose(series["time_s"], seconds, rtol=0, atol=1e-9)]
    if len(matched) != 1 or not np.isclose(matched[0]["inlet_pressure_pa"], pressure):
        raise ValueError("Field time and imposed pressure must match a single published flow sample")
    sample = matched[0]
    lines = mesh.decode().splitlines()
    tets = read_block(lines, "Tetrahedra", int)[:, :4] - 1
    vertices = read_block(lines, "Vertices", float)[:, :3]
    boundary = read_block(lines, "Triangles", int)
    cell_nodes = field.GetCells()
    if (field.GetNumberOfCells() != len(tets)
            or not np.all(np.diff(vtk_to_numpy(cell_nodes.GetOffsetsArray())) == 10)
            or any(field.GetCellType(i) != VTK_QUADRATIC_TETRA for i in range(len(tets)))):
        raise ValueError("Expected the ten quadratic velocity samples of every refined tetrahedron")
    local_nodes = vtk_to_numpy(cell_nodes.GetConnectivityArray()).reshape(-1, 10)
    coords = vtk_to_numpy(field.GetPoints().GetData())[local_nodes] / 1e-4
    velocity = vtk_to_numpy(field.GetPointData().GetVectors())[local_nodes]
    if not np.isfinite(velocity).all():
        raise ValueError("Nonfinite velocity field")
    ids = {tuple(point): i for i, point in enumerate(vertices.round(8))}
    connectivity = np.array([[ids[tuple(p)] for p in cell] for cell in coords[:, :4].round(8)])
    if set(map(tuple, np.sort(connectivity))) != set(map(tuple, np.sort(tets))):
        raise ValueError("Field tetrahedra differ from the proposed figure mesh")
    edge_velocity = {}
    for k, (i, j) in enumerate(EDGE_NODES, start=4):
        if not np.allclose(coords[:, k], (coords[:, i] + coords[:, j]) / 2, atol=1e-9, rtol=0):
            raise ValueError("Field is missing the quadratic edge samples")
        for edge, value in zip(np.sort(connectivity[:, [i, j]]), velocity[:, k]):
            edge_velocity[tuple(edge)] = value
    walls = boundary[boundary[:, 3] == 2, :3] - 1
    wall_nodes = np.isin(connectivity, np.unique(walls))
    wall_edges = np.sort(walls[:, [(0,1), (1,2), (2,0)]].reshape(-1, 2), axis=1)
    wall_speed = max(np.linalg.norm(velocity[:, :4][wall_nodes], axis=1).max(),
                     max(np.linalg.norm(edge_velocity[tuple(edge)]) for edge in np.unique(wall_edges, axis=0)))
    if wall_speed > 1e-8 * np.linalg.norm(velocity, axis=2).max():
        raise ValueError(f"Nonzero wall velocity in the exported field: {wall_speed} m/s")
    flux_errors = []
    wide, narrow = (5, 6) if np.any(boundary[:, 3] == 5) else (3, 4)
    for tag, quantity, angle in [(1, "inlet_m3_s", 0), (wide, "wide_m3_s", 30), (narrow, "narrow_m3_s", -38)]:
        faces = boundary[boundary[:, 3] == tag, :3] - 1
        xyz = 1e-4 * vertices[faces]
        areas = np.linalg.norm(np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1) / 2
        # For a quadratic triangle the exact area mean is the mean of the three edge values.
        means = np.mean([[edge_velocity[tuple(sorted((face[i], face[(i+1) % 3])))]
                          for i in range(3)] for face in faces], axis=1)
        direction = np.array([np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle)), 0.0])
        flux = areas @ (means @ direction)
        flux_errors.append(abs(flux / sample[quantity] - 1))
    if max(flux_errors) > 1e-7:
        raise ValueError(f"Exported velocities do not reproduce the three solver fluxes: {flux_errors}")
    print(f"Field time {seconds:g} s, head {pressure:g} Pa, maximum relative cap flux error {max(flux_errors):.3g}")
    return seconds, pressure, int(sample["step"])


def sample_volume(field, mesh_digest, seconds, pressure, step, output):
    # The image plane uses the same orthographic basis as the exterior mesh.
    vtkSMPTools.Initialize(1)  # Stable cell selection for samples near shared faces.
    eye = np.array([0.36, -0.04, 0.93])
    eye /= np.linalg.norm(eye)
    horizontal = np.array([1., 0., 0.]) - eye[0] * eye
    horizontal /= np.linalg.norm(horizontal)
    basis = np.column_stack((horizontal, np.cross(eye, horizontal), eye))
    original = vtk_to_numpy(field.GetPoints().GetData())
    coords = (original / 1e-4) @ basis
    field.GetPoints().SetData(numpy_to_vtk(np.ascontiguousarray(coords), deep=1))
    field.Modified()
    lo, hi = coords.min(axis=0) - 0.05, coords.max(axis=0) + 0.05
    dims = np.ceil((hi - lo) / [0.018, 0.018, 0.04]).astype(int) + 1
    resample = vtkResampleToImage()
    resample.SetInputDataObject(field)
    resample.SetUseInputBounds(False)
    resample.SetSamplingBounds(*np.column_stack((lo, hi)).reshape(-1))
    resample.SetSamplingDimensions(*dims)
    resample.Update()
    point_data = resample.GetOutput().GetPointData()
    vector = vtk_to_numpy(point_data.GetArray("velocity_m_s"))
    mask = vtk_to_numpy(point_data.GetArray("vtkValidPointMask")).astype(bool)
    speed = np.full(len(mask), np.nan, dtype=np.float32)
    speed[mask] = np.linalg.norm(vector[mask], axis=1) * 1000
    if not mask.any() or not np.isfinite(speed[mask]).all() or np.max(speed[mask]) <= 0:
        raise ValueError("VTK volume sampling did not produce a physical velocity field")
    np.savez_compressed(output, speed_mm_s=speed.reshape(tuple(dims[::-1])), lower=lo, upper=hi,
                        mesh_sha256=mesh_digest, time_s=seconds, inlet_pressure_pa=pressure,
                        source_step=step, basis=basis)
    print(f"Sampled {mask.sum()} interior points; maximum speed {np.max(speed[mask]):.5g} mm/s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field", required=True, type=Path)
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument("--volume", required=True, type=Path)
    parser.add_argument("--vtu", type=Path, help="Write a compressed ParaView field in original SI coordinates")
    args = parser.parse_args()
    reader = vtkXMLUnstructuredGridReader() if args.field.suffix == ".vtu" else vtkUnstructuredGridReader()
    reader.SetFileName(str(args.field))
    if isinstance(reader, vtkUnstructuredGridReader):
        reader.ReadAllVectorsOn()
    reader.Update()
    field = reader.GetOutput()
    mesh = args.mesh.read_bytes()
    trace = np.genfromtxt(args.trace, delimiter=",", names=True, ndmin=1, dtype=None, encoding="utf-8")
    if "closure" in trace.dtype.names:
        trace = trace[(trace["full_mesh_sha256"] == sha256(mesh).hexdigest()) & (trace["closure"] == "full")]
        if not len(trace):
            raise ValueError("No resolved reference for the supplied mesh in the combined CSV")
        trace = trace[(trace["steps_per_cycle"] == max(trace["steps_per_cycle"]))
                      & (trace["mean_head_pa"] == min(trace["mean_head_pa"]))]
    seconds, pressure, step = validate(field, mesh, trace)
    if args.vtu:
        # Merge equivalent corner and edge nodes for a connected ParaView mesh.
        clean = vtkStaticCleanUnstructuredGrid()
        clean.SetInputData(field)
        clean.ToleranceIsAbsoluteOn()
        clean.SetAbsoluteTolerance(1e-14)
        clean.Update()
        connected = clean.GetOutput()
        connected.GetFieldData().ShallowCopy(field.GetFieldData())
        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(str(args.vtu))
        writer.SetInputData(connected)
        if writer.Write() != 1:
            raise OSError(f"Could not export {args.vtu}")
        field = connected
    sample_volume(field, sha256(mesh).hexdigest(), seconds, pressure, step, args.volume)
