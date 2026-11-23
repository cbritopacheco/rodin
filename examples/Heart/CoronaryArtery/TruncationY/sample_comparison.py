"""Co-register three peak-pressure P2 velocity solutions on the retained Y."""

import argparse
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy
from vtkmodules.util.vtkConstants import VTK_QUADRATIC_TETRA
from vtkmodules.vtkCommonCore import vtkPoints, vtkSMPTools
from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkUnstructuredGrid
from vtkmodules.vtkFiltersCore import vtkResampleToImage
from vtkmodules.vtkIOLegacy import vtkUnstructuredGridReader
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

from sample_velocity import EDGE_NODES, read_block, validate


def retained_geometry(full_mesh, cut_mesh):
    full = full_mesh.decode().splitlines()
    cut = cut_mesh.decode().splitlines()
    nodes = read_block(cut, "Vertices", float)[:, :3]
    tets = read_block(cut, "Tetrahedra", int)[:, :4] - 1
    full_nodes = read_block(full, "Vertices", float)[:, :3]
    full_tets = read_block(full, "Tetrahedra", int)[:, :4] - 1
    distance, ids = cKDTree(full_nodes).query(nodes)
    if distance.max() > 1e-8 or not set(map(tuple, np.sort(ids[tets]))).issubset(set(map(tuple, np.sort(full_tets)))):
        raise ValueError("The retained cells must belong to the refined reference mesh")
    edges, numbering = np.unique(np.sort(tets[:, EDGE_NODES], axis=2).reshape(-1, 2),
                                 axis=0, return_inverse=True)
    coordinates = np.vstack((nodes, nodes[edges].mean(axis=1)))
    connectivity = np.column_stack((tets, len(nodes) + numbering.reshape(-1, 6)))
    # This is the same camera as draw_geometry in the diagnostics figure.
    eye = np.array([.36, -.04, .93])
    eye /= np.linalg.norm(eye)
    right = np.array([1., 0., 0.]) - eye[0] * eye
    right /= np.linalg.norm(right)
    basis = np.column_stack((right, np.cross(eye, right), eye))
    screen = coordinates @ basis
    points = vtkPoints()
    points.SetData(numpy_to_vtk(np.ascontiguousarray(screen), deep=1))
    cells = vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(np.arange(len(tets) + 1, dtype=np.int64) * 10, deep=1),
                  numpy_to_vtkIdTypeArray(connectivity.astype(np.int64).ravel(), deep=1))
    grid = vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.SetCells(VTK_QUADRATIC_TETRA, cells)
    return grid, coordinates * 1e-4, screen.min(axis=0) - .05, screen.max(axis=0) + .05, basis


def sample(args):
    full_mesh = (args.mesh_dir / "full.mesh").read_bytes()
    cut_mesh = (args.mesh_dir / "cut.mesh").read_bytes()
    trace = np.genfromtxt(args.trace, delimiter=",", names=True, dtype=None, encoding="utf-8")
    trace = trace[trace["full_mesh_sha256"] == sha256(full_mesh).hexdigest()]
    if not len(trace):
        raise ValueError("The reference mesh differs from the one recorded in the campaign")
    steps = max(trace["steps_per_cycle"])
    head = min(trace["mean_head_pa"])
    trace = trace[(trace["steps_per_cycle"] == steps) & (trace["mean_head_pa"] == head)]
    grid, coords, lower, upper, basis = retained_geometry(full_mesh, cut_mesh)
    cases = ("full", "variable", "constant")
    narrow_flux = []
    for case in cases:
        path = getattr(args, case)
        reader = vtkXMLUnstructuredGridReader() if path.suffix == ".vtu" else vtkUnstructuredGridReader()
        reader.SetFileName(str(path))
        if isinstance(reader, vtkUnstructuredGridReader):
            reader.ReadAllVectorsOn()
        reader.Update()
        field = reader.GetOutput()
        run = trace[trace["closure"] == case]
        seconds, pressure, step = validate(field, full_mesh if case == "full" else cut_mesh, run)
        if step != steps + steps//4:
            raise ValueError("Expected peak inlet pressure in the second cycle")
        narrow_flux.append(run["narrow_m3_s"][run["step"] == step].item())
        # Copy velocities at coincident P2 nodes. All three solutions are then
        # interpolated on exactly the same retained tetrahedra (no extrapolation).
        distance, ids = cKDTree(vtk_to_numpy(field.GetPoints().GetData())).query(coords)
        if distance.max() > 1e-12:
            raise ValueError(f"{case} is missing velocity nodes on the refined retained mesh")
        vectors = numpy_to_vtk(np.ascontiguousarray(vtk_to_numpy(field.GetPointData().GetVectors())[ids]), deep=1)
        vectors.SetName(case)
        grid.GetPointData().AddArray(vectors)
    vtkSMPTools.Initialize(1)
    dims = np.ceil((upper-lower) / [.018, .018, .04]).astype(int) + 1
    resample = vtkResampleToImage()
    resample.SetInputDataObject(grid)
    resample.SetUseInputBounds(False)
    resample.SetSamplingBounds(*np.column_stack((lower, upper)).ravel())
    resample.SetSamplingDimensions(*dims)
    resample.Update()
    values = resample.GetOutput().GetPointData()
    valid = vtk_to_numpy(values.GetArray("vtkValidPointMask")).astype(bool)
    if not valid.any():
        raise ValueError("No volume samples in retained domain")
    result = {}
    for case in cases:
        vectors = vtk_to_numpy(values.GetArray(case))
        speed = np.linalg.norm(vectors[valid], axis=1) * 1000
        if not np.isfinite(speed).all() or np.max(speed) <= 0:
            raise ValueError(f"Invalid sampled velocity in {case}")
        result[case] = speed.astype(np.float32)
        print(f"{case}: maximum speed {speed.max():.6f} mm/s")
    np.savez_compressed(args.output, **result, valid=valid.reshape(tuple(dims[::-1])),
                        lower=lower, upper=upper, basis=basis, cases=cases,
                        cut_mesh_sha256=sha256(cut_mesh).hexdigest(), full_mesh_sha256=sha256(full_mesh).hexdigest(),
                        time_s=seconds, inlet_pressure_pa=pressure, source_step=step,
                        narrow_m3_s=narrow_flux, mean_head_pa=head)
    print(f"Compared {valid.sum()} common interior positions at {seconds:g} s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-dir", required=True, type=Path)
    parser.add_argument("--trace", required=True, type=Path)
    for name in ("full", "variable", "constant"):
        parser.add_argument("--"+name, required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    sample(parser.parse_args())
