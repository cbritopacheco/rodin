"""Generate conforming full and truncated Y meshes in MEDIT format."""

import argparse
from pathlib import Path

import gmsh
import numpy as np


ANGLES = np.radians([30.0, -38.0])
RADII = (0.72, 0.46)
CUT = 3.0
ENDS = (6.6, 8.2)
INLET = -4.5
DIRECTIONS = np.column_stack((np.cos(ANGLES), np.sin(ANGLES), [0.0, 0.0]))


def geometry(dim):
    occ = gmsh.model.occ

    def tube(direction, start, end, radius):
        if dim == 3:
            return dim, occ.addCylinder(*(start * direction), *((end - start) * direction), radius)
        tag = occ.addRectangle(start, -radius, 0, end - start, 2 * radius)
        occ.rotate([(dim, tag)], 0, 0, 0, 0, 0, 1, np.arctan2(direction[1], direction[0]))
        return dim, tag

    parent = tube(np.array([1.0, 0.0, 0.0]), INLET, 0.0, 1.0)
    junction = occ.addSphere(0, 0, 0, 1.02) if dim == 3 else occ.addDisk(0, 0, 0, 1.02, 1.02)
    daughters = [tube(direction, 0.37, CUT, radius) for direction, radius in zip(DIRECTIONS, RADII)]
    retained, _ = occ.fuse([parent], [(dim, junction), *daughters])
    if len(retained) != 1:
        raise ValueError("The retained Y must form a single domain")
    extensions = [tube(direction, CUT, end, radius) for direction, end, radius in zip(DIRECTIONS, ENDS, RADII)]
    domains, mapping = occ.fragment(retained, extensions)
    occ.synchronize()
    if len(domains) != 3 or any(len(parts) != 1 for parts in mapping):
        raise ValueError("The Y needs a retained domain and two separate extensions")
    return domains, mapping[0][0]


def write_mesh(path, dim, vertices, cells, faces):
    used = np.unique(cells)
    numbering = np.zeros(len(vertices), dtype=int)
    numbering[used] = np.arange(1, len(used) + 1)
    if any(not np.all(numbering[nodes]) for nodes, _ in faces):
        raise ValueError("A tagged face is missing from the selected cells")
    with path.open("w") as stream:
        stream.write(f"MeshVersionFormatted 2\n\nDimension {dim}\n\nVertices\n{len(used)}\n")
        np.savetxt(stream, np.column_stack((vertices[used, :dim], np.zeros(len(used)))), fmt=["%.16g"] * dim + ["%d"])
        sections = ("Tetrahedra", "Triangles") if dim == 3 else ("Triangles", "Edges")
        stream.write(f"\n{sections[0]}\n{len(cells)}\n")
        np.savetxt(stream, np.column_stack((numbering[cells], np.ones(len(cells), int))), fmt="%d")
        stream.write(f"\n{sections[1]}\n{sum(len(nodes) for nodes, _ in faces)}\n")
        for nodes, label in faces:
            np.savetxt(stream, np.column_stack((numbering[nodes], np.full(len(nodes), label))), fmt="%d")
        stream.write("\nEnd\n")


def generate(output, dim, mesh_size):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.model.add("truncation-y")
        domains, retained = geometry(dim)
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.4 * mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)
        gmsh.model.mesh.generate(dim)
        node_tags, positions, _ = gmsh.model.mesh.getNodes()
        vertices = positions.reshape(-1, 3)
        node_index = np.empty(int(node_tags.max()) + 1, dtype=int)
        node_index[node_tags] = np.arange(len(node_tags))
        element_types = {1: 1, 2: 2, 3: 4}  # linear segment, triangle, tetrahedron

        def connectivity(entity):
            dimension, tag = entity
            _, nodes = gmsh.model.mesh.getElementsByType(element_types[dimension], tag)
            return node_index[nodes.reshape(-1, dimension + 1)]

        caps = np.vstack(([INLET, 0, 0], CUT * DIRECTIONS, np.array(ENDS)[:, None] * DIRECTIONS))
        cap_labels = (1, 3, 4, 5, 6)
        face_labels = {}
        for entity in gmsh.model.getEntities(dim - 1):
            center = gmsh.model.occ.getCenterOfMass(*entity)
            distances = np.linalg.norm(caps - center, axis=1)
            nearest = np.argmin(distances)
            face_labels[entity] = cap_labels[nearest] if distances[nearest] < 1e-7 else 2
        if [sum(label == cap for label in face_labels.values()) for cap in cap_labels] != [1] * len(cap_labels):
            raise ValueError(f"Expected five unique cap entities, found {face_labels}")

        all_cells = np.concatenate([connectivity(entity) for entity in domains])
        retained_cells = connectivity(retained)
        faces = [(connectivity(entity), label) for entity, label in face_labels.items()]
        retained_faces = [(connectivity(entity), face_labels[entity]) for entity in gmsh.model.getBoundary([retained], oriented=False)]
        output.mkdir(parents=True, exist_ok=True)
        write_mesh(output / "full.mesh", dim, vertices, all_cells, faces)
        write_mesh(output / "cut.mesh", dim, vertices, retained_cells, retained_faces)
        print(f"{dim}D Y, h_max={mesh_size:g}: {len(all_cells)} full cells, {len(retained_cells)} retained cells")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--dim", type=int, choices=(2, 3), default=3)
    parser.add_argument("--mesh-size", type=float, default=0.46)
    args = parser.parse_args()
    if args.mesh_size <= 0:
        parser.error("--mesh-size must be positive")
    generate(args.output, args.dim, args.mesh_size)
