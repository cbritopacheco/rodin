# Integrations, I/O, and bindings

Everything outside the native Eigen path. All integrations are additive:
they live in their own directories, mirror core interfaces, and core never
includes them.

## PETSc

Full mirror tree (`PETSc/Math`, `PETSc/Solver`, `PETSc/Assembly`,
`PETSc/Variational`, `PETSc/IO`) over PETSc Mat/Vec/KSP
via `Object.h` RAII wrappers. All specifics — error-handling idiom, the
no-resize LinearSystem contract, CI version skew, known traps — in
[petsc.md](petsc.md). 2D examples generally do NOT touch PETSc (native
Eigen); 3D paths do.

## MMG

Wrapper over the MMG remesher: `MMG::Mesh` (Rodin mesh + MMG metadata),
`Adapt` (metric-based adaptation), `MeshOptimizer` (quality optimization),
`LevelSetDiscretizer` (implicit-domain discretization from a level set),
`MMG5.h` conversion utilities, MEDIT loader/printer specializations. On
this branch MMG is the remesher and implicit-domain discretizer; native
alternatives exist only on other branches.

## Scotch

`Scotch::MeshPartitioner` — graph-partitioner implementation of the
`Geometry::MeshPartitioner` interface (used for MPI distribution).

## MPI

Mirrors the local stack per module: `MPI/Context`, `MPI/Geometry`
(distributed mesh/sharding), `MPI/Assembly`, `MPI/Variational`, `MPI/IO`.
The context template parameter (`Context::MPI` vs `Context::Local`) selects
it; no `#ifdef`s in core classes.

## IO (src/Rodin/IO/)

Loader/Printer class-template pattern: `MeshLoader`/`MeshPrinter`,
`GridFunctionLoader`/`GridFunctionPrinter`, specialized per
`IO::FileFormat`:

- **MEDIT** (`.mesh`/`.sol`) — primary mesh interchange (MMG ecosystem).
- **MFEM** (`.mesh`/`.gf`) — meshes and grid functions.
- **XDMF + HDF5** — visualization/time series: `IO::XDMF` writes a domain
  with `xdmf.grid().setMesh(mesh).add("name", gf)` per step; heavy data in
  HDF5 side files; auto-detects curved (order>1) transformations and
  samples them. This is the standard example output (open in ParaView).
- **HDF5.h** — direct mesh/grid-function persistence.

Examples write outputs into the CWD — run from scratch directories.

## Serialization (src/Rodin/Serialization/)

boost::serialization adapters for Rodin/Eigen/boost containers
(Connectivity, Array, FlatMap/FlatSet, EigenMatrix, Optional, ...) —
binary checkpointing of meshes and fields.

## Python bindings (py/)

pybind11 module (`py/rodin.cpp`, `type_cast.h`), built via scikit-build
(`pyproject.toml`, `setup.py`, `RODIN_WITH_PY`). Currently a thin surface —
check `rodin.cpp` for what is actually exposed before promising Python
functionality.

## Blender plugin (plugins/rodin_blender/)

Blender integration for visualization/mesh interchange. Independent of the
C++ build.
