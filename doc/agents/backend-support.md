# Backend and support matrix

Use this page to decide which mirror implementations and tests a change needs.
Rodin keeps one mathematical surface while varying storage, assembly, and
parallel execution through template specializations.

## Backend vocabulary

| Dimension | Local Eigen path | OpenMP path | MPI path | PETSc path |
| --- | --- | --- | --- | --- |
| Context | `Context::Local` | `Context::Local` plus OpenMP assembly | `Context::MPI` | Local or MPI PETSc objects |
| Mesh | `Geometry::Mesh<Context::Local>` | same mesh | `Geometry::Mesh<Context::MPI>` | same mesh contexts |
| Assembly | `Assembly::Sequential` | `Assembly::OpenMP` | `Assembly::MPI` | `PETSc::Assembly::*` mirrors |
| Linear algebra | `Math::Matrix`, `Math::SparseMatrix`, `Math::Vector` | same objects | usually PETSc for distributed solves | `PETSc::Math::Matrix`, `Vector`, `LinearSystem` |
| Solvers | Eigen/SuiteSparse wrappers | same solvers after assembly | usually PETSc KSP/SNES | `Solver::KSP`, PETSc `CG`/`GMRES`, `SNES` |
| IO | MFEM, MEDIT, HDF5, XDMF | same | HDF5/XDMF shard-aware paths | PETSc vector/grid-function paths |

## What "complete support" means

| Feature family | Local | OpenMP | MPI | PETSc |
| --- | --- | --- | --- | --- |
| Mesh topology/geometry | required | same as Local | mirror specialization and reconciliation | consumed by PETSc forms |
| Finite element space | required | same as Local | mirror specialization when distributed use is claimed | must work with PETSc trial/test if claimed |
| Form-language node | required | same expression node | works if the FES/backend traits support it | works if PETSc form assembly supports it |
| Named bilinear forms | Eigen Sequential | Eigen OpenMP | PETSc MPI | PETSc Sequential/OpenMP/MPI for Mass, Diffusion, Helmholtz, and LinearElasticity |
| Linear/bilinear assembly | `Sequential` | `OpenMP` if enabled | `MPI` if distributed | PETSc Sequential/OpenMP/MPI if PETSc support is claimed |
| Solver | sparse/dense specialization as appropriate | same solver | usually not standalone | PETSc `LinearSystem` specialization if claimed |
| IO format | local loader/printer | same | shard-aware loader/printer if distributed | PETSc data specialization if PETSc fields are claimed |

A feature can deliberately support only one column. Say that explicitly in the
class docs, guides, tests, and PR notes. Do not imply backend parity merely
because the primary template compiles.

## Kokkos local mesh support

| Backend | Context | Current support |
| --- | --- | --- |
| Kokkos | `Context::Kokkos<Context::Local, Backend>` | MeshBase and ConnectivityBase specializations with resident arrays, coordinate/attribute/map updates, device traversal and geometry, quadrature measure reductions, indexed assembly region traversal, host compatibility adapters |

The optional `Rodin::Kokkos` target requires `RODIN_USE_KOKKOS=ON` and an
external Kokkos installation (4.2 or newer). `Context::CUDA` requires Kokkos built with CUDA;
there is no implicit host fallback. Native execution spaces can exercise the
same specialization on the host without establishing CUDA correctness.

The mesh derives directly from MeshBase and retains no LocalMesh. Temporary
host builders/loaders create independently owned Kokkos arrays. Reference-
returning host queries materialize compatibility entries; those are never
authoritative. Copies deep-copy resident storage. Updates execute through
Kokkos before returning; geometry edits invalidate caches with LocalMesh
flush semantics. Default real P1 maps use resident vertices. Independent P1
and identity maps are supported; custom and higher-order maps are rejected.

Connectivity compute/discover/restrict operate directly on resident arrays,
using device sorting and scans. Existing IDs, reference orientation and CPU
traversal order are preserved. Missing edges/faces are discovered from cell
keys; relation arrays are independently owned CSR tables. Mesh attribute/map
slots grow lazily in resident memory while retaining existing entity data.
Only small sizes and reference-cell constants cross the host boundary during
computation. Device views
borrow resident arrays and require completed kernels before mutation, owner
replacement or destruction. Existing host iterators and caches remain host
interfaces, with transfers as necessary. KokkosIteration follows the existing
indexed assembly region contract and keeps device handles private. Boundary
and Interface require explicit face-to-cell incidence; their counts include
all face candidates, and filtering uses exactly one or more than one incident
cell respectively. Launching kernels remains the assembly backend's concern.
Device spaces, fields, operator assembly,
solvers, SubMesh and MPI distribution are not provided.
See [the Kokkos guide](../Guides/Kokkos.dox) for lifecycle and build details.

In CUDA builds, Kokkos's launcher compiles targets linked to `Kokkos::kokkos`
with `nvcc_wrapper`, together with every Rodin header they include. NVCC
imposes constraints that host Kokkos backends do not detect:

- Functors launched by class methods must be accessible at namespace scope:
  NVCC's generated host stubs name them outside the class, so nested kernel
  functors and the value types they carry cannot be private or protected.
- Extended lambdas (`KOKKOS_LAMBDA`) cannot appear in private or protected
  member functions, including GoogleTest's `TestBody`; tests launch device
  work from free functions.
- A Kokkos change is CUDA-verified only after building and running the Kokkos
  tests with a CUDA-enabled Kokkos installation.

## Common extension decisions

- If the math is backend-independent, implement the core expression or
  integrator once and let assembly backends consume it.
- If the storage changes, specialize at the linear algebra or grid-function
  layer, not in the mathematical API.
- If the mesh is distributed, reason about ownership, shared entities, ghosts,
  and reconciliation before writing assembly code.
- If PETSc is involved, reason against PETSc 3.19 behavior even if the local
  PETSc is newer.
- If OpenMP is involved, test determinism and accumulation semantics. Most
  speedups should come from per-cell caching, not a duplicate fast path.

