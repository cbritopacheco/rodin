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
| Linear/bilinear assembly | `Sequential` | `OpenMP` if enabled | `MPI` if distributed | PETSc Sequential/OpenMP/MPI if PETSc support is claimed |
| Solver | sparse/dense specialization as appropriate | same solver | usually not standalone | PETSc `LinearSystem` specialization if claimed |
| IO format | local loader/printer | same | shard-aware loader/printer if distributed | PETSc data specialization if PETSc fields are claimed |

A feature can deliberately support only one column. Say that explicitly in the
class docs, guides, tests, and PR notes. Do not imply backend parity merely
because the primary template compiles.

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

