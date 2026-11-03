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
external Kokkos installation (4.2 or newer; validated with 5.2.2). `Context::CUDA` requires
Kokkos built with CUDA and `Context::HIP` Kokkos built with HIP; there is no
implicit host fallback. Native execution spaces can exercise the same
specialization on the host without establishing device correctness.

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
Device kernels assemble variational problems on this mesh (P1, H1<2> and P0g
spaces, integrals, constraints, linear elasticity and the Neo-Hookean, Saint
Venant-Kirchhoff, Mooney-Rivlin and Yeoh hyperelastic laws) into Eigen or PETSc systems solved by Rodin solvers. A
Problem named with a `Math::LinearSystem<Kokkos::SparseMatrix<Real, Backend>,
Math::Vector<Real>>` keeps its matrix in device memory as a CSR, equal
coefficient by coefficient to the host matrix of the same assembly, and copies
only the right-hand side to the host; `CuDSS`, `RocSolverRF` and the Kokkos-Kernels `CG`
and `GMRES` (`Rodin/Kokkos/Solver.h`) read it in place, and `BoomerAMG`
(`Rodin/Kokkos/Solver/BoomerAMG.h`) hands it to hypre on the device. Such a
problem takes no preassembled bilinear form. SubMesh
is not provided. `Mesh<Context::Kokkos<Context::MPI, Backend>>` makes the
shard of a `Mesh<Context::MPI>` resident and marks owned entities in its
device view; single field `P1` PETSc problems are assembled on it into
distributed Kokkos systems (`Rodin::KokkosMPI`). Periodic conditions,
preassembled forms, mixed problems and Eigen systems are not distributed.
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
- A `using namespace Rodin;` directive at namespace scope makes `Kokkos`
  ambiguous with `Rodin::Kokkos`; NVCC resolves the `using namespace Kokkos;`
  directives of the Kokkos templates at instantiation and fails. Examples and
  tests place the directive in a function.
- NVCC parenthesizes host expressions before handing them to GCC, and GCC
  then rejects `Integral(Grad(u), Grad(v)).over(1) + ...` with "template
  placeholder not permitted in parameter": the deduced `Grad(u)` is parsed as
  a parameter type. Such a restricted term goes after another term of the sum,
  or in a named variable.
- NVCC defines `EIGEN_GPUCC`, which turns off Eigen's vectorization in host
  code too. A sparse factorization or Krylov solve instantiated in an NVCC
  translation unit is about 20% slower than under GCC. The sparse `CG` and
  `SparseLU` are therefore instantiated in `RodinSolver`, which is compiled by
  the host compiler (solvers-assembly.md). Their layout does not depend on
  the compiler: `EIGEN_MAX_ALIGN_BYTES` is 16 under both.
- The NVCC 12 front end checks a requires clause of a member template when
  its class is instantiated, so a clause naming a member of an incomplete
  specialization (such as `DeviceProjection<Context::Local>`) fails on
  host spaces; such a clause goes through a named concept. It also finds a
  braced hinted `insert` of a Boost 1.74 `flat_map` ambiguous;
  `emplace_hint` is used instead.
- Kokkos checks kernel launches only when its bounds checks are enabled. A
  kernel refused at launch, for instance when a nearly full device has no
  memory left for its local arrays, leaves its outputs unwritten without an
  error. The device assemblers call `Kokkos::Internal::checkLaunches()` after
  each integrator fill (`Rodin/Kokkos/Utility.h`); a new path launching
  integrator kernels does the same.
- A Kokkos change is CUDA-verified only after building and running the Kokkos
  tests with a CUDA-enabled Kokkos installation.

In HIP builds, the ROCm `clang++` compiles Rodin and the Kokkos target adds
`-xhip` to the targets linking it, as the Kokkos launcher hands only those to
`nvcc_wrapper` under CUDA. Configuring with `hipcc` compiles the host-only core
as HIP too, where `constexpr` functions calling host code and Eigen
`long double` math fail. Clang enforces rules that GCC and NVCC let pass:

- A partial specialization that refines a core specialization through a
  `requires` clause, such as the device `Integral<Dot<...>>`, declares its
  template parameters in the order of the core one. Clang orders partial
  specializations by their constraints only when their template parameter
  lists are equivalent, and otherwise reports them ambiguous.
- `-Wunused-lambda-capture` is an error. A capture used only in one branch of
  an `if constexpr` is discarded explicitly in the other branch, and a member
  called from a generic lambda capturing `this` is named `this->template f<T>`.
- With the ROCm `clang++`, `FindOpenMP` resolves the ROCm libomp. Under
  `hipcc` it returns `-fopenmp=libomp` without a library and the links miss
  libomp.
- `Kokkos::Internal::checkLaunches()` reads `hipGetLastError()` as it reads
  `cudaGetLastError()` under CUDA.
- hypre needs `--with-gpu-warp-size=32` on RDNA GPUs, since ROCm 10 lacks
  the `rocm-smi` it detects the warp size with, and MPI runs need
  `HWLOC_COMPONENTS=-pci` against the zlib that ROCm loads (see the HIP
  build section of `doc/Guides/Kokkos.dox`).
- A Kokkos change is HIP-verified only after building and running the Kokkos
  tests with a HIP-enabled Kokkos installation.

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

