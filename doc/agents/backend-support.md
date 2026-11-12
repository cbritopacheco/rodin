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
Kokkos built with CUDA, `Context::HIP` Kokkos built with HIP and
`Context::SYCL` Kokkos built with SYCL; there is no
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
interfaces, with transfers as necessary. Each per-entity host query costs a
copy and a device synchronization, so host code visiting a whole mesh goes
through `toHost()` or `Connectivity::materialize(d)`, which copy each resident
table once. KokkosIteration follows the existing
indexed assembly region contract and keeps device handles private. Boundary
and Interface require explicit face-to-cell incidence; their counts include
all face candidates, and filtering uses exactly one or more than one incident
cell respectively. Launching kernels remains the assembly backend's concern.
Device kernels assemble variational problems on this mesh (P1, H1<K> and P0g
spaces, integrals, constraints, linear elasticity and the Neo-Hookean, Saint
Venant-Kirchhoff, Mooney-Rivlin and Yeoh hyperelastic laws) into Eigen or PETSc systems solved by Rodin solvers. A
Problem named with a `Math::LinearSystem<Kokkos::Math::SparseMatrix<Real, Backend>,
Math::Vector<Real>>` keeps its matrix in device memory as a CSR, equal
coefficient by coefficient to the host matrix of the same assembly, and copies
only the right-hand side to the host; `CuDSS`, `RocSolverRF` and the Kokkos
`CG`, `BiCGSTAB` and `GMRES` (`Rodin/Kokkos/Solver.h`) read it in place, and
`BoomerAMG`
(`Rodin/Kokkos/Solver/BoomerAMG.h`) hands it to hypre on the device. Such a
problem takes no preassembled bilinear form. Ordinary local Kokkos assembly
policies retain ordered structural groups and validate their coordinates,
source positions and ranks exactly on the backend before reuse. Numerical
zero filtering, constraint values and lifting are recomputed; copies start
with empty plans. PETSc retains its separate COO cache. SubMesh
is not provided. `Mesh<Context::Kokkos<Context::MPI, Backend>>` makes the
shard of a `Mesh<Context::MPI>` resident and marks owned entities in its
device view; `P1` and `P0g` PETSc problems, single or multi-field, are
assembled on it into distributed Kokkos systems (`Rodin::Kokkos::MPI`). Each
field keeps its global numbering shifted by its offset, as on the host, for
the constraints, the preassembled forms and the problem, but the PETSc rows of
a process hold its owned degrees of freedom of each field in turn, so that the
halo of a product stays that of its part of the mesh. With the offsets alone
the owned rows of a process are not its degrees of freedom, and the halo is
about as large as the local part. The host MPI assembly of several fields
(`PETSc::Assembly::MPI`) lays out its rows the same way through
`PETSc::Assembly::FieldLayout`, which maps an index of the system to its row
only where it enters PETSc; the constraints stay in the numbering of the
system. `LinearSystem::getFieldStarts()` gives the
first row of each field on the process, and `Problem` reads and writes the
fields and builds its field splits from it; code that slices the system by
the offsets must use it too. Each test
function must share the space of the trial function of its block. PETSc grid functions enter
device expressions, local or distributed: their coefficients, ghosts included,
are read from the local form of the vector without communication. `H1<K>`
spaces, scalar or vector, are distributed alike: the H1 space of the
resident shard numbers the local degrees of freedom as the host H1 space of
the shard does. Standalone device `BilinearForm` and `LinearForm` on
distributed spaces assemble distributed Kokkos `Mat` and `Vec`, rectangular
included, and enter distributed problems as preassembled forms: their owned
rows are complete and are renumbered locally, so each of their entries must
couple degrees of freedom the process holds. Periodic and identification
Dirichlet conditions are distributed: their pairs relate global indices of
their fields, each process may give any of them (an identification Dirichlet
condition gives those of the boundary it holds), and they are resolved in a
global constraint map after the fixed values are synchronized; masters the
process does not hold become virtual local degrees of freedom. Eigen systems
are not distributed.
See [the Kokkos guide](../Guides/Kokkos.dox) for lifecycle and build details.

### Transfer audit

The following is a source audit, not a measured ranking. Keep the public
assembly, mesh and solver contracts when reducing transfers:

- `Kokkos/Solver/CuDSS.h` delegates exact graph validation to
  `CuDSSOperator::Pattern`. Kokkos snapshots and comparisons stay in the
  matrix memory space and return only a scalar mismatch result. Compare the
  final CSR graph: reusing assembly groups does not certify it, because
  numerical zero filtering can change the stored graph.
  cuDSS requires CUDA-enabled Kokkos and the `Rodin::Kokkos` target. Include
  `Rodin/Kokkos/Solver.h` or `Rodin/Kokkos/Solver/CuDSS.h`, and keep Kokkos
  initialized during use. The former explicit `Rodin/Solver/CuDSS.h` include
  forwards to it; the core `Rodin/Solver.h` umbrella no longer includes the
  vendor integration or carries its CUDA dependencies.
- `Kokkos/Solver/RocSolverRF.h` shares the exact backend graph snapshot with
  cuDSS. Its gathered pivots are checked by a Kokkos reduction on the
  rocSPARSE HIP stream: only a scalar invalid count returns, independently of
  the solution vector storage type. The same reduction runs on host execution
  spaces in tests. Initial host UMFPACK/AMD analysis still needs its input
  arrays. rocSOLVER requires a HIP-enabled Kokkos installation, the
  `Rodin::Kokkos` target and an initialized Kokkos runtime. Its implementation
  is included by `Rodin/Kokkos/Solver.h`; the former explicit
  `Rodin/Solver/RocSolverRF.h` include forwards to it. The core
  `Rodin/Solver.h` umbrella no longer includes this backend integration.
  No standalone HIP source or CMake HIP language setup is needed.
- Kokkos boundary assembly compacts prescribed values on the backend, then
  returns them through `IndexMap` on the host. `SystemAssembly::build` uploads
  fixed masks and values again. Keep a backend constraint representation
  internally to avoid this round trip; retain host maps where the public
  API requires them. Distributed fixed-DOF reconciliation in `Assembly/PETSc.h`
  also uses host masks and a host-accessed ghost exchange. Any replacement
  must preserve shared-owner averaging, identification and ghost semantics.
- PETSc coefficient lowering reads backend arrays through
  `VecGetArrayReadAndMemType` and caches copied/gathered coefficients by
  PETSc state. Its distributed fill rebuilds and uploads the local-position
  map whenever coefficient state changes. Cache the position map by the
  finite element layout independently of numeric coefficient state.
- PETSc matrix COO preallocation downloads coordinates only when its cached
  graph/preallocation is rebuilt. Numeric `MatSetValuesCOO` and
  `VecSetValuesCOO` already receive device pointers for supported GPU types;
  host objects deliberately use a host compatibility path.
- cuDSS RHS/solution staging and ordinary `Problem::solve` field updates
  include device-to-device copies with backend vectors. Descriptor rebinding
  or shared storage would require lifetime, reallocation and failure-contract
  changes; these are distinct from host-transfer removal. A failed cuDSS
  solve currently leaves the public solution untouched.
- Explicit `toHost`, file output, host coefficient containers and legacy
  mesh/connectivity queries intentionally cross to host memory. GPU loops
  should use views instead of per-entity host queries. Prefer bulk
  materialization for host algorithms. Geometry counts, allocation sizes,
  convergence norms and comparison flags are small control results; reducing
  their synchronization frequency is a separate optimization.

In CUDA builds, Kokkos's launcher compiles targets linked to `Kokkos::kokkos`
with `nvcc_wrapper`, together with every Rodin header they include. NVCC
imposes constraints that host Kokkos backends do not detect:

- Functors launched by class methods must be accessible at namespace scope:
  NVCC's generated host stubs name them outside the class, so nested kernel
  functors and the value types they carry cannot be private or protected.
- Extended lambdas (`KOKKOS_LAMBDA`) cannot appear in private or protected
  member functions, including GoogleTest's `TestBody`; tests launch device
  work from free functions.
- Library headers launch no extended lambda: every kernel is a named functor,
  a public nested struct template or a namespace scope struct, never a class
  local to the launching function. NVCC wraps an extended lambda in
  `__nv_hdl_wrapper_t`, whose copy calls a function pointer stored once per
  translation unit and set by the wrapper constructor of that unit. When two
  sources of an executable instantiate the same header kernel, the linker
  keeps one copy of the out of line launch, which can read the still null
  pointer of the other source and jump to address zero (Kokkos 5.2, CUDA 12.5
  and 13.1). A lambda of a test or example source is unique to it and is
  safe; a user lambda carried as a functor member is too.
  A nested functor whose static members use private constants of its class
  is declared after them, in a public section at the end of the class.
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
  specialization (such as `ContextProjection<Context::Local>`) fails on
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
- PETSc builds STRUMPACK and SLATE with the ROCm `clang++`, both of them:
  SLATE headers embed `omp_nest_lock_t`, whose size differs between libgomp
  and libomp, and a mixed pair faults in every multi-rank factorization.
  SLATE requires `MPI_THREAD_MULTIPLE` from the application.
- A Kokkos change is HIP-verified only after building and running the Kokkos
  tests with a HIP-enabled Kokkos installation.

In SYCL builds, `icpx` of oneAPI 2026.2 compiles Rodin and the device code
of Intel GPUs ahead of time. Its toolchain and the GPU runtime carry defects
that CUDA and HIP do not:

- The OpenMP lowering of `-fiopenmp` miscompiles the host assembly from -O1
  on: in an `omp for`, a test on an `Optional` under a loop-invariant branch
  skips every iteration, so `.over(attribute)` assembles to zero. The top
  CMakeLists.txt and RodinConfig.cmake compile OpenMP with `-fopenmp=libomp`
  (the front end of Clang) and keep libiomp5 at link; the default
  `-fp-model=fast` of icpx is overridden by `-fp-model=precise`.
- IGC 2.40 and 2.41 drop the stores of some lanes in irreducible loops, which
  LLVM builds when it threads a branch on a value of the thread into a loop.
  A test on the geometry inside a loop over the dimension did it in the
  device projections: keep such tests out of device loops, as
  `Geometry::getReferenceVertexCoordinate` does with one mask per geometry,
  and compare a new device kernel with the host on a mixed mesh.
- Double atomic additions in team scratch memory hang about one launch in two
  on the Arc Pro B70: sum each entry in the registers of one thread instead,
  as the team kernels of `InternalVirtualWork.h` do for the displacement
  gradients, and run a new team kernel several times.
- IGC gives a thread 128 registers by default, and the bilinear integration
  kernels spill hundreds of them at SIMD16 (the warnings of ocloc in the
  build log say how many). `src/Rodin/Kokkos/CMakeLists.txt` passes
  `-ze-intel-enable-auto-large-GRF-mode` to ocloc for INTEL_BMG, so the
  kernels that would spill get 256 registers; this took the Oseen atrium
  assembly from 36.0 to 16.5 ms on the renumbered mesh (104 to 72 ms raw).
  Device code is compiled ahead of time by each compile command, not at
  link, so a backend option must reach the compile line. On the raw mesh the
  COO merge stays bound by the locality of its gather (53 ms against 3.3 ms
  renumbered): renumber the mesh.
- IGC 2.41 and 2.42 crash on a kernel of oneDPL that hypre compiles just in
  time, so the image stays on IGC 2.40.
- The `nrm2` of oneMKL 2026.0 to 2026.2 on the GPU reads past the end of a
  vector longer than about 65536 entries, up to the next multiple of 2048.
  Kokkos-Kernels routes `KokkosBlas::nrm2` there and PETSc `VecNorm` of its
  Kokkos vectors; since `VecDuplicateVecs` puts the Krylov work vectors in
  one allocation, CG and BiCGStab norms take in the next vector and a field
  split reads garbage. The SYCL image patches Kokkos-Kernels
  (`kokkos-kernels-sycl-nrm2.patch`) to keep its native kernel on SYCL; a
  PETSc built against an unpatched Kokkos-Kernels with oneMKL fails the
  same way.
- Open MPI does not know SYCL memory, so PETSc runs with
  `-use_gpu_aware_mpi 0`; `HWLOC_COMPONENTS=-levelzero` keeps mpirun from
  crashing and `ONEAPI_DEVICE_SELECTOR=level_zero:gpu` keeps PETSc off the
  OpenCL devices, on which Kokkos refuses to run.
- No direct solver runs on an Intel GPU: CuDSS needs NVIDIA and RocSolverRF
  AMD. Saddle point systems go to PETSc with a Schur field split or to a
  host factorization.
- A Kokkos change is SYCL-verified only after building and running the
  Kokkos tests with a SYCL-enabled Kokkos installation.

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

