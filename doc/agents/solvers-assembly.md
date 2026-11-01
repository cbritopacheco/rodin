# Assembly and solvers

For new assembly behavior, solver support, or backend claims, read
[workflows.md](workflows.md), [backend-support.md](backend-support.md), and
[numerical-contracts.md](numerical-contracts.md) with this file. They define
the extension checklist, supported-configuration matrix, and numerical
contracts that assembly/solver changes are expected to preserve.

## Assembly (src/Rodin/Assembly/)

Named kernels return entries with test rows and trial columns, already using
the trial-linear/test-conjugate-linear convention. `ScatterMap` and the PETSc
named insertion paths insert these values unchanged on both first assembly
and reassembly. Conjugating at insertion changes complex coefficients and is
not a storage conversion. See the complex named-form regressions for matched
Integral operators and the PETSc local/distributed variants.

- Backends are mesh-iteration strategies: `Sequential` and `OpenMP`
  (per-thread local evaluation, then merge). `Default.h` selects by build
  configuration. The MPI and PETSc trees provide their own mirrored
  backends.
- `AssemblyBase::execute` is non-const across all backends. Assemblers may
  retain sparsity maps and working storage between calls. Input wrappers
  remain const; identification inputs explicitly borrow a non-const
  shape-function expression for updating its integration-point cursor.
  Local named-form kernels also use non-const `compute`, with
  one copy per assembly or OpenMP worker. Concurrent executions must use
  separate assembler instances. Custom backends must override non-const
  `execute`; const assembler callers must hold a non-const instance instead.
- Inputs are typed bundles: `BilinearFormAssemblyInput<TrialFES, TestFES>`
  (and the linear counterpart) carry the spaces plus **two integrator
  lists: local (per-polytope quadrature) and global** (whole-mesh
  couplings, e.g. P0g global-constant terms). New assembly features must
  handle both lists.
- `ConstraintMap` — read-only expansion map for essential constraints,
  covering both DirichletBC flavours: *value* constraints (DOF = number)
  and *identification* constraints (DOF = linear combination of other
  DOFs). Elimination happens at assembly via this map and
  `Math::LinearSystem`'s eliminate support — constraints are structural,
  not penalties.
- Performance model (measured, repeatedly): assembly cost is dominated by
  **integrand evaluation**, which is memory-bandwidth-bound — parallel
  speedup plateaus around 4–5 threads; the serial merge/setFromTriplets
  tail is ~1%. The effective lever is per-cell caching of reference-basis
  values/Jacobians in the integrator's bind step, never more threads and
  never a second "fast path" (conventions.md).

## Linear solvers (src/Rodin/Solver/)

All are thin `LinearSolver` wrappers presenting
`Solver(problem-or-system).solve()` with chainable setters
(`setTolerance(...).setMaxIterations(...)`); CTAD picks the specialization
from the system type.

- Eigen iterative: `CG`, `BiCGSTAB`, `GMRES`, `DGMRES`, `MINRES`, `IDRS`,
  `IDRSTABL`, `LeastSquaresCG`.
- Eigen direct: `SparseLU`, `SparseQR`, `SimplicialLLT`, `SimplicialLDLT`,
  `LDLT`, `HouseholderQR`, `PartialPivLU`.
- SuiteSparse (configure-gated): `UMFPack`, `ParU`, `CHOLMOD`, `SPQR`.
- MUMPS (configure-gated): `MUMPS`, the only direct solver that exploits a
  symmetric matrix (`setSymmetric`), halving the factor it stores.
- cuDSS (`RODIN_USE_CUDSS`): `CuDSS`, a direct solver on an NVIDIA GPU for
  real systems. On a `Math::SparseMatrix` it copies the CSC arrays to the
  device as the CSR arrays of the transpose and solves with the transpose,
  so a general matrix needs no conversion. On a `Kokkos::SparseMatrix` in
  device memory it reads the CSR arrays in place and only copies the
  right-hand side and the solution; `Internal::CuDSSOperator` describes
  both layouts. Like `MUMPS`, it retains the analysis
  and the factors: `solve` reuses them, `factorize(axb)` refactorizes after
  the values change and keeps the analysis, and a changed pattern requires
  `clear(Factorization::Symbolic)`. Only the row and nonzero counts are
  checked against the analysis. `setSymmetric` selects `LDL^T` or `LL^T`,
  `setRefinementSteps` adds iterative refinement. Configure with
  `-Dcudss_DIR=<cudss>/lib/cmake/cudss`.
- Kokkos-Kernels (`RODIN_USE_KOKKOSKERNELS`, target `Rodin::KokkosKernels`,
  header `Rodin/Kokkos/Solver.h`): `CG` and `GMRES` specialized for a system
  whose matrix is a `Kokkos::SparseMatrix`, so a device assembly is solved
  without PETSc. Both are Jacobi preconditioned and read the CSR arrays in
  place; only the right-hand side and the solution cross to the host. `CG`
  is written on Kokkos with one SpMV handle per solve and keeps the
  iteration, defaults and stopping test of the sparse Eigen CG. `GMRES` is
  the restarted GMRES of Kokkos-Kernels (`setRestart`, 30 by default), with
  the iteration limit rounded up to whole cycles. Configure with
  `-DKokkosKernels_DIR=<prefix>/lib/cmake/KokkosKernels`; Kokkos-Kernels must
  be built on the Kokkos Rodin uses, with `int` ordinals and offsets.
- hypre (`RODIN_USE_HYPRE`, target `Rodin::KokkosHypre`, header
  `Rodin/Kokkos/Solver/BoomerAMG.h`): `BoomerAMG`, the conjugate gradient of
  hypre preconditioned by one BoomerAMG V-cycle, for an SPD system whose
  matrix is a `Kokkos::SparseMatrix`. hypre copies the CSR arrays in the
  memory of the backend (`HYPRE_IJMatrixSetValues2`), so nothing but the
  right-hand side and the solution crosses to the host. The hierarchy uses
  GPU settings (PMIS, ext+i, l1 Jacobi, Gaussian elimination on the coarsest
  grid); `setStrongThreshold`, `setAggressiveLevels` and
  `setMaxInterpolationElements` tune it and release the retained hierarchy.
  The hierarchy and the hypre copy of the matrix it comes from are retained
  across solves (`setup` rebuilds them from the current matrix, `clear`
  releases them, a copy of the solver does not share them), like
  `KSPSetReusePreconditioner`: the CG runs on a fresh copy of the current
  matrix and the V-cycle on the retained one (`HYPRE_PCGSetPrecondMatrix`,
  with a no-op preconditioner setup), so a stale hierarchy only costs
  iterations. A size change raises. The stop test is the unpreconditioned
  relative residual, 1e-10 by default. MPI and hypre are initialized on the
  first setup if needed, their finalization is left to the application,
  after the last solver is destroyed. hypre must use 32-bit indices and doubles.
  Configure with `-DHYPRE_DIR=<prefix>`.

Factorization-based solvers share no base class, but report uniformly through
`getInfo()` and `success()` (`Solver/Info.h`): the success of the most recent
operation, the retained `Factorization` stage, and the backend's own `status`.
A factorization that fails is an outcome, not a defect: the solve is skipped,
the solution vector is left untouched, and nothing is raised. A violated
contract, such as a matrix that is not square, still raises. `HouseholderQR`
and `PartialPivLU` report nothing, because Eigen exposes no status for them:
a caller that needs to know checks the residual, or picks another solver.
The dense solvers act on `Math::Matrix` systems. A `Problem` assembles into
whichever `LinearSystem` it is given, and only its deduction guide assumes a
sparse one, so a dense problem is written by naming the template arguments:
`Problem<Math::LinearSystem<Math::Matrix<Real>, Math::Vector<Real>>, U, V>`.
- Platform: `AppleAccelerate`.
- PETSc KSP wrappers live under `PETSc/Solver` (petsc.md).

Choosing: SPD → CG (iterative) or SimplicialLDLT/CHOLMOD (direct);
nonsymmetric → BiCGSTAB/GMRES or SparseLU/UMFPack/ParU; symmetric indefinite
(saddle point) → MUMPS with `setSymmetric`; small dense → LDLT/
PartialPivLU. Examples default to CG for Poisson-like and SparseLU/LDLT
for Newton tangents. A sequence of systems with one sparsity pattern, such
as the time steps of a flow, keeps one `CuDSS` or `MUMPS` instance and calls
`factorize` at each step. On a 2D Poisson system with 1M unknowns, `CuDSS`
takes 6.6 s for the first solve against 38 s for `SparseLU`, then 234 ms per
refactorization and 15 ms per solve (RTX 2080 Ti, 32 cores). The Kokkos
`NavierStokes` example does this on a device matrix: 160 steps take 3.3 s
with `CuDSS` (3.6 s through the host matrix) against 60 s with a new
`SparseLU` at each step.

PETSc on the device, measured on the Kokkos examples (RTX 2080 Ti,
32 cores, `examples/Kokkos/PDEs/Seq_*.cpp`):

- Vector elasticity: CG with `-pc_type gamg` and the rigid body modes given
  by `MatSetNearNullSpace` takes 28 iterations and 3.1 s at 822k unknowns,
  where hypre takes 14 s, ILU 27 s and SuperLU_DIST 19 s, and Jacobi does not
  converge. The vector P1 numbering is component-major, so the modes are
  built by hand: a block size cannot describe them.
- Scalar 3D Poisson, Q1: Jacobi CG solves 1M unknowns in 208 ms (159
  iterations at the default `rtol` of the wrapper), GAMG in 6.1 s. The
  Kokkos-Kernels `CG` on the device matrix (`Krylov_Poisson.cpp`) takes
  171 ms for 158 iterations at 1e-10 and 121 ms at 1e-5. The default BoomerAMG settings run out of device memory at 1M
  unknowns (5 GB at 262k); `-pc_hypre_boomeramg_strong_threshold 0.7
  -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_P_max 4` solves it in
  1.9 s and is also faster at 262k. `BoomerAMG` with the same settings on
  the device matrix (`BoomerAMG_Poisson.cpp`) takes 644 ms at 1M (42
  iterations, against 32 for PETSc with the unpreconditioned norm) and
  404 ms at 262k, against 592 ms. Of these, the setup takes 463 ms at 1M
  (343 ms at 262k, MPI and hypre initialization included) and a solve with
  the retained hierarchy 260 to 280 ms (90 to 113 ms).
- Taylor-Hood Oseen steps: FGMRES with a Schur field split (`selfp`, GAMG
  on the velocity block, a few GMRES iterations on the pressure) beats ILU
  and SuperLU_DIST, 6.1 s against 29 s and 25 s for 5 steps at 149k
  unknowns. `Problem::setFieldSplits` names the blocks after the trial
  functions.
- SuperLU_DIST on one process uses the GPU for a handful of large Schur
  updates only (about 12 ms of DGEMM in a 3 s factorization).
  `-mat_superlu_dist_3d -mat_superlu_dist_d 1` divides the factorization
  time by three; `SUPERLU_ACC_SOLVE=1` makes it ten times slower.

`RodinSolver` is a compiled library. It instantiates the sparse `CG` and
`SparseLU` for `Real` and `Complex` (`CG.cpp`, `SparseLU.cpp`), and its
`RODIN_SOLVER_EXTERN_TEMPLATES` definition turns the matching
`extern template` declarations on in every target that links it. Their
`solve` is defined out of class so that it is not inline. This keeps the
Eigen kernels out of NVCC translation units, which run Eigen without
vectorization (backend-support.md). A target that includes these headers
without linking `Rodin::Solver` still instantiates them itself. Another
wrapper used by device code gets the same treatment: an out-of-class `solve`,
the `extern template` lines, and a `.cpp` added to `RodinSolver_SRCS`.

## Nonlinear: NewtonSolver

`Solver::NewtonSolver` drives residual/tangent pairs assembled from the
form language. **It has constant damping only — no line search, no
validity gating.** Consequences, confirmed by measurement:

- Barrier-type energies (log barriers, det>0 metrics) cannot be minimized
  reliably with it; solvers that need Armijo backtracking + validity
  rejection hand-roll the loop today (Adaptation examples/benchmarks).
- A global line-searched Newton is a known, wanted-but-unbuilt feature. If
  you build it, it belongs in NewtonSolver as an opt-in policy, defaulted
  off, with existing behavior bit-preserved.

Related but distinct: `Math::NewtonRaphson` is a small scalar root finder
(local constitutive solves), not the FEM Newton driver.

## The LinearSystem contract

`Math::LinearSystem` (and its PETSc mirror) binds to fixed spaces: sizes
are immutable for its lifetime; DOF elimination implements essential
constraints; a new mesh/space means a new system object (this is a hard
rule on the PETSc side — petsc.md). Warm starts: solution vectors are
reused, not zeroed, across assemblies.

When a feature claims backend support, update the expectation in
[backend-support.md](backend-support.md) and cover every listed assembly path
that can exercise it.
