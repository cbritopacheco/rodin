# Assembly and solvers

## Assembly (src/Rodin/Assembly/)

- Backends are mesh-iteration strategies: `Sequential` and `OpenMP`
  (per-thread local evaluation, then merge). `Default.h` selects by build
  configuration. The MPI and PETSc trees provide their own mirrored
  backends.
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
- SuiteSparse (configure-gated): `UMFPack`, `CHOLMOD`, `SPQR`.
- Platform: `AppleAccelerate`.
- PETSc KSP wrappers live under `PETSc/Solver` (petsc.md).

Choosing: SPD → CG (iterative) or SimplicialLDLT/CHOLMOD (direct);
nonsymmetric → BiCGSTAB/GMRES or SparseLU/UMFPack; small dense → LDLT/
PartialPivLU. Examples default to CG for Poisson-like and SparseLU/LDLT
for Newton tangents.

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
