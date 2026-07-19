# Conventions and hard rules

These rules exist because violating them has already cost real time. Each one
states the rule, then why.

## Scope discipline

**Keep changes minimal and behavior-preserving. Never bundle a speculative
optimization or "fast path" into a requested change.** If a performance idea
looks worthwhile, propose it separately with a measurement plan. A past
"make the assembly backends uniform" task grew a MatAXPY preassembly fast
path that regressed Dirichlet elimination on CI's PETSc version; the whole PR
was abandoned. When you do performance work, the acceptance bar is *identical
numerics* (iteration counts, final energies/fits) against a baseline run —
not just "tests pass".

## Error handling in PETSc-facing code

**Use `assert(ierr == PETSC_SUCCESS)` after PETSc calls.** This is the
uniform house idiom. Do not introduce a checking macro, and do not convert
existing asserts to one. Better still: backend-independent code (form
language, `Problem`, `LinearSolverBase`) should not touch PETSc at all —
keep PETSc calls inside `src/Rodin/PETSc/`.

## Design preferences

- **Internal variables are first-class DOFs.** For constitutive models with
  internal state (active extension, Hill–Maxwell γ/β, etc.), expose the state
  as `GridFunction`s on a discontinuous space and solve the coupled residual
  with the global Newton. Do not do per-quadrature-point Schur condensation /
  hidden local Newton solves — the mixed approach produced confirmed-wrong
  tangents. Reserve condensation for a measured performance need.
- **Constraints inside optimization solves are smooth penalties, never hard
  projections** (after Knupp–Kolev–Mittal–Tomov 2021, arXiv:2105.12165).
  Snapping nodes to a manifold inside a line search creates α-invariant
  discontinuities that make Armijo reject every step. Projection is fine
  *outside* the solve (initialization, coarsener feature projection).
- **Dimension-generic code** uses `Math::SpatialVector`/`SpatialMatrix`
  (runtime dims, stack storage) rather than fixed 2×2 types, so 2D/3D share
  one path. Intrinsically-2D constructs (e.g. triangle-specific targets) are
  fine as such — don't fake genericity.
- **FES-independence.** Code that maps between meshes, displacement fields,
  and local DOFs must not assume a DOF layout (`node*vdim+component`), a
  nodal-interpolatory basis, or P1. Rodin's high-order H1 uses Fekete/Dubiner
  bases, so "coefficient = nodal coordinate" is false at P2. Get geometry
  from `Geometry::Point(cell, rc)` (physical coords / Jacobian) and field
  values/Jacobians from the `GridFunction` itself.

## Testing conventions

- Derivative code (residual/tangent pairs) is gated by finite-difference
  consistency tests; keep them green and add one for any new term.
- Mesh rebuilds/compaction (remeshing, SubMesh extraction) re-index
  polytopes; tests must compare coordinate/key *sets*, not per-index values.
- Numerical acceptance in the adaptation pipelines is measured, not eyeballed:
  report fit rms/max, qmin, inverted count, coverage per stage.

## Output artifacts

Examples write `*.h5`/`*.xdmf`/`*.log` into the current directory. Run them
from a scratch directory. Never commit these files.
