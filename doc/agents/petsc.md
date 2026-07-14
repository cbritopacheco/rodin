# PETSc backend notes

## Version skew: CI is the real gate

CI (Ubuntu gcc-14 "Unit"/"Manufactured" jobs) builds against **PETSc 3.19.6**;
local dev machines are typically newer (e.g. MacPorts 3.22). Semantics differ.
Concrete example: `MatZeroRowsColumns(..., diag=1, ...)` on 3.19 **aborts** if
a zeroed row has no preallocated diagonal entry ("Matrix is missing diagonal
entry"), while 3.22 tolerates it. Any change under `src/Rodin/PETSc/Assembly/`
must be reasoned against 3.19 behavior, not just pass locally.

## Error-handling idiom

`assert(ierr == PETSC_SUCCESS)` after each call — see
doc/agents/conventions.md. Remember asserts compile out under `-DNDEBUG`
(Release/RelWithDebInfo), so a swallowed failure surfaces later as a
confusing downstream error (historically: "Row too large" long after a failed
resize). When debugging such symptoms, suspect an earlier silently-failed
PETSc call.

## No in-place resize — LinearSystem lifetime contract

A PETSc `LinearSystem` is bound to fixed finite element spaces; its global
sizes never change over its lifetime. PETSc has no in-place resize for an
assembled `Mat`/`Vec` (`MatSetSizes`/`VecSetSizes` are pre-layout only).
`MatrixSetup` therefore: reuses on size match, sets up virgin objects, and
**raises** on an assembled-but-different-size request. `VectorSetup.h`
mirrors this for `b`/`x`/`res`, with `zeroOnReuse=true` for `b`/`res` and
`false` for `x` (preserves the solver's warm start). A different mesh/space
means a new `LinearSystem` — do not "fix" size mismatches by resizing.

Pattern-stability regression tests:
`tests/manufactured/Rodin/PETSc/TargetedAssemblyTest.cpp` (reuse keeps the
sparsity pattern with 0 mallocs; structural zeros stay allocated).

## Solver factor reuse depends on the sparsity pattern

PETSc-backed direct solvers reuse factorizations only while the matrix
sparsity pattern is unchanged (`MatGetNonzeroState` /
`SAME_NONZERO_PATTERN`). Assembly changes that perturb the nonzero
pattern silently invalidate factor reuse — treat pattern stability as
part of the assembly contract (the regression tests above pin it).

Note: further PETSc-backed adaptation solvers (e.g. for WNGIR in 3D)
live on the `module/Adaptation` lineage, not on `develop` — verify
presence before citing them.
