# Testing and CI

## Suite layout

```
tests/unit/Rodin/<Module>/...      GoogleTest, mirrors src/Rodin exactly
tests/manufactured/Rodin/...       Convergence + targeted regressions
tests/benchmarks/                  Google Benchmark (target RodinBenchmarks)
tests/installation/                Install-tree consumption checks
src/Rodin/Test/                    Library-side helpers (Random functions, Utility)
```

- ctest is **not** wired up (`ctest -N` → 0 tests). Run the gtest
  executables directly: `build/tests/unit/Rodin/<Module>/
  Rodin<Module><Component>Test` (`--gtest_filter`, `--gtest_list_tests`).
- Manufactured tests are organized by space and cross-space combination
  (`H1/`, `P0/`, `P1/`, `P0_P1/`, `P1_H1/`, `PreassembledMixed/`,
  `Assembly/`, `Solver/`, `Models/`, `MPI/`, `PETSc/`) — the place for
  "does this converge at the right rate" and backend-contract regressions
  (e.g. `PETSc/TargetedAssemblyTest.cpp` pins sparsity-pattern reuse:
  0 mallocs on reuse, structural zeros preserved).

## House test patterns

- **FD-consistency gates** for every residual/tangent pair: assemble R and
  J, compare J·du against finite differences of R — required < 1e-6, at
  P1 AND at higher order if the term claims order-genericity (a term can
  be FD-self-consistent yet wrong, and P1-consistent yet P2-broken; both
  have happened).
- **Manufactured solutions** for operators/solvers: pick exact u, derive
  f, assert convergence order.
- **Identity checks** for algebraic structure: R = M·u, symmetry, energy
  identities (R·u = 2E) — these hold exactly for consistent quadrature,
  independent of resolution.
- **Rebuild-safe comparisons**: mesh optimizers re-index on rebuild;
  compare coordinate/key *sets*, never per-index values.
- **Randomized inputs** via `Rodin::Test::Random` where a property should
  hold for arbitrary data.
- **Benchmarks as oracles**: performance claims cite
  `RodinBenchmarks --benchmark_filter=...` numbers plus identical-numerics
  evidence (iteration counts, final energies) — never "it feels faster".
- New capability defaults off; its tests prove both the feature and the
  default-off no-regression.

## CI (.github/workflows/)

`Build.yml`, `Tests.yml` (Unit + Manufactured jobs), `Benchmarks.yml`
(published to the benchmarks site), `Documentation.yml` (Doxygen site),
`Installation.yml`, `StaticAnalysis.yml`, `Labeler.yml`. Branches
`master` and `develop` both gate.

CI facts that bite:

- Ubuntu gcc-14 builds against **PETSc 3.19** (local machines are usually
  newer) — assembly changes can pass locally and fail CI; reason against
  3.19 semantics (petsc.md).
- Asserts vanish in optimized CI configs — a silently swallowed error can
  look like an unrelated downstream failure.
- Doxygen documentation is published per-branch; malformed doc comments
  can fail Documentation.yml.
