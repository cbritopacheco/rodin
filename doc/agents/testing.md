# Testing and CI

## Suite layout

```
tests/unit/Rodin/<Module>/...      GoogleTest, mirrors src/Rodin exactly
tests/manufactured/Rodin/...       Manufactured solutions + fixed-mesh regressions
tests/convergence/{h,p,hp,...}/... Multi-resolution convergence-rate validation
tests/benchmarks/Rodin/<Module>/   Module/component Google Benchmark executables
tests/benchmarks/Workloads/        Shared benchmark forms and numerical oracles
tests/installation/                Install-tree consumption checks
src/Rodin/Test/                    Library-side helpers (Random functions, Utility)
```

- ctest **is** wired up: suites register through `gtest_discover_tests` with
  the labels `unit`, `manufactured`, `convergence`, `slow`, `distributed` and
  `petsc`. Scope a run to what changed:
  `ctest --test-dir build/tests -L unit -LE "slow|distributed"
  --output-on-failure`, or `-R <pattern>` for one suite. Running a gtest
  executable directly still works when you want its raw output
  (`build/tests/unit/Rodin/<Module>/Rodin<Module><Component>Test`,
  `--gtest_filter`, `--gtest_list_tests`).
- Build type is part of correctness here: unit tests may run `Debug`
  (assertions, sanitizers), but manufactured and convergence tests solve PDEs,
  so they run `Release`/`RelWithDebInfo` — `Debug` is impractically slow.
- Manufactured tests are organized by space and cross-space combination
  (`H1/`, `P0/`, `P1/`, `P0_P1/`, `P1_H1/`, `PreassembledMixed/`,
  `Assembly/`, `Solver/`, `Models/`, `MPI/`, `PETSc/`) — the place for
  "does this converge at the right rate" and backend-contract regressions
  (e.g. `PETSc/TargetedAssemblyTest.cpp` pins sparsity-pattern reuse:
  0 mallocs on reuse, structural zeros preserved).
- Convergence tests are organized by refinement mechanism (`h/`, `p/`, `hp/`,
  `isoparametric/`). They use enough discretizations to establish the stated
  rate interval, integrate error norms independently of the solved field, and
  assert every observed rate.
- Library benchmarks mirror source modules, including MPI/PETSc mirrors.
  Module and component executables are opt-in and are not CTest timing gates.
  `RodinBenchmarks` retains the historical publication registration set;
  verification-support profiling lives separately in `benchmarks/convergence`.

A change is complete after the affected unit test executable passes. Assembly
and solver changes also require the relevant manufactured tests to pass.

## Required coverage

Whenever applicable, a change ships every one of these that makes sense for
it. "It passes the happy path" is not coverage; a test that cannot fail
proves nothing. If a category does not apply, that is a judgement to state,
not one to skip silently.

| Category | What it checks |
| --- | --- |
| Happy-path / positive (success-path, expected-success) | Valid input gives the expected result |
| Negative / error-path (failure-path, expected-failure) | Invalid input is rejected correctly |
| Boundary / edge-case | Values at or near limits |
| Exceptional behavior | Correct exception, error code, or failure state |
| State-transition | Object moves between valid states correctly |
| Invariant / property | A rule remains true across many inputs |
| Regression | A previously fixed bug stays fixed |
| Interaction / collaboration | Calls to dependencies occur correctly |
| Side-effect | Files, logs, events, database writes, etc. are correct |
| Idempotency | Repeating an operation has the intended effect |
| Determinism | Same input gives the same result when required |
| Numerical | Accuracy, tolerance, stability, convergence |
| Resource / lifecycle | Allocation, cleanup, ownership, leaks |
| Concurrency | Races, ordering, deadlocks, thread safety |
| Timeout / cancellation | Operations stop or cancel correctly |
| Serialization / round-trip | Encode–decode preserves information |
| Compatibility | Behavior remains valid across versions or formats |

Several of these have a house form, given under *House test patterns*
below: FD-consistency gates and manufactured solutions are the numerical
category; identity checks (R = M·u, symmetry, R·u = 2E) are the invariant
category; rebuild-safe set comparisons are how determinism is stated for
re-indexing code. Others map onto standing hazards rather than patterns:
PETSc handle ownership is resource/lifecycle, the OpenMP and MPI assembly
backends are concurrency, and CI building against PETSc 3.19 while local
builds are newer is compatibility.

Use [backend-support.md](backend-support.md) to decide which configurations a
feature must test, and [numerical-contracts.md](numerical-contracts.md) to
choose the invariant or manufactured check that proves the numerical behavior.

Two rules that make the difference between coverage and theatre:

- **A regression test must be shown to catch the bug.** Revert the fix,
  watch the test fail, restore the fix. A regression test never observed
  red is an assumption.
- **A tolerance-based test must fail when the quantity is wrong.**
  Perturb the code (drop a term, change a factor) and confirm the test
  rejects it. Otherwise the tolerance, not the code, is what passed.

## House test patterns

- **FD-consistency gates** for every residual/tangent pair: assemble R and
  J, compare J·du against finite differences of R — required < 1e-6, at
  P1 AND at higher order if the term claims order-genericity (a term can
  be FD-self-consistent yet wrong, and P1-consistent yet P2-broken; both
  have happened).
- **Manufactured solutions** for operators/solvers: pick exact u, derive
  f, assert convergence order.
- **All-geometry coverage** for every new solver backend: exercise an
  assembled problem on Segment, Triangle, Quadrilateral, Tetrahedron,
  Hexahedron, Pyramid, and Wedge meshes.
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

`Build.yml`, `Tests.yml` (unit, manufactured, convergence, installation and
coverage jobs), `TestGroup.yml` (the reusable directory-group runner),
`Benchmarks.yml` (published to the benchmarks site), `Documentation.yml`
(Doxygen site), `StaticAnalysis.yml`, `Labeler.yml`. Branches
`master` and `develop` both gate.

All test checks use the `Tests /` umbrella. Test group paths follow their
CMake registration directories, with named component leaves for oversized
directories. Short sibling suites share their common parent job while retaining
their full directory labels. MPI retains its module position; `slow` remains an independent
execution label. Compiler, optimization, threading, sanitizer, scalar and MPI
rank configurations must survive scheduling changes. Splitting jobs does not
justify weakening a numerical check or its timeout.

`dev/ci_test_plan.py` generates the current directory matrices. The configured
`RodinTestTargets-<configuration>.tsv` manifest, emitted by CMake, is authoritative
for enabled executable targets. `rodin_discover_tests` forwards GoogleTest
registrations unchanged while recording their execution profiles. Explicit
CTest commands supply the MPI target/profile associations. New test targets
in an existing registration directory are included automatically.

Use the scoped runner to avoid discovering unbuilt sibling executables:

```sh
python3 dev/ci_tests.py build --build build --group unit/Rodin/Geometry --profile unit
python3 dev/ci_tests.py run --build build --group unit/Rodin/Geometry --profile unit
```

The build must already be configured for the desired backend and sanitizer
settings. The runner retains generated commands and properties, restores the
original working directories, adds hierarchical CTest labels, verifies its
selection, and runs CTest serially. Each distributed invocation keeps its
registered MPI rank count. Build-time and discovery-only runtime profiles are
redirected outside the measured coverage directory.

Unsanitized Unit, Manufactured and Slow lanes upload production-only LCOV
independently as they finish. Native/MPI/PETSc convergence and restriction lanes
also capture and upload their own reports; sanitizer and installation checks do
not contribute coverage. There is no separate coverage test execution lane.
The former Coverage backend configuration is retained under ordinary Unit and
Manufactured jobs as `Backends / Baseline`; their slow tests run only in Slow.
Native MPI H1 restriction jobs partition the complete 420-case matrix by MPI
rank, retaining all geometries, degrees, mappings and per-test watchdogs. This
limits each instrumented job to one rank under the existing 90-minute limit.
These jobs use passive OpenMP waiting, no Open MPI CPU binding and yielding
inside blocking MPI calls to avoid starving compute ranks on shared runner
CPUs. They record the runner CPU topology and rank affinity before execution;
rank counts, thread counts and per-test watchdogs remain unchanged.
Original optional-backend Unit and default-backend Manufactured configurations
remain intact. Compilation contexts distinguish backend, optimization and
threading; final line coverage is their union, never an average of percentages.

`dev/ci_test_plan.py` reads the explicit workflow matrices with pinned PyYAML
(provisioned by preflight) to inventory every targeted coverage partition.
Long PETSc convergence lanes partition independent PDEs and communicator
sizes rather than sharing one job deadline. Curved elasticity invokes each
GoogleTest study separately; each study retains its complete refinement
history. Pyramid nonlinear rate studies have separate scheduling partitions,
with complementary controls executed once. The CI registration guards verify
the concrete CMake targets, MPI rank registrations, matrix-size limits, and
selection of every curved elasticity and nonlinear study exactly once.
Zero-hit baseline jobs retain previously compiled test/template instantiations
and production libraries even when those targets are absent from a runtime
selection. Tests and benchmarks supply execution hits, never denominator
entries. `Tests / Coverage / Finalize` requires every current-commit producer,
successful uploads, report digests, and the complete configuration-specific
compilation inventory before
triggering Codecov notifications and enforcing the 95% production line floor.
It neither rebuilds nor reruns tests. Missing producers cannot be replaced by
carryforward coverage. The local guard tests run serially in `Tests / CMake / CI`.

CI facts that bite:

- Ubuntu gcc-14 builds against **PETSc 3.19** (local machines are usually
  newer) — assembly changes can pass locally and fail CI; reason against
  3.19 semantics (petsc.md).
- Asserts vanish in optimized CI configs — a silently swallowed error can
  look like an unrelated downstream failure.
- Doxygen documentation is published per-branch; malformed doc comments
  can fail Documentation.yml.
- Documentation.yml is also the review gate for specialization-family class
  pages: if a class has multiple supported specializations, its Doxygen
  description must include a complete `Specialization` / `Description` table
  with each specialization linked.
- Documentation review also treats unlinked public references as failures:
  `@see` blocks and prose lists should use explicit `@ref` references or HTML
  links for public classes, templates, specialization groups, and generated
  file pages instead of bare comma-separated names.
