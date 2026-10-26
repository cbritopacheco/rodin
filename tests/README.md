# How to build and run the tests

## Suite layout

| Directory | Purpose | CMake enable option |
| --- | --- | --- |
| `unit/Rodin/<Module>/` | Component behavior, contracts and regression tests. | `RODIN_BUILD_UNIT_TESTS` |
| `manufactured/Rodin/` | Fixed-mesh numerical checks with known solutions. | `RODIN_BUILD_MANUFACTURED_TESTS` |
| `convergence/{h,p,hp,isoparametric}/` | Multi-resolution numerical verification. | `RODIN_BUILD_CONVERGENCE_TESTS` |
| `benchmarks/Rodin/<Module>/` | Performance measurements. | `RODIN_BUILD_BENCHMARKS` |
| `installation/` | Installed-package consumption checks. | Dedicated CI configuration; see its [guide](installation/README.md). |
| `nafems/` | Published numerical verification cases and their qualification records. | See the [case catalog](nafems/README.md). |

Optional backends require their corresponding library configuration as well.
CTest labels describe a test's category and execution requirements. `slow` is
an additional label within its existing directory; it does not introduce a
parallel directory tree. MPI tests retain their module location and explicit
rank registrations.

## Add a test

For another case in an existing executable, add the GoogleTest case and add any
new source file to that executable's CMake source list. For a new executable,
follow the neighboring target's links and register it with
`rodin_discover_tests`, defined in [CMakeLists.txt](CMakeLists.txt). Include new
directories from their parent's `CMakeLists.txt`. MPI/PETSc tests must retain
their backend initialization and launcher pattern.

Ordinary grouped lanes discover registered targets through the configured
CMake manifest. They do not need a workflow entry for each executable. Dedicated
convergence and rank partitions have explicit execution selectors, which must
be extended when a new case falls outside them. The
[registration and lane editing guide](../dev/ci-testing.md) gives the procedures,
file ownership and local checks for each kind of change.

## Build and inspect a selection locally

Configure a build with the required suite and backend options first, using the
repository's [build procedures](../doc/agents/workflows.md). From the repository
root, a configured unit build can select and run Geometry as follows:

```sh
python3 dev/ci_tests.py build --build build --group unit/Rodin/Geometry --profile unit
python3 dev/ci_tests.py list --build build --group unit/Rodin/Geometry --profile unit
python3 dev/ci_tests.py run --build build --group unit/Rodin/Geometry --profile unit
```

Use an actual group from the generated plan, as described in the registration
guide. An optional-backend target is absent when that backend is disabled.
The runner executes separate CTest invocations serially; an MPI invocation
still uses its registered ranks. For direct CTest runs, use
`ctest --test-dir build -j1 --output-on-failure` with a suitably scoped label
or name filter.

## How CI builds and covers these tests

Every execution lane restores compatible shared builds. Common library objects
are compiled once per compiler/backend/build-type/threading/sanitizer
configuration. Test executables are then compiled by their physical registration
directory, with bounded numbered subdivisions for large directories. Regular,
Slow, geometry and MPI-rank selections can reuse an executable while retaining
their separate runtime checks and counters.

Each eligible lane saves its own LCOV report and receipt. Family publishers
validate those contributions and submit incremental unions to Codecov. The final
gate requires every expected current-commit contribution and preserves zero-hit
compiled production lines before enforcing 95% line coverage. Test sources,
benchmarks and `src/Rodin/Test` are excluded from the coverage denominator;
their execution contributes hits in production code. Sanitizer and installation
checks remain required but do not contribute coverage.

The [workflow architecture](../.github/workflows/README.md) explains the job
graph, artifacts, failure diagnosis and publication checks.

## Numerical convergence tests

[`convergence/`](convergence/README.md) contains multi-resolution validation of
theoretical finite-element rates. Enable it with `RODIN_BUILD_CONVERGENCE_TESTS`;
tests carry the `convergence` CTest label. See its guide for the refinement axes
and acceptance rules.

## Benchmarks

[`nafems/`](nafems/README.md) is reserved for published numerical verification
cases and the Rodin Verification Manual. It is distinct from `benchmarks/`,
which measures performance; see the [benchmark guide](benchmarks/README.md).
The NAFEMS case catalog currently lists candidates only; no NAFEMS result is
claimed yet.
