---
name: build-and-test
description: Build Rodin and run the right test executables. Use when compiling the library, running unit/manufactured/benchmark tests, or verifying a change.
---

# Build and test Rodin

## Build

```sh
git submodule update --init --recursive   # first clone only
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
```

A configured `build/` tree usually already exists — prefer building just the
target you need:

```sh
cmake --build build -j --target RodinVariationalH1Test
```

Target names equal executable names. An ASAN tree may exist at `build-asan/`
(use it for memory bugs, e.g. the flaky PETSc metricDot crash — see
doc/agents/petsc.md).

## Test

ctest is NOT wired up (`ctest -N` lists 0 tests). Run GoogleTest executables
directly; the path mirrors the source tree:

```sh
build/tests/unit/Rodin/<Module>/Rodin<Module><Component>Test
# e.g.
build/tests/unit/Rodin/Variational/RodinVariationalH1Test
build/tests/unit/Rodin/Solver/RodinNewtonSolverTest
build/tests/unit/Rodin/Math/RodinMathMatrixTest
```

Use `--gtest_filter=<Suite>.<Name>` to narrow, `--gtest_list_tests` to
enumerate.

Suites:
- `tests/unit` — fast, per-module; always run the executable(s) covering the
  code you touched.
- `tests/manufactured` — convergence and targeted regressions. Mandatory for
  assembly/solver changes (PETSc assembly regressions live in
  `tests/manufactured/Rodin/PETSc/TargetedAssemblyTest.cpp`).
- `tests/benchmarks` — Google Benchmark, target `RodinBenchmarks`; use
  `--benchmark_filter=<regex>`. Used as the measurement oracle for
  performance claims.

## Pass criteria

- The affected unit test executable exits 0.
- Derivative/term changes: the FD-consistency tests still pass (and at P2,
  not only P1, for Adaptation terms).
- Performance changes: identical numerics (iteration counts, final
  energies/fits) vs a baseline run, plus the benchmark delta.
- PETSc assembly changes: remember CI runs PETSc 3.19 — local green is not
  the gate (doc/agents/petsc.md).
