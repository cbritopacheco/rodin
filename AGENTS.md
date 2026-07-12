# Agent Guide for Rodin

Rodin is a modular C++20 finite element framework oriented toward shape and
topology optimization (variational form language, meshing/remeshing, level
sets, PETSc/MMG integrations). This file is the entry point for AI coding
agents (Claude, Codex, Cursor, ...). Read it fully; follow the pointers only
when the task touches that area.

## Knowledge base

The knowledge base lives in `doc/agents/` and is hierarchical — index at
`doc/agents/README.md`. Read top-down, stop at the depth your task needs:

- Level 1 (always): `doc/agents/philosophy.md` — **read before writing any
  code**: the typed graded expression algebra the library models, node
  patterns, minimality rules, surface style. And
  `doc/agents/conventions.md` — the hard rules.
- Level 2 (orientation): `doc/agents/architecture.md` — the full module map.
- Level 3 (per-domain, open what the task touches): `core.md` (vocabulary,
  Alert, Threads, Context, Math), `geometry.md` (mesh model),
  `variational.md` (form language, FES, Problem), `solvers-assembly.md`,
  `physics.md` (Solid, Heart, level-set toolkit, Adaptation),
  `integrations.md` (MMG, MPI, IO, bindings), `petsc.md`, `testing.md`.
- Level 4 (foundations, `doc/agents/theory/`): the mathematics tied to the
  code that embodies it — `variational-formulation.md`,
  `finite-elements.md`, `meshes-and-geometry.md`, `shape-optimization.md`,
  `continuum-mechanics.md`. Open these when writing new formulations or
  when the *why* of a design matters.

Step-by-step procedures (build, run examples, debug) live in
`.claude/skills/*/SKILL.md`. Claude Code loads them automatically; other
agents should read them as plain markdown when the task matches.

If `graphify-out/` exists, `graphify-out/GRAPH_REPORT.md` and
`graphify query "<question>"` can answer cross-module structure questions
faster than grep.

## Build

```sh
git submodule update --init --recursive   # first time only
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
```

Key CMake options (all default ON): `RODIN_BUILD_EXAMPLES`,
`RODIN_USE_PETSC`. An existing configured `build/` tree is usually present;
prefer incremental builds of the target you need
(`cmake --build build -j --target <name>`) over full rebuilds.

## Test

Tests are GoogleTest executables (ctest is NOT wired up — `ctest -N` shows 0
tests). Run the binary for the module you touched:

```sh
build/tests/unit/Rodin/Adaptation/RodinAdaptationTargetMatrixOptimizationTest
build/tests/unit/Rodin/Variational/RodinVariationalH1Test
```

Naming: `build/tests/unit/Rodin/<Module>/Rodin<Module><Component>Test`.
Suites: `tests/unit` (fast), `tests/manufactured` (convergence/regression,
includes PETSc assembly regressions), `tests/benchmarks` (Google Benchmark,
target `RodinBenchmarks`).

A change is not done until the affected unit test executable passes and, for
assembly/solver changes, the relevant manufactured tests pass too.

## Non-negotiable rules (summary — details in doc/agents/conventions.md)

1. **PETSc error handling:** `assert(ierr == PETSC_SUCCESS)` after each call
   is the house idiom. Do not introduce checking macros or convert existing
   asserts.
2. **Minimal, behavior-preserving changes.** Do not bundle speculative
   optimizations or fast paths into a requested change; propose them
   separately. For performance work, verify identical numerics (iteration
   counts, final energies) against a baseline run.
3. **CI builds against PETSc 3.19** (local dev is typically newer). PETSc
   assembly changes can pass locally and fail CI — check 3.19 semantics.
4. **No hard geometric projections inside mesh-optimization solves** — surface
   fitting is always a smooth penalty term (see `doc/agents/conventions.md`).
5. **Internal variables are first-class DOFs** in Solid — no per-quadrature
   Schur condensation (see `doc/agents/conventions.md`).

## Housekeeping

- Example runs dump `*.h5` / `*.xdmf` / `*.log` output into the CWD. Run
  examples from a scratch directory, and never commit these artifacts.
- Branches: `master` is the default; active development happens on
  `module/*`, `model/*` topic branches off `develop`. Do not commit or push
  unless asked.
