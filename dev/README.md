# Developer Tooling

The `dev/` directory contains repository-maintenance tools: ratchets,
linters, aggregate translation units, baselines, and helper scripts used by
developers and CI. These files do not define Rodin's numerical API; they keep
the tree reviewable while large historical cleanups are handled incrementally.

Use this directory when you need to answer one of these questions:

- Did my change introduce a new documentation warning?
- Did my change add a new naming or house-style violation?
- Are the lines I touched formatted?
- What baseline should shrink after I fixed old warnings?
- What aggregate translation unit should clang-tidy analyze?

## Ratchets and Baselines

### Production Coverage and Template Audit

`RodinCoverage` merges a zero-hit baseline from compiler notes with runtime
hits, then extracts only `src/Rodin/` records. Compiled source units that no
test executes remain in the denominator. Tests, benchmarks,
examples, library-side test helpers (`src/Rodin/Test/`), and dependencies
may execute production code but never enter the
coverage denominator. The regular and slow CI reports are merged by Codecov.
Optimized optional-backend unit jobs also save audited production coverage
artifacts in GitHub Actions for inspection and offline merging.

```sh
cmake --build build --target RodinCoverage
python3 dev/audit_coverage.py build/RodinCoverage.info --output build/coverage-audit.json
```

The audit unions line hits across input reports, rejects nonproduction
records, lists uncovered production lines, and inventories source files absent
from the report, with template declaration locations. CI uploads this JSON
for each coverage configuration. Absence is a review candidate: declarations,
static tables, disabled optional modules, and never-instantiated templates can
all be absent. The reported percentage alone cannot distinguish them.

Numerical specialization tests must instantiate supported scalar, storage,
geometry, and backend variants explicitly. Current audit findings include:

- Spatial matrices: real and complex arithmetic for all shapes, including empty matrices, and
  compatible products through dimension three, plus rectangular pseudoinverses.
- Iterative solvers: dense/sparse real MINRES, and dense/sparse real/complex
  least-squares CG and IDRS, with cold and nonzero initial guesses.
- Optional direct solvers: locally configured UMFPACK, CHOLMOD and SPQR.
  CHOLMOD/SPQR tests instantiate real and complex Hermitian systems and verify
  fresh solves after copy, move and polymorphic clone. MUMPS and ParU also
  pass their real sparse numerical and lifecycle suites, including assembled
  problems on all seven supported cell geometries. MUMPS checks retained
  symbolic/numeric stages and symmetric factorization; ParU checks index
  conversion and repeated right-hand sides. Both also rebuild independent
  resources after copy, move and polymorphic clone while preserving settings.
- P1 scalar derivatives: real/complex affine fields on every geometry,
  including face traces. Instantiation exposed obsolete basis, coefficient,
  and optional-attribute APIs, now covered by numerical regressions.
- Periodic constraint assembly: real/complex partition-of-unity weights and
  preservation of assembled maps through copy and move. These tests check the
  constraint API; they do not certify periodic enforcement in a PDE solver.
- Face normals: boundary orientation and both interior traces in two and three
  dimensions; missing and nonmatching trace attributes raise exceptions.
- Solid mechanics: Yeoh energy/stress/tangent checks in two and three
  dimensions, and P1 follower-pressure residual/tangent checks on segment and
  triangle faces. Higher-order pressure spaces are not claimed by this API.
- Nonlocal multi-field Eigen problems: real/complex dense/sparse storage,
  checked against exact separated double integrals with independent trial/test
  restrictions, repeated assembly, copy/move/assignment and polymorphic clone.
  Four manufactured mass-plus-nonlocal solves recover prescribed constant fields.
  Two additional dense real/complex solves combine preassembled diagonal and
  off-diagonal blocks, preassembled loads, nonlocal couplings and exact boundary
  elimination; both reject an injected preassembled-matrix coefficient error.
  These checks run through sequential and OpenMP backends. Constructor
  instantiation exposed unused real-only tuple assembler allocations; removing
  those obsolete private caches enables the existing generic assembler paths.
  Lifecycle regressions also exposed nonassignable backend state and loss of
  the formulation on copy/move, now retained explicitly. This does not extend
  the separate direct tuple-assembly API beyond its existing specializations.
- Standalone cardiac physics evaluators: exact midpoint kinematics, passive
  stress chain-rule differentiation with respect to the new displacement,
  viscous stress and its midpoint derivative, and every current/previous valve
  branch including transition equalities. Independent finite differences and
  analytic identities reject velocity, tangent and valve-slope perturbations.
- Distance wrappers: real local P1 signed/unsigned Eikonal intervals, missing
  interfaces, exact P1 Poisson nodal values and quadratic refinement, P2 H1
  exact quadratic recovery, region-restricted signed Poisson overloads, and
  analytic affine Rvachev/Spalding-Tucker normalization. Rvachev instantiation
  exposed an obsolete gradient-output API; the repaired call passes physical
  value checks. Perturbed signs, loads and normalization factors are rejected.
- Affine cell geometry cache: exact skew-triangle Jacobians, inverses, metrics
  and P1 gradients in both orientations, independent per-cell measures and
  index maps, singular-cell rejection, and enforcement of the documented
  planar-triangle scope. Including this previously absent header exposed
  silent acceptance of quadrilateral cells; a failing regression preceded
  the input validation fix.
- Hilbert H1 extension: scalar/vector P1 cosine manufactured solutions with
  homogeneous Neumann data, scalar refinement, and exact affine Dirichlet
  data through repeated solves. Instantiation exposed an obsolete problem
  signature, solver construction and vector range; the repair uses the
  stored trial/test functions and preserves added boundary conditions.
- Scotch partitioning: valid assignments for cells and faces, repeated
  repartitioning, empty/unused lifetimes, invalid count/dimension checks,
  and independently owned custom strategies. An unused partitioner was
  observed aborting before the lifecycle repair. Graph/strategy resources
  now belong to each call; a saved strategy string preserves value semantics
  rather than aliasing the caller's opaque C allocation.
- Eigen 3.4 MINRES does not compile for complex scalars. The Rodin header's
  generic scalar documentation does not establish dependency support.
- IDRstab's wrapper names `Eigen::IDRSTABL`, which is absent in Eigen 3.4's
  unsupported iterative solver module. This remains an API/dependency gap,
  rather than a variant certified by the runtime coverage percentage.
- Nonlocal `Potential` quadrature currently specializes matching P1 trial/test
  families. Instantiating a P0 trial with a P1 test fails because its
  `QuadratureRule` is undefined; mixed P0/P1 nonlocal forms are an API gap,
  rather than a variant certified by the existing P1 nonlocal tests.

These findings are not exclusions from production coverage. A complete
template audit requires checking the supported variants, not merely adding
header includes or explicit instantiations without numerical assertions.

Several checks use a ratchet pattern. A ratchet records the current known
backlog in a committed baseline file, compares the current run against that
baseline, and fails only on new findings. When old findings are fixed, rerun
the check with `--update-baseline` so the baseline shrinks.

Never grow a baseline to silence a new finding. Fix the finding instead, or
discuss why the rule should change.

### Doxygen Warnings

`check_doxygen_warnings.py` collects Doxygen warnings and audits generated XML
for parameter and return coverage, then compares the combined findings against the
normalized entries in `doxygen_warnings.baseline`.

The normal pass preserves Doxygen's warning policy. A separate XML pass enables
full extraction, including private and static helpers, and requires descriptions
for every function parameter and non-void return. This catches unnamed inputs,
deduction-guide parameters, and documentation hidden by `@internal`. Constructors,
destructors, deleted functions' return types, and unnamed variadic ellipses do not
require return or parameter descriptions where no such value exists. `@retval`
descriptions also satisfy return coverage. The `--log` option checks an existing
warning log only; run Doxygen through the script for the complete XML audit.

```sh
python3 dev/check_doxygen_warnings.py
python3 dev/check_doxygen_warnings.py --update-baseline
python3 dev/check_doxygen_warnings.py --doxygen /path/to/doxygen
python3 dev/check_doxygen_warnings.py --log /path/to/warnings.log
```

The baseline records the Doxygen version that produced it. Warning sets are
not compared across different versions because Doxygen changes wording and
diagnostics between releases.

### clang-tidy Naming

`check_clang_tidy.py` runs the naming rules from `.clang-tidy` over
`TidyTU.cpp` and compares findings against `clang_tidy.baseline`.

```sh
cmake -S . -B build -DRODIN_BUILD_TIDY=ON
python3 dev/check_clang_tidy.py --build-dir build
python3 dev/check_clang_tidy.py --flags "<compiler flags>"
python3 dev/check_clang_tidy.py --flags "<compiler flags>" --update-baseline
```

`TidyTU.cpp` is an aggregate translation unit that includes Rodin headers so
clang-tidy can inspect the library consistently. The CMake-backed invocation
also checks optional modules enabled in that build with their actual dependency
include paths and feature definitions.

### House Style

`style_lint.py` checks conventions that clang-format and clang-tidy do not
express:

- Boost Software License header blocks
- path-derived include guards
- no `#pragma once`
- Doxygen `@brief` coverage in headers
- `/** ... */` for multiline Doxygen documentation; `///` for single-line
  documentation (`///<` remains available for trailing member documentation)
- PETSc includes only under `src/Rodin/PETSc/`

```sh
python3 -m unittest discover -s dev -p test_style_lint.py -v
python3 dev/style_lint.py
python3 dev/style_lint.py src/Rodin/Variational
python3 dev/style_lint.py --update-baseline
```

Existing findings live in `style_lint.baseline`. As with every ratchet,
update it only to shrink the backlog.

### Changed-Line Formatting

`check_format.py` runs `git-clang-format` against the merge base and only
checks lines touched by the current change.

```sh
python3 dev/check_format.py
python3 dev/check_format.py --base origin/develop
python3 dev/check_format.py --binary clang-format
```

To fix reported hunks locally:

```sh
git clang-format <base>
```

## Dependency Graphs

`include_dependency_graph.py` builds an include graph for the codebase. It is
useful when checking whether a header-level change introduces an unwanted
dependency or when planning a refactor across modules.

Install the Python packages in `Requirements.txt` before running graph tools:

```sh
python3 -m pip install -r dev/Requirements.txt
python3 dev/include_dependency_graph.py --help
```

## CI Status

Static-analysis reports are published for the main branches:

| Branch  | Static Analysis |
|:-------:|:---------------:|
| master  | [![Static Analysis](https://github.com/cbritopacheco/rodin/actions/workflows/StaticAnalysis.yml/badge.svg?branch=master)](https://cbritopacheco.github.io/rodin/static/dev/refs/heads/master/) |
| develop | [![Static Analysis](https://github.com/cbritopacheco/rodin/actions/workflows/StaticAnalysis.yml/badge.svg?branch=develop)](https://cbritopacheco.github.io/rodin/static/dev/refs/heads/develop/) |

## What Belongs Here

Keep repository-wide developer utilities in `dev/` when they are meant to be
run directly by maintainers or CI. Prefer small scripts with explicit command
line interfaces, deterministic output, and committed baselines when the check
is intentionally incremental.

Do not place generated reports, experiment outputs, meshes, solver logs, or
one-off scratch files here. Those belong in a scratch directory and should not
be committed.

Every resource file belongs under `resources/` and must be tracked through Git
LFS. Before committing a new resource file, check `git check-attr filter --
<path>` and confirm the intended state with `git lfs status`.
