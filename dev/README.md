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
python3 dev/audit_coverage.py build/RodinCoverage.info --functions --output build/coverage-functions.json
```

The audit unions line hits across input reports, rejects nonproduction
records and records for deleted production files, lists uncovered production
lines, and inventories source files absent from the report, with template
declaration locations. CI uploads this JSON
for each coverage configuration. Absence is a review candidate: declarations,
static tables, disabled optional modules, and never-instantiated templates can
all be absent. Constant evaluation and optimization can also remove runtime
records from tested helpers. Check test registrations and instantiations
before interpreting an absent record as missing tests. The reported percentage
alone cannot distinguish these cases.

The optional `--functions` inventory merges function execution counters by
production file and symbol, supporting both `FN`/`FNDA` and `FNL`/`FNA` LCOV
records. Zero-count functions are candidates for inspection when shared line
coverage hides a compiled variant. Optimization and inlining affect these
counters; they do not establish whether a never-instantiated variant works.
This inventory does not change the measured line denominator.

Numerical specialization tests must instantiate supported scalar, storage,
geometry, and backend variants explicitly. Current audit findings include:

- Preassembled value constraints: dense/sparse P1 and quadratic H1 vector
  mass/load forms satisfy independent projector/lifting matrix and RHS
  identities with nonzero affine boundary values. Full solves recover the
  interpolated affine field, and targeted assembly preserves the unrequested
  side according to each backend's contract. Sequential and OpenMP one/four
  thread runs use explicit polynomial-exact quadrature; callable coefficients
  have unknown order, so the unconstrained mass/load identity is checked too.
  Reversing the dense Dirichlet lifting sign fails the RHS oracle.
- Single-field affine identification: dense/sparse P1 and quadratic H1 vector
  problems recover a prescribed constant field subject to an affine component
  tie on triangle and quadrilateral meshes. An admissible affine perturbation
  gives an independent exact mass-energy integral. Repeated assembly and
  targeted operations preserve the system; OpenMP checks one/four threads and
  supports targeted assembly, while sequential assembly rejects it. The latter
  originally cleared the RHS before rejection; vector initialization now
  follows constraint validation.
- Parent/submesh derivative transport: real/complex P1 and quadratic H1
  gradients, coordinate derivatives, Jacobians and divergences satisfy physical
  polynomial oracles on all seven cell geometries. Both inclusion and
  restriction are evaluated through geometric points, generic integration
  points and quadrature-associated integration points. Nonconstant quadratic
  derivatives distinguish physical evaluation locations. Doubling the
  integration-point gradient fails the complex H1 oracle.
  Local evaluation also exercises tabulated and generic integration points.
  The generic case crashed before repair: H1 gradient, coordinate derivative
  and divergence lacked the null-quadrature geometric fallback already present
  in the H1 Jacobian operator.
- Linear elasticity: distinct trial/test space instances and shared spaces
  reproduce analytic affine strain work for real/complex P1 and quadratic H1
  fields on all seven cell geometries. Field interpolation avoids coefficient
  layout assumptions, and repeated assembly preserves the result. Complex
  instantiation exposed eight unsupported mixed-scalar in-place spatial-matrix
  products; assigning the promoted product repairs both integrator families.
  Doubling the generic distinct-space shear term fails the strain-work oracle.
- Cardiac recruitment: all five linear/constant segments and the two inactive
  tails satisfy independent active-stiffness and active-stress residual and
  tangent formulas. The instantaneous relaxation limit checks its algebraic
  residual and both derivatives. Individual active-state entries agree with
  central finite differences; removing the relaxation strain derivative fails.
- Euclidean primitives: `Point2D` owns a two-component `Math::SpatialVector`,
  and displacements and segment directions use the same fixed-capacity type.
  The previous `Eigen::Vector2` storage was already fixed-size; this change
  standardizes the spatial API rather than claiming an allocation reduction.
  Float and real tests cover value semantics, invalid coordinate dimensions
  through native Alert exceptions, affine arithmetic, segment interpolation,
  barycenters, projections, circle intersections and rectangle corners.
  A doubled segment direction is rejected by the interpolation test.
  Direct Eigen inheritance and expression APIs are replaced by the point's
  coordinate accessors, spatial-vector view and owned arithmetic results.
- Form ownership: real/complex P1 and quadratic H1 forms reject foreign trial
  and test contributions without changing the assembled mass operator or load.
  Addition/subtraction and linear assignment exercise named-trial and unnamed
  test UUID diagnostics; P1 also checks global nonlocal contributions. Constant
  field interpolation supplies the independent unit-interval action oracle.
- Contextual Alert exceptions: class, member-function and namespace payloads
  retain exact unformatted messages and captured output. Literal and owned
  function-name variants propagate through the native exception base.
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
- P1/quadratic H1 scalar coordinate derivatives: real/complex affine fields
  on every geometry, including face traces. Direct and scalar-base expression
  views report the geometry-dependent total derivative degree. Instantiation
  exposed obsolete basis/coefficient APIs, missing P1 order metadata and
  missing shared-base order forwarding; compiler failures preceded the repairs.
  Manufactured mass solves with coordinate-derivative loads recover exact
  constants on two meshes for every geometry, including complex coefficients
  and OpenMP with one/four threads. Doubled derivatives fail both solves.
- P1/quadratic H1 one-sided gradient/Jacobian/divergence traces: a continuous piecewise
  affine field has different prescribed derivatives on either side of an
  aligned interface. Real and complex scalar/vector cases recover both traces
  on all seven geometries, including boundary faces, and reject missing trace
  attributes. Always choosing the first neighbor fails all three numerical
  oracles in both scalar fields and element families. All four cases in the
  dedicated trace executable pass. Coordinate derivatives, including copied,
  moved and cloned expressions, also recover the selected one-sided trace and
  reject first-neighbor mutations. The affine/order tests remain separate.
- Boundary-normal point transport: exact outward normals on retained half-box
  submeshes across all six 2D/3D geometries survive parent inclusion and
  restriction, copying, moving and polymorphic cloning. Excluded parent faces
  and unrelated coincident meshes are rejected. A halved normal fails the
  independent coordinate-plane oracle. All seven boundary-normal tests pass.
- Global integrator negation: real/complex constant kernels produce the exact
  negative outer product of independently integrated nodal basis functions.
  Direct expressions, lists, copy/move/clone and bilinear-form subtraction
  preserve this matrix and trial/test attribute restrictions. Instantiation
  exposed missing global negation wrappers; OpenMP execution then exposed
  omitted global contributions in standalone sparse forms. Both repairs are
  covered by failing regressions. Constant manufactured solutions on two
  meshes of every geometry pass in serial and with one/four OpenMP threads.
  Lost negation signs and doubled OpenMP accumulation fail the numerical tests.
- Single-pair PETSc MPI assembly: an isolated communicator recovers exact
  restricted nonlocal nodal moments and a preassembled restricted load.
  Repeated assembly replaces prior output; targeted operator/vector updates
  preserve the opposite side. All 13 MPI form cases pass on two ranks, and
  tripled nonlocal contributions fail the independent matrix oracle. This
  does not establish full-domain nonlocal integration across arbitrary shards.
- Standalone PETSc MPI global integrators: an independent restricted-moment
  test exposed omitted nonlocal contributions. The repaired path retains test
  DOFs across trial lookups and respects both attribute filters. Exact matrices
  survive repeated assembly and copied/moved/cloned assemblers. All 14 MPI
  form cases pass on two ranks; a Release manufactured constant solution also
  passes on all seven geometries and two mesh sizes, with repeated and
  copied/moved forms. Tripled contributions fail both oracles. These nonlocal
  tests use isolated communicators and do not certify arbitrary-shard quadrature.
- Mixed affine identification: real/complex dense/sparse mass systems recover
  prescribed constant fields with nonzero affine defects and fixed master
  data on all seven geometries and two mesh sizes. Repeated, copied and moved
  problems agree. All four manufactured cases pass in serial and with one/four
  OpenMP threads; discarding the defect fails all four cases.
- Sequential assembler lifecycle: real/complex dense, sparse, triplet and
  vector outputs retain exact two-cell mass/load values after copy, move and
  polymorphic cloning. Repeated execution replaces prior output. Both lifecycle
  cases and all 72 existing sequential tests pass; doubled assembly fails each
  storage and lifecycle oracle.
- HDF5 convenience exporters: scalar/vector/matrix values and dataset shapes
  agree with independent affine data on linear and curved meshes. Failed
  persistence and visualization opens preserve fields and leave missing
  directories absent. Substituting node output for cell output fails all three
  numerical exporters. All 99 HDF5 tests pass.
- PETSc HDF5 checkpoints: matrix-valued P1 fields retain independent coefficient
  values and range metadata on all seven geometries. Equal-size transposed
  ranges, incompatible coefficient counts and unsupported complex payloads
  are rejected before changing destination data. Scalar/vector/matrix fields
  also retain values after failed opens and unsupported streams, and accept
  legacy real-only files. All 18 tests pass with real PETSc; removing matrix
  shape validation fails on every geometry. Complex-specific assertions are
  conditional and require a complex PETSc configuration for execution.
- Shard-builder contracts: uninitialized/incompatible operations and missing
  vertices are rejected. Invalid direct connectivity preserves pending entities,
  insertion indices and topology; duplicate insertion preserves original state.
  Owner/halo views and attribute overrides survive finalization without mutating
  the parent. All three contract cases and all 31 existing shard tests pass;
  compiled corrupted attributes fail both direct and parent metadata oracles.
- Foreign form functions: real/complex local and global contributions with
  mismatched trial/test identities are rejected before changing integrator
  lists or assembled operators/vectors. Named and unnamed foreign trials
  exercise the diagnostic metadata paths.
- Grid-function failure contracts: scalar/vector MFEM and MEDIT failed opens
  preserve every coefficient and create no output directory. Geometrically
  coincident but unrelated mesh points are rejected through both geometric
  and integration-point evaluation; valid local evaluation and coefficients
  remain intact. All 38 grid-function tests pass.
- MEDIT malformed-input validation: numeric headers, entity counts, invalid
  records and truncation are rejected for all eight entity sections. A
  partially read vertex section is also rejected. Coordinates, topology and
  attributes of the destination remain unchanged; premature replacement is
  rejected by both new tests. All 30 mesh-loader tests pass.
- Taylor boundary-shift policy: real affine P1 fields recover exact boundary
  values after inward shifts on every boundary face of all seven geometries.
  Physical offsets, signed time-step caps, accumulated correction and repeated
  completed hits are checked independently. All nine Lagrangian tests pass;
  outward nudges, reversed correction and enlarged time caps are rejected.
- Distinct trial/test spaces: 21 quadrature cases check every mass and
  derivative entry against the independent generic evaluator for real/complex
  scalar/vector P1 and quadratic H1 spaces on all seven geometries. Quadratic
  H1 trials with linear tests also exercise rectangular operators. Bound-cell
  accessors retain the supplied cell. Targeted coefficient mutations fail all
  seven P1 and seven mixed-order H1 cases; seven unaffected same-order H1 cases
  remain green. All 51 existing quadrature-expression tests also pass.
  Real coefficients are checked in both scalar fields; complex coefficients
  use complex-valued spaces. Complex coefficients on real-valued spaces are
  rejected by existing static checks and are not certified by these tests.
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
- PETSc MPI mixed nonlocal problems: exact four-block P1 matrices with
  independent trial/test attribute restrictions, missing attributes and
  repeated assembly. Coupled mass-plus-nonlocal solves recover two prescribed
  pairs of constant fields with affine cross-field boundary constraints.
  Both tests reject an injected nonlocal coefficient error. These checks use
  isolated self communicators, including concurrent execution on two ranks;
  they certify that specialization, not complete-domain nonlocal integration
  across arbitrary distributed shards.
- PETSc mixed problem lifecycle: constrained mass-plus-nonlocal systems
  recover exact constant fields after copy, move, assignment and polymorphic
  cloning in sequential, OpenMP and isolated MPI contexts. Affine manufactured
  fields are recovered on three interval meshes. Instantiation exposed deleted
  base/backend assignment calls; mixed problems now independently own and copy
  their assembly strategy. Sequential and OpenMP numerical tests reject an
  injected nonlocal coefficient error.
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
- Retired affine cell geometry cache: earlier checks exercised exact triangle
  geometry and rejected unsupported input. PR #367 removes the unused cache
  and analytic adapters; the coverage branch preserves those deletions and
  removes cache-only tests. Active cell-deformation checks remain. Coverage
  of the merged source requires fresh compiler notes and runtime profiles.
- GLL nodes and angle wrappers: an absent optimized runtime record does not
  imply missing tests. `Variational/H1/GLLTest.cpp` checks node locations,
  ordering, symmetry and mapping; `Math/RadTest.cpp` checks angle arithmetic.
  Declaration-only headers, type aliases and compile-time traits likewise
  require inspection rather than fabricated runtime coverage.
- Attribute-restricted mesh measures: analytic half-unit-box area, volume and
  exterior measures for all six 2D/3D cell geometries, including combined,
  empty and missing attribute sets. Local and MPI suites pass; the distributed
  checks include an empty root at two/four ranks. Both tests reject doubled
  production measure accumulation.
- MPI uniform-grid validation: inconsistent shape dimensions and zero/one
  vertex extents are rejected on every rank. A uniform point grid retains one
  global owner and zero-dimensional ghost replicas. The nine-test mesh suite
  passes at one, two and four ranks.
- XDMF reference-node vocabulary: linear nodes agree with the geometry module's
  reference vertices for all positive-dimensional cells; point nodes have no
  reference coordinates. Point topology metadata and invalid-geometry rejection
  are checked explicitly. All 99 HDF5 tests pass, and displaced reference nodes
  fail the new oracle.
- Boundary attribute tracing: selected exterior faces inherit mapped cell
  attributes while unselected and interior faces retain theirs. Repeated
  tracing and missing mappings preserve the result; an incorrect mapped
  attribute is rejected by the exact oracle.
- MMG mesh-copy policies: sizing parameters, level-set base references,
  material split rules and their inverse lookup, and input/output names retain
  independent ownership after source mutation and destruction. Native material
  queries recover the prescribed parent and split references. All 110 MMG
  tests pass; copied sizing and material-reference mutations fail the new test.
- Flow boundary starts: signed inward characteristics on every boundary face
  of all seven positive-dimensional cell geometries agree with exact constant
  velocity trajectories. Outward starts stop immediately at the same physical
  point with the full remaining time. Zero-time interior-face starts also
  preserve the input point.
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
