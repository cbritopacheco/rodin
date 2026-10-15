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
coverage denominator. The regular, slow and optional-backend CI reports are merged by Codecov.
Unsanitized Debug and Release optional-backend unit jobs save audited production
coverage artifacts in GitHub Actions and upload only the filtered LCOV report.
Sanitizer configurations never contribute coverage profiles. Debug compiler
notes retain assertion-enabled implementations that optimized builds can erase.
CI audits include function inventories for supported-template review.

`audit_source.py` separately inventories source function bodies and exception
sites without preprocessing or instantiation. Install its optional developer
dependencies (`tree-sitter`, `tree-sitter-cpp`) outside the project, then run:

```sh
python3 dev/audit_source.py --output /tmp/rodin-source-inventory.json
```

The optional `--line-inventory` accepts JSON mapping repository-relative source
paths to line-number/hit dictionaries. This permits comparing bodies with the
validated LCOV union, including bodies absent from otherwise measured headers.
Absence is an investigation candidate: inlining and constant evaluation can
remove attribution. Parser recovery locations and unexpanded macro bodies
remain explicit review limitations. This inventory does not certify all template
specializations or produce an alternative coverage percentage.

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

- Single-field PETSc targets: scalar, vector and rectangular matrix P1 and
  quadratic H1 fields match exact interval mass matrices, load moments and
  nonzero boundary values. Coordinate/component labels recover the nodal map
  without assuming coefficient ordering. Complex cases use non-real mass
  coefficients and boundary values. Direct and preassembled contributions,
  cloned backends, repeated full assembly and separate matrix/RHS updates
  preserve their exact contracts. All six scalar/backend/thread combinations
  pass in Release, assertion-enabled and AddressSanitizer builds; fresh CMake
  checks also pass individually. Every range/degree variant rejects doubled
  partial boundary diagonals in both scalar and backend modes. The GCC/AArch64
  P1 vector-load warning reflects its unconstrained byte-sized length bound
  in Release; the configured dimension two is within capacity three, and the
  focused path also passes O3/NDEBUG AddressSanitizer checks in both scalar modes.
- Distributed mixed PETSc targets: the local exact block oracle now gathers
  owned matrix rows and vector entries across the system communicator. Real
  and complex runs at one, two, three, four and eight ranks check full,
  matrix-only and RHS-only updates, including an empty root shard and
  additional empty ranks. All 20 Release/assertion-enabled MPI runs and all
  twelve retained local backend/thread runs pass individually. Wrong partial
  boundary diagonals fail both scalar variants at two and four ranks. Fresh
  CMake builds register the five MPI configurations and the local target;
  all eight fresh-binary checks also pass serially.
- Mixed targeted PETSc assembly: an exact two-field interval matrix combines
  local mass blocks, constant-kernel nonlocal blocks, preassembled cross-field
  terms and loads, and distinct nonzero Dirichlet values. Full assembly checks
  column lifting and the recovered constant fields; matrix-only and RHS-only
  updates check every entry and preserve the other target, through cloned
  backends and repeated assembly. All six real/complex sequential/OpenMP runs
  pass in both Release and assertion-enabled builds, with OpenMP at one and
  four threads. Doubling the partial boundary diagonal fails each of the four
  scalar/backend variants. The target is selected by explicit real/complex CI.
  A fresh coherent CMake PETSc/OpenMP build also passes its sequential and
  one-/four-thread OpenMP cases individually, with registration checked.
- Index-generator ownership: bounded, vector and set generators preserve every
  remaining suffix through polymorphic copies and moves, including empty and
  exhausted sequences. Bounded copies previously discarded the cursor;
  vector/set copies retained iterators into the source storage and were
  observed using freed memory under AddressSanitizer. Copies now preserve
  position within their own containers. Borrowed pointer, vector-iterator and
  forward-list-iterator variants retain independent cursors while their
  external range stays alive. All 34 cases pass individually in Release and
  with assertion-enabled, leak-enabled AddressSanitizer.
- Distributed runtime space interfaces: real/complex scalar, vector and
  rectangular matrix P0, P0g, P1 and quadratic H1 spaces retain independent
  topology counts and complete DOF maps through copies and moves. Runtime
  metadata, disjoint contiguous ownership, local/global round trips,
  supported traces and complete global index coverage are checked on
  triangular and tetrahedral meshes, including an empty root shard. All 40
  new Release/Debug MPI runs pass at one, two, three, four and eight ranks;
  20 existing interpolation runs pass at their registered rank counts.
  Shifting production P0g entity indices fails all four new interface cases.
- Native MEDIT loading: scalar, vector and rectangular matrix P1/quadratic
  H1 fields are checked against analytic values on all eight cell families
  in real and complex arithmetic. Malformed metadata and truncated payloads
  preserve existing coefficients; repeated failures have stable member
  diagnostics and valid input recovers on the same loader. P1 now validates
  and stages its complete payload before committing. H1 uses native
  interpolation of the parsed P1 field, replacing obsolete scalar-only DOF
  injection and covering hexahedra. Truncated P1 overwrite and the Debug H1
  vector node-count assertion were reproduced before repair. All twelve
  cases pass separately in Release and Debug; scaled production payloads fail
  every analytic oracle.
  Fresh CMake builds pass all twelve new cases and all 22 existing P1 and
  spatial-matrix I/O cases, each executed separately. Explicit real/complex
  PETSc CI also builds and selects the I/O, affine-constraint and nonlocal
  regression targets.
- PETSc MEDIT vertex output: scalar, vector and rectangular matrix P1 and
  quadratic H1 fields retain their analytic affine values in real and complex
  builds. Public printer dispatch emits the expected version, dimension,
  vertex count and field types, with exactly one complete component row per
  vertex and a terminal marker. Repeated printing is checked. All twelve
  cases pass serially; doubling printed production values fails every oracle.
  The registered CMake target also builds and passes against freshly compiled
  current production libraries.
- Zero-dimensional projection: a point reference vertex now has zero active
  coordinates, matching its centroid and reference dimension. The CI Debug
  assertion in `EveryReferenceVertexIsFixed` was reproduced before repair;
  the test now explicitly checks every geometry's vertex dimensions. All 22
  projection cases pass serially with assertions enabled and in Release.
- SWIFT termination: a rejected line-search trial and an insufficient active
  hinge Newton budget preserve the accepted displacement and cannot report
  geometric success. The inner limit retains its final stationarity residual
  and typed stopping reason. Selecting unavailable MUMPS now identifies
  `setParameters` and preserves existing controls; its diagnostic regression
  fails before repair. The complete configured solver executable passes
  with every case run separately.
- Native single-field nonlocal assembly: sparse/dense real/complex systems
  combine constant-kernel double integrals with closed-form interval mass
  matrices and load moments. All trial/test region combinations include
  unrestricted, selected, unmatched and absent cell attributes. Repeated
  execution through runtime assembler clones survives destruction of their
  original owners; targeted updates preserve the other side of the system.
  Sequential and OpenMP variants pass serially, with one and four OpenMP
  threads. Bypassing production trial-region filters fails the new oracles.
- Matrix derivative traces: real/complex matrix P0, P0g, P1 and quadratic
  H1 fields exercise boundaries, independently selected interior sides,
  copies, moves, clones and parent/submesh transport on all seven cell
  geometries. Analytic affine gradients and zero constant-space gradients
  agree entry by entry. Missing trace domains and unrelated meshes are
  rejected with member-function context; all eight diagnostic regressions
  fail before the contextual exception change. Perturbing a returned
  production gradient entry fails every matrix oracle. All twelve registered
  scalar/vector/matrix trace cases pass serially.
- Rank-one Jacobian forms: real and complex P1 and quadratic H1 retain
  analytic contractions through copies, moves, rebinding and polymorphic
  clones on intervals and triangles. Both shared and distinct finite-element
  spaces are checked. Complex coefficient regressions failed before repair:
  the test action now uses conjugation, and the shared-coefficient matrix
  uses an adjoint when completing its upper triangle. All four lifecycle
  cases and eleven existing real-valued rank-one cases pass serially. A real
  matrix coefficient with a complex Jacobian is rejected earlier by the
  current `Dot` range-type contract; this unsupported combination was audited
  separately and did not motivate a broader expression-language change.
- Runtime-polymorphic local quadrature: independently owned linear and
  bilinear integrator collections retain their cell bindings, explicit orders
  and analytic field contractions after their original concrete owners are
  destroyed. Rebinding and cloning are checked through the abstract assembly
  interfaces. Real/complex P1 and quadratic H1 include twelve form families
  on intervals and triangles, including weighted divergence. Doubling the
  generic local matrix entries used by that family fails both real and complex
  oracles; all four registered cases pass serially.
- Polymorphic finite-element spaces: runtime collections contain scalar,
  vector and matrix P0, P0g, P1 and quadratic H1 variants in both real and
  complex configurations. Public virtual metadata agrees with independent
  topology-based dimensions; every local/global DOF map agrees, contains
  distinct in-range entries and collectively covers the global index set.
  Supported traces and object-identity comparisons are checked through the
  base interface without assuming a high-order coefficient layout. Shifting
  the production global index by one fails both scalar configurations; all
  six registered lifecycle/interface cases pass serially.
- Local quadrature ownership: real/complex P1 and quadratic H1 rules retain
  explicit order and current polytope bindings through copies, moves and
  polymorphic clones. Eleven scalar/vector mass, gradient, Jacobian,
  divergence and load families are checked on interval and triangle meshes.
  Local coefficient contractions agree with independent integration of
  analytic affine fields; rebinding on successive cells checks cached state.
  One-dimensional divergence couplings and distinct-space vector Jacobians
  are instantiated explicitly. Doubling production local matrix entries is
  rejected in real and complex cases. All four registered cases pass serially.
- PETSc single-field nonlocal restrictions: real and complex P1 kernels
  reproduce the outer product of independently assembled restricted moments
  and the exact action on two different affine fields. All combinations of
  unrestricted, selected and unmatched trial/test attributes are checked,
  including zero operators and repeated full assembly. Sequential and OpenMP
  variants are instantiated explicitly; OpenMP runs at one and four threads.
- Higher-order nonlocal audit: `Integral(Potential(kernel, u), v)` does not
  compile for quadratic H1 because the existing nonlocal `QuadratureRule`
  specialization accepts only P1 spaces. This is an absent capability,
  not runtime coverage of a supported H1 variant; no new quadrature algorithm
  or coverage exclusion was introduced for the finding.
- Locator reference clipping: near-boundary queries on all seven
  positive-dimensional cell families return contained reference coordinates
  whose remapped physical residual satisfies the configured tolerance.
  Lower boundaries, simplex sum constraints, tensor-product bounds, wedge
  heights and pyramid side/apex bounds are checked, together with rejection
  outside the tolerance band. All 55 locator cases pass serially; skipping
  both production clipping calls fails the new oracle. Production sources
  and coverage exclusions remain unchanged for this batch.
- Matrix-space lifecycle: real/complex P0, P0g, P1 and quadratic H1 copies
  retain owned DOF and element maps after source destruction. Copy/move
  assignment replaces a different mesh and matrix shape consistently; later
  source reassignment leaves the copy independent. Interpolated constant
  matrix fields recover every entry, and public local/global DOF mappings
  agree. Zero/oversized row and column counts now identify the constructor
  and class through native member-function exceptions; the diagnostic checks
  failed before this change. Doubling production matrix-field evaluation is
  rejected in both scalar configurations.
- Embedded-surface face normals: affine, vertical-edge and folded patches
  select the outward conormal in each incident cell's tangent plane. A curved
  quadratic patch checks unit length and orthogonality to independently known
  surface normals and edge tangents at multiple points. The vertical-edge
  regression failed before repair; the previous edge-only construction did
  not determine the selected cell's tangent plane. The repaired path maps the
  face point into that cell before evaluating its Jacobian. A centroid-plane
  mutation fails the curved oracle. All four face-normal cases pass serially.
- MMG affine cuts: a planar surface retains unit area and a unit-length
  isovalue contour. A unit-cube cut checks both material volumes, every
  interface coordinate and unit interface area while changing labels and
  switching them back on the same discretizer. Reuse retained stale generated
  label mappings before repair; the derived mapping is now cleared per call.
  All 112 MMG cases pass serially; shifting the production cutting level
  rejects the physical-coordinate and measure oracles.
- Unresolved solver-template audit: standalone inclusion of
  `Solver/IDRSTABL.h` fails with the configured Eigen 3.4 dependency because
  that dependency does not define `Eigen::IDRSTABL`. Its absent LCOV record
  does not establish tested support. No solver substitution or coverage
  exclusion was introduced for this finding.
- PETSc single-field affine identifications: linear and quadratic vector spaces
  on triangles and quadrilaterals recover an exact constant mass solution and
  the independent admissible-direction energy. Real sequential and OpenMP
  assembly and complex sequential assembly are instantiated explicitly.
  Targeted assembly remains unsupported for these ties; rejection now precedes
  matrix/vector setup so the previously assembled system survives unchanged.
  The preservation regression failed before this repair. Doubling the
  production reconstruction offset is rejected by the exact-field and recovery
  oracles. Tests run serially, including OpenMP at one and four threads.
- MFEM mesh validation: the dimension keyword and bounded topological/spatial
  dimensions are validated before mesh replacement. Straight meshes require
  the separate spatial-dimension record specified by the
  [MFEM format](https://mfem.org/mesh-format-v1.0/); the former fallback for an
  omitted record was unsupported and corrupted memory. All 94 MFEM resource
  files inspected include that record. An embedded triangle checks distinct
  topological/spatial dimensions, exact coordinates and physical area.
  Physical-line diagnostics count comments and blank lines and restart when a
  loader is reused. Wrong-keyword and diagnostic tests failed before repair;
  the malformed-coordinate reproducer also aborted with an invalid free.
- Constraint maps: real/complex affine chains carry both coefficients and
  offsets, including elimination of fixed masters and cancellation of repeated
  masters. Finalization is idempotent, and copied maps remain independent.
  The affine-offset and fixed-master oracles failed before their repair.
  Invalid indices, conflicting and singular rows, and cycles are rejected;
  reset restores usable state. Region iteration matches indexed selection for
  cells, faces, boundary facets and interior interfaces on a two-dimensional mesh.
- Preassembled value constraints: real/complex dense/sparse P1 and quadratic H1 vector
  mass/load forms satisfy independent projector/lifting matrix and RHS
  identities with nonzero affine boundary values. Full solves recover the
  interpolated affine field, and targeted assembly preserves the unrequested
  side according to each backend's contract. Sequential and OpenMP one/four
  thread runs use explicit polynomial-exact quadrature; callable coefficients
  have unknown order, so the unconstrained mass/load identity is checked too.
  Reversing the dense Dirichlet lifting sign fails the RHS oracle.
  Complex variants use nonzero imaginary boundary and load data and verify
  complete complex coefficients through the same projection identities.
  Discarding imaginary prescribed values fails the complex RHS oracle.
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
  Debug CI exposed a constructor assertion that rejected the implemented
  distinct-space path. Both real/complex aborts were reproduced locally;
  the constructor now permits distinct spaces as the assembly code requires.
- Cardiac recruitment: all five linear/constant segments and the two inactive
  tails satisfy independent active-stiffness and active-stress residual and
  tangent formulas. The instantaneous relaxation limit checks its algebraic
  residual and both derivatives. Individual active-state entries agree with
  central finite differences; removing the relaxation strain derivative fails.
- Boundary extraction and submesh ownership: segment skins and zero-dimensional
  keep/trim selections reject unsigned-loop underflow regressions reproduced
  as crashes before the repair. Boundary copies retain coordinates, parent
  maps and attributes across all supported cell geometries, including vertex
  labels introduced implicitly through closure. Nested surface selections
  retain ancestry and point restriction after the original submesh dies.
  Counting measure is one per point. Five new cases and 28 existing mesh and
  submesh cases pass separately; the new cases also pass with assertions,
  ASan and leak checks. Doubling production boundary measure fails the oracle.
- Geometry metadata ownership: assigning sparse attributes over populated
  storage replaces both growing and shrinking extents, preserves engaged
  zero and maximum labels, and remains independent after source mutation
  and destruction. Interface tracing respects current face selection and
  unordered material pairs, preserves exterior faces, and skips missing
  material attributes. Embedded line and surface samples own cached
  Jacobians and reference coordinates after their inputs and formulas die;
  physical maps, pseudoinverses and Gram metrics agree with affine formulas.
  All three cases pass individually with assertions and sanitizer checks.
  Doubling the production embedded determinant fails the metric oracle.
- Matrix derivative and Jacobian lifecycle: real/complex P0, P0g, P1 and
  quadratic H1 owning clones preserve constant/affine derivatives on all
  seven cell geometries. Every coordinate direction agrees with its analytic
  derivative and with reconstruction from cloned basis expressions after
  their originals are destroyed. The first excluded and largest unsigned
  directions reject both field and basis evaluation through member-function
  exceptions while preserving the field. Diagnostic tests failed before the
  contextual repair. All cases pass individually with assertions and
  AddressSanitizer/leak detection; doubled production Jacobian and partial
  derivative values are rejected in both scalar configurations.
- Matrix divergence lifecycle: real/complex P0, P0g, P1 and quadratic H1
  expressions preserve their field or leaf identity through copy, move and
  owning clones. Constant and affine fields match analytic row divergence
  and independent basis reconstruction on all seven cell geometries. A valid
  three-column field in two dimensions rejects divergence evaluation through
  member-function exceptions without changing the field. Both diagnostic
  tests failed before the contextual exception repair. All cases pass
  individually with assertions, AddressSanitizer and leak detection enabled.
  Doubling production divergence is rejected in both scalar configurations.
- Polytope iterator ownership: cold, cached and dirty releases retain entity
  identity, ordered vertex keys and attributes after iterator destruction.
  Empty and exhausted iterators release null. Advanced cell, face and vertex
  iterators replace a generic iterator cache while preserving the remaining
  cursor. All seven cell geometries pass normal and assertion-enabled
  AddressSanitizer runs with leak detection, one case at a time.
- Physical point ordering: exhaustive finite grids in embedding dimensions
  one through three agree with an independent lexicographic oracle. Pairwise
  checks precede sorting and duplicate-key insertion. The original comparison
  failed in two dimensions because both early-exit conditions tested the same
  inequality; the second condition now rejects a larger leading coordinate.
- Euclidean primitives: `Point<Dimension, T>`, `LineSegment`, `Line`, `Sphere`,
  `Hyperplane` and `Box` use `Math::FixedSizeVector` and `FixedSizeMatrix`
  storage with compile-time dimensions. Existing planar names remain aliases.
  Scalar/determinant/cross-product paths handle dimensions one through three;
  higher-dimensional line calculations use column-pivoted QR with normalized
  direction columns. Object symbols contain QR matrices only for the tested
  higher dimensions, four and six. Zero-dimensional points and boxes have
  their unique-point and counting-measure contracts; primitives requiring a
  nonzero direction, distinct endpoints or nonempty sphere reject construction
  through native member-function exceptions.
  Float and real variants in dimensions zero, one, two, three, four and six
  cover affine arithmetic, projections, skew distances, unique intersections,
  tangent hyperplanes, closed-box operations and all box corners. Expressions
  are evaluated into owned fixed storage before assignment, including reversal
  of the point's own coordinates. Dynamic sizes and orientations are validated
  before Eigen conversion. Small and large direction scales preserve geometry.
  Release counterexamples exposed spurious coincident-line intersections in
  the Gram determinant calculation; direct small-dimensional formulas and QR
  rank checks repair them. A separate planar counterexample projected $(1,1)$
  to $(0,1)$ with strongly unequal coefficients; anchoring on the largest normal
  coefficient preserves the geometry. All 74 cases pass separately in Debug,
  fresh Release CTest, and assertion-enabled ASan/LSan runs. Doubling production
  projection and removing alias-safe assignment fail the respective float and
  real manufactured oracles.
  Production coverage retires all previous Euclidean implementation records
  and captures matching fresh compiler notes. The Euclidean batch measurement is
  32,001 of 33,867 lines (94.49%); Euclidean contributes 480 of 482 lines
  (99.59%). Compared with the previous mesh batch, Euclidean adds 359 measured
  lines and 357 hits; normal Debug instantiation also reveals 79 Alert formatter
  lines (30 hits) and four Unit lines (all hit). These changes explain the larger
  denominator and 51 additional uncovered lines. Tests, benchmarks, sanitizer
  profiles and mutation profiles remain excluded. This measurement covers the
  compiled float/real and dimension variants, not every possible instantiation.
- Debug instrumentation audit: existing independent Alert formatting, named-form,
  HDF5 and periodic-boundary oracles are also compiled without optimization.
  All four Alert cases, 35 applicable named-form cases, 99 HDF5 cases and both
  periodic cases pass individually; segment linear elasticity is explicitly
  skipped because its contract requires dimension two or three. The complete
  merged measurement becomes 34,213 of 36,606 lines (\(93.46\%\)). Relative to
  the Euclidean batch, 2,739 additional production lines become measurable,
  with 2,212 additional hits. These compiler notes remain in the denominator,
  including newly exposed unexecuted bodies and failure paths. Reusing the
  older denominator would conceal this instrumentation gap.
  The curved HDF5 node-attribute oracle now uses the exact affine reference
  field \(4r+9s\) on the triangle mapped by \((r,s)\mapsto(2r,3s)\), instead
  of reevaluating the quadratic transformation under test. The previous
  expected value differed by approximately \(1.07\times10^{-14}\) in Debug
  because that second evaluation introduced additional rounding. The original
  absolute tolerance \(10^{-14}\) is unchanged; doubled production exports
  fail the corrected oracle. Failed and mutation profiles are not merged.
  Subsequent normal Debug checks also pass individually: 47 P1-element cases,
  121 H1-element cases, 39 matrix cases, 35 mesh-format cases, six promoted
  algebra cases and 69 high-order field IO cases. The resource-dependent
  private builds use the existing hydrated resource directory, leaving the
  frozen production snapshot untouched. The combined audit now retains
  37,161 measured lines with 34,826 hits (\(93.72\%\)); another 477 hits are
  required for \(95\%\) at this denominator. The corrected HDF5 affine oracle
  also passes separately with Release optimization.
- File-format diagnostics: all known enum values retain their format names;
  unknown values are streamed numerically instead of passing a null string
  to the stream. Both regressions fail before the fix: the stream enters a
  failed state and grid-function exceptions lose the rejected format value.
  Scalar/vector fields retain their coefficients and create no output directory
  after unsupported load/save requests. The native member-function context is
  checked. Both cases pass individually in normal Debug, fresh Release CTest
  and assertion-enabled ASan/LSan builds; all 38 existing grid-function cases
  also pass individually in normal Debug. Every previous coverage record for
  the changed formatter header is retired before fresh compiler notes are
  captured. The matching-source merge covers 34,935 of 37,233 production lines
  (\(93.83\%\)); the formatter header contributes 14 of 14 lines. Tests,
  pre-fix failures and sanitizer profiles remain outside the measured sources.
- Region interpolation: real/complex P1 and P2 scalar/vector fields reproduce
  independent affine traces on a selected face or interface, preserve every
  unselected vertex, ignore false predicates and preserve the result on repeated
  interpolation. Selected submeshes accept retained parent-cell points and
  reject other parent cells while preserving coefficients. The fixtures compute
  the incidence relations required by H1 explicitly. All four typed cases pass
  separately in normal Debug, fresh Release CTest and assertion-enabled
  ASan/LSan builds. Routing faces to cells or interfaces to the boundary fails
  both real and complex numerical cases for each mutation. This batch adds 11
  production line hits without enlarging the denominator: 34,946 of 37,233
  lines (\(93.86\%\)). Mutation and sanitizer profiles remain excluded.
- MEDIT parser variants: raw-character ranges and string const iterators cover
  all 28 keyword spellings, signed/unsigned integer limits, empty and trailing
  input, structured declarations, vertex dimensions zero through three, and
  entity records with one through eight vertices. Four pre-fix checks reproduce
  the whitespace skipper joining `Triangles tail` into one word and accepting
  split declarations such as `Dim ension 3`. Keyword lexemes now require
  contiguous alphabetic characters while preserving surrounding whitespace
  and whitespace before numeric values. The native unexpected-keyword
  exception type is retained.
  All nine parser cases pass individually in normal Debug, fresh Release
  CTest and assertion-enabled ASan/LSan builds. All 35 public mesh-loader cases
  and 12 real/complex P1/P2 scalar/vector/matrix field-loader cases also pass
  individually in Release; the 12 field-loader cases additionally pass in
  normal Debug. The coverage merge discards every old line and function record
  for the three changed parser classes. Other records in the header are retained
  only after verifying byte-identical source lines and unchanged positions.
  The combined audit covers 35,079 of 37,311 production lines (\(94.02\%\));
  MEDIT contributes 464 of 473 lines. At this denominator, another 367 hits
  are required for \(95\%\). Tests, benchmarks, failed runs and sanitizer
  profiles do not enter the measured sources or normal coverage evidence.
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
- MMG cached adjacency and solution ownership: planar triangles, embedded
  surface triangles and tetrahedra retain independent cached neighbors after
  source destruction, checked against Rodin face incidence and exact geometric
  measures. The tetrahedron-only adjacency allocation was observed overrunning
  its buffer under AddressSanitizer before the repair. Named scalar/vector/tensor
  solutions preserve independent payloads and names; empty solutions copy and
  swap ownership correctly. All 115 MMG cases pass serially, the three new
  lifecycle cases also pass under AddressSanitizer, and deliberate adjacency
  and coordinate corruption fail the new oracle. Sanitizer and mutation
  profiles are excluded from the production coverage merge.
  Leak-enabled CI also exposed omitted cleanup of level-set base-reference
  arrays in the bundled MMG dependency. The material-copy and repeated volume
  cut leaks were reproduced locally; Rodin now releases the array before
  native mesh cleanup. All 115 cases pass individually with AddressSanitizer
  and leak detection enabled.
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
- The legacy `Variational/TraceOperator.h` wrapper cannot be included against
  the current templated `FunctionBase`. A fresh standalone compile, including
  `Variational.h` first, fails at its template/type mismatch. Its documentation
  example does not establish an available API; the absent-source audit records
  this wrapper separately from the runtime line percentage.
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
- PETSc KSP ownership: six pre-fix real-scalar copy/move cases across KSP,
  CG and GMRES fail with missing reference increments, nonempty move sources,
  and premature context destruction. Explicit copy and move constructors now
  retain or transfer counted context ownership, preserving shared handle
  identity and independent stored settings. The borrowed problem and explicit
  preconditioner matrix keep their existing lifetime contracts. A frozen
  pre-fix compile probe also rejects the const communicator query; casting
  the handle value fixes the qualifier error.
  All nine lifecycle/settings/field-split cases pass individually with real
  and complex PETSc 3.19.6, and with assertion-enabled real Debug builds.
  Four existing factor-reuse and two SNES state cases pass with real PETSc;
  the two SNES cases also pass with complex PETSc. Nine real lifecycle cases
  pass under AddressSanitizer with the repository's external MPI/PETSc leak
  detection policy (`detect_leaks=0`); explicit native reference-count oracles
  check ownership. This is not an LSan certification.
  Poisson manufactured fixtures use the native PETSc scalar field, including
  complex amplitude \(1-2i\). The original fixed-grid error bounds remain;
  an independent physical \(L^2\) check verifies order in \([1.8,2.2]\)
  over three refinements for triangles and quadrilaterals, with order-eight
  and order-ten error quadratures agreeing to relative \(10^{-8}\).
  All fourteen manufactured cases pass separately in real and complex Release
  builds. An isolated doubled-solution mutation fails all nine lifecycle
  cases and both physical convergence cases; mutation and sanitizer profiles
  never enter the coverage merge.
  The source-matched audit retires all old KSP.cpp, KSP.h and Object.h records
  and merges fresh normal Release/Debug counters against production snapshot
  run66. The combined 143-report result is 35,296 of 37,555 production lines
  (\(93.9848\%\)), with no nonproduction or missing-source records. Relative
  to run64, the batch adds 244 measured lines and 217 covered lines, leaving
  382 additional hits needed for \(95\%\) at this denominator. The expanded
  denominator is retained rather than hiding supported template variants.
- MFEM native coefficient loading: independent real/complex scalar, vector and
  rectangular-matrix P0/P1 inputs exercise both component orderings and check
  printer coefficients directly against exact values. The three complex P1
  variants fail before the fix because the first coefficient passes through
  `std::stod`; malformed first/later and truncated inputs also demonstrate
  standard exceptions or silent acceptance before the fix. The P1 reader now
  extracts its native scalar type and validates every coefficient before
  modifying the destination, raising native member-function diagnostics.
  All fourteen cases pass individually in Release, normal Debug and ASan/LSan.
  Eighteen existing P1 and nine P0 IO cases pass individually in fresh Release
  builds, including all supported cell geometries. The change does not alter
  MFEM headers, node mappings, component ordering or matrix shape conventions.
  Source snapshot run73 verifies all 570 production files. The coverage merge
  retires the complete changed P1 MFEM loader class and its function records;
  retained MFEM records map only to byte-identical source before and after the
  changed class. All other production sources match the prior audited snapshot.
  The combined 147-report measurement is 35,372 of 37,574 production lines
  (\(94.1396\%\)); MFEM has 1,480 of 1,540 lines covered. The batch adds
  19 measured production lines and 76 hits, leaving 324 further hits needed
  for \(95\%\) at the current denominator. Tests, benchmarks, sanitizer
  counters and pre-fix counters remain excluded.
- OpenMP lifecycle and spatial boundaries: exact two-cell mass plus constant-kernel
  nonlocal matrices test real/complex dense, sparse and triplet assemblers across
  copying, moving, polymorphic cloning and repeated assembly. Vector assemblers
  obey independent nodal-integral values. Explicit two-thread settings survive
  the lifecycle transitions. The ten-case lifecycle executable passes one case
  at a time in Release, normal Debug and ASan/LSan. A private mutation doubles
  all tested OpenMP matrix/vector contributions; both new scalar configurations
  reject it with 24 numerical failures each, covering all four operator families,
  three lifecycle instances and both repeated assemblies.
  Nine additional float/real/complex cases check empty spatial vectors, streaming,
  empty trace/determinant/inverse identities, rectangular row/column extraction,
  and zero-input matrix-vector products across every supported boundary size.
  They pass separately in Release, normal Debug and ASan/LSan. Existing exact
  PETSc mixed-target assembly oracles also pass individually in normal Debug
  with both real and complex PETSc, for Sequential and OpenMP.
  Production sources are unchanged and byte-match all 570 files in snapshot
  run73. Source-path normalization precedes the accepted run79 audit; the
  unnormalized run78 diagnostic is not a certified measurement. The combined
  159-report result is 36,075 of 38,255 production lines (\(94.3014\%\)).
  Compared with run74, normal Debug instrumentation exposes 681 more production
  lines and adds 703 hits. OpenMP has 631 of 658 lines covered, PETSc OpenMP
  698 of 775, SpatialMatrix 1,167 of 1,221, and SpatialVector 327 of 355.
  The expanded denominator is retained; 268 more hits are needed for \(95\%\)
  at this denominator. Tests, benchmarks, sanitizer profiles and mutation
  profiles remain excluded. The final lifecycle file comment was updated
  after compilation; its executable test code is identical to the archived
  validation sources.
- Complex affine identification: independent matrix and load oracles now check
  \(T^{*}AT\) and \(T^{*}(b-Ad)\), including native dense/sparse,
  Sequential/OpenMP and real/complex PETSc single/mixed systems. Previous
  transpose-only row expansion fails the complex projection and physical-energy
  checks, although the exact constant solution can still be recovered. Row
  factors now use conjugation; trial columns and reconstruction retain their
  original coefficients. Local PETSc and MPI partial identification assembly
  rejects before resetting system outputs. MPI rejection is collective, so an
  empty shard does not continue into PETSc collectives after another rank throws.
  Pre-fix real mixed tests demonstrate cleared outputs; complex tests also
  reject incorrect full matrix/load values. The three-rank empty-shard case
  times out before the collective rejection fix and passes afterward.
  Eight native cases, sixteen local PETSc cases, twenty distributed groups
  (real/complex, one/two/three/four/eight ranks), twenty-eight Release Poisson
  cases and ten native lifecycle cases pass separately. Final normal Debug
  validation passes twenty local cases and eight distributed groups. Native
  affine tests pass ASan/LSan; six complex PETSc cases pass ASan with external
  leak detection disabled under the PETSc test policy. Sanitizer and pre-fix
  counters remain isolated from the coverage merge.
  Snapshot run84 contains all 570 production files. The audit retires old
  hits and function hits throughout every changed assembler specialization,
  retaining byte-identical mapped lines only as zero-count denominator
  candidates. Corrected Release counters and normal Debug counters produce
  179 reports: 36,632 of 39,138 production lines (\(93.5970\%\)). The
  expanded denominator and retired hits explain the decrease from run79;
  tests, benchmarks, deleted sources and private diagnostic runs are excluded.
  Further matching-source regression refreshes are pending; this is not a
  \(95\%\) completion claim.
- CI distributed timeout: Tests run 38518314738, job 115702442764
  (`Unit (gcc-14, Debug, OFF, ON, ON)`) times out the three-rank whole
  `RodinPETScMPIFormTest` suite after 120 seconds. The first named-form case
  passes in about 111 seconds, leaving insufficient aggregate time for the
  remaining thirteen cases; no numerical assertion fails. That case now has
  its own invocation for each existing rank count, while the complementary
  filter retains all remaining and future cases. Both registrations keep the
  original 120-second timeout, mesh, tolerances, labels and MPI leak policy.
  All fourteen cases pass individually in a fresh normal Debug build on three
  ranks; CTest JSON confirms the six complementary registrations. This
  confirmed CI failure satisfies the user's condition for resuming pushes
  after batch validation. At the last check, 270 of 274 CI jobs passed, one
  was cancelled, the timeout job failed, and two remained in progress.
  Fresh normal Debug builds of the existing native Sequential/OpenMP suites
  also pass all 139 cases individually. Their matching-source counters produce
  the 183-report run89 audit: 36,730 of 39,193 production lines
  (\(93.7157\%\)), with no nonproduction or missing-source records.
  This adds 55 measured lines and 98 hits relative to run87; 504 further hits
  are required for \(95\%\) at this denominator.
- Matching-source nonlocal refresh: all sixteen native nonlocal cases and both
  complex PETSc restriction cases pass individually in normal Debug. The
  187-report run91 audit measures 36,928 of 39,237 production lines
  (\(94.1152\%\)), with no nonproduction or missing-source records. This
  restores fresh observations for changed assembler specializations instead
  of retaining their retired historic hits. A further 348 hits are needed
  for \(95\%\) at this denominator. The configured OpenMP default means
  mixed `Problem` fixtures do not certify explicit Sequential mixed assembly;
  direct backend oracles are being added rather than crediting those variants
  from the default backend's line hits.
  Four added explicit Sequential mixed cases cover real/complex dense/sparse
  storage against separated moment products, independent trial/test attribute
  restrictions and unmatched attributes. All twenty native nonlocal cases pass
  individually in normal Debug, Release and ASan/LSan. A private mutation doubles
  mixed Sequential local/global contributions: all four new cases fail and
  the sixteen unaffected cases pass. Its profiles are isolated and excluded.
  The final formatting change only wraps the assembly-input constructor call;
  executable tokens match the archived Debug/sanitizer test sources. Production
  remains byte-identical to snapshot run84. The 191-report run93 audit measures
  36,977 of 39,237 production lines (\(94.2401\%\)), adding 49 hits without
  reducing the denominator; 299 additional hits are required for \(95\%\).
- Explicit Sequential preassembled targets: four real/complex dense/sparse cases
  extend the projection-and-lifting oracle to the Sequential backend regardless
  of the configured default. P1 and quadratic H1 vector fields on triangles and
  quadrilaterals check full assembly, repeated assembly, row-only targeted
  elimination and preservation of the opposite system output. All fourteen
  lifecycle cases pass individually in Release, normal Debug and ASan/LSan.
  Doubling only Sequential preassembled insertion fails all four new cases
  while the ten unaffected cases pass. The first private mutation compile was
  killed under memory pressure; its isolated retry completes with one heavy
  compiler. Sanitizer and mutation counters are excluded.
  Before the subsequent root/quadrature repairs, production still matched
  snapshot run84. The 195-report run95 audit measures 37,001 of 39,237 lines
  (\(94.3013\%\)), adding 24 hits without reducing the denominator. There
  are 275 remaining hits needed for \(95\%\) at this denominator.
- Bracketed root endpoints: float and real regressions expose a valid exact left
  endpoint being discarded by the first interior bracket update, returning no result for
  a linear function. Both new typed endpoint cases fail on the preserved
  pre-fix binary, while its fourteen other cases pass. Endpoint roots are now
  accepted before any interior iteration; the header also owns its required
  standard includes. Sixteen cases pass individually in Release, normal Debug
  and ASan/LSan, covering invalid initial guesses, zero/absent derivative
  safeguards and depleted work budgets. All 21 heart-model caller regressions
  pass individually in Release, including branch-flow derivatives and dynamic
  Jacobian finite differences. A doubled returned-root mutation fails all four
  new initial-guess/derivative cases; endpoint and zero-budget controls pass.
- Quadrature ownership and point dimension: six rule-family cases check copies,
  moves and polymorphic clones against closed-form moments, finite data and
  reference-domain membership. Clones outlive their source object and receive
  independent cache identities. The direct Gaussian and Lobatto point rules
  fail the pre-fix dimension oracle: both stored a one-component coordinate for
  a zero-dimensional reference element. They now store an empty coordinate
  with unit weight. All 94 quadrature/Lobatto normal Debug and Release cases
  pass individually. All twenty exactness/ownership cases pass ASan/LSan using
  seven consistently instrumented QF implementation objects.
  The initial mixed sanitizer/library binary had incompatible Eigen allocation
  conventions: normal allocation used plain malloc while the sanitizer-built
  destructor expected handmade alignment metadata. A consistent private library
  resolves that diagnostic. The counter-instrumented Debug sanitizer sweep then
  times out at 600 seconds. An optimized, assertion-enabled sanitizer build
  without private coverage counters passes the unchanged sweep in 0.745 seconds;
  matching normal Debug/Release timings are 12.573 and 1.440 seconds. Numerical
  fixtures, tolerances and timeout are unchanged. Doubling measured quadrature
  weights fails all six new ownership oracles. All private diagnostic, sanitizer
  and mutation profiles remain excluded from coverage.
- Euclidean CI repair: Style run 38601562645 reports three naming violations
  in the operation-detection traits and nine undocumented aliases. Detection
  structs now use `HasIntersectMethod`, `HasDistanceMethod` and
  `HasConnectMethod`; positive/negative static assertions check their public
  signatures. The five planar aliases, three fixed-coordinate vector aliases
  and line/sphere intersection result alias now have documentation. Local
  clang-tidy 19 naming checks pass for the affected header. Targeted Doxygen
  1.9.1 XML confirms descriptions for all nine reported aliases; the pinned
  1.14.0 remote ratchet remains required. All 148 existing Euclidean Debug and
  Release cases pass individually. Changed-line formatting and whitespace
  checks pass. Remote Format and StyleLint passed before this repair.
  Final snapshot run102 verifies all 570 production files. The audit retires
  all prior NewtonRaphson and GaussLobatto hits/function hits and the changed
  GaussLegendre point builder. Other retained executable lines are byte-identical;
  Euclidean comments map without changing code, and renamed trait declarations
  carry no previous executable records. The 215-report run104 measurement is
  37,137 of 39,338 production lines (\(94.4049\%\)), with no nonproduction
  or missing-source records. Compared with run93, it adds 101 measured lines
  and 160 hits; 235 further hits are required for \(95\%\). Tests, benchmarks,
  sanitizer, pre-fix, failed-run and mutation profiles remain excluded.
- PETSc target refresh: all four existing real/complex single-field target cases
  pass separately in fresh normal Debug builds against snapshot run102. The
  219-report run107 audit measures 37,234 of 39,386 production lines
  (\(94.5361\%\)), restoring 97 matching-source hits while exposing 48 more
  measured lines. The denominator is retained; 183 further hits are needed
  for \(95\%\) at this point.
- Public PETSc problem lifecycle: new single/mixed interval mass-system cases
  check copies, moves, assignments, polymorphic clones, destruction of one
  owner, const access, partial target updates and lazy full reassembly before
  solving an assigned body. Independent closed-form mass matrices and loads
  verify nonzero value constraints and the solution written to borrowed fields.
  PETSc LinearSystem copies explicitly share reference-counted matrix/vector
  handles: updates are therefore checked through both owners. The first private
  diagnostic incorrectly assumed deep system copies; it is excluded and no
  production behavior is changed to satisfy that assumption.
  All eight normal real/complex Debug/Release runs pass individually. Both
  complex ASan cases pass with external leak detection disabled under the PETSc
  policy. An isolated mutation retains the stale assembled flag on body
  assignment; both cases reject the resulting wrong solution. CTest JSON
  verifies the real/complex registrations, labels and 120-second timeout.
  The 227-report run109 audit measures 37,300 of 39,425 production lines
  (\(94.6100\%\)), with no nonproduction or missing-source records. It adds
  39 measured lines and 66 hits relative to run107, leaving 154 further hits
  needed for \(95\%\) at that denominator. Production sources remain
  byte-identical to run102 for this batch; failed-oracle, sanitizer and mutation
  profiles are excluded. The subsequent Flow repair requires separate fresh
  observations before its coverage can be certified.

- Flow boundary and scalar regressions: a continuing boundary callback that
  changes neither time, cell, coordinates nor correction now stops at a
  zero-time boundary hit, using the existing finite-time progress criterion.
  The unchanged implementation fails the boundary oracle in Debug and Release;
  its interior-start counterpart already passes. Complex scalar composition
  previously fails to compile because the exit value fixes the deduced return
  type to Real. The zero value now has the type of operand plus real correction,
  retaining Boolean-to-Real promotion. Real, complex, Boolean, copy and correction
  checks pass in all fourteen Debug, Release and native ASan/LSan cases.
  Private guard removal fails the boundary case; doubled field evaluation fails
  both scalar-composition cases. The three unaffected mutation counterparts pass.
  Seven further oblique-characteristic cases cross internal faces and corners on
  all supported positive-dimensional cell geometries. Exact straight trajectories
  and affine field values agree in both time directions in Debug, Release and
  native ASan/LSan, one case at a time. Halving the returned reference
  coordinates fails every geometry's physical-position oracle in isolation.
  All old Flow hits and function hits are retired before collecting matching
  observations against snapshot run114. The 235-report run118 audit measures
  37,497 of 39,637 production lines (\(94.6010\%\)); no nonproduction or
  missing-source records remain. Newly exposed executable lines are retained,
  leaving 159 further hits needed for \(95\%\). Sanitizer, pre-fix, diagnostic
  and mutation profiles remain excluded. Remote Style, Documentation and
  Installation checks pass for head 518e2b172; Build and Tests remain queued.

- Native solver ownership: all forty real/complex Debug/Release case-runs pass
  for CG (dense/sparse), BiCGSTAB (sparse), GMRES and DGMRES (dense/sparse),
  LDLT, HouseholderQR and PartialPivLU (dense), and SimplicialLDLT,
  SimplicialLLT and SparseQR (sparse). Explicit solutions and relative residuals
  check changed operators/loads, cold and warm guesses, copies, moves,
  polymorphic clones, associated problem references and original destruction.
  Fresh compilation exposes invalid Eigen copy/move initializers in both GMRES
  move constructors and DGMRES's dense copy and both move constructors. The
  repair follows the existing sparse-copy fix in e15528c52: retain the problem
  and construct fresh workspace. Tests reapply iteration options; they do not
  certify copied option values or factors. Pre-fix compile diagnostics are
  excluded. BiCGSTAB has no dense specialization; its failed dense instantiation
  records an unavailable API rather than a passed or substituted algorithm.
  Fourteen Release manufactured cases check copied GMRES/DGMRES solvers on
  affine harmonic Dirichlet fields across all seven positive-dimensional cell
  geometries, including every cell vertex and centroid after repeated solves.
  All sixteen affine chart cases pass in Debug, Release and ASan/LSan for
  tabulated P1 and pointwise-only elements across eight cell geometries and
  admissible embeddings. The first private fixture assigned to PointCloud's
  value-returning accessor and left its controls zero; the corrected fixture
  writes matrix entries. That diagnostic is excluded without changing geometry.
  Attribute initialization, sparse chart clones, connectivity-key lookup,
  rectangular spatial identities, binary-exact scaling and contextual warnings
  have independent value oracles. The warning's plain context and styled output
  are checked separately; a private diagnostic incorrectly expected styled
  output to contain the entire unstyled string contiguously and is excluded.
  All 178 scoped normal case-runs pass. Twenty-five isolated numerical mutations
  reject doubled algebraic solutions, mapped coordinates and manufactured fields.
  All four changed-solver, sixteen chart and five metadata cases pass ASan/LSan.
  All six Alert and fifteen spatial cases also pass ASan/LSan.
  Named iteration budgets retain the tested values; four final helper smoke
  case-runs pass. Generated CTest scripts verify fifty unique registrations,
  with 36 unit cases, fourteen manufactured cases and 120-second timeouts.
  Final snapshot run132 verifies all 570 production files. Old GMRES/DGMRES and
  SubMeshBuilder hits/function hits are retired. Fresh normal observations
  include matching out-of-line library notes and isolated successful runtime
  profiles, retaining their complete compiler-note baselines. The 303-report
  run137 audit measures 37,971 of 40,059 production lines (\(94.7877\%\)),
  with no nonproduction or missing-source records. Relative to run118, it adds
  422 measured lines and 474 hits; 86 further hits are needed for \(95\%\).
  Failed fixtures, discovery, sanitizer, pre-fix and mutation counters remain
  excluded. Test/benchmark/helper sources never contribute measured lines.
- CI head 518e2b172: Build, Style, Documentation and Installation pass.
  Tests run 38609137345 is cancelled overall, but all sixteen Unit jobs fail.
  Representative Debug jobs 116008410521 and 116008410503 fail four cases:
  GaussLegendreTest.PointGeometry plus the segment skin and point-cloud
  keep/trim cases. Release job 116008410416 fails only PointGeometry.
  The quadrature fixture expected one reference coordinate and read it; the
  point contract has zero coordinates. The fixture now checks empty coordinates
  and unit weight. Debug skin cases assert in Connectivity::polytope because
  SubMesh::Builder inserts a point through the positive-dimensional API.
  Vertex emission now belongs to finalize(), which already creates the vertex
  map and copies labels. Initialization sets the builder dimension before its
  first read. All forty relevant normal Debug/Release checks and all five skin
  ASan/LSan checks pass separately. Commit fe96b663a pushes these validated
  repairs after the observed CI failures; fresh CI is required for that head.

- Cache and optional ParU contracts: reference-point extensions exercise owned
  samples, supplied Jacobians and rebinding across eight geometries and every
  supported embedding. Quadrature cache growth preserves slots and dimension
  buckets; an empty ambient point retains unit metric density. ParU rejects
  nonsquare operators and invalid loads without modifying output or existing
  factors, then recovers from a singular factorization after a same-pattern
  value change. All six cases pass separately in Debug, Release and ASan/LSan
  with leak detection enabled. Private mutations doubling mapped coordinates
  or computed solutions fail the corresponding numerical oracles; unaffected
  empty-coordinate and pre-solve validation cases remain green.
- Fresh coverage run150 retains 315 normal initial/runtime reports and the
  complete production-source baseline. It measures 38,097 of 40,133 production
  lines (\(94.9269\%\)); test, benchmark, helper, failed, sanitizer and mutation
  profiles are excluded. Current production hashes match source snapshot132.
  Run155 adds the passing format-keyword and NaN component checks, measuring
  38,158 of 40,184 lines (\(94.9582\%\)) across 317 reports. New observations can
  expose additional measured template lines as well as cover old misses.

- Formulation clones preserve raw and base expressions, partially assembled
  operators and partially assembled loads after the originals are destroyed or
  zeroed. Real/complex and dense/sparse systems agree with an independently
  integrated interval mass matrix and load. Builder pointer coordinates remain
  owned after caller arrays expire; repeated indexed lvalue insertion deduplicates
  connectivity. The six builder cases pass in Debug and Release.
- PETSc regional interpolation checks Cells, Faces, Boundary and Interface on a
  nine-vertex mesh against an affine field. Selected vertices equal their exact
  values and unselected vertices retain a distinct sentinel; a false predicate
  preserves the entire field. All nineteen affected real Debug cases pass.
- Fresh run159 includes 325 initial/runtime reports, current production hashes
  and the retained complete baseline: 38,209 of 40,211 production lines are hit
  (\(95.0213\%\)). No nonproduction or missing-source records occur. Sanitizer,
  mutation, failed and discovery counters remain excluded. This is a local
  measured threshold; pending mode checks and remote CI remain separate gates.

- Final local validation: the two formulation-clone cases pass in Release and
  ASan/LSan with leak detection enabled. Omitting copied local bilinear terms
  in a private production-header mutation fails both real and complex numerical
  tests. All six builder cases pass ASan/LSan; the 38 format/math cases pass in
  Release. All nineteen real PETSc field cases pass separately in Debug,
  Release and ASan (external leak detection disabled under the PETSc policy).
  The eight new cache/ParU/formulation CTest registrations preserve their unit
  labels and 120-second per-case timeout.
- Run166 merges 331 normal reports and retains the run159 measurement:
  38,209 / 40,211 production lines (\(95.0213\%\)). All 570 production hashes
  match the worktree. Their sorted compact JSON SHA-256 is
  `f50e45a02451890ca46eb69e3d999d55ecc64b37b28a54387f29dee1f36ba271`.
  This audit combines normal Debug and Release observations; the main remote
  coverage job uses Release. A local threshold does not establish the current
  remote Codecov percentage. The absent-source inventory contains 172 audit
  candidates, including declarations, aliases, constant data, optimized
  functions and the unsupported variants discussed above; none is silently
  removed from the production-source scope.

- Final complex PETSc validation uses freshly rebuilt current Geometry, IO,
  QF and PETSc libraries. All nineteen field cases pass separately in Debug
  and Release, including preservation and non-real analytic values on every
  regional projection path. Run171 adds both normal producers to the merge:
  335 reports retain 38,209 / 40,211 production lines (\(95.0213\%\)).
  Function observations expand to 91,852, with 67,833 hit; line coverage does
  not conceal the separate zero-hit-function inventory. All 570 production
  hashes also match tracked committed files, with no production worktree diff.
- The nineteen complex field cases also pass individually under ASan, with
  external leak detection disabled under the PETSc policy. These profiles are
  excluded from run171. All local validation for the final coverage batch is
  complete; fresh remote CI/Codecov remains required after pushing it.

- CI head fe96b663a, Slow job 116290790575: GitHub's failure annotation states
  that the hosted runner lost communication with the server. Build remains the
  last in-progress step; no slow test started. Logs stop while the instrumented
  DistinctSpaceQuadrature and WeightedComplexQuadrature units compile. This is
  not evidence of a numerical failure or a 60-minute workflow timeout.
  Separate current-source ARM64 GCC14 Release/coverage recompilations succeed:
  compiler peak RSS is 3,095,936 KiB and 3,049,824 KiB, respectively, with wall
  times 127.73 s and 159.88 s. All 21 distinct-space and seven weighted-complex
  cases pass individually. These local measurements do not establish the cause
  of an x86 hosted runner's communication loss. Existing CI build parallelism
  remains unchanged; the new head supplies a fresh CI attempt.
- XDMF topology vocabulary: an external protocol oracle checks linear mixed
  identifiers, quadratic mixed identifiers and quadratic topology names for
  all eight cell geometries, independently of the writer's helper outputs.
  Values are taken from the published
  [XDMF topology table](https://www.xdmf.org/index.php/XDMF_Model_and_Format#Topology).
  All four format-contract cases pass separately in Debug and Release. A private
  mutation mapping Segment to Polygon's mixed identifier fails the new test,
  while the three unrelated format cases remain green. No production code or
  measured-source exclusion changes in this batch.
- Fresh run179 retains 343 normal initial/runtime reports and all production
  hashes: 38,214 / 40,211 production lines (\(95.0337\%\)). The five additional
  line hits come from the format checks; the measured denominator is unchanged.
  Mutation, sanitizer, failed/discovery and compiler-resource-only profiles are
  excluded. Successful normal runtime captures from the two measured quadrature
  executables are included. This follow-up batch is local while CI validates
  pushed head 16e6a576b, whose documented audit remains run171.

- Spatial expression assignment regression: assigning an Eigen matrix or array
  expression that references the destination cleared its coordinates before
  evaluating the right-hand side. Float, Real and Complex regressions fail on
  the prior implementation, covering identity, scaling/translation, reversed
  coordinates and shrinking views at every active size from zero to three.
  Both assignment overloads now construct an independent fixed-capacity spatial
  vector before copying it into the destination. The original zeroed inactive
  storage contract is preserved, including after shortening and regrowing.
  No dynamic allocation is introduced. Empty initializer-list construction and
  in-place addition also receive explicit boundary checks.
- All 21 spatial-boundary cases pass separately in Debug, Release and ASan/LSan
  with leak detection enabled. The 74 main vector/matrix cases pass in Debug;
  the six scalar-product cases and eight tensor cases pass in Debug, and the
  seven enabled tensor cases pass in Release. The unsupported-access tensor
  assertion test is intentionally Debug-only.
- Source snapshot182 freezes all 570 current production hashes. The 343 prepared
  reports retire old observations of both changed SpatialVector assignment
  overloads and retain only byte-identical mapped lines elsewhere. Fresh
  normal math profiles replace those observations. Run187 (353 reports) measures
  38,233 / 40,223 production lines (\(95.0526\%\)); the increased template
  denominator is retained. Pre-fix, sanitizer, failed and discovery observations
  remain excluded. Only root/Alert archives and current header-only math code
  are linked by these producers, so stale geometry/backend assignment objects
  cannot enter their fresh observations.
- Final run190 includes the 74 main vector/matrix cases in both Debug and
  Release and all six scalar-product cases in both modes. Together with the
  boundary/tensor checks, this batch has 217 successful normal case-runs and
  21 successful ASan/LSan case-runs. The 359-report audit measures 38,265 /
  40,266 production lines (\(95.0305\%\)), with no nonproduction or missing
  production records. Fresh vector/matrix compiler notes expose 43 more
  measured lines than run187, and that denominator is retained. The regression
  repair and its coverage evidence remain local while CI checks 16e6a576b.

- SpatialMatrix aliasing regression: both Eigen expression assignment overloads
  shared the vector defect, clearing the destination before reading aliased
  entries. Float, Real and Complex regressions fail on the prior matrix code
  across every shape with row/column extents zero through three. They cover
  self-assignment, scaling/translation, rectangular transpose, shrinking views
  and zeroed inactive entries after regrowth. Matrix expression assignment now
  writes into an independent fixed-capacity SpatialMatrix and only then replaces
  the destination. The existing unrolled one/two/three-entry copy kernels remain;
  array expressions forward to the same canonical matrix expression path.
  The historical zeroStorage contract is preserved rather than removed.
- All 24 spatial-boundary cases pass separately in Debug, Release and ASan/LSan
  with leak detection enabled. Source snapshot192 retains all 570 production
  hashes. The 359 prepared reports retire both changed matrix expression
  assignment methods; unchanged lines retain only byte-identical mapped hits.
  This continues the earlier vector-method retirement rather than reinstating
  stale observations. Fresh normal math profiles will replace the retired hits.
- Final matrix batch: all 80 main matrix/vector/scalar-product cases pass in
  both Debug and Release. Including boundary cases, this gives 208 successful
  normal case-runs and 24 ASan/LSan case-runs. Run196 retains 367 normal reports
  and measures 38,267 / 40,268 production lines (\(95.0308\%\)), with no
  nonproduction or missing-source records. All production hashes match
  snapshot192. Pre-fix and sanitizer counters are excluded. The matrix repair
  remains local alongside the earlier vector and XDMF follow-up batches while
  CI validates pushed head 16e6a576b.

- Root-solver numerical contracts: multiplying finite residuals to compare
  signs underflows for small function scales, loses the bracket and returns a
  coordinate near one for the analytic root three quarters. A zero function
  tolerance also misses an exact interior zero, while summing large finite
  bounds overflows midpoint construction. Float and Real tests reproduce all
  three defects (four failed cases and two finite-value assertion failures).
  Direct sign comparisons now preserve the bracket without residual products;
  exact interior zeros terminate independently of tolerance; all three midpoint
  sites use overflow-safe `std::midpoint`. The former endpoint, safeguard and
  iteration-budget contracts remain tested.
- Scaling checks use a linear function with independently known root, both
  slope signs, initial guesses on both sides and ordinary/small/large finite
  scales. Function tolerances scale with the function, and coordinate accuracy
  is checked separately. Large positive and negative bounds exercise invalid
  initial guesses, safeguarded iterations and interval-tolerance exits.
  A distinct narrow-bracket test verifies the returned midpoint, its error
  bound and early termination. A private wrong-midpoint mutation fails all four
  midpoint-sensitive cases; exact-zero cases remain green.
- All 24 root cases pass separately in Debug, Release and ASan/LSan with leak
  detection enabled. All 21 cardiac-model caller cases pass separately in
  Release. Snapshot199 retains all 570 production hashes; the 367 prepared
  reports retire all prior root-header hits and function hits before fresh
  normal root/caller counters replace them. Pre-fix, sanitizer and mutation
  observations are excluded from the merge.
- Final root batch: 69 normal case-runs and 24 ASan/LSan case-runs pass.
  Test work budgets allow all binary64 bisection levels plus safeguarded startup
  steps; coordinate tolerances include rounding margin for scaled residuals.
  Run203 includes 373 normal initial/runtime reports and measures 38,271 /
  40,271 production lines (\(95.0336\%\)). No nonproduction or missing-source
  records occur, and all 570 production hashes match snapshot199. The larger
  denominator remains visible. The follow-up repair remains local while CI
  validates pushed head 16e6a576b.

- XDMF file-reference and vector-output contracts: two new tests check empty
  filename placeholders with repeated path separators, normalized relative XML
  references, static scalar data, and static vector data at both nodes and cell
  centers. Independent HDF5 reads verify rank, extents and exact affine samples;
  the second write changes the field to prove first-snapshot retention. XML
  checks verify Vector metadata, centers, shapes and both snapshot references.
- All six XDMF cases pass separately in Debug, Release and ASan/LSan with leak
  detection enabled. A normal Debug compilation of the unchanged XDMF source
  also passes all six cases, giving 18 normal and six sanitizer case-runs.
  Native IO/QF dependencies were rebuilt from current source. The private
  component source is hash-identical to snapshot199 and its coverage source
  mapping is recorded. The first helper compilation used an incorrect header
  search path and failed; it contributes no observations. Quoted-header lookup
  fixes the helper without changing production code.
- A private Scalar-for-Vector metadata mutation fails the new vector test while
  the other five cases pass. Disabling filename separator normalization fails
  the new filename test while the attribute-registration control passes.
  Mutation and sanitizer profiles are excluded from the normal coverage union.
- Run207 retains 379 normal initial/runtime reports and measures 38,330 /
  40,338 production lines (\(95.0221\%\)), with no nonproduction or missing-source
  records. All 570 production hashes still match snapshot199. Compared with
  run203, XDMF.cpp exposes 65 additional measurable lines and gains 53 covered
  lines; XDMF.h exposes two additional lines and gains six covered lines.
  Thus 59 more covered lines accompany 67 more measured lines, and the expanded
  denominator is retained. The batch remains local while CI checks pushed head
  16e6a576b; the latest live job read reports 72 successes, 40 running and 162
  queued, with no failures. Remote CI and Release-only Codecov confirmation
  remain pending.

- CI baseline-convergence timeout on pushed head16e6a576b: job116321006178
  reports "The job has exceeded the maximum execution time of 45m0s".
  Setup/build occupied about 24 minutes. The serial convergence run reached
  397 of 427 cases, with 396 completed passes and no reported failed case,
  before cancellation during HP Stokes on Hexahedron. Completed case timings
  total 1,285.11 seconds. Increase the baseline job budget from 45 to 90 minutes
  to accommodate setup/build plus the complete numerical gate. All selected
  tests, refinement levels, assertions and ctest -j1 execution are unchanged.
  The revised workflow parses as YAML; completion still requires fresh CI.
  This observed CI timeout satisfies the user's deferred-push condition.

- Common file-interface template audit: the generic Loader and Printer file-path
  overloads had compiler notes but zero runtime hits. New public-interface
  tests use base references to concrete MFEM/MEDIT mesh and P1 field codecs.
  Mesh checks independently verify coordinates, dimensions, cell attributes,
  replacement of old output, and missing-file rejection without changing mesh
  counts. Field checks cover scalar, vector and rectangular-matrix ranges and
  verify each coefficient against independently assigned integer values.
- All six format-contract cases pass separately in Debug, Release and native
  ASan/LSan with leak detection enabled: 12 final normal case-runs and six final
  sanitizer case-runs. An append-instead-of-truncate mutation fails both new
  cases cleanly while the four existing cases pass. Fatal stale-output checks
  prevent the deliberately malformed output from reaching format parsers.
  Superseded compile attempts, mutation runs and sanitizer counters are excluded.
- Run211 retains 383 normal initial/runtime reports and measures 38,344 /
  40,338 production lines (\(95.0568\%\)). Four previously unhit Loader lines,
  four Printer lines and six MEDIT lines now execute. No additional measured
  lines, nonproduction records or missing-source records occur. All 570
  production hashes match snapshot199. Formatting, style lint and all nine
  coverage-auditor tests pass. Fresh remote CI/Codecov remain required after
  publishing the local batches and the separately committed timeout fix.

- Unit/angle template audit: Unit.h previously had only six measured lines,
  despite documented arithmetic that returns concrete unit types. Binary/unary
  arithmetic actually returned the base Unit and compound assignments returned
  base references, preventing ordinary Rad/Deg expression composition. Tests
  cover all six arithmetic operators and all four compound operators for a
  custom Real unit, Rad, Deg, Float storage and Integer storage, checking exact
  values, concrete result types and original-object references during chaining.
  The operators now preserve Derived values/references as documented.
- Degree conversion multiplied by pi before dividing by 180, overflowing for
  large finite inputs whose radian result is finite. Forming the conversion
  scale before multiplying preserves these results. The independent oracle
  divides before multiplying and uses standard-library pi in long double;
  its ordering also stays finite where long double has only double precision.
  Tests include cardinal, fractional, signed, near-minimum-normal and large
  finite degree values with rounding-scaled tolerances.
- The final source reproduces eleven failed cases with frozen pre-fix Unit/Deg
  headers, while the ten existing unit cases pass. All 21 final unit cases pass
  separately in Debug, Release and native ASan/LSan with leak detection enabled.
  A private half-angle mutation fails the conversion oracle while the other
  twenty cases pass. Both 75-case Euclidean runs pass, including a concrete
  angle-expression consumer checking independently known sphere points and
  tangents; both eleven-case Rad runs pass. This batch contributes 214 successful
  normal case-runs and 21 sanitizer case-runs. No PDE formulation changed.
- Snapshot213 retains all 570 production hashes. The 383 prepared reports retire
  all previous Unit/Deg header and function observations; other production
  hashes remain identical to snapshot199. Run218 includes 395 normal reports
  and measures 38,386 / 40,380 production lines (\(95.0619\%\)). Unit.h now
  measures 45 covered lines instead of six, and Deg.h adds three covered lines:
  42 additional measured lines and 42 additional hits. Nonproduction and
  missing-source records are empty. Pre-fix, mutation, sanitizer and superseded
  setup attempts are excluded. Formatting, style lint and all nine auditor
  tests pass. The repair remains local while CI validates f6d558fe5; the latest
  live job read has four successes, 30 running and 240 queued, with no failures.

- Unit comparison boundary audit: implementing <= as !(>) and >= as !(<)
  incorrectly made floating NaN values ordered. New exact comparisons against
  the stored scalar fail for the custom Real unit, Rad, Deg and Float storage
  on the prior code; Integer storage and the other 21 unit cases pass.
  Both operators now delegate directly to the scalar, preserving unordered
  floating semantics and ordinary integer ordering. The class contract states
  this behavior explicitly.
- Each of the five unit types checks all six comparison operators over every
  pair of representative values. Floating values include signed zero,
  subnormals, adjacent representable values around one, finite extremes,
  infinities and quiet NaN; integer values include both finite extrema. A
  private <=-to-< mutation fails all five new cases and the existing comparison
  case while twenty unrelated cases pass, demonstrating equality boundaries
  are checked even for the integer control.
- All 26 cases pass separately in Debug, Release and native ASan/LSan with leak
  detection enabled: 52 normal case-runs and 26 sanitizer case-runs. Snapshot221
  retains all 570 production hashes. The 395 prepared reports retire all prior
  Unit header/function hits, retaining earlier Unit/Deg retirements in history.
  Run223 contains 399 normal reports and measures 38,386 / 40,380 production
  lines (\(95.0619\%\)); the line counts are unchanged although the boundary
  semantics now have regression checks. Nonproduction and missing-source
  records are empty; pre-fix, mutation and sanitizer observations are excluded.
  Style lint and all nine auditor tests pass. The batch remains local while
  CI validates f6d558fe5: six successful test jobs, 34 running and 234 queued,
  with no failures at the latest live read. Style, Installation and
  Documentation have passed for that published head; Build and Tests continue.

- CI coverage-collection audit: normal optional-backend Release reports were
  saved only as artifacts, while assertion-enabled Debug implementations never
  contributed CI coverage. The Unit job now collects both unsanitized build
  modes and uploads only the filtered RodinCoverage.info file to Codecov.
  Artifact/report names include the build mode to keep all four threading/mode
  combinations distinct. Sanitizer profiles remain excluded. Unit, Slow and
  Coverage audits now include function inventories. The 95% floor and production
  source filters are unchanged; additional Debug notes can expand the measured
  denominator and must remain visible.
- Workflow validation confirms exactly four normal Unit coverage combinations,
  excludes sanitizer uploads, verifies unique artifact names, restricts uploads
  to the filtered LCOV file with discovery disabled, and checks function-audit
  collection in all three reporting jobs. YAML parsing, style lint and all nine
  coverage-auditor tests pass. This is configuration verification, with actual
  optional-backend CI/Codecov validation still pending after publication.
- An end-to-end native CMake Debug smoke build explicitly matches GCC14/gcov14,
  builds Unit, Rad and OptionalReference tests, and passes all forty registered
  cases one at a time. Build-time discovery counters use a separate excluded
  GCOV_PREFIX; runtime commands come from CTest JSON without executing binary
  discovery. The actual RodinCoverage target captures initial plus runtime
  notes and emits only production records; its function audit remains enabled
  in the merge. An earlier local smoke used the default GCC13 with gcov14;
  collection correctly rejected that mismatch, and none of it is merged.
- Run228 adds the valid CMake Debug report to the 399-report source221 union,
  retaining 400 reports and all 570 matching production hashes. It measures
  38,390 / 40,384 production lines (\(95.0624\%\)), with no nonproduction or
  missing-source records. OptionalReference contributes four newly measured,
  covered lines; existing identity/mutation/const-referent assertions exercise
  its previously absent runtime methods. The partial native smoke validates the
  collector, not the complete backend CI matrix. At the latest live read,
  f6d558fe5 has fifteen successful test jobs, 36 running and 223 queued, with
  no failures. This collection update stays local under the deferred-push rule.

- XDMF publication/recovery tests: a directory placed at the XML output path
  causes close() to report contextual MemberFunctionException diagnostics while
  preserving the open state and snapshot count. After removing the obstruction,
  flush() publishes both existing and later snapshots without closing the
  writer. Explicit close finalizes; a later close does no IO, even if the
  generated XML has been removed, and writes remain rejected after closure.
- A separate analytic tetrahedron test independently reads HDF5 coordinate,
  topology and nodal-value datasets. Exact shapes are checked before reading
  fixed buffers; coordinate sets, mixed-topology identifier/index sets, XYZ
  metadata and affine field samples are verified against independently known
  values. The tetrahedron identifier follows the external XDMF format table:
  https://www.xdmf.org/index.php/XDMF_Model_and_Format#Topology . The 3D branch
  already had runtime hits; this test strengthens its data/metadata oracle.
- All eight XDMF cases pass separately in Debug, Release and native ASan/LSan
  with leak detection enabled: 16 normal and eight sanitizer case-runs. Normal
  runner discovery profiles are isolated from actual case profiles. Native
  IO/QF dependencies were rebuilt from current headers. Closing before flush
  succeeds fails only the new recovery test; returning XY for dimension three
  fails only the new 3D test, with the other seven cases passing in each
  private mutation. Neither mutation nor sanitizer counters enter coverage.
- History confirms commit0046901a9 deliberately removed destructor close to
  avoid MPI collective lifecycle hazards. Correct the stale auto-close wording
  without changing executable behavior. Snapshot232 preserves all 570 production
  hashes; only XDMF.h lifecycle comments differ from snapshot221. Physical line
  counts stay unchanged, and every old/fresh executable coverage location maps
  to byte-identical code at the same line. All previous implementation-hit
  retirements remain in provenance history.
- Run233 retains 404 normal reports and measures 38,392 / 40,384 production
  lines (\(95.0674\%\)), with two additional publication-error diagnostic hits
  and no denominator change. Nonproduction and missing-source records are empty.
  Formatting, style lint and all nine auditor tests pass. The latest live CI
  read for published f6d558fe5 has 26 successful test jobs, 40 running and 208
  queued, with no failures; this batch remains local under the deferred-push
  condition.
- Additional template-audit boundary: OptionalReference currently exposes only
  non-const pointer/dereference member functions. Its existing tests exercise
  mutable wrappers borrowing mutable and const referents. Const-wrapper access
  is an API gap, rather than a supported specialization certified by its four
  covered lines; no new interface or source exclusion is introduced for it.

### Euclidean planar entry-point cleanup

- Retain the planar aliases in ForwardDecls.h: Point2D, LineSegment2D and
  Line2D bind the corresponding primitive to dimension two; Circle binds
  Sphere and Rectangle binds Box. Compile-time contracts verify each alias
  for both Float and Real in the existing fixed-storage tests.
- Remove the five include-only planar .hpp wrappers and their CMake source
  entries. Public .h entry points remain. Point.hpp, LineSegment.hpp and
  Line.hpp retain actual cross-primitive template definitions.
- Run235 passes all 75 Euclidean cases separately in Debug and all 75 in
  Release, with one case process at a time and discovery profiles isolated.
  All five public planar .h entry points compile in a separate consumer.
  CMake regeneration, formatting, style lint and all nine auditor tests pass.
- Snapshot236 contains 565 production files. Every retained file matches
  snapshot232 byte for byte; the five removed headers contain no measured
  executable lines, functions or branches. Run236 preserves the 404-report
  coverage union at 38,392 / 40,384 lines (\(95.0674\%\)); nonproduction and
  missing-source records remain empty. No implementation hits are retired.
  These repeated validation runs add no new coverage claims. Historical
  source snapshots and implementation-hit retirements remain preserved.
- The live CI read for published f6d558fe5 reports 29 successful test jobs,
  40 running and 205 queued, with no failures. Keep the cleanup local under
  the deferred-push condition.

### PETSc scalar-variant collection and broader compiler measurement

- Enable production coverage in the existing four PETScMatrixRanges jobs:
  real/complex storage and OpenMP on/off. Keep RelWithDebInfo, the same target
  and test selection, MPI ranks, two compiler jobs, sequential CTest execution
  and 45-minute budget. Explicitly disable ASan/UBSan; install lcov and select
  gcov-14. Capture through RodinCoverage, audit function inventories, save
  uniquely named scalar/threading artifacts and upload only filtered LCOV
  with automatic discovery disabled. The production source scope is unchanged.
- Strengthen the existing coefficient-norm test in complex builds: coefficients
  3+4i and -12i give analytic 2-, 1- and infinity-norms 13, 17 and 12. The real
  control retains its previous coefficients and expected values. Run238 builds
  fresh complex PETSc 3.19 dependencies and the test in RelWithDebInfo using
  GCC/gcov14. All nineteen complex cases and nineteen real controls pass as
  separate sequential processes. A private coefficient mutation stripping
  imaginary parts fails only the norm case; eighteen controls pass. All
  discovery and mutation profiles remain outside the normal collector.
- Run239 exercises the actual CMake RodinCoverage target on the clean complex
  build and verifies 69 hit complex GridFunction symbols. All 565 retained
  production hashes match snapshot236; no nonproduction or missing-source
  records occur. Workflow checks verify four unique unsanitized uploads,
  unchanged test selection/ranks/concurrency/budget, and function inventories.
- Run241 merges this verified combined initial/runtime report into the global
  union, now 405 normal reports. It measures 38,572 / 42,433 production lines
  (\(90.9009\%\)). Hits increase by 180, but the fresh RelWithDebInfo compiler
  notes expose 2,049 additional measured lines. Largest additions are
  Polytope.cpp (482), MPI Mesh.cpp (325), local Mesh.cpp (154), MEDIT.cpp (86),
  Shard.cpp (72) and MPI HDF5.h (66). The scalar smoke does not exercise many
  of these library paths. This broader union supersedes the narrower run236
  percentage; retain its denominator and run broader suites against the fresh
  libraries. No production implementation has changed or been excluded.
- Current goal remains incomplete: recover at least 1,740 additional hits if
  the denominator stays fixed, then verify the published CI/Codecov result.
  Latest published f6d558fe5 has 34 successful test jobs, forty running and
  two hundred queued, with no failures. Other workflow gates have passed;
  defer pushing until the Tests run fails or completes.

### Exhaustive polytope identity and fresh geometry coverage

- Add one independent multiset-identity case for every supported Polytope::Key
  size, including empty keys, unique vertices and repeated vertices. Enumerate
  all permutations, require equality in both directions and require invariant
  hashes. Replacing a vertex with a new value must break equality; changing
  multiplicity must also break equality. No inequality requirement is imposed
  on hashes of unequal keys, since collisions are permitted.
- Run242/243 builds the final key, polytope and distributed mesh tests against
  the fresh RelWithDebInfo libraries. Run244 passes 29 key and 55 polytope
  cases separately, plus nine distributed mesh cases at each of 1, 2, 3, 4
  and 8 ranks: 45 sequential MPI case groups. Each participating rank reports
  success. Every group finishes before the next starts; native cases run only
  after the distributed groups finish. Discovery profiles remain isolated.
- A private Polytope.cpp mutation that makes all five-vertex keys equal fails
  only the new exhaustive case; the other 28 key cases pass. Run247 keeps
  mutation profiles outside the normal build. The initial private compile
  attempt lacked the original local include directory and produced no normal
  coverage; the corrected private compile uses the original include path.
- Run245 merges the cumulative normal CMake coverage report, preserving all
  565 current production hashes and all previous source/retirement history.
  The union now contains 406 reports and measures 39,477 / 42,437 production
  lines (\(93.0250\%\)): 905 recovered hits and four newly measured lines.
  Nonproduction and missing-source records are empty. Keep the expanded
  denominator; the remaining target is at least 839 more hits if it stays
  fixed. Largest newly exposed gaps remain in Polytope/Point/Mesh methods,
  MEDIT/MFEM IO and PETSc SNES. More fresh broad-suite runs are required.
- Changed-line formatting, style lint and all nine auditor tests pass. The
  latest live CI read for f6d558fe5 reports 43 successful test jobs, forty
  running and 191 queued, without failures. Keep the batch local pending CI
  failure or completion; the goal is still incomplete.

### Broader optimized coverage and non-real manufactured SNES root

- Run248/251 builds fresh RelWithDebInfo geometry, native IO, quadrature and
  complex PETSc lifecycle suites. Run250/252 passes 324 additional cases as
  separate sequential processes: 141 geometry, 35 mesh-loader, twelve MEDIT,
  99 HDF5, fifteen Gauss-Legendre, twenty exactness/ownership and two SNES
  cases. Compiler-note capture starts only after all compilers finish.
  Successful counters are preserved before subsequent batches; discovery
  profiles stay outside the normal build.
- Run253 adds the cumulative filtered report to the global union, now 407
  normal reports. It measures 40,111 / 42,476 production lines
  (\(94.4321\%\)): 634 recovered hits and 39 new measured lines. All 565
  retained production files match snapshot236; nonproduction and missing-source
  records remain empty. At this denominator, 242 more hits are required.
- Add a non-real P1 manufactured cubic root to the SNES state suite:
  \(F(u)=u+u^3-(z+z^3)\), \(z=1+i/2\), with tangent \(1+3u^2\) and an initial
  iterate near that root. Check the full complex infinity-norm error and
  zero-iteration warm-start convergence. The solver uses the same bounded
  iteration/evaluation budgets as the existing real-root control. The error
  tolerance allows for assembly-scaled residual stopping, and is checked by
  deliberate imaginary-output deletion rather than assumed sufficient.
- Run254 passes all three complex cases against freshly rebuilt PETSc 3.19
  dependencies in RelWithDebInfo, plus both unchanged real Release controls.
  Run255 privately inserts VecRealPart after SNESSolve; only the new non-real
  case fails, while both existing controls pass. Every mutation profile is
  excluded. No production implementation changes are needed.
- Include the SNES state executable in the existing four PETScMatrixRanges
  jobs, so both real and complex controls execute remotely and the non-real
  root runs in complex builds. Preserve all previous tests/ranks, serial
  execution, compiler limits, unsanitized coverage filtering and job budget.
  Workflow validation confirms SNES is the only added selection.
- Run257 adds four source-matched initial/runtime reports for the successful
  real/complex cases, retaining the run253 line counts and percentage across
  411 reports. Observed function symbols increase to 92,303, with 68,312 hit;
  numerical template coverage is stronger despite unchanged line coverage.
  Formatting, style lint, workflow checks and all nine auditor tests pass.
- Latest f6d558fe5 CI read: 56 test jobs successful, forty running, 178 queued,
  no failures. Keep this batch local under the deferred-push condition. The
  next mesh/index/quadrature and distributed IO build remains active; the
  95% goal and published CI/Codecov verification remain incomplete.

### Fresh distributed and lifecycle validation retains the full denominator

- Run256/260 builds mesh/indexing/partitioning, Gauss-Lobatto, MPI geometry
  and HDF5, PETSc linear-system/KSP lifecycle, and parser/file-interface
  executables against the same fresh RelWithDebInfo libraries. All compiler
  and runtime phases preserve the two-compiler and one-case-group limits.
- Run258 passes 142 native cases. Run259 completes 325 sharder groups, eighty
  submesh groups and 230 HDF5 groups at 1, 2, 3, 4 and 8 ranks. The existing
  sharder rank-limit guards skip 137 complete groups, and the existing
  submesh empty-rank guard skips one complete group at one rank. Every other
  group passes on every participating rank. Do not count these 138 skipped
  groups as executed numerical tests. Run261 passes 26 further linear-system,
  KSP and file-interface/parser controls individually.
- Run262 merges the cumulative filtered normal report, retaining 412 reports,
  565 current production hashes and every previous source/retirement record.
  It measures 40,281 / 42,506 production lines (\(94.7654\%\)), recovering
  170 hits and exposing thirty further measured lines. No failed, mutation,
  sanitizer or discovery profiles enter the report. Successful counter
  snapshots are preserved before subsequent batches.
- Run263/264 builds and passes one hundred more cases individually: six
  metadata lifecycle, two iterator lifecycle, 58 connectivity, seven centroid,
  eleven parametric-transformation and sixteen chart-ownership cases.
  Capture begins only after all compilation and runtime work is terminal.
- Run265 retains 413 normal reports and measures 40,323 / 42,523 production
  lines (\(94.8263\%\)): another 42 hits and seventeen newly measured lines.
  Nonproduction and missing-source records remain empty. All 565 production
  files still match snapshot236. At this measured denominator, 74 additional
  hits are required; further supported instantiations may expand it again.
- Remaining-source inspection distinguishes legitimate gaps from invariant
  guards and unreachable continuations after Alert::Raise. P1Element.hpp's
  many default local-index branches reject invalid basis indices; do not
  fabricate such states for coverage. Supported MPI H1 face-interior indexing
  remains a useful next template-audit target. Retain all measured lines.
- Style lint and all nine coverage-auditor tests pass. The run265 CI job read
  has 73 successful test jobs, forty running and 161 queued, with no failures.
  A separate all-check-run and commit-status inspection finds no failed gate
  and no published Codecov result. Keep the commits local until CI failure or
  completion, and keep the goal active pending both 95% coverage and remote
  verification.

### Full complex H1 restriction matrix and CI coverage publication

- Run267 compiles the complete PETSc MPI H1 restriction target against the
  fresh complex PETSc 3.19 RelWithDebInfo libraries, using one compiler for
  this large template translation unit. Run269 executes its full numerical
  case at 1, 2, 3, 4 and 8 ranks, one distributed group at a time. All five
  groups pass without skips and within the existing 300-second gate: 162.12,
  133.54, 125.17, 128.77 and 204.02 seconds, respectively.
- The case retains degrees one through six, scalar/vector/matrix backend
  ranges, every supported positive-dimensional geometry, affine and curved
  parents, and full/boundary/sparse/nested submesh restriction. It checks
  every held functional and independent polynomial samples; the paired
  getGlobalIndex overload is also checked against the returned DOF array.
  Discovery profiles are isolated, and prior successful counters are saved
  before execution. No assertions, combinations or rank counts are removed.
- Run270 retains 414 source-matched normal reports. It measures
  40,375 / 42,546 production lines (\(94.8973\%\)): 52 additional hits and
  23 more measured lines. Observed function symbols increase from 92,317 to
  93,818, and hit symbols from 68,327 to 69,598. This is runtime template
  evidence, not include-only instantiation. All 565 retained production
  hashes match snapshot236; nonproduction and missing-source records are empty.
  The remaining line target is 44 additional hits at this denominator.
- The two existing PETScComplexConvergence submesh jobs now capture and
  publish production coverage rather than discard it. Select gcov-14,
  explicitly disable sanitizers, enable coverage only for suite=submesh,
  audit function inventories and save unique threading-specific artifacts.
  Upload only filtered RodinCoverage.info with automatic discovery disabled.
  Every existing build/run command, suite/rank matrix, thread setting, test
  timeout and job budget remains unchanged. Other complex-convergence suites
  keep coverage disabled. Local actual CMake capture validates the collector.
- Workflow validation confirms that exact scope; existing Unit and matrix-range
  collector checks, style lint and all nine coverage-auditor tests pass.
  Latest run270 CI read: 128 successful test jobs, forty running and 106 queued,
  no failures. Keep the changes local until CI failure or completion and keep
  the goal active pending both 95% and published CI/Codecov verification.

### Face-projection regressions, correction and fresh implementation coverage

- Strengthen projection oracles with exact fixed-point checks on every planar
  edge and a parameterized exterior-point sweep on all seven supported
  positive-dimensional geometries. Independently check parent halfspace
  inequalities and the requested face equality. The tolerance allows sixteen
  rounding units for short coordinate projections/halfspace products.
- Run272 observes five genuine failures before changing production code:
  triangle/quadrilateral projections reflect parent coordinates of points
  already on their face; tetrahedron/pyramid/wedge clipping can return points
  outside the requested simplex face. All 26 unrelated controls pass. These
  pre-fix profiles are excluded. History and Flow/Advection callers confirm
  that output coordinates are consumed in the parent reference cell.
- Preserve those parent coordinates and delegate faulty triangular clipping
  to the existing triangle-cell projector instead of duplicating its
  algorithm. Clarify the reference-coordinate contract. No input preconditions
  are bypassed, and no numerical test thresholds are loosened.
- Archive all old normal compiler notes/counters in run273 before rebuilding;
  reset the 74 working counter files to avoid stamp/implementation mixing.
  Snapshot274 retains 565 production files. Only Polytope.cpp and its .h
  contract change. Retire every old Project::face DA/branch/function record
  across all 414 prior reports (old range 1595-1910, new range 1595-1898),
  including nested lambdas. Retain only byte-identical CPP prefix/suffix
  locations; header documentation preserves every executable location and
  physical line count. All other production hashes match snapshot236.
- Run275 rebuilds current dependencies. Run276 passes 31 projector, fourteen
  Flow, seven Flow-crossing and nine Advection cases one process at a time.
  Run277 passes all 31 projector cases with assertions and again under native
  ASan/LSan, with instrumentation in the changed CPP and leak detection enabled.
  The private-copy Debug build is validation-only; its copied CPP source path
  is not silently omitted from a claimed production measurement.
- Run279 recompiles Debug CPP from its original production path and passes
  all 31 cases separately. Run280 asserts that production CPP appears in both
  its initial and runtime LCOV reports before merging; sanitizer, pre-fix and
  discovery profiles remain excluded. It measures 40,476 / 42,665 lines
  (\(94.8693\%\)) across 417 reports. New normal measurements add 101 hits
  and 119 measured lines relative to the prior source audit, including a
  legitimate complete assertion-enabled CPP baseline.
- Run281 links existing controls to that Debug CPP object and passes 129 more
  cases individually: 29 key, 55 polytope, 28 traits and seventeen extended
  polytope cases. Run282 retains 419 reports and measures 40,501 / 42,665
  production lines (\(94.9279\%\)), recovering 25 further hits without a
  denominator change. Nonproduction and missing-source records remain empty;
  all 565 current production hashes match snapshot274. The remaining target is
  31 additional hits if this denominator stays fixed.
- Changed-line formatting, style lint and all nine auditor tests pass. Latest
  run281 CI read for f6d558fe5: 203 successful jobs, forty running and 31 queued,
  with no failures. The goal remains active and this fix remains local until
  the deferred-push condition is met and published CI/Codecov is verified.

### Public mesh file and absent-face contracts

- Exercise public Mesh::save/load for MFEM, MEDIT and HDF5, replacing an
  existing two-dimensional mesh with a one-dimensional mesh. Compare physical
  coordinate and segment sets independently of vertex numbering, and verify
  the association between cell attributes and physical segments. Permit valid
  format sidecars and remove all scratch output after the checks.
- Reject unsupported format values with native member-function diagnostics
  while preserving geometry, topology, attributes and the absent filesystem
  destination. A closed triangular surface and a single quadrilateral also
  verify empty boundary/interface exceptions, independent edge valences,
  preserved topology and the complementary nonempty face selection.
- The first file fixture assumed normalized UniformGrid spacing. Actual grid
  spacing is one, and all formats preserved it correctly. Scale the source
  explicitly by one half; exclude every failed run285 profile. No production
  change was necessary for these file or face-selection contracts.
- Run286 passes eleven mesh cases individually in both normal Debug and
  Release. Run288 preserves their initial and runtime notes, including the
  production Mesh.cpp records, and measures 40,515 / 42,667 lines
  (\(94.9563\%\)) across 423 reports. The two extra measured lines remain
  in the denominator. Run289 freshly rebuilds and passes twelve SubMesh,
  five skin-lifecycle and three UniformGrid cases one process at a time.
  Run290 retains their cumulative report and reaches 40,530 hits.
- Run292 passes all twelve final mesh cases individually in each normal
  build mode. Run293 retains both initial/runtime report pairs and all prior
  valid reports: 40,534 / 42,667 production lines (\(95.0008203\%\)),
  428 reports, 93,850 observed and 69,649 hit function symbols. All 565
  production hashes match snapshot274. Nonproduction and missing-source
  records are empty; failed, sanitizer and discovery counters stay excluded.
- This reaches the local line target with a narrow margin. It does not prove
  every uninstantiated template, and published CI/Codecov confirmation remains
  required. Keep the goal active and defer the batch push until the current
  Tests run fails or completes. Changed-line formatting, house-style lint and
  all nine coverage-auditor tests pass.

### GLL absent-template audit in normal Debug

- Run294 recompiles the existing GLL mathematical suite in normal Debug and
  passes all 51 cases individually. These check analytic nodes, ordering,
  symmetry, intervals and the affine mapping. Compiler notes and runtime
  counters expose eight GLL.h lines, all executed. Symbol inspection retains
  scalar access/count methods and node arrays, with no emitted compute method;
  constant initialization is not counted as runtime execution of its algorithm.
- Run295 retains both initial and runtime reports together with run293's
  complete inventory. It measures 40,542 / 42,675 production lines
  (\(95.0017575\%\)) across 430 reports and 565 matching production hashes.
  Observed function symbols are 93,864, with 69,663 hit. The absent-source
  inventory drops from 163 to 162; declaration-only headers, constants and
  previously documented API gaps remain audit candidates rather than source
  exclusions. Nonproduction and missing-source records remain empty.
- Run296 privately changes the GLL01 affine factor from one half to one
  quarter. MappingConsistency_K15 fails for the incorrect finite nodes;
  KnownValues_K3 remains green. Neither the shadow header nor any mutation
  counters enter the production measurement. The real production tree is
  unchanged. The local target is reached; remote CI/Codecov is still pending.
- A completed current-head optional-backend CI artifact contains 31,614 hits
  across 33,725 measured lines, with no nonproduction or missing-source records.
  This partial report is not merged with newer local source snapshots. The
  regular CI coverage jobs depend on Unit/Manufactured/Slow completion, so
  their final aggregate cannot yet be inferred from that artifact.

### High-order quadrilateral face orientation audit

- Audit unexecuted valid reflection branches in H1.hpp rather than forcing
  invalid-index assertions. Add 24 parameterized cases spanning all eight
  square symmetries for hexahedron, pyramid and wedge faces, each at degrees
  three and four. Explicitly construct the face before its volume cell and
  verify that the requested ordering survives mesh construction.
- Map face and cell reference nodes independently through Geometry::Point.
  Every physical face node must match exactly one physical cell node and
  their global DOFs must agree. The coordinate matching tolerance is far
  below distinct node separation and allows the existing GLL node accuracy.
  No node match is inferred from the DOF map being checked.
- Run299 passes the preliminary 24 Debug cases; run300 passes all 24 final
  cases individually in both normal Debug and Release. Run301 retains both
  initial/runtime report pairs and all prior valid observations: 40,552 /
  42,675 production lines (\(95.0251904\%\)) across 434 reports. Ten valid
  orientation lines are recovered without expanding the denominator. Observed
  function symbols remain 93,864, with 69,723 hit. All 565 production hashes
  match snapshot274; nonproduction and missing-source records remain empty.
- Run302 privately ignores the chosen bottom-face transform in the hexahedron
  DOF assignment. The identity-orientation control passes and the reflected
  case fails on incorrect finite global DOFs. Mutation profiles and shadow
  headers never enter coverage, and production code remains unchanged.
- The target is registered with the unit label so the existing unit and
  coverage CI jobs include it. Formatting and house-style lint pass. The
  current remote run still has no failures; keep these batches local until
  the deferred-push condition is met and require published coverage validation.

### Retained CI aggregate inputs

- Save filtered RodinCoverage.info alongside coverage-audit.json in the
  existing regular and slow coverage artifacts, matching the optional/PETSc
  collectors. This preserves the line records needed to reproduce the final
  same-head union independently of the external coverage service.
- Parse the workflow and compare its complete structure before/after: only
  the two artifact path fields change. Existing matrices, commands, test
  selections, serial execution, dependencies and time budgets are identical.
  All five coverage-producing job families still upload only the filtered
  production LCOV file with automatic discovery disabled.
- Never combine raw reports from different heads without source provenance.
  The latest local measurement stays 40,552 / 42,675 lines
  (\(95.0251904\%\)); this artifact-retention change does not alter it.

### Dense preassembled mixed-block manufactured checks

- Add four manufactured cases: Sequential/OpenMP assembly on triangle and
  quadrilateral meshes, coupling an affine P1 field and a constant P0 field
  through four dense blocks with coefficients two, three, five and seven.
  Compare the complete matrix/vector against integral-only assembly, then
  independently check the manufactured coefficient vector's residual and the
  solved coefficients. Explicit tuples, offsets and UUID maps follow the
  mixed assembly input contract; each bilinear form uses its trial's solution
  type. This exercises rectangular off-diagonal block placement.
- The first runtime fixture (run314) shows matching assembly paths but an
  incorrect manufactured residual: its arbitrary callable supplies no known
  polynomial order, leaving the automatic load rule underintegrated. Declare
  integration order two for affine forcing times P1 and one for affine
  forcing times P0. Preserve the tight numerical tolerances; no production
  change is needed. Exclude the entire failed runtime producer.
- Run317's first corrected compile exhausts the container filesystem without
  a compiler diagnostic. Hash-verified recoverable archives run319/run321 move
  only old validation executables and object files to the host build mount.
  Preserve all compiler notes, counters, reports, source snapshots and logs.
  Retry the exact command only after matching the fixture bytes, verifying
  usable space and confirming that no runtime counters exist for this build.
- Run317 then passes all four corrected cases one process at a time in normal
  RelWithDebInfo with two OpenMP workers. Run318 retains initial and runtime
  compiler-note reports alongside all earlier valid evidence: 40,595 / 42,690
  production lines (\(95.0925275\%\)), 436 reports and 565 matching hashes.
  The broader denominator gains fifteen lines, while hits increase by 43.
  Function inventories record 94,064 observed and 69,877 hit symbols.
  Nonproduction and missing-source records remain empty.
- Run316 privately doubles dense preassembled coefficients in Sequential
  insertion. The numerical test fails on incorrect finite results, while the
  unaffected OpenMP control passes. Mutation and failed-fixture counters stay
  excluded. Run322 builds the actual CMake target and registers all four cases
  with the manufactured label. Formatting and house-style lint pass.
- The remaining current-head CI jobs are the two slow-test configurations;
  272 jobs have succeeded, with no observed failure. Published coverage and
  the deferred push remain pending.


### Source-body and exception inventory independent of instantiation

- Run341 syntactically scans all 565 production files plus the umbrella header,
  retaining inactive preprocessor branches and uninstantiated template bodies.
  The inventory contains 6,188 function bodies, 4,933 enclosed by template
  declarations. Comparing against the validated run318 line union identifies
  1,200 nonempty bodies with no attributed lines and 155 measured bodies with
  no hits. These are review candidates, not a new coverage denominator:
  constant evaluation and inlining can erase attribution.
- The scanner reports 94 parser recovery locations across 32 files. None spans
  more than two source lines after blanking the declaration-generating
  BOOST_SERIALIZATION_SPLIT_MEMBER spelling. Macro-generated bodies are not
  expanded, and parser recovery sites still require review. Empty/defaulted
  bodies, lambdas nested in ordinary bodies, and type-only declarations are
  not independently counted as semantic runtime functions.
- The exception scan confirms five direct standard runtime-error throws and
  one string-literal throw, plus 58 generic Alert exception sites in instance
  members that should carry member-function context. The MEDIT, MMG, Scotch
  and variational custom exception classes already derive from Alert.
  Static/free functions remain distinct from members with an instance.
  Complete-definition counters also include constructor initializers: 1,394
  definitions have initializer lists; 1,491 nonempty or initialized definitions
  have no attributed lines, and 152 measured definitions have no hits.
  These counts refine the body-only candidates above without treating
  initializer execution as an empty constructor body.
- A focused synthetic fixture verifies that uninstantiated and inactive bodies
  remain visible, that helper files are excluded, and that instance/static/free
  exception contexts and direct throws are distinguished. Keep parser tooling
  optional; it is not part of Rodin's build/runtime dependencies.


### CI timeout repair with retained test selections

- Tests run 38716344479 at f6d558fe5 ends cancelled after both Slow jobs exceed
  GitHub's six-hour limit. The failure annotations explicitly state the limit.
  Logs also show identification watchdog failures at one/three ranks, five
  PETSc H1 restriction watchdog failures, and accumulated manufactured-group
  timeouts. No completed assertion-failure groups were observed. This is the
  confirmed CI failure that releases the requested push deferral.
- Split Slow into four partitions per thread configuration. A read-only CTest
  guard verifies that the partitions are nonempty, disjoint and exhaustive.
  A minimal CMake fixture validates CTest index/label selection without running
  tests. Preserve every existing exclusion, the serial runner, and all degrees,
  ranks, numerical assertions and watchdog values.
- Register manufactured elasticity and Poisson groups by geometry. Generated
  CMake filters retain all original 450 Poisson and 64 elasticity rank/case
  combinations exactly once. Register PETSc H1 restriction by geometry at
  ranks one, two, three, four and eight. A token comparison verifies that all
  numerical body statements and tolerances are unchanged.
- Run342 passes all 35 PETSc geometry/rank cases, one MPI group at a time,
  using optimized GCC14 and complex PETSc 3.19.6. Every case stays below the
  original 300-second watchdog. Pyramid takes approximately 131/114/85/76/109
  seconds at the five rank counts. No geometry, order or range is dropped.
- The identification tetrahedron case passes separately at ranks one and three
  in approximately 154 and 64 seconds. Reduced compiler debug metadata avoids
  a local full-debug-metadata compiler OOM while preserving optimization.
  These timings do not establish the cause of the old CI timing differences.
  New CI confirmation remains necessary.
- Isolate build-time and PRE_TEST discovery profiles outside the capture tree
  in all five coverage-producing job families. Controlled threading and MPI
  binding settings match existing dedicated backend jobs. All other job
  families and every production-only Codecov upload remain unchanged.
- Regenerate actual CMake registrations, parse the complete workflow, run
  formatting and house-style checks, and verify all 565 production hashes
  still match run318: 40,595 / 42,690 lines (\(95.0925275\%\)). These changes
  do not alter production source or the local measured line inventory.


### Implicit exception candidates

- Extend the optional source inventory with calls to checked accessors,
  optional-value access, numeric string conversions, `std::get` and
  `std::any_cast`. Resolve callees and preconditions before treating a
  candidate as a violation: compile-time tuple access and sparse iterator
  values do not have the optional/variant failure modes suggested by their
  spelling. Exclude comments and avoid counting an outer member call merely
  because its receiver contains a nested `std::get` call.
- Run360 records 227 candidates against the local Alert-correction source.
  Public checked list/quadrature access and MMG numeric text conversion are
  concrete boundary follow-ups. Internal topology-map lookups, bound-element
  state and guarded optionals require their documented invariants to be
  reviewed; do not fabricate inconsistent internal states for coverage.
- Four tooling regression cases pass in separate Python processes: inactive
  and uninstantiated definitions, instance/static/free exception context,
  constructor-initializer attribution, implicit-call extraction, and explicit
  parser-recovery reporting. The tooling tests report an explicit skip if
  their optional parser dependencies are unavailable.


### Contextual Alert exception correction

- Replace five explicit standard-runtime-error throws and one string-literal
  throw with Alert variants. Upgrade all 58 identified generic Alert sites in
  instance members to MemberFunctionException, retaining specialized domain
  exceptions already within the Alert hierarchy. Static/free functions remain
  independent of an instance. Distributed downcast failures explicitly throw
  the contextual Alert object, preserving non-returning behavior.
- Convert SpatialTensor axis bounds failures from std::array::at to an explicit
  member-context Alert check. Invalid resize and axis requests preserve active
  extents and storage size. Correct Alert documentation to describe throwing
  after diagnostic output rather than unconditional program termination.
- Observe all three new PETSc/tensor regressions fail against the original
  implementations, then pass individually against the correction. The new
  MPI global-ID-cap regression also fails at both ranks against the original
  source-hash-matched Mesh.cpp, privately linked ahead of the current archive;
  it passes at ranks two and three against the correction. All red profiles
  and the private original implementation are excluded from coverage.
- Native MPI/MMG library builds and dependencies pass with one compiler job.
  All three PETSc linear-system cases and nine Release tensor cases pass in
  separate processes. All 330 MPI geometry case/rank invocations pass at ranks
  one, two, three, four and eight; 140 invocations explicitly skip unsupported
  rank counts. No runtime test groups overlap. Formatting and production style
  checks pass.
- Run354 hash-verifies and archives 129 old note/counter files before removing
  obsolete runtime counters. Run362 freezes all 565 production files plus the
  umbrella header. Normal post-correction profiles remain isolated in runs
  364-366. The historical run318 percentage belongs to the previously
  published production source: the changed-source batch still needs fresh
  aggregate CI/coverage confirmation and does not inherit that percentage.


### Degree/mapping isolation for the remaining complex restriction timeout

- New-head Tests run 38775961397 at 7c2c18013 reports one complex ON/submesh
  failure: one-rank Pyramid exceeds the 300-second watchdog. The two-rank
  Pyramid case passes at 298 seconds; the other 34 geometry/rank cases pass.
  The confirmed job failure releases the next deferred push.
- Split the existing matrix into degree, mapping and geometry parameters at
  every existing rank. All 420 CTest/GTest cases retain exactly the original
  5,040 scalar/vector/matrix and full/boundary/sparse/nested configurations.
  Numerical oracle statements, coefficient/physical samples, tolerances and
  the rank-zero local-read protocol are byte-identical ignoring whitespace
  and timer declarations. The watchdog remains 300 seconds.
- Optional diagnostics measure space construction, field setup, projection
  and verification, with local held-cell/DOF-slot/sample counts and cached
  global DOF sizes. Diagnostics add no collectives. Retain CTest logs in the
  complex collector even on failure; production coverage filtering is unchanged.
- Run370's initial test compile exposes structured-binding commas within a
  GoogleTest macro argument; parenthesize the name-generator lambda. No
  runtime counters exist for that failed compile. Run373 then exposes a
  diagnostic bug: a rank-zero count call is collective, causing the two-rank
  case to hang. Replace all diagnostic count calls with local iterator counts.
  Preserve the entire failed producer under its explicit excluded label and
  confirm that no live MPI children remain. The corrected two-rank control
  completes in about 0.23 seconds with all 24 range records.
- Run376 passes all 420 cases, one MPI group at a time, and verifies every
  range/submesh configuration actually executes. It records all 18,144
  rank-local range records. Observed maximum case durations at ranks one,
  two, three, four and eight are approximately 58, 45, 36, 36 and 36 seconds.
  These are validation timings, including permitted compilation during part
  of the run, rather than a production performance-improvement claim.
- Single-rank diagnostic data locate the large Pyramid cost in projection.
  Cell/DOF-slot counts alone do not explain its ratio to tensor-product
  hexahedra. Preserve phase/work data for further investigation; no numerical
  kernel optimization or tolerance change is introduced by this repair.
- Regenerate and verify all 420 exact CTest filters, timing environment and
  watchdog properties. The workflow changes only the complex collector's
  selection and diagnostic retention; all existing test families remain.
- Additional native tensor validation passes all ten assertion-enabled cases
  and all ten AddressSanitizer cases individually with leak detection enabled.
  Tensor/header object instrumentation links the existing optimized native
  support libraries. Sanitizer, failed-timing and original-implementation red
  profiles are excluded from coverage. Current-source aggregate coverage
  remains pending fresh CI confirmation.


### Point-element contracts from the source audit

- Add six real/complex P1 point-element cases covering constant reproduction,
  zero derivative extensions, vector component functionals at every supported
  extent, rectangular matrix component ordering, and zero polynomial order.
  Non-real coefficients prevent complex component loss from going unnoticed.
  The reference point and nodal functional both have zero coordinates.
- All 53 cases in the existing P1 executable pass separately in normal Release
  (run384) and assertion-enabled mode (run391). No production source changes
  are needed. Run385 privately changes the point first derivative from zero
  to one: all six new cases reject that wrong finite value, while the segment
  derivative control passes. Mutation headers and counters stay excluded.
- Run396 verifies the current P1 source hashes against snapshot362, captures
  the two normal test-object profiles, and confirms execution of the previously
  unhit point derivative assertion/zero return and vector point-order return.
  This is a file-specific check, not a changed-source aggregate percentage.


### Public checked-list bounds use contextual Alerts

- The implicit exception audit identifies FormLanguage::List::at as a public
  checked accessor explicitly forwarding std::out_of_range. Replace mutable
  and const bounds failures with MemberFunctionException and preserve the
  selected owned object for valid indices. Empty, one-past-end and maximum
  size_t requests are rejected without changing size, value or object identity.
- Both new regressions fail against the original behavior (run390), then all
  22 list cases pass individually in normal Release, assertion-enabled and
  AddressSanitizer validation (runs394/399/400). Leak detection is enabled for
  sanitizer validation; those profiles and all red profiles remain excluded.
- Declare Rodin::Alert as a public FormLanguage dependency so consumers receive
  the necessary termcolor include path and exception implementation linkage.
  The actual CMake target builds with one compiler job and discovery profiles
  isolated outside capture. Formatting, production style and the source
  exception-context scan pass. Aggregate current-source coverage remains
  pending; no historical percentage is carried onto the changed header.


### Quadrature bounds and instantiation-independent audit follow-up

- Replace PolytopeQuadrature::getPoint's standard checked-container exception
  with MemberFunctionException. The new invalid-index regression fails against
  the original implementation (run404), then all 12 quadrature-index cases
  pass individually in Release, assertion-enabled and AddressSanitizer modes
  (runs411/412/414). Leak detection is enabled. The sanitizer test object links
  optimized native support libraries; its profiles and the red profiles are
  excluded from coverage. Native Geometry dependencies are rebuilt serially,
  with their previous coverage notes retained in run407 before replacement.
- Run410 inventories all 566 production files and 6,188 function definitions
  without preprocessing or template instantiation. No explicit standard/string
  throws or generic Alert::Exception calls in instance members remain in the
  parsed inventory. Custom exceptions also derive from Alert; Scotch's variant
  already derives from MemberFunctionException. Specialized domain exceptions
  need a separate context review rather than being classified solely by name.
- The inventory still reports 224 potential implicit exception calls and 94
  parser recoveries. These are review candidates, not 224 proven violations.
  Remaining public checked-container forwarders include Scotch::Partitioner's
  getPartition and Geometry::SubMesh's getPolytopeMap. AttributeIndex dimension
  access needs contract review. MMG::MeshLoader has six std::stoul parsing sites
  that can propagate invalid_argument/out_of_range. Internal maps, guarded
  optionals, tuple get and typed any_cast must be reviewed by their actual
  preconditions; invalid internal states are not manufactured for coverage.
- Run409 uses four sequential syntax-only probes, producing no coverage notes
  or runtime counters. Including IDRSTABL fails because configured Eigen lacks
  Eigen::IDRSTABL. Including legacy TraceOperator fails on its obsolete
  FunctionBase template argument. Derived Solid Input/Output examples with
  explicit default constructors fail conversion to their std::function aliases:
  populate/extract do not supply a call operator. These findings concern usable
  source/API support, independently of whether coverage tools emit lines.
- Fresh CI run38790622528 artifacts are audited against an exact archive of
  published head ba57395, not the subsequently modified working tree. Run408's
  union of three PETSc matrix-range collectors measures 23,576 lines and hits
  14,730 (62.4788 percent), with no nonproduction records or missing sources.
  This is a partial collector union: regular, slow and complex-H1 results are
  still needed. It establishes neither final coverage nor a coverage regression.
  Current-source aggregate confirmation remains pending.


### Scotch checked-assignment exception contract

- Replace Partitioner::getPartition's std::out_of_range forwarding with
  MemberFunctionException. Document that an entity must have a stored assignment.
  Exercise access before the first partition, one-past-end and maximum indices,
  member context, and preservation of every valid assignment after rejection.
  The previous InvalidArguments test now expects the public Alert contract.
- The new case rejects the original implementation (run416); all eight Scotch
  cases pass separately in normal Release (run420) and assertion-enabled ASan
  with leak detection (run422). Sanitizer and red profiles remain excluded.
  Both the Scotch wrapper and test object are sanitizer-instrumented, linking
  the existing optimized support libraries and external Scotch library.
- The complex local build has Scotch disabled, so private commands compile the
  current wrapper and test explicitly and link Scotch. The inherited test flags
  promote GCC's existing fclose deleter attribute warning to an error; retain
  the initial diagnostic and demote that warning for private validation only.
  This does not establish a Scotch-enabled CMake configuration check.
- Run419 retains current production/configuration source hashes. Run423 checks
  the normal producer's eight passing cases and matching Scotch source hash,
  then confirms positive coverage counts on the bounds guard, Alert path and
  valid return. No overall percentage is inferred from this file-specific check.
  Run418 reduces the remaining implicit-call inventory to 223 candidates;
  public local/MPI SubMesh dimension access and MMG numeric parsing remain open.


### Local and distributed SubMesh checked map dimensions

- Replace both SubMesh getPolytopeMap checked-container forwarders with
  MemberFunctionException. Preserve stored map capacity: builders retain empty
  slots above the selected submesh dimension, so rejection is based on vector
  size rather than the submesh's topological dimension. The local builder uses
  parent topological dimension; the MPI builder uses ambient dimension.
- The new local regression fails against the original implementation (run427).
  All 13 local cases pass individually in normal Release (run436), and all 85
  MPI case/rank invocations pass individually at ranks one, two, three, four
  and eight (run438). The new MPI regression runs at every rank count without
  skipping and performs its maximum-index query only on rank zero, followed by
  an explicit test barrier. Guard rejection introduces no collective operation.
- All 13 local assertion-enabled AddressSanitizer cases pass separately with
  leak detection (run441). The test object is instrumented; optimized native
  support libraries remain linked. Red and sanitizer profiles stay excluded.
  Native dependencies are rebuilt serially against the final formatted headers;
  runs428/431 retain hash-verified copies of previous coverage notes first.
- Run439 records source/configuration hashes. Run442 checks the two normal
  producers and confirms positive guard, Alert and valid-return counts against
  matching local/MPI SubMesh header hashes. This is file-specific evidence;
  aggregate current-source coverage remains pending. Run435 inventories 221
  remaining implicit-call candidates, including MMG numeric parsing.
- The PR367 integration audit identifies one newer commit, 5c2a9697a, after
  previously merged add865281. Its independent Jacobian/distortion strengths
  preserve the old default products. Git's merge preflight succeeds without
  conflicts; integration and validation of that newer commit remain pending.


### Integrate the current PR367 hinge controls

- Merge current PR367 head 5c2a9697a after validating the independent SubMesh
  batch. The earlier cleanup head add865281 remains an ancestor; its removals
  are preserved. The new commit replaces a global hinge coefficient and row
  weights with independent Jacobian/distortion strengths. Defaults preserve
  the old products: 100 times 0.1 equals 10, and 100 times 1 equals 100.
- The merge is conflict-free, including the expanded solver coverage tests and
  contextual Alert changes. Scan source, tests, examples and documentation for
  retired weight fields and flags; the obsolete global CLI flag remains only
  in its deliberate rejection test. The integrated develop base remains
  cf9ff43f9, matching the current GitHub PR base.
- Build all affected tests serially with optimized flags. Nine hinge cases,
  48 solver cases, 13 parameter/CLI cases and six Adapt entry-point cases pass
  individually (runs448/449/450/452). These include product-strength energy and
  Newton-system equivalence, independent 2D/3D finite-difference tangents,
  higher-order assembly, invalid controls and physical mesh-application checks.
- Run447 retains current source/configuration hashes. Run453 validates all
  four normal producers, captures their matching test-object notes/counters,
  and inventories six changed SWIFT headers. Hinge.h's 40 measured lines and
  Parameters.h's 18 measured lines are all hit in that scoped union. Other
  headers retain investigation candidates; this partial file-specific union
  establishes no aggregate goal percentage. No stale-source or sanitizer
  profiles are included. Source inventory, production style and formatting pass.


### MMG metadata conversion errors use member-context Alerts

- Centralize the six metadata numeric conversions in MeshLoader::readIndex,
  which owns entry reading, conversion and section-specific diagnostics.
  Translate std::stoul's invalid_argument/out_of_range categories to a local
  MemberFunctionException and throw it explicitly on the non-void error path.
  Preserve existing numeric conversion and one-based assignment semantics.
- Twelve malformed-input checks (invalid text and overflow in each section)
  reject the original implementation, while the original valid-input control
  passes (run457). All 117 optimized MMG cases then pass individually against
  the corrected native library (run462). Valid controls check all six metadata
  sets, zero-based assignment, surrounding whitespace and retained topology.
- Both focused sanitizer cases pass separately with leak detection (run464).
  The two fixtures are copied verbatim and both they and MeshLoader.cpp are
  instrumented; optimized native support libraries remain linked. Fixture and
  loader hashes are retained. Red, sanitizer and failed compile profiles remain
  excluded. Run454 is compile-only: string loop copies are corrected to const
  references before the successful original-implementation fixture build.
- Native MMG and its dependencies rebuild with one compiler job. Run458 retains
  previous MMG coverage notes; earlier run431 retains the IO/QF notes before
  their refresh. Run461 records source/configuration hashes. Run468 checks the
  matching native loader and normal full-suite test-object profiles: all 49
  measured MeshLoader.cpp lines and all five MeshLoader.h lines are hit. Its
  eight required positive checks cover conversion, contextual throw and all
  six metadata setters. GCC attributes the wrapped tetrahedron call to the
  setter statement, leaving its continuation line unmeasured. No overall goal
  percentage is inferred from this file-specific result.
- Run459 inventories 216 implicit-call candidates, including the now-caught
  std::stoul call. No generic Alert::Exception calls in parsed instance members
  remain. Production style and changed-line formatting pass.


### Retain completed private builds and resolve workspace output I/O wait

- A full container root filesystem is relieved by archiving 15 completed
  private builds to /build/archived-validation-run456. Every file is SHA-256
  verified before removing its original copy, and original directory paths
  become symlinks to the preserved archives. Normal, sanitizer and failed
  producer classifications stay distinct; never capture the archive wholesale.
- A private MMG compile writing to the workspace volume is observed at 356
  seconds wall time with 31 seconds CPU and request_wait_answer wait state.
  Stop only its verified compiler and driver, observe terminal status, retain
  its partial output and mark it excluded with no runtime counters (run460).
  The unchanged-source local-output rerun completes in about 73 seconds with
  about 35 seconds child CPU (run466). These are build/storage observations,
  not numerical-kernel performance claims. The successful normal and focused
  sanitizer compilations use local outputs and retain their coverage notes.


### Instantiate Solid CRTP input/output forwarding contracts

- The instantiation-independent audit finds Input::populate and Output::extract
  absent from runtime attribution. Their documented direct conversion of CRTP
  objects to std::function fails the earlier run409 probes. Correct both examples
  to capture a concrete object in a lambda and call through its const CRTP base;
  provide an explicit derived default constructor to avoid protected-base
  aggregate initialization. Production tokens outside comments are unchanged.
- Add the native RodinSolidLocalCallbacksTest target. Three cases cover owned
  input payloads surviving callback copies, mutation of the original point,
  preservation of unrelated tags and kinematics, and output forwarding for
  NeoHookean and Yeoh cache types in two and three dimensions. Independent
  diagonal-deformation Jacobian/invariant values and original object addresses
  check the forwarded cache and point, rather than only callback invocation.
- All three cases pass separately in normal Release, assertion-enabled and
  AddressSanitizer validation (runs470/473/474). Leak detection is enabled;
  header/test-object instrumentation links optimized native support libraries.
  The first fixture compile's linked-list indexing error is corrected to
  iteration before any runtime validation. Discovery profiles stay excluded.
- Private input mutation overwrites the populated activation with zero: only
  the input case fails. Private output mutation forwards a cache copy with its
  invariant increased by one: both cache-type cases fail on wrong finite values
  and reference identity, while the input control passes (runs481/482). Preserve
  private headers, commands and logs; all mutation/sanitizer profiles are excluded.
- Run476 successfully compiles both corrected documentation adapter forms with
  syntax-only probes and verifies no coverage notes/counters are emitted.
  Run475 retains source hashes. Run479 checks matching normal Release/Debug
  notes and profiles: Input's forwarding line79 has four hits, Output's line89
  has sixteen. All three measured lines in each header are hit in this scoped
  union. This closes the two forwarding-body audit candidates without making
  an aggregate coverage claim. Source inventory, style and formatting pass.
