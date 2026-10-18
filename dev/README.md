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
- Function-valued `Potential::getValue` is an incomplete runtime path: it
  contains an unconditional assertion and dereferences a quadrature pointer
  that its constructor leaves null, with no public quadrature setter. Its
  missing counters must not be treated as a supported but merely untested
  evaluation. Shape-function nonlocal potentials are a separate specialization
  with existing numerical coverage.
- The P1 and H1 quadrature classes for `Mult<FunctionBase, Dot<Trial, Test>>`
  are dormant under the current public expression route. Multiplication on
  either side of a bilinear `Dot` normalizes to `Dot(coefficient * trial, test)`,
  which selects a different specialization. Their absent counters are source
  audit findings, not evidence that the existing outer-weighted public form
  lacks tests; the dispatch and numerical specialization suites cover that form.
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

`check_format.py` compares complete-file canonical formatting against lines
touched by the current change. It only reads Git state. `format_cpp.py` uses
clang-format followed by the plain-struct indentation rule: access labels are
inset two spaces, explicit-section members four, and plain-struct members two.

```sh
python3 dev/check_format.py
python3 dev/check_format.py --base origin/develop
python3 dev/check_format.py --binary clang-format
```

To fix reported hunks locally:

```sh
python3 dev/check_format.py --base origin/develop --fix [paths...]
```

This applies only complete canonical edits intersecting changed lines; unrelated
legacy layout is left untouched. `--fix` cannot be combined with `--head`, and
never writes Git's index. For new files or an intentional complete-file pass,
use `python3 dev/format_cpp.py -i <file>` instead.
Indentation-insensitive anchors prevent distant matches in repetitive code.
Before writing a changed-region result, its code/literal token sequence and
preprocessor directives are checked against the original; an unsafe partial
edit is rejected rather than applied.
Changed ranges are recomputed until no further edits are needed, including when
Git realigns neighbouring repeated braces after an indentation change.

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


### Rectangular SWIFT hinge integrators and bound-copy contracts

- Add a mixed P1 trial/P2 test-space regression in two and three dimensions,
  exercising previously uninstantiated integrator variants and HingeMetric's
  distinct-space test-Jacobian branch. Explicit common witness subdivisions
  define one frozen discrete energy for both spaces; no degree-exact continuum
  quadrature or convergence claim is made. Field coefficients are obtained
  through public interpolation, without assumptions about high-order DOF layout.
- Check the rectangular assembled tangent against the mixed second directional
  finite difference of pointwise hinge energy, and the assembled Newton residual
  against its first directional derivative. Frozen current states have positive
  Jacobians and satisfy their distortion budget; both hinge curvatures remain
  active throughout the perturbations. These energy evaluations do not call
  the assembled metric or load being checked.
- Bind, query and clone both integrators. Check all cached local entries, then
  change the referenced Jacobian strength and rebind the copies. Fresh bindings
  match the copies, cached originals remain unchanged, and both copied local
  systems change. This checks retained parameter references and independent
  matrix/vector storage rather than only clone allocation.
- All ten hinge cases pass individually in optimized normal Release, optimized
  assertion-enabled and optimized AddressSanitizer modes (runs485/489/490).
  Leak detection is enabled; instrumented header/test objects link optimized
  native support libraries. The first mixed fixture compile's inconsistent
  multi-declarator CTAD is corrected by separating different space types; no
  runtime profiles exist for that failed compile.
- A private mutation scales both distinct-space test rows by 1.5. Only the new
  case fails, with finite nonzero tangent values disagreeing with the independent
  energy curvature in both dimensions; all nine prior controls pass (run491).
  Preserve private headers, commands and logs; sanitizer and mutation profiles
  stay excluded from coverage.
- Run492 validates matching source hashes and two normal ten-case producers.
  Hinge.h hits all 40 measured lines; HingeForce.h hits 33 of 34, HingeMetric.h
  57 of 58. All distinct-space assembly statements and bound-polytope assertion/
  return lines are hit. The two remaining misses are the inadmissible frozen-
  deformation skips; no invalid internal state is manufactured to hit them.
  These scoped counts establish no aggregate percentage. Formatting and style
  pass; production sources are unchanged in this batch.
- Run494 preserves nine further completed private builds with checksum-verified
  archives and original-path symlinks, retaining normal/excluded classifications.


### Expand the independently audited CI collector union

- Run495 downloads seven newly completed coverage artifacts, validating every
  artifact's workflow-run identity and published head ba57395. Run496 adds them
  to the three retained PETSc collectors and verifies the frozen source snapshot
  byte-for-byte against git archive of that head. Raw artifact IDs/checksums and
  normalized reports are retained separately from local, newer-source profiles.
- The ten-collector partial union measures 39,367 production lines and hits
  35,581 (90.3828 percent), with no nonproduction records or missing sources.
  Main and other collectors are unfinished; these counts are neither the final
  CI aggregate nor certification of the newer local commits.
- Run499 inventories all 6,188 definitions on that exact published source against
  the partial line union, retaining inactive/uninstantiated bodies. It identifies
  1,636 nonempty/initialized definitions without attribution and 257 measured
  definitions with no hits, plus 94 parser recoveries. These are partial-report
  investigation candidates, not a semantic uncovered-function count. Current
  local forwarding, guard and parsing closures retain their separate evidence.


### MPI coordinate and transformation setter contracts

- Add dimension-aware affine coordinate updates through the scalar setter for
  every held vertex, including ghost copies. Check each coordinate, fresh cell
  physical coordinates/Jacobians after warming caches, and distributed area/
  volume scaling by six/twenty-four in two/three dimensions. All shared copies
  receive the same physical transformation before collective measure queries.
- Register a geometry-preserving tracked transformation only on the final rank.
  Check local map identity, unchanged physical/Jacobian values, and exactly one
  destruction after the mesh releases its ownership. Other ranks proceed to an
  explicit test barrier: registration and mapping queries remain non-collective.
- The first fixture incorrectly calls the dimension-three volume API on a
  dimension-two mesh and crashes before setter evaluation. Symbolized frames
  locate this in getVolume/isOwned. Correct the fixture to select area or volume
  by dimension, confirm no live MPI launcher remains, and retain the entire
  failed producer/profiles under excluded labels (run498).
- Run501 completes all 340 case/rank invocations, one group at a time at ranks
  one, two, three, four and eight. Existing pass/skip contracts are preserved,
  including 793 rank-local skips in older cases. Both new cases execute at every
  rank count without skips. No production source changes are needed.
- A private native Mesh.cpp mutation adds 0.125 to scalar-assigned coordinates.
  At two ranks the new geometry case rejects finite physical-coordinate errors,
  while measure scaling stays correct. The transformation-ownership and original
  vector-setter controls pass (run504). Mutated source, commands and profiles
  remain excluded; the original native library is unchanged.
- Run506 checks the matching native source hash and successful producer, then
  confirms positive counts at both setter entries and their forwarding calls.
  This is file-specific evidence, not a new aggregate percentage. Formatting
  and style pass. Run505 archives three completed private builds with verified
  checksums and original-path symlinks, retaining producer classifications.


### Starling-envelope branches and active-state equilibrium

- Add RodinSolidActiveFiberLawTest with a bounded, continuous normalized-envelope
  check and five segment/plateau equilibrium cases. Expected calibration gains
  initialize independent stationary states satisfying $\gamma^2=n_0k_0$ and
  $\gamma\beta=n_0\sigma_0$. Positive constant activation and zero extension change
  must preserve these states across four time steps and two activation levels.
  This couples the gain branches to physical state evolution rather than merely
  comparing the piecewise implementation with itself.
- All six cases pass separately in normal Release, assertion-enabled and
  AddressSanitizer modes (runs508/512/513), with leak detection enabled and
  optimized native support libraries linked. Discovery and sanitizer profiles
  remain excluded. The first fixture compile's dangling-else warnings are
  corrected with explicit braces around GoogleTest assertions before validation.
- A private header mutation scales only the 0.3-to-0.73 segment by 0.9. The
  continuity case and that segment's stationary-state case reject the defect;
  all four other segment/plateau controls pass. Failed state values remain
  finite and positive (run514). Private mutation artifacts remain excluded.
- Run515 captures the two current normal producers and confirms positive hits
  on all eight previously missed Starling lines. Historical run318 file evidence
  is used only after verifying current header bytes equal retained snapshot274
  and checking the saved run338 line inventory agrees exactly with run318's
  97-line/89-hit record. The scoped union now hits all 117 measured header lines,
  including newly emitted normal Debug lines. No historical aggregate percentage
  is carried onto current source. Production code is unchanged; formatting and
  style pass.


### PETSc field-split copy ownership

- Add three public-contract tests for FieldSplits: independent names and IS
  indices after source mutation/destruction, copy assignment releasing the old
  owned reference and preserving self-assignment, and copying empty/unbound
  default descriptors. The old destination IS retains an external observer
  reference so its reference count can be checked without accessing freed data.
- All six LinearSystem cases pass individually in optimized and
  assertion-enabled builds (runs517/520), linked to optimized native support
  libraries. No sanitizer run is claimed for this batch. Production is unchanged.
- A private, safely reference-counted mutation shares source handles in copyFrom.
  Both copy/assignment cases reject it; the other four controls pass (run524).
  All mutation profiles remain excluded. Run523 records the current header hash,
  passing producers and positive counts for every copyFrom statement, including
  both nonnull and null branches. The scoped header capture hits 81 of 82 emitted
  lines; the remaining line403 follows an always-raising unsupported-merge
  diagnostic. No aggregate percentage is inferred.
- The first coverage verifier (run522) used historical line positions, including
  current braces/blank lines, and correctly failed its positive-count assertions.
  Run523 instead checks current executable statement positions; neither verifier
  changes runtime profiles or source. Discovery profiles remain excluded.


### Source-wide exception idiom audit

- Run518 inventories 566 production files and 6,189 definitions directly from
  source, without preprocessing or requiring template instantiation. It retains
  file hashes, 390 exception records, 216 possible implicit-exception calls and
  94 parser recovery locations. Four audit-tool regressions pass separately.
  The raw inventory and manual review are retained as
  run518-current-exception-source-audit.json and run518-exception-review.json
  under build-coverage-linux.
- There are 305 MemberFunctionException references and no generic
  Alert::Exception references classified inside instance functions. All seven
  explicit throw statements resolve to Alert types: Alert's own self-throw,
  MEDIT's derived UnexpectedKeywordException twice, translated MMG parsing,
  both MPI asSubMesh overloads, and ProblemBase::getBody. The eight MPI
  UniformGrid NamespacedException references are in a static factory, confirmed
  against its header declaration; syntax-only out-of-line classification does
  not recognize the static declaration elsewhere.
- Standard checked-access calls remain a separate API-contract review. Priority
  candidates include AttributeIndex dimension access, local/MPI finite-element
  DOF access, MPIMesh::getGlobalIndex, PETSc ghost lookup and unbound integrator
  getPolytope access. Their presence prevents claiming that every possible
  escaping exception is an Alert. Previously validated List/quadrature/Scotch/
  SubMesh checks and MMG numeric parsing already translate their public errors.
- Several candidates are safe or unrelated: OpenMP optional thread counts and
  QuadratureRule cached values are engagement-guarded; Point's variant access
  follows a check of its other alternative; tuple get and Eigen iterator value
  calls are not throwing optional/variant access. ConstitutivePoint's private
  any values use the same Tag::Type for insertion/access, with a missing-tag
  assertion contract. No blanket catch is added around internal invariants,
  allocation failures, user callbacks or third-party backends. Parser recoveries
  and macro expansion limits remain visible; this is not an exhaustive proof of
  all transitive exception behavior.


### Individual MPI identification geometry watchdogs

- Published-head Tests run38790622528 fails in Slow multithreaded partition2:
  166 of 167 entries pass; RodinMPIAssemblyIdentificationSlowTest_np3 exhausts
  its shared 600-second timer. The verbose rerun passes tetrahedron in about
  462 seconds, pyramid in about43 seconds and hexahedron in about42 seconds,
  then expires during wedge. The job logs and monitor trigger are retained.
- Replace each four-geometry identification group with four individual entries
  at ranks one, two, three, four and eight. Preserve every filter/input, numerical
  check and 600-second watchdog. Set each entry's processor count explicitly.
  Run526 proves that the original geometry-filter union is unchanged and that
  the four Slow partitions remain nonempty, exhaustive and disjoint. Other
  registration command multisets are unchanged. Eight existing connectivity
  parameter names contain ASLR-dependent raw pointer bytes, so names alone are
  unsuitable for the unchanged-command comparison.
- The first local launcher (run527) finds no MPI assembly executable in the
  partial native build; no test case executes. Retain this failed launch under
  excluded labels. Run528 builds the current test source with one optimized
  compiler job, and run529 passes all twenty geometry/rank cases separately,
  one MPI group at a time. The maximum local case time is about155 seconds;
  this is local execution evidence, not a cross-machine performance comparison.
  The dedicated normal producer is retained, with discovery/failed-launch
  profiles isolated. Support libraries remain optimized; no sanitizer claim
  is made for this registration repair.
- Runs531/532 archive completed private validation builds only, with verified
  checksums and original-path symlinks. Producer classifications are preserved;
  these mixed archives must never be captured wholesale. The current-head CI
  failure satisfies the user's condition for pushing the accumulated batches.


### Euclidean source and specialization attribution audit

- Revisit the two historical Box misses against byte-identical current source.
  All 75 existing cases pass individually in a fresh assertion-enabled build
  (runs538/539). Run540 captures 472 of 475 measured Euclidean lines, including
  the same Box misses and one LineSegment continuation-line miss. These scoped
  numbers are not a current global aggregate.
- GCC's per-specialization block graph explains all three lines: Box.h81/125
  and LineSegment.hpp128 contain blocks whose incoming edges are exclusively
  zero-count exception edges in otherwise executed functions. Run541 verifies
  30 such records. Run551 repeats the proof against the final normal producer.
  Keep the raw misses and denominator; no exclusion or invalid-state test is
  introduced to force compiler-generated exceptional cleanup.
- Raw zero-count functions also need interpretation: the initial capture has
  106 zero-count symbols, but grouping identical source locations/C++ signatures
  to expose ABI duplicates leaves 32 zero-count source-signature groups. This
  grouping is an audit aid, not a semantic coverage proof (run542).
- Add Float/Real zero-dimensional ownership checks: mutable/const coordinate
  access borrows the same fixed storage; valid empty affine expressions construct
  and assign the unique point; box corners and copied corners are separately
  owned subobjects. All 77 cases pass individually with assertions and under
  AddressSanitizer with leak detection enabled (runs546/547), linking optimized
  support libraries. Sanitizer/discovery counters remain excluded. Run543 is a
  superseded test-object compile with no runtime producer, and is not captured.
- A private mutation aliases lower() to the upper corner only in dimension zero.
  Both new scalar cases reject the alias, while four original zero-dimensional
  controls pass (run553). The entire mutation producer remains excluded.
- Run549 confirms ten previously zero-count source-signature groups now execute.
  It classifies all 22 remaining groups: four index accesses have no valid index
  in dimension zero, six private tolerance helpers are unused by their dimension
  dispatch paths, and twelve operations require zero-dimensional lines/spheres
  rejected by construction. No rejected internal object is fabricated. Raw
  zero-count symbols fall to 96; ABI duplicates remain explicitly retained.
- The source-only audit reports 105 Euclidean definitions, zero parser recoveries
  in this module, and no nonempty/initialized definition without line attribution
  in the fresh normal capture (run550). Compiler-absent bodies are checked from
  source rather than silently omitted. The scoped line result stays 472/475 even
  though ten specialization groups are closed, illustrating why source-line
  coverage alone does not certify template coverage. Production is unchanged.
- Run537 preserves the completed MPI build with checksum-verified archival and
  original-path links. This test/audit batch stays local until the current
  published-head CI fails or completes.


### CellDeformation alias-safe gradient replacement

- The historical source-attribution guide is restricted to byte-identical files
  and excludes explicitly superseded historical ranges. Run557 initially matches
  515 files; after this production repair, run567 excludes CellDeformation and
  matches 514 files. Its 1,261 unattributed definitions, including 891 without
  constexpr in their signatures, are investigation candidates rather than a
  semantic uncovered-body count. Inlining, constant evaluation, parser recovery
  and unsupported APIs still require individual review. No historical percentage
  is reused to certify current source.
- Inspect the setter's history before changing its established formula. A valid
  public call passes getDeformationGradient() back as the displacement-gradient
  input. The old implementation overwrites that borrowed input with identity
  before adding it, silently producing $2I$ rather than $I+F_{old}$. Run559 observes
  finite wrong entries/Jacobians in dimensions one, two and three; the original
  nonaliasing displacement-gradient control passes. The entire pre-fix producer
  and profiles remain excluded, with the old header/hash retained.
- Evaluate Identity + H as an owned SpatialMatrix value before assigning m_F.
  The existing operator+ returns owned stack storage, so input aliasing is safe
  and nonaliasing arithmetic remains the same. Document that the input may alias
  the owned gradient. Preserve cache invalidation and the setter's return value.
- The regression checks each updated matrix entry, analytic triangular Jacobians,
  stable borrowed-gradient identity, inverse consistency, and independence of a
  copied state whose caches were populated before mutation. All thirteen cases
  pass separately in optimized and assertion-enabled modes (runs563/564), and
  in the final AddressSanitizer build with leak detection enabled (run569).
  Adaptation is header-only; support libraries remain optimized.
- The style check catches missing braces around the new nested assertion loops.
  Add Allman braces without changing their statements before the final sanitizer
  compile. Run562 is a superseded sanitizer compile with no runtime producer;
  run568 builds the final braced fixture. Sanitizer/discovery profiles are excluded.
- Run566 captures only corrected normal producers with matching current header
  bytes: all 74 emitted header lines are hit, including the getter and all four
  cache resets. This is file-specific evidence, not a new global aggregate.
  Published-head CI artifacts for the earlier header are superseded for this
  file. The local batch remains unpushed until that CI run fails or completes.


### Public tensor conjugation and contraction instantiations

- The unchanged-source guide identifies Math::conj(SpatialTensor) without source
  attribution. Existing exhaustive tensor tests call the member conjugate()
  directly, leaving the public free overload uninstantiated in that evidence.
  Route existing numerical contracts through Math::conj and Math::dot, retaining
  member coverage transitively and all axis-extent combinations, including zeros.
- Compute expected complex conjugates/contractions from explicit real/imaginary
  components; expected squared norms use component squares. The oracle no longer
  invokes the same scalar conjugation helper used by the tensor implementation.
  Check conjugate shape/size and involution in addition to every active entry.
- All nine optimized cases and all ten assertion-enabled cases pass individually
  (runs573/574). The extra case is assertion-only. All ten AddressSanitizer cases
  pass separately with leak detection enabled (run575), linking optimized support
  libraries. Sanitizer/discovery counters remain excluded.
- Run576 captures only matching normal producers and confirms positive return-line
  attribution plus entry counts for all six public conjugation specializations:
  Real/Complex tensors at ranks three, four and five. Raw optimized and Debug
  notes, reports and per-specialization data are retained. This is evidence for
  the public template overload, not a new global aggregate. Production is unchanged;
  formatting and style pass. The local batch stays unpushed under the current
  CI deferral condition.


### Spatial-vector Eigen-expression compound updates

- Existing compound-update tests use spatial-vector operands. Add Real/Complex
  Eigen-expression cases covering column and row views in every active size from
  zero through three. Reverse and scale an external fixed-size operand, check
  each expected entry independently, and subtract the same expression to recover
  the initial vector. Verify returned object identity, unchanged logical size,
  preserved inactive storage and an unchanged external operand.
- Inactive values are populated through a valid three-component vector followed
  by public shrinking, and observed through the documented full-storage accessor.
  No private or invalid internal state is constructed. The expressions borrow
  external operands that remain alive throughout each synchronous check.
- All 37 cases pass individually in optimized, assertion-enabled and
  AddressSanitizer modes (runs580/581/582), with leak detection enabled and
  optimized support libraries linked. Sanitizer/discovery counters are excluded.
- Run583 validates positive normal entry counts for all eight compound-expression
  specializations: addition/subtraction, row/column expression forms and both
  scalar types. All four legal sizes execute for each expression form. Internal
  invariant defaults are retained; no artificial state or coverage exclusion is
  introduced. This is specialization-specific evidence, not a global percentage.
  Production is unchanged. Formatting/style pass, and the local batch remains
  unpushed until the current published-head CI fails or completes.


### Spatial-matrix dynamic-vector product scalar pairings

- Confirm the candidate overloads are active code, then exercise direct products
  with dynamic vectors. Existing all-shape checks use a spatial-vector operand,
  leaving these separate overloads without attribution in the historical guide.
  Add all four Real/Complex scalar pairings and all sixteen row/column shapes
  from zero through three, including zero-row and zero-inner-dimension products.
- Check coefficient scalar promotion at compile time and output rows/columns at
  runtime. Expected real/imaginary components are computed explicitly from the
  input-generating constants, including the unconjugated complex-product signs.
  Verify both inputs and full matrix storage remain unchanged after evaluation.
- The first fixture assumes a size() accessor shared by both result types and
  fails to compile. Same-scalar products return a dynamic vector; mixed products
  return a spatial matrix with rows/cols accessors. Correct the fixture to check
  their common shape interface. Run586 is retained under a failed/excluded label
  with its original path preserved; no test executes from that build.
- All 43 cases pass individually in optimized, assertion-enabled and
  AddressSanitizer modes (runs590/591/592), with leak detection enabled and
  optimized support libraries linked. Sanitizer/discovery counters are excluded.
- Run593 verifies positive normal entry counts for four product specializations,
  spanning the same-scalar and mixed-scalar overload families. Numerical checks
  cover all legal shapes without forcing internal invariant defaults. The raw
  reports and matching source hash are retained. This is specialization-specific
  evidence, not a global aggregate. Production is unchanged; formatting/style
  pass. The local batch remains unpushed pending current CI failure/completion.


### CI prepared-report upload repair (run645)

- At published head a8ee3733e, Unit Debug with multithreading and sanitizers
  disabled passes its build, unit, distributed and all 35 slow MPI geometry
  cases. Production coverage capture and artifact upload also pass. The final
  Codecov step receives a runner shutdown signal and exits 143. The origin of
  that shutdown is not established; no test failure is inferred from it.
- Its upload log shows Codecov CLI 11.3.1 running gcov over build notes,
  including third-party and test objects, despite receiving RodinCoverage.info
  with coverage-file search disabled. RodinCoverage has already captured both
  compiler-note baselines and successful normal runtime counters with gcov-14,
  restricted to production sources and excluding src/Rodin/Test.
- Set the official Codecov action's plugins input to noop in all five upload
  steps, so the prepared LCOV report is uploaded without plugin recollection.
  Keep explicit report paths, disabled file search and fail-on-error behavior.
  This does not weaken the 95% production threshold or alter test execution.
- Run645 compares parsed workflows before and after removing the five new
  plugin inputs: all remaining content is identical, including every matrix,
  test command and capture. Official action input documentation and source are
  retained with blob provenance. Production style and diff whitespace pass.
  The new-head CI upload still requires confirmation. Failed test, discovery,
  sanitizer, mutation and pre-fix profiles remain excluded from coverage.

### Normal Debug unit coverage partition repair (run661)

The normal Debug, multithreaded Unit job at `6eec68ed4` passed its unit,
distributed and slow MPI geometry tests, then reached the six-hour runner
limit during coverage capture. The job lasted 21,623 seconds; its unit phase
alone took 12,391 seconds. This is a job budget failure, not a failing numerical
case. No cancelled capture is included in coverage evidence.

Split only that configuration's serial unit selection into two alternating
CTest index partitions. Keep all other configurations as one partition, and
retain the complete distributed and slow geometry selections in partition 1.
Every runtime command remains sequential. Each job checks that its discovered
unit partitions are nonempty, disjoint and exhaustive before running tests.
Coverage artifacts include the partition number; each retains the compiled
baseline, including unexecuted production lines. The collector count becomes
21 when all producer jobs finish.

The configured local optimized complex inventory contains 794 selected unit
cases; its two partitions each contain 397 cases with exact disjoint union.
The workflow matrix grows from 16 to 17 entries, only the affected normal Debug
configuration is duplicated, and every other workflow job is unchanged.
The selection audit ran discovery only, with profiles redirected away from
runtime coverage. Fresh CI must still confirm the affected job's elapsed time
and upload completion.

Recovering the completed run's available reports gives 44,965 hits among
47,737 measured production lines (94.1932%), from 17 of the previous 20
collectors. Tests, benchmarks and library test helpers are excluded. The
missing Debug multithreaded collector and two dependent regular aggregates
make this incomplete evidence, not the final project percentage. The next
coverage target is 96%; prioritize gaps only after accounting for these
missing producers and the source-level template audit.


### Eigen matrix-view product template batch (run672)

Add exact-value tests for both operand orders between Eigen reverse views and
spatial matrices, for Real/Real, Real/Complex, Complex/Real and Complex/Complex.
Each variant checks all 64 legal dimension triples, including empty inner
products and zero result extents. Together they check 512 products. Independent
real/imaginary scalar sums verify unconjugated multiplication and promotion;
poisoned inactive storage and complete input snapshots verify shape, active
bounds and operand preservation. No unsupported dimensions or invalid internal
states are introduced, and production code is unchanged.

All 47 matrix unit cases passed individually in optimized run667 and
assertion-enabled run669. Source-matched normal captures in run672 verify eight
product specializations at SpatialMatrix.h entries 2488 and 2511, each with
128 positive calls across those producers. These are finite supported template
variants, not a claim that every possible Eigen expression has been covered.
The optimized snapshot differs from the formatted test only in whitespace.
Sanitizer and discovery profiles are excluded; no project-wide percentage is
inferred from the focused captures.

Run671 also passed all 47 cases individually under AddressSanitizer with
`detect_leaks=1`. The test translation unit is instrumented, while the reused
Rodin support archive is optimized; this does not establish whole-library
sanitizer coverage. Its profiles remain excluded from measured coverage.


### MFEM matrix diagnostics and loader recovery batch (run680)

Add real/complex P0/P1 matrix-field rejection tests using the existing
independent MFEM coefficient fixture. For each of four supported variants,
eight malformed dimension/ordering headers and twelve empty, truncated or
nonnumeric coefficient payloads must raise a native member-function diagnostic
and leave distinct destination coefficients unchanged. After every rejection,
the same loader accepts independent valid data in both component orderings and
matches every expected coefficient. This checks 80 invalid payloads and 160
valid recovery loads per build. Production code is unchanged.

All 18 MFEM coefficient unit cases passed individually in optimized run675
and assertion-enabled run676. The final optimized build is run679; the earlier
run673 compile-only snapshot precedes loop-style/include/warning cleanup and is
excluded. Only successful final normal profiles enter the run680 source-matched
capture. It verifies all four matrix loader specializations at MFEM.h entry
3922 and positive counters for both native header and coefficient rejection
statements. Compiler line counters combine block/template attributions and
are used as execution evidence, not as unique exception-event counts.
Sanitizer and discovery profiles are excluded; the capture does not establish
a new project-wide percentage.

The initial private ASan run678 failed in the pre-existing scalar round-trip
control before reaching the new tests: an Eigen allocation using the handmade
aligned convention was freed by optimized support code using the direct-malloc
convention. Preprocessor checks showed `EIGEN_MALLOC_ALREADY_ALIGNED=0` for the
ASan translation unit and `=1` for the optimized support archive. Run681 changes
only that macro to `1` in the private ASan compile command; the test source is
byte-identical. Run682 then passed all 18 cases individually with leak detection
still enabled. The controlled command/source comparison and result are retained
in `run682-mfem-asan-compatibility-audit.json`. This corrects the private mixed
build's allocation ABI, not production code or CI flags. The test translation
unit is instrumented and the support archive remains optimized, so no
whole-library sanitizer claim is made. Failed, sanitizer and discovery profiles
remain excluded from coverage.

### Spatial row-vector product template batch (run689)

Exercise spatial row-vector times spatial matrix multiplication for Float,
Real and Complex over every legal row/column shape, including zero inner
extent and the single-row empty result. The 48 products use independent exact
real/imaginary scalar sums, with nonzero complex parts to distinguish ordinary
multiplication from conjugation. Shape and complete operand storage snapshots
check active bounds and input preservation. No invalid internal states or
production changes are introduced.

All 27 spatial-boundary cases passed individually in optimized run685 and
assertion-enabled run686. Source-matched run689 captures verify three product
specializations at SpatialMatrix.h entry 2399, each with 32 positive calls
across those two normal producers. This verifies the named scalar variants;
it does not claim every possible scalar template instantiation. Sanitizer and
discovery profiles are excluded, and no project-wide percentage is inferred.

Run688 also passed all 27 cases individually under ASan with leak detection
enabled. Its private compile matches Eigen's allocation ABI to the reused
optimized support archive, following the controlled run682 diagnosis. The test
translation unit is instrumented; the support archive is optimized, so this is
not whole-library sanitizer coverage. These profiles remain excluded.

### Native HDF5 dataset helper template batch (run696)

Check all four documented native storage types (U64, I32, F64 and U8) using
numeric lowest/maximum values and independent raw HDF5 type, shape and data
reads. Each type writes nine datasets: empty/single/vector/scalar datasets,
three rank-two shapes including zero rows or columns, and two rank-three
shapes including a zero extent. Duplicate vector and matrix writes attempt
different data, must raise native diagnostics, and preserve the existing
schema and values. File-associated object counts return to their baseline
after both successful helper calls and rejected duplicates; this count does
not measure unrelated global HDF5 resources. Scalar readers reject empty and
multi-element datasets. Production code is unchanged.

All 100 HDF5 unit cases passed individually in final optimized run698 and
assertion-enabled run693. Run696 captures verify twenty helper specializations:
four native types for vector/scalar reads and vector/scalar/matrix writes, at
HDF5.h entries 1093, 1140, 1160, 1206 and 1222. These captures include only
successful final normal producers. Their test snapshots are byte-identical;
the earlier review run690/run692, sanitizer and discovery profiles are excluded.
The fixture directly checks 36 datasets; it does not claim all possible shapes
or numeric values. No project-wide percentage is inferred from this audit.

Run695 also passed all 100 cases individually under ASan with leak detection
enabled. Its private compile matches Eigen's allocation ABI to the reused
optimized support archive; the test translation unit is instrumented while
that archive remains optimized. No whole-library sanitizer claim is made,
and sanitizer profiles remain excluded from measured coverage.

### Marginal coverage and next priorities (run699)

Verify that all production sources remain unchanged from `6eec68ed4`, then
merge the eight final normal reports from runs672, 680, 689 and 696 with the
17 available CI collectors for that source snapshot. Check every local raw
report hash and source-header hash before normalizing paths. This diagnostic
hybrid contains no tests, benchmarks or library test helpers; review, failed,
sanitizer and discovery producers remain excluded. It is not current-head CI
confirmation and the three missing CI collectors still prevent completeness.

The incomplete baseline has 44,965 hits among 47,737 measured lines. The local
batches add 58 hit lines and 56 measured lines, giving 45,023 / 47,793 (94.2042%).
The Eigen matrix-view batch adds eight measured/hit lines; the MFEM matrix
rejection batch covers four previously uncovered lines; row-vector products
add 48 measured lines and 46 hits. The native HDF5 batch adds no line hits to
this union while independently verifying twenty template/helper variants and
schema/resource contracts. Template validation and line coverage gains must
therefore be reported separately. Adding tests can expand the measured code
through newly instantiated bodies, including unexecuted fallback statements.

Prioritize reachable line gaps next. The 306 lines in wholly-zero files include
Scotch and manufactured distance/Solid code whose existing tests belong to the
missing regular collectors; do not add duplicate tests just to compensate for
an incomplete report. Suitable independent candidates are native rejection of
MPI-backed fields by HDF5's convenience attribute exporters and duplicate-group
rejection by its direct open-file visualization helpers. Do not exercise
invalid internal geometry indices, synthesize impossible dataspace failures,
or treat compiler cleanup attribution as an untested numerical algorithm.

The detailed marginal/source provenance is retained in
`run699-marginal-coverage.json`, with normalized diagnostic inputs under the
run699 scratch audit. The fixed-denominator shortfall to 96% is only a guide;
new collectors or instantiated code can change both numerator and denominator.

### MPI convenience-export guards and datatype-cache lifecycle (run715)

Use a valid WORLD-distributed segment mesh with 17 grid points (16 cells) to
check HDF5's
convenience node/cell attribute exporters on real scalar, two-component vector
and rectangular matrix P1 fields. Each convenience overload must reject the
MPI-backed FES with native guidance to pass a shard explicitly, preserve every
coefficient, and create no file. The full existing MPI IO suite runs one case
per process group at one rank; the new case additionally runs at ranks 2, 3, 4
and 8. No groups overlap, and the full CI rank registration remains unchanged.

The first leak-enabled ASan run706 passed all one-rank cases, then reported a
cached derived MPI datatype after the new two-rank case's assertions passed.
P1's nonblocking requests already have matching receives and wait_all; the
leak is not an outstanding Rodin request. The test main manually initialized
and finalized MPI, so no owning Boost environment cleared its datatype cache
before shutdown. The installed Boost destructor explicitly clears that cache
before MPI_Finalize when it owns initialization. A minimal two-rank Boost-only
control in run708 reproduces the derived-type leak with manual initialization
and passes with Boost-owned initialization; it does not link Rodin.

Change only this test executable's main to use boost::mpi::environment ownership,
following existing distributed test harnesses. Run710 then passes all 51 MPI
groups with leak detection enabled and the unchanged external-runtime patterns
in cmake/lsan.supp. The controlled old/new test-source comparison differs only
in main's lifecycle handling. No production FES, library code, CI sanitizer
settings or leak suppression patterns are changed. The private test translation
unit is instrumented with the matching Eigen allocation ABI; reused support
archives remain optimized, so this is not whole-library sanitizer coverage.

Final optimized run712 and assertion-enabled run714 also pass all 51 MPI groups.
Run715 captures only their post-repair normal profiles and verifies six exporter
specializations at HDF5.h entries 1726/1780, plus both native rejection guards
at lines 1733/1738 and 1798/1803. Final optimized, Debug and ASan test snapshots
are byte-identical. Failed, sanitizer, discovery, standalone-control and all
pre-repair producers are excluded from measured coverage. No project-wide
percentage is inferred from this focused capture.

### Direct HDF5 duplicate-group rejection batch (run722)

Write topology, coordinates and region attributes from a valid triangle mesh
to an open writable file with its mesh parent group. Attempt each direct writer
again with a quadrilateral mesh having different vertex count, coordinates,
topology and labels. Require native path-specific duplicate-group diagnostics,
then compare the original topology and its size, coordinate shape and values,
attributes in every dimension, an unrelated sentinel dataset, and the baseline
file-associated object count. This exercises public duplicate-write errors;
no invalid handles, corrupted internal mesh state or allocation failure
injection are used. Production code is unchanged.

All 101 HDF5 unit cases passed individually in optimized run718 and
assertion-enabled run719. Source-matched run722 captures verify the three
writer entries at HDF5.h 1374, 1504 and 1582, and positive counters on all nine
targeted rejection lines (1382–1384, 1511–1513 and 1589–1591). Compiler line
counters prove execution, not unique exception-event counts. Only successful
normal producers enter the capture; sanitizer and discovery profiles are
excluded. No project-wide percentage is inferred from this focused audit.

Run721 also passed all 101 cases individually under ASan with leak detection
enabled. Optimized, Debug and ASan test snapshots are byte-identical. The private
ASan translation unit uses the matching Eigen allocation ABI with the reused
optimized support archive; this is not whole-library sanitizer coverage and
its profiles remain excluded.

### Full normal translation-unit coverage refresh (run725)

Expand the six batches' diagnostic captures from their focused headers to all
production headers emitted by their twelve final normal test translation units.
Verify every report hash and every recorded production source against both the
worktree and the frozen `6eec68ed4` source snapshot. Exclude test/benchmark code,
Rodin/Test helpers, sanitizer, discovery, review, failed and pre-repair producers.
Out-of-line support-library counters are not included in these private TU
captures; the incomplete seventeen-collector CI baseline remains separate in
provenance. This hybrid is not current-head CI confirmation.

The result is 45,039 hits among 47,795 measured production lines (94.2337%).
Against that incomplete CI baseline, the six batches add 74 hits and 58 measured
lines. The MPI guard batch adds six hits and two measured lines; the direct
HDF5 group batch adds ten hits, closing its nine rejection statements and the
previously unhit attribute-assignment line 1627. The fixed-denominator shortfall
to 96% is 845 hits, but the missing CI collectors can still change the result.

All raw normal inputs and per-source hashes are retained in
`run723-full-normal-inputs`; the diagnostic marginal audit is
`run725-full-normal-union.json`. Completed sanitizer build directories were
copied, hash-verified and replaced by archive symlinks in run724 to preserve
failure/repair evidence while recovering Docker capacity.

Further Flow audit distinguishes the existing all-geometry cell/face tracking
checks from an untested lower-dimensional start fallback: valid vertices in
2D/3D and edges in 3D are below the supported starting cell/face dimensions.
A focused contract test can check failure/zero results, unchanged input-point
identity and absence of operand/velocity evaluation, without invalid internal
geometry indices or a tolerance chosen only to hit a branch. The tangential
face-selection rounding fallback requires a separate numerical justification;
it should not be forced by inconsistent cached physical/reference coordinates.

### Flow lower-dimensional start contract batch (run733)

Use valid vertices in 2D/3D and edges in 3D to exercise the existing trace
fallback below the supported cell/face starting dimensions. Across six cell
geometries and both nonzero time signs, twenty combinations must abort with
zero result/correction and preserve the starting mesh, polytope and reference
coordinates. A nonzero operand and nonzero velocity track callback counts;
neither may be evaluated on this fallback. Empty reference points are checked
by size and active components rather than reading inactive spatial storage.
No invalid internal indices or inconsistent physical/reference caches are used.

All 15 Flow cases passed individually in optimized run729 and assertion-enabled
run730. Source-matched run733 captures verify the new Flow specialization at
entry 338 and 80 normal executions of the previously unhit fallback line 518.
The initial run726 fixture compile used an unsupported spatial-point equality
operator; it produced no runtime evidence and remains excluded. Final normal
captures use only the corrected component-wise fixture. Sanitizer, failed and
discovery profiles are excluded, and no global percentage is inferred here.
Production code is unchanged.

Run732 also passed all 15 cases individually under ASan with leak detection
enabled. Final optimized, Debug and ASan test snapshots are byte-identical.
The private ASan translation unit uses the matching Eigen allocation ABI with
an optimized support archive; its profiles remain excluded and it does not
establish whole-library sanitizer coverage.

### Advection stop-inside boundary contract batch (run741)

Exercise the previously untested `StopInsideBoundaryPolicy` directly on every
boundary face of Segment, Triangle, Quadrilateral, Tetrahedron, Pyramid,
Hexahedron and Wedge grids. An independent axis-aligned physical normal fixes
the expected displacement, checking the requested inward offset and its time-step
cap for positive and negative time steps. The policy must consume the remaining
time, retain the cell and additive value correction, and become idempotent.
Separate exhausted-time and documented invalid-local-face checks preserve the
reference coordinates, cell and correction. No degenerate mesh is fabricated.

All 11 advection cases passed individually in optimized run735 and
assertion-enabled run737. Source-matched run741 captures from these two normal
producers show 16 previously uncovered lines executed relative to the run725
diagnostic. Raw reports and source hashes are retained in
`build-coverage-linux/run741-advection-inputs/`. This is marginal execution
evidence against an incomplete baseline, not final CI coverage confirmation.
Only production sources are captured; discovery and mutation profiles are
excluded. Production code is unchanged.

In isolated run738, replacing the inward displacement with zero makes the new
physical-offset test fail with 126 rejected offset assertions. The shadow header
and its profiles remain separate from normal producers; no production-header
edit was needed to establish sensitivity.

Run740 passed all 11 cases individually under ASan with leak detection enabled.
The optimized, Debug and ASan test snapshots are byte-identical. The private
ASan translation unit uses the matching Eigen allocation ABI with optimized
support archives; this does not establish whole-library sanitizer coverage.
Its profiles are excluded. Existing manufactured advection tests also exercise
this policy indirectly, so the 16-line diagnostic delta is not a claim of 16
additional lines in the eventual complete CI union.

### Distributed XDMF shape and time-series contract batch (run750)

Use MPI WORLD with independently sized triangular pieces and the last rank as
publisher. A two-snapshot transient mesh series carries a static node tensor,
a transient rectangular cell matrix and a transient node vector. Explicit
rank-aware patterns keep every piece distinct. Check each rank's XML dimensions
and heavy-data references, raw HDF5 dataset ranks/shapes and all nonsymmetric
matrix entries. Changing all fields after the first snapshot must retain the
original static tensor while updating the other values. Rank-dependent values
and piece cardinalities distinguish swapped, duplicated or omitted pieces.

Optimized run744 and assertion-enabled run747 each passed 52 sequential MPI
groups: the full 48-case I/O executable at one rank, followed by the new WORLD
case at ranks 2, 3, 4 and 8. The original CI rank matrix is unchanged. Run749
also passed the same 52 groups under ASan with leak detection enabled and the
unchanged external-library suppression file. The ASan test translation unit
uses the matching Eigen allocation ABI with optimized support archives; this
is not whole-library sanitizer coverage.

Run743 freshly compiled the current XDMF implementation and linked it into
both normal test binaries, so run750 captures include attributable out-of-line
XML execution as well as header templates. Raw production-only reports and
source hashes are preserved in `build-coverage-linux/run750-mpi-series-inputs/`.
Four executed attribute helper specializations include local P1 scalar, vector
and matrix node output and P0 matrix cell output. Five previously uncovered
XDMF lines are positive relative to the incomplete run725 diagnostic: the
non-publisher early return and distributed static-file, tensor/matrix and
vector-shape output. This is not final CI coverage confirmation. Shared header
lines add no marginal line coverage, despite the additional supported types.

The first run742 fixture failed compilation because a GTest macro needed
explicit braces; it produced no runtime evidence. The run750 audit initially
asserted counters at a declaration without a counter and an empty-name helper
fallback. Inspection established that `flush()` replaces an empty grid name
with `mesh` before calling `makeGridPieceName`, so its empty-name return is
unreachable through this caller. The corrected audit verifies actual recorded
XML lines without modifying coverage filters or production code. Failed
fixtures, discovery, sanitizer and mutation profiles remain excluded.

In isolated run752, forcing every snapshot to overwrite the static attribute
makes the new two-rank case fail with 180 rejected value assertions. The
production source remains unchanged and all mutation counters are excluded.
The initial run751 mutation harness failed compilation because an ordinary
include directory could shadow system `hdf5.h` with Rodin `HDF5.h`; a quote-only
include directory fixes the isolated harness without changing normal builds.

### CI unit-partition identity repair (run754)

Published d3aa5f0 failed job 117371174442 after a successful build and discovery,
before runtime tests, at `Unit-test partitions overlap`. The verifier compared
only GTest names. Distinct component executables share names: the P1 executable
and the separate LinearForm, TestFunction and TrialFunction executables each
register matching named cases. CTest partitions invocation registrations, so
name-only sets can report overlap even when executable selections are disjoint.

Compare the full invocation/property record, excluding the report-local
backtrace index, and reject duplicate identical invocation records explicitly.
The runtime commands, labels, index stride, MPI ranks and partition placement
are unchanged. Run754 extracts the exact published and repaired workflow
snippets: a six-case CTest fixture with a name repeated across distinct
component invocations makes the published verifier fail and the repaired
verifier pass with two sets of three. Deliberately overlapping partitions and
omitting a case still fail. The configured optimized complex inventory also
passes with 397 + 397 distinct invocations. Slow/distributed fixture cases stay
excluded. These are discovery checks, with no accepted runtime coverage.
The durable audit is `build-coverage-linux/run754-partition-verifier-audit.json`.
Fresh CI must confirm the repair against its complete configured inventory.

### Matrix reference-element guards and source priorities (run767/run768)

Directly exercise the public matrix-range guards of P0Element, P0gElement,
P1Element and H1Element through order three for Real and Complex, on Point,
Segment, Triangle, Quadrilateral, Tetrahedron, Pyramid, Wedge and Hexahedron.
Both geometry and scalar-element constructor overloads reject zero or excessive
extents with native member-function diagnostics. Every local functional rejects
empty and incompatible values, including a transposed shape with the same flat
component count, before reading a component. The same functional then recovers
with a valid nonsymmetric real/complex matrix, checked against its scalar
functional. Element dimensions and order remain unchanged. No internal indices
or mesh incidence tables are fabricated.

All 22 MatrixRange cases passed individually in optimized run761; the new
case also passed alone with assertions in run764. Source-matched run767 captures
verify 12 constructor variants, 24 rejection/recovery functional variants and
all 20 targeted exception lines. Those 20 lines were unhit in the incomplete
run725 diagnostic. The captured production-only files include four matrix
reference headers and their two included scalar implementation headers;
LCOV's header include patterns also select the latter. Exact allowed-source
hashes and raw reports are retained in
`build-coverage-linux/run767-matrix-element-inputs/`. No project-wide percentage
is inferred from this marginal result. Failed fixtures and discovery profiles
are excluded. Run757 failed a range-loop-copy warning and was stopped after
that fatal diagnostic; run758 revealed a noncopyable scalar P1 functional.
The final fixture uses references for both without changing production code.

Source-priority audit run768 classifies 249 gaps from the incomplete run725
baseline with current/frozen source hashes matched. All 138 P1Element.hpp gaps
are invalid-index/coordinate/geometry fallback markers, asserts or NaN returns.
Of 61 SpatialMatrix.h gaps, 57 are incompatible-size/capacity or unsupported
operation fallbacks, three zero-row subcases are bypassed by the earlier
zero-row return, and one closing-brace counter needs compiler-block attribution.
The 49 H1.hpp gaps contain 46 invalid-index/transform/incidence fallbacks, one
reachable pyramid face-order skip candidate, and two skips incompatible with
valid pyramid/hex face sets. XDMF's empty piece-name return is bypassed by its
sole caller normalizing empty grid names to `mesh`. This is a prioritization
audit, not a blanket claim that invalid public inputs are untestable. Counters,
filters and the coverage denominator are unchanged; further valid-path tests
should follow the completed CI union rather than manufacture invariant failures.
The durable source classification is
`build-coverage-linux/run768-source-gap-priorities.json`.

Run766 passed the new case alone under ASan with leak detection enabled. The
final optimized, Debug and ASan translation-unit snapshots are byte-identical.
The private ASan translation unit uses the matching Eigen allocation ABI with
optimized support archives; it does not establish whole-library sanitizer
coverage, and its profiles are excluded. The completed normal builds and raw
capture directory are checksum-archived with retained symlinks in run769;
the original immutable normal-capture manifest remains separate from the later
sanitizer validation metadata.

### Existing-workload calibration and current source identity (run778/run779)

Rebuild and run the existing Distance (7 cases), Hilbert (3), Scotch (8) and
P1 HyperElasticity (11) workloads from current source, one runtime case at a
time, with optimized numerical builds. All 29 pass. Scotch's current wrapper
and tests are compiled explicitly with retained snapshots and the established
private ignored-attributes warning demotion; this is not a Scotch-enabled
CMake configuration check. No test or production source changes are needed.
Run778 preserves four normal captures, source hashes and complete successful
runtime records in `build-coverage-linux/run778-model-gap-inputs/`.

Run779 re-audits the 17 incomplete old-head CI collectors plus 24 successful
source-equivalent local normal reports. Every local source hash matches both
the current production tree and the frozen 6eec68ed4 snapshot. The earlier
run725 checkpoint reproduces exactly before adding the newer reports. The
combined diagnostic measures 47,799 production lines with 45,331 hit
(94.836712%). Existing workloads contribute 250 additional hit lines and four
additional measured lines; the newer Flow/advection/XDMF/reference-element
batches contribute another 42 hit lines without increasing that checkpoint's
denominator. Thus whole-file zero-hit gaps in these models were substantially
missing-producer evidence, not a reason to duplicate their existing tests.

At this diagnostic's fixed denominator, 79 additional hits would reach 95%
and 557 would reach 96%. Neither number establishes the final CI shortfall:
normal Debug and regular aggregate reports are still missing from the old CI
baseline. The current published-head run must supply all 21 collectors before
its union can confirm the target. Tests, benchmarks, library test helpers,
failed, mutation, discovery and sanitizer producers remain excluded. Raw input
provenance, the detailed delta and remaining gaps are retained in
`build-coverage-linux/run779-current-normal-union.json` and its audit companion.

Current-run source normalization is independently checked: completed job
117392874622 checked out synthetic merge da9b5c8364253286c7fb61698497a7ffb7526f96.
Its GitHub Git tree exactly equals published head c65b9848c, and current develop
cf9ff43f9 is an ancestor of that head. The coverage snapshot therefore represents
the tested tree rather than relying solely on the workflow's reported head SHA.
The durable attestation is
`build-coverage-linux/run779-ci-source-identity.json`. This proves source
identity, not coverage completeness or CI completion.

### SWIFT representable geometric-scale guard batch (run797)

Exercise finite positive reference sizes whose FE-order geometric scale
underflows to zero or overflows to infinity. Both automatic and explicit
positive geometric tolerances must reject those unusable reference scales
with native solve-member diagnostics. Reuse the existing translated-line
fixture at P1 and H1 orders two and three for twelve rejection combinations.
Strengthen the separate existing geometric-constant control case to check the
initial-target reason, explicit target, no linear solve and zero displacement.
Its three positive solves are not repeated in the new rejection case. No
production code or solver template body changes are needed.

All final 49 SWIFTSolver cases passed individually in optimized run792; both
focused cases passed with assertions in run794. Focused validation contains
verbatim includes, fixture and positive/negative test blocks from the full test
source, with hashes and complete original snapshots retained. This does not
claim that full and focused translation units are byte-identical. The controls
perform geometry checks and initial-target exits, not unoptimized iterative
PDE solves.

Run797 verifies the three existing default-fixture solver specializations and
positive counts on exception lines 892 and 893. Raw final-source normal reports
and successful runtime records are retained in
`build-coverage-linux/run797-swift-scale-inputs/`. Pre-review, failed, discovery,
mutation, diagnostic and sanitizer profiles are excluded. Three additional
MUMPS-availability guard lines are also positive from an existing case in the
full suite; they are not attributed to the new test. No project-wide percentage
is inferred from this file-specific capture.

The initial fixture duplicated the existing positive controls after its
rejections. Its default ASan run passed but took about 117 seconds; the same
immutable binary reproduced about 119 seconds wall/119 seconds user time,
zero major faults and no swap. Two stopped GDB samples landed in ASan's
FakeStack::Allocate through Eigen evaluator construction, SpatialVector zero
storage and H1 Jacobian quality sampling. Changing only the fake-stack option
in an excluded diagnostic reduced the same case to about 0.76 seconds, with
leak detection still enabled. Retain the extra positive assertions in their
existing case and isolate the new rejections so the PR adds no duplicate
quality-sampling work. Production code and normal sanitizer settings remain
unchanged. Original default-ASan acceptance passed before diagnosis; the
changed-option control and intentionally killed GDB run are not acceptance or
coverage evidence. The causal record and stack log are retained in
`build-coverage-linux/run798-swift-scale-asan-cause-audit.json` and its referenced
excluded diagnostic records.

Run796 passed both final focused cases individually under the normal
`ASAN_OPTIONS=detect_leaks=1` setting. The final focused Debug/ASan manifests
match the complete optimized source snapshot and each verbatim fixture/control/
rejection block. The private ASan TU retains matching Eigen allocation ABI and
optimized support archives; it does not establish whole-library sanitizer
coverage. Only final normal runs 792/794 feed the run797 capture; all pre-review,
ASan, changed-option timing controls and GDB profiles remain excluded.

### Finite non-strict MPI reconciliation caps (run806)

Check that repeated reconciliation with finite non-strict owner/global-ID
round caps preserves a fully converged mesh. Cover triangle, quadrilateral,
tetrahedron, hexahedron, pyramid and wedge intermediate entities. Compare
global entity counts against an independent local mesh, and compare global
IDs, owners, halo maps and owned/shared/ghost flags against the converged
distributed baseline. Compute the required connectivity explicitly before
reconciliation. SELF uses zero rounds; WORLD uses its communicator size with
a larger check period to exercise finite epoch clipping. These checks concern
preservation of converged metadata, not convergence from insufficient rounds.

Optimized run801 and assertion-enabled run803 each passed fourteen sequential
MPI process groups: all ten MeshTest cases at one rank and the new case at
two, three, four and eight ranks. Both the test TU and MPI Mesh.cpp were freshly
compiled in each configuration; their snapshots match the current sources.
Run806 captures only these successful normal producers, excludes test/helper
sources, and verifies positive counts on six previously uncovered finite-cap
lines: 853, 857, 859, 975, 979 and 981. Raw reports, source hashes and runtime
records are retained in `build-coverage-linux/run806-mpi-cap-inputs/` and the
companion file audit. No production change or coverage-filter change is made.

ASan run805 also passed all fourteen sequential groups with leak detection
enabled and the unchanged external MPI-library suppressions. Both the private
test and MPI Mesh.cpp objects are instrumented, using matching Eigen allocation
ABI; the remaining support archives are optimized. This does not establish
whole-library sanitizer coverage. The final optimized, Debug and ASan test and
implementation snapshots all match the current files. Sanitizer profiles are
excluded from the coverage capture.

The run807 diagnostic union adds the final SWIFT normal reports and these MPI
reports to run779: 47,799 measured lines and 45,342 hit lines, approximately
94.8597 percent. The denominator is unchanged; the two batches add five and
six hit lines respectively. At this fixed denominator, another 68 hits would
reach 95 percent and another 546 would reach 96 percent. This combines seventeen
incomplete old-run CI collectors with thirty source-matched local normal
reports; it is not current-head CI confirmation. Failed, review, discovery,
mutation and sanitizer profiles remain excluded, as do tests, benchmarks and
library test helpers. Detailed deltas and remaining gaps are retained in
`build-coverage-linux/run807-current-normal-union.json` and its audit companion.

### Targeted dense mixed manufactured assembly (run816)

Extend the existing affine P1/P0 dense mixed manufactured workload to call
both LHS-only and RHS-only assembly, twice each, starting from unrelated
matrix/vector values. Compare the selected side with the full system already
verified against the manufactured field; require exact preservation of the
opposite side. Cover Sequential and OpenMP on triangles and quadrilaterals.
The existing full-system residual and recovered-solution checks remain in place.

All four cases passed individually in optimized run811 and optimized ASan
run815 with leak detection enabled. Run818 also passed the two OpenMP cases
individually with two worker threads. ASan uses the matching Eigen allocation
ABI and optimized support archives; this does not claim whole-library
sanitizer coverage. Final normal, mutation and sanitizer test snapshots match
the current source. The earlier run809 compile-only snapshot precedes a
formatting adjustment and has no acceptance or coverage role.

Run812 force-includes private copies of the two assembler headers with only
the final mixed target branch changed to always select the matrix. Both
backend mutation cases failed, with four expectations each, in run813.
Workspace production sources remain unchanged. Mutation profiles are excluded.

Run816 records eight calls to each supported targeted dense mixed backend
variant. Optimized Eigen moves have no separate assignment-line counters;
the function, full-assembly call, target branch and return line are positive.
Selected-side values and opposite-side preservation are established by the
numerical assertions and the mutation check. Source hashes, successful normal
runtime records and the production-only report are retained in
`build-coverage-linux/run816-dense-mixed-target-inputs/` and its file audit.
Only the final normal run811 contributes to this capture; discovery, sanitizer,
mutation and pre-review profiles are excluded, together with tests, benchmarks
and library test helpers.

Run817 adds that report to the previous diagnostic union: 47,803 measured
lines and 45,378 hits, approximately 94.9271 percent. The additional 36 hits
comprise eight previously missed targeted-overload lines, twenty-four other
previously missed lines from the existing full-assembly workload, and four
newly measured hit lines. Do not attribute all 36 to the new assertions.
At this fixed denominator the remaining planning shortfalls are 35 hits to
95 percent and 513 to 96 percent. The union still contains seventeen incomplete
old-run CI collectors; it is not current-head CI confirmation. Detailed deltas
and remaining gaps are retained in
`build-coverage-linux/run817-current-normal-union.json` and its audit companion.

### Boundary-normal error diagnostics and geometry recovery (run830)

Test collapsed triangle edges and tetrahedron faces through both stack and
dynamic output overloads. Require the native normalize-member diagnostic,
repair the shared mesh coordinates, then reuse the same normal object and
check the independently known outward coordinate vector. Separately check
both endpoints of a one-dimensional mesh: this implementation rejects their
one-row Jacobians with a native member-function diagnostic.

The endpoint test was observed red in optimized run822: SpatialMatrix::rows()
returns an eight-bit integer, which streamed as a control character instead
of the number one. Cast that diagnostic value to size_t before insertion.
No numerical algorithm or exception family changes. The exact row-count
expectation now passes. The retained pre-fix log is
`build-coverage-linux/run832-boundary-normal-before-fix.log`.

All nine cases passed individually in optimized run825, assertion-enabled
run827 and ASan run829. ASan retains leak detection and the matching Eigen
allocation ABI; optimized support archives remain outside its instrumentation.
An isolated private-header reversal of the normalized vector was rejected
by four new recovery expectations in run832. No mutation touches workspace
production sources; the excluded failure log is retained alongside the
pre-fix log. Final optimized, Debug, ASan and mutation test snapshots match.
Run821's failed pre-fix producer and run823's pre-fix compile-only snapshot
have no acceptance or coverage role.

Run830 verifies positive counts on all eight error-diagnostic lines in the
two guards. Raw successful fixed-source normal reports, source hashes and
runtime records are retained in `build-coverage-linux/run830-boundary-normal-inputs/`
and its file audit. Failed, pre-fix, discovery, mutation and sanitizer profiles
are excluded, as are tests, benchmarks and library test helpers.

The production-header change invalidates its older counters. Run833 discards
all nineteen old BoundaryNormal records before adding only fresh successful
normal runs for that header. Every other production source remains hash-matched
to the frozen CI source; the changed header has an explicit new source snapshot.
The replacement header has 130 measured/112 hit lines versus the old 125/99.
The resulting diagnostic union is 47,808 measured and 45,391 hit lines,
approximately 94.9444 percent: another 27 hits would reach 95 percent and 505
would reach 96 percent at that fixed denominator. This is still an incomplete
old-CI/local diagnostic, not current-head CI confirmation. Detailed replacement
provenance and gaps are retained in
`build-coverage-linux/run833-current-normal-union.json` and its audit companion.

Refresh the preprocessing-independent source inventory after the header fix:
566 files, 6,189 definitions, 390 exception references and 216 implicit exception
candidates. The 94 parser recovery locations and other limitations remain
explicit in `build-coverage-linux/run833-current-source-audit.json`. This keeps
inactive configurations and uninstantiated definitions visible; it does not
claim execution of every specialization or convert the inventory into a
source-coverage percentage.

### Power-degree metadata and automatic quadrature (run841)

Exercise the public Pow degree query with signed, unsigned, boolean and
floating exponent types on real P1 fields. Cover zero, negative and fractional
exponents, unknown base degree, maximum representable degree, and the exact
boundary before degree multiplication overflows. Construct large degrees by
nested expressions rather than a fabricated metadata provider. Huge powers
are queried only for degree metadata; they are never evaluated or passed to
quadrature. Check the conservative unknown-order policy explicitly.

Use automatic quadrature to integrate powers of the affine field 1+x over
the unit interval for exponents zero, one, two, three, six and nine. Sum the
P0 load entries and compare against the independent antiderivative formula,
computed with ldexp rather than Pow. No quadrature-order override is used.
All nineteen cases passed individually in optimized run836, assertion-enabled
run838 and ASan run840 with leak detection enabled. ASan retains the matching
Eigen allocation ABI and optimized support archives. The optimized test
snapshot differs only in one verified EXPECT_EQ argument line wrap; its C++
tokens are unchanged. Final Debug, ASan and mutation test snapshots match.

An isolated private-header mutation returning degree zero for every ordinary
nonzero integral power was rejected by the new automatic-quadrature case in
run843: fourteen expectations failed, including integral-value comparisons.
The workspace production header is unchanged; the mutation log is retained
in `build-coverage-linux/run843-pow-order-mutation-excluded.log`.

Run841 confirms seven executed degree-query specializations spanning int,
unsigned int, size_t, bool and Real exponents, with positive unknown-base,
negative, zero, overflow and product-return counters. Raw successful normal
reports and source hashes are retained in `build-coverage-linux/run841-pow-order-inputs/`
and its file audit. Failed, discovery, mutation and sanitizer producers are
excluded, together with tests, benchmarks and library test helpers.

Run844 adds seventeen hit lines and two measured lines to run833. All sixteen
previously missed Pow lines become positive; its floating-exponent unknown-order
return is newly measured and hit. The other new measured line is the empty
OpenMP linear-form copy-constructor body, emitted but unhit for this workload.
Both attributions are retained in `build-coverage-linux/run844-new-line-attribution.json`.
The diagnostic union is 47,810 measured and 45,408 hit lines, approximately
94.9759 percent. At this fixed denominator the planning shortfalls are twelve
hits to 95 percent and 490 to 96 percent. This retains the changed-header
counter reset from run833 and remains an incomplete old-CI/local diagnostic,
not current-head CI confirmation. Detailed provenance and gaps are retained
in `build-coverage-linux/run844-current-normal-union.json` and its audit companion.

### Public mixed targeted assembly with independent block orders (run853)

Extend the dense P1/P0 manufactured fixture through the public virtual
ProblemBase::assemble(target) entry point. Keep trial blocks ordered P1/P0
but order test blocks P0/P1, whose sizes differ on both meshes. Assert trial
and test offsets independently and compare against an explicit row permutation
of the already verified manufactured system. Start a fresh problem for each
target, without prior full assembly, seed unrelated matrix/vector values,
repeat the call, and require exact opposite-side preservation.

Public Problem chooses its backend at configuration time, so this block runs
only in the fixture case matching that default: OpenMP locally, Sequential
in a configuration without OpenMP. The direct assembler checks continue to
exercise both enabled backends; this does not claim a public Sequential run
under the local OpenMP configuration.

All four manufactured cases passed individually in optimized run847 and ASan
run849 with leak detection enabled. Run852 also passed the two OpenMP cases
individually with two worker threads. ASan uses matching Eigen allocation ABI
and optimized support archives; it does not establish whole-library sanitizer
coverage. Final optimized, ASan and mutation test snapshots match exactly.

Run850 force-includes a private Problem header that reverses only the public
mixed target dispatch, with valid indices. The new public checks rejected it
with eight failed expectations in run851; both selected-side numerical values
and opposite-side preservation are checked. Workspace production sources stay
unchanged. The excluded mutation log is retained in
`build-coverage-linux/run851-public-mixed-target-mutation-excluded.log`.

Run853 confirms eight executed calls to the permuted mixed public target
specialization. Only successful normal run847 contributes to its production
capture; tests, benchmarks, library test helpers, failed, discovery, mutation
and sanitizer producers are excluded. Raw reports, source hashes, function
execution counts and runtime records are retained in
`build-coverage-linux/run853-public-mixed-target-inputs/` and its file audit.

Run854 adds twenty-three previously missed Problem.h lines without increasing
the measured denominator. The diagnostic union is 47,810 measured and 45,431
hit lines, approximately 95.0241 percent. It crosses the local 95 percent
planning threshold; another 467 hits would reach 96 percent at this fixed
denominator. The changed-header reset from run833 remains in effect. This is
still an incomplete old-CI/local diagnostic, not current-head CI confirmation.
Detailed provenance and gaps are retained in
`build-coverage-linux/run854-current-normal-union.json` and its audit companion.

### Nonmatrix P0 MFEM coefficient-count errors (run862)

Cover real/complex scalar and vector P0 loaders with empty, truncated,
overlong, first-invalid and partially invalid coefficient payloads in both
MFEM layouts. Require the native load-member diagnostic with exact expected
and observed counts, preserve distinct existing field coefficients, then reuse
the loader for both independent valid layouts. The new typed cases perform
forty rejections and eighty recovery loads per build.

All twenty MFEM coefficient cases passed individually in optimized run857,
assertion-enabled run859 and ASan run861 with leak detection enabled. ASan
retains matching Eigen allocation ABI and optimized support archives; it does
not establish whole-library sanitizer coverage. All final test snapshots match.
Run862 verifies the four nonmatrix P0 loader variants and positive counts on
all three previously missed count-diagnostic lines. Successful normal inputs,
source hashes and runtime records are retained in
`build-coverage-linux/run862-mfem-p0-count-inputs/` and its file audit.

Run866 uses a private-header mutation that clears the destination immediately
before raising the native count exception. Both typed cases reject it in
run867, with twenty failed preservation checks each. The audit is retained in
`build-coverage-linux/run867-mfem-p0-count-mutation-audit.json`. The first
mutation compile in run863 used a normal include path for Rodin/IO, causing
the system hdf5.h lookup to select Rodin's HDF5.h on the mounted source tree.
Use quote-only lookup for that private dependency directory. The failed
compile has no runtime or coverage role. No workspace production code changes.

Run865 adds three hit lines without changing the measured denominator:
47,810 measured and 45,434 hit lines, approximately 95.0303 percent. Another
464 hits would reach 96 percent at that fixed denominator. The changed-header
reset from run833 remains in effect; this is still an incomplete old-CI/local
diagnostic, not current-head CI confirmation. Detailed provenance and gaps are
retained in `build-coverage-linux/run865-current-normal-union.json` and its audit
companion. Failed, mutation, discovery and sanitizer producers remain excluded,
together with tests, benchmarks and library test helpers.

### Current CI calibration and PETSc signed body snapshots (run881)

Run870 incorporates fifteen successful collectors from current published-head
CI alongside the seventeen older collectors and source-matched local reports.
The additional CI data expose 249 measured lines and 212 hits: the prior
95.0303 percent diagnostic becomes 94.9791 percent at 48,059 measured/45,646
hit lines. This is denominator expansion, not evidence that the tests stopped
executing previously covered code. Source hashes match; all thirty-two older
records for the changed BoundaryNormal header are discarded. The detailed
per-file delta is retained in `build-coverage-linux/run870-current-ci-delta.json`.

Run873 associates that incomplete line inventory with the preprocessing-
independent source audit. Twelve fully unhit definitions have at least four
measured body lines, excluding parser recoveries. Candidates include the two
PETSc body operations exercised below; definitions and lambdas overlap, so
their counts must not be summed as project gaps. The absent IDRSTABL and
legacy TraceOperator findings already documented above remain unchanged.
Run869 independently reproduces the configured Eigen IDRSTABL include failure
with a syntax-only probe and no runtime profiles.

Add an exact P1 mass-system check for PETSc body-minus-bilinear-form and
body-plus-linear-form operations. The expected operator is twice a prescribed
scalar times the analytic two-element mass matrix; the residual convention
gives the negative prescribed load moments. Change the original preassembled
matrix/vector after forming the body to verify snapshot ownership, and assemble
twice. The complex case uses nonzero imaginary coefficients.

The first private named-form compile selected complex PETSc for the existing
real-valued-space component and failed the intended scalar compatibility
checks; run871 has no runtime or coverage role. Run872 freshly compiles all
54 current native support units needed by the real component, excluding the
unneeded MMG/Scotch backends, with one compiler, optimization and assertions.
All twenty-two named-form cases pass individually against that support in
run873. Only the 24 linked fresh native objects with runtime profiles, plus
the test object, feed the real capture.

Run875 passes the new complex case using verbatim includes, test block and
main from the complete source; full snapshots and block hashes are retained.
Its optimized configured complex support archive counters are not captured.
Run877 also passes the focused complex case under ASan with leak detection,
matching Eigen allocation ABI and unchanged external MPI suppressions. The
support archives remain optimized; this is not whole-library sanitizer coverage.
An isolated private-header mutation dropping both negative scaling signs is
rejected by twenty numerical expectations in run880. No production source
changes. The excluded mutation log is retained in
`build-coverage-linux/run880-petsc-body-sign-mutation-excluded.log`.

Run878 confirms four executed body-operator variants (two operations, real
and complex) and sixteen named-form OpenMP execution variants. Successful
normal reports and source/producers are retained in
`build-coverage-linux/run878-petsc-body-inputs/` and its file audit. Discovery,
failed, mutation and sanitizer producers are excluded, together with tests,
benchmarks and library test helpers.

Run881 freezes the same fifteen current CI inputs for an attributable delta:
33 additional hits and two measured lines. Sixteen previously missed PETSc
body-operation lines and sixteen PETSc OpenMP assertion lines become positive.
The remaining hit and new miss are compiler-generated class-declaration
attributions in AssemblyBase and TrialFunctionMismatchException respectively.
Detailed deltas and attribution are retained in
`build-coverage-linux/run881-petsc-delta.json` and its new-line companion.
The diagnostic is 48,061 measured and 45,679 hit lines, approximately 95.0438
percent; another 460 hits would reach 96 percent at this fixed denominator.
This remains incomplete CI/local evidence, not final CI confirmation. Detailed
provenance and gaps are retained in
`build-coverage-linux/run881-current-normal-union.json` and its audit companion.
Run884 refreshes the source audit against these counters, retaining 566 files,
6,189 definitions, 390 exception references and the explicit parser/reachability
limitations in its audit and prioritization artifacts.

### Embedded MMGS surface optimization batch (run893)

Exercise the previously unhit surface optimizer on a distorted four-triangle
patch in the plane $z=x+2y$. An independent triangle calculation checks positive
orientation, total area $\sqrt{6}$, strictly improved minimum triangle quality,
material references and preservation of the four corner coordinates. Compare
coordinate sets after rebuilding rather than assuming stable mesh indices.

All 118 MMG cases pass individually in optimized assertion-enabled run888.
The verbatim focused case passes Debug run890 and ASan run892 with leak
detection and matching Eigen allocation ABI. All nine native MMG C++ wrappers
are rebuilt for each configuration; other support archives and the MMG C
backend remain optimized. The initial optimized smoke is excluded from capture.
The isolated run895 mutation skips MMGS while reporting success; precisely the
quality-improvement assertion rejects the unchanged patch. Production sources
remain unchanged. Its log is retained in
`build-coverage-linux/run895-mmg-surface-mutation-excluded.log`.

Run893 captures only the two successful normal producers, retaining full test
snapshots, focused-block hashes, current native-source hashes and paired
objects/profiles. It confirms two executions of optimizeMMGS. Tests, benchmarks,
library test helpers, discovery, initial smoke, sanitizer and mutation profiles
are excluded from measured sources. Inputs and provenance are retained in
`build-coverage-linux/run893-mmg-surface-inputs/` with its file audit companion.

Run896 freezes the same fifteen current CI inputs used by run881. The batch
adds twelve hits: eleven native optimizer lines and one dispatcher line. A
newly measured, uncovered MMG5.h line increases the denominator by one. The
result is 48,062 measured and 45,691 hit lines, approximately 95.0668 percent;
another 449 hits would reach 96 percent at this fixed denominator. The union
still combines incomplete CI with successful local runs, not final CI proof.
Detailed provenance and the per-file delta are retained in the run896 audit
and delta artifacts. The current CI watcher reports no test-job failures;
normal Debug partition two passes while partition one remains in progress.

### MMG builder move-assignment batch (run908)

Move a three-dimensional two-tetrahedron builder into a populated planar
builder, destroy the source, then finalize the destination. Check replacement
of geometry, dimension, coordinates, cell count and all six MMG constraint
sets against explicit expected values. The assignment returns the destination
reference. Different old and new indices expose stale or missing metadata.

All 119 MMG cases pass individually in optimized run901; the focused case
passes Debug run903 and ASan run905 with leak detection. Native MMG objects
are reused from the preceding source-hashed current builds with matching
configuration and Eigen allocation ABI; test objects are rebuilt from exact
current snapshots. The C backend and remaining support archives stay optimized.
Run907 omits the required-triangle move in an isolated private implementation;
the new test rejects the stale index set. Its excluded log is retained in
`build-coverage-linux/run907-mmg-builder-mutation-excluded.log`.

Run908 pairs the successful normal profiles with their actual test/native
objects and verifies two builder move-assignment executions. No sanitizer,
mutation or discovery counters enter coverage. Tests, benchmarks and library
test helpers remain excluded. Retain raw inputs and source/producers in
`build-coverage-linux/run908-mmg-builder-inputs/` and its file audit companion.
Run909 adds nine hits in MeshBuilder.cpp, with no denominator expansion:
48,062 measured and 45,700 hit lines, approximately 95.0855 percent. Another
440 hits would reach 96 percent at that denominator. The same fifteen current
CI inputs are frozen for this comparison; this remains incomplete diagnostic
CI/local evidence. Detailed provenance and the delta are retained in run909
artifacts. CI has 164 successful jobs and no test-job failures at this check.

### Cosine polynomial-order template batch (run918)

Exercise cosine order metadata for literal constants, P0 fields, positive-degree
P1 fields and unknown-order callbacks. Equal-valued constant and unknown-order
operands must retain distinct metadata. Verify GetOrderIfConstant on the same
public function variants, compare point values with the standard cosine, and
check the analytic integral of the cosine of two cellwise constant values under
automatic quadrature.

All 23 cosine cases pass individually in optimized run913, Debug run915 and
ASan run917. Leak detection is enabled and the Eigen allocation ABI matches the
optimized support archives; their counters are not captured. Both new cases
reject an isolated sign-reversed cosine implementation in run920. Production
headers remain unchanged; the excluded mutation evidence is retained in
`build-coverage-linux/run920-cosine-sign-mutation-excluded.json`.

Run918 captures only successful final optimized and Debug Test objects and
confirms thirteen executed order-method variants across Cos and
GetOrderIfConstant. A null-state check and retry-directory handling in the
capture script were repaired before its successful capture; runtime producers
and their profiles were unchanged. Raw reports and source/producers are retained
in `build-coverage-linux/run918-cosine-order-inputs/` with the file audit
companion. Tests, benchmarks, library test helpers, discovery, sanitizer and
mutation counters are excluded.

Run922 freezes the preceding fifteen current CI reports and adds ten hits,
five each in Cosine.h and Function.h, without denominator expansion. The
diagnostic is 48,062 measured and 45,710 hit lines, approximately 95.1063 percent;
another 430 hits would reach 96 percent at this denominator. It remains
incomplete CI/local evidence. Provenance and per-file delta are retained in
run922 artifacts. The latest CI snapshot has 167 successful jobs and no
test-job failures; CI confirmation remains pending.

### Native hybrid MMG copy ownership batch (run932)

Construct disjoint tetrahedral and prismatic cells through the native MMG API,
with triangle, quadrilateral and edge boundary records. Supply the public native
boundary-cache records without running a remesher. Copy the mesh, verify distinct
allocations, change the source references and boundary data, destroy the source,
then check retained coordinates, connectivity, references and cache indices.
This tests native snapshot ownership, not remeshing support for hybrid cells.

All 120 MMG cases pass individually in optimized run925. The verbatim focused
case passes Debug run927 and ASan run929 with leak detection and matching Eigen
allocation ABI. Current native wrapper objects are reused only after their
source hashes match; matching ASan wrapper objects instrument allocation and
copy operations. The C backend and other support archives remain optimized.
Run931 omits quadrilateral copying in a private implementation; the new case
rejects the null destination allocation. Production sources remain unchanged.
The excluded mutation log is retained in
`build-coverage-linux/run931-mmg-hybrid-mutation-excluded.log`.

Run932 confirms execution of all six formerly missed native buffer-copy
branches and eight total copyMesh calls across its successful normal producers.
Inputs and source/producers are retained in
`build-coverage-linux/run932-mmg-hybrid-inputs/` with the file audit companion.
Tests, benchmarks, library helpers, discovery, sanitizer and mutation counters
are excluded. Run933 adds twelve hits in MMG5.cpp without adding measured lines:
48,062 measured and 45,722 hit lines, approximately 95.1313 percent. Another
418 hits would reach 96 percent at this denominator. Run936 separately adds the
sixteenth current CI report, slow-mtON-part4, with no additional lines or hits.
Neither union is final CI confirmation.

Run937 refreshes the preprocessing-independent source inventory against those
counters: 566 files, 6,189 definitions and 390 exception references, retaining
94 parser recoveries. Six fully unhit measured definitions have at least four
body lines after excluding parser recoveries. A separate list retains 706
unmeasured template definitions for classification; this is not a count of
missing runtime tests, because constant evaluation, inlining and inactive or
unsupported variants can suppress compiler attribution. Definitions and lambdas
overlap. Inventory, limitations and priorities are retained in the run937
artifacts. The default host Python lacked the parser modules; the audit succeeds
with the existing source-audit dependency directory on PYTHONPATH. CI has 183
successful jobs and no test-job failures at this check; the first normal Debug
partition has passed its serial unit stage and is running distributed cases.

### Analytic rotating Flow characteristics and template audit (run946)

Trace the autonomous field $v=(1/2-y,x-1/2)$ through triangle and quadrilateral
meshes. Compare signed quarter turns and full revolutions against the analytic
rotation about $(1/2,1/2)$, including the transported affine operand, zero
unconsumed time/correction, no boundary callbacks and preservation of the start
point. Fixed RK4 substeps and tight face-bisection settings isolate curved
characteristics from the existing constant-velocity tests.

All sixteen Flow cases pass individually in optimized run939, Debug run941 and
ASan run943. Leak detection is enabled and Eigen allocation ABI matches the
optimized support archives; their counters are not captured. Run945 reverses
the signed velocity in a private header; the analytic coordinate and operand
checks reject it. The excluded log is retained in
`build-coverage-linux/run945-flow-rotation-mutation-excluded.log`.

Run946 captures only the successful final normal Test profiles, retaining exact
current test snapshots and production-source hashes. The actual new trace
method executes 32 times across the two configurations. Its audit distinguishes
that method from 65 related emitted allocator/lambda functions whose demangled
names also contain the outer trace signature. Inputs and provenance are retained
in `build-coverage-linux/run946-flow-rotation-inputs/` with the file audit
companion. Discovery, sanitizer, mutation, tests, benchmarks and library helpers
are excluded from measured sources.

Run947 adds two hits at Flow.h lines 848 and 872: bracket-expansion exhaustion
and unsuccessful face-bracketing return. Other tracing paths were already
covered. The diagnostic is 48,062 measured and 45,724 hit lines, approximately
95.1355 percent, using the same sixteen current CI reports. Another 416 hits
would reach 96 percent at this fixed denominator; this is not final CI proof.
The per-file delta is retained alongside its union/provenance artifacts.

Run950 independently instantiates the public dense and sparse LinearSystem
merge templates with an offset argument. Both compile with the configured
compiler/dependencies; this is syntax-only evidence, with no runtime or profiles,
not numerical validation. Full probe sources, commands and diagnostics are
retained in its excluded audit artifact for the next runtime audit.

Run951 classifies all seven explicit production throw expressions as Alert
exception construction or forwarding, including MEDIT's Alert-derived
UnexpectedKeywordException. Eight NamespacedException sites flagged by the
syntactic instance-member heuristic belong to static MPI UniformGrid methods;
the header declarations establish their static context. No production exception
change follows from those flags. The 216 potential indirect exception sites
remain callee/precondition candidates, not findings or a whole-program proof.
CI continues without test-job failures at this check.

### Dense/sparse merge offsets and lifecycle coverage (run973)

The source audit's unmeasured public merge templates pass syntax probes but
fail independent numerical checks for their documented nonzero offset. Four
real/complex dense/sparse cases agree at zero offset and reject offset two in
run953. The oracle forms $T^*AT$ and $T^*b$, imposes the constraint equation,
and checks the full operator, right-hand side, solved vector and projected
residual. Constraint weights are real, including positive and negative values.
The stored solution vector is preserved. Empty maps are no-ops.

Both implementations ignored the offset. Preserve local indices in the graph
queue/map and apply the offset to every constrained/supporting matrix and
vector index, including final constraint rows. No algorithm or exception
policy changes. The regression body observed failing before the correction is
byte-identical in the final expanded unit source; run976 retains its proof.
The excluded four-case failure evidence is retained in
`build-coverage-linux/run953-merge-prefix-failures-excluded.json`.

All nineteen initial cases pass optimized run955, Debug run957 and ASan run959;
all four existing manufactured mixed-assembly cases pass individually in
optimized run961 against the corrected header. These preliminary unit profiles
are superseded in the final capture. Run962/963 initially replace the header
with 242 measured/hit lines, but run964 identifies 31 previously measured,
unchanged sparse lifecycle/const-access lines absent from that attribution.
Its final mapping uses the native header's full path to exclude PETSc.

Add real/complex dense/sparse copy construction, copy assignment, move
construction, move assignment and self-assignment checks, with all three arrays
verified through const access. Clear source arrays after copying to establish
independent ownership. All 23 final cases pass individually in optimized
run968, Debug run970 and ASan run972 with leak detection and matching Eigen
allocation ABI. Remaining support archives stay optimized; their counters are
not captured. The four manufactured cases need no repeat: production header
and manufactured source are unchanged after their successful run.

Run973 captures only final normal unit producers 968/970 and manufactured
producer 961. It verifies four merge instantiations and 32 entries. Exact test
snapshots and fixed-header hashes/timestamps are retained. All 39 older native
LinearSystem header records are discarded rather than reused across the change.
Run975 confirms that fresh profiles recover all 31 attribution gaps. Final
inputs/provenance are retained in
`build-coverage-linux/run973-merge-lifecycle-inputs/` and its file audit companion.
Tests, benchmarks, library helpers, discovery, pre-fix, preliminary unit and
sanitizer profiles are excluded from the final capture.

Run974 retains the existing BoundaryNormal reset and freezes the same sixteen
current CI inputs. It reports 48,207 measured and 45,872 hit lines, approximately
95.1563 percent. The fixed native header has 273 measured/hit lines. Relative
to run947 the union adds 145 measured lines and 148 hits; the other three hits
are two OpenMP assembly lines and one Sum line. Another 407 hits would reach
96 percent at this denominator. This is incomplete diagnostic CI/local evidence,
not final exact-head confirmation. Run977 refreshes the preprocessing-independent
source audit against final counters; uninstantiated and inactive variants still
require classification even when a measured header reports complete coverage.
The latest CI snapshot has 245 successful jobs and no test-job failures.

### Public linear-system base forwarding and sparse storage (run994)

Instantiate mutable and const accessors, merge and elimination through each
real/complex dense/sparse system's public Parent interface. Verify object
identity, independent transformed operators/right-hand sides, an offset weighted
constraint followed by an offset essential condition, solved residuals and
preservation of the stored solution vector.

Run981 fails to compile the complete new test source because the CRTP derived
argument names LinearSystem of the scalar types rather than the concrete
matrix/vector system. Correct the two inheritance declarations and two Parent
aliases. The dense case then passes, but the sparse case crashes after merge
in run983. GDB run984 identifies raw outer-column traversal in eliminate.
Diagnostic run985 shows uncompressed columns with unused reserved capacity:
capacities 4/5/6/4 versus active counts 3/3/4/2. Raw compressed-array traversal
reads those unused slots as row indices. Establish compressed storage before
obtaining elimination's raw pointers. No constraint algebra changes.

Retain the excluded pre-fix compile log, intermediate crash audit, GDB trace and
storage diagnostic in the run981/run983/run984/run985 artifacts. The same final
test source observed failing before repair passes afterward. All 27 unit cases
pass individually in optimized run987, Debug run989 and ASan run991, with leak
detection and matching Eigen allocation ABI. Both original public-base syntax
probes also pass in run998. The four existing manufactured mixed-assembly cases
pass individually in optimized run993 against the final header. Support archives
remain optimized; their counters are not captured.

Run994 captures only successful final current-header normal unit/manufactured
producers. It confirms 32 base-forwarder instantiations (eight methods across
four system types) and four merge variants with forty entries. Retain exact test
snapshots, header hashes/timestamps, inputs and provenance in
`build-coverage-linux/run994-base-storage-fixed-inputs/` and its file audit
companion. Failed compile/crash runs, GDB, diagnostics, discovery, sanitizer,
tests, benchmarks and library test helpers are excluded from measured sources.

Run995 discards all 42 preceding native LinearSystem header records. Run1000
separately adds the seventeenth current CI report, optional-Release-mtOFF-part1,
with no additional measured lines or hits. Its old-header resets discard 43
LinearSystem records and 34 BoundaryNormal records. Fresh profiles supply 290
measured/hit native header lines. Run1001 verifies that all 273 previously
measured unchanged lines remain measured and hit, using source correspondence
without historical counter reuse.

The diagnostic is 48,224 measured and 45,889 hit lines, approximately 95.1580
percent. This batch adds seventeen measured lines and seventeen hits; another
407 hits would reach 96 percent at this denominator. The percentage remains
incomplete CI/local evidence, not final exact-head confirmation. Run1002 refreshes
the preprocessing-independent inventory: 566 files, 6,189 definitions, 390
exception references and 94 retained parser recoveries. CI continues without
test-job failures; coverage capture for the first normal Debug partition remains
in progress at this check.

### Unsupported MFEM collection diagnostics (run1010)

Exercise the existing fallback loader and printer for scalar-real, vector-real
and scalar-complex global-constant spaces. Reject MFEM with the specific
member-qualified Alert messages, without changing coefficient data, input
position/state or an existing output prefix. Save and reload the same fields
through HDF5 afterward to verify supported persistence still works.

The initial private compile used a nonexistent vector P0g deduction guide;
run1003 is excluded. Use the existing explicit SpatialVector range instead.
All seven IO format cases pass individually in optimized run1005, Debug run1007
and ASan run1009. Leak detection is enabled with matching Eigen allocation ABI;
support archives remain optimized and their counters are not captured. No
production changes or numerical tolerance changes are involved.

Run1010 confirms six executed rejection variants, with two entries each across
the normal configurations. Successful final normal profiles and exact source/
producer hashes are retained in `build-coverage-linux/run1010-mfem-fallback-inputs/`
and its file audit companion. Tests, benchmarks, library helpers, failed compile,
discovery and sanitizer profiles are excluded.

Run1011 verifies that production sources still match the frozen run1000 source
snapshot, preserves both changed-header resets and appends only the two new
normal reports. Seven prior misses in GridFunctionLoader.h and four in
GridFunctionPrinter.h become hits, with no denominator expansion. The diagnostic
is 48,224 measured and 45,900 hit lines, approximately 95.1808 percent; another
396 hits would reach 96 percent at that denominator. It remains incomplete
CI/local evidence, not final exact-head confirmation.

Run1016 separately confirms that both normal Debug partitions on published head
c65b9848cefbdae9059f4b3576ccd95a918d6d1f completed successfully, including their
coverage stages. This validates the earlier partition repair, not the later
local commits. Their exact CI jobs and step conclusions are retained in
`build-coverage-linux/run1016-ci-debug-split-confirmation.json`.

### Integrator-list and boundary-condition ownership (run1025)

The previously unexecuted list-plus-Dirichlet ProblemBody overload now has
independent constrained-mass checks for real and complex scalars in dense and
sparse storage. Destroy the source list and condition before assembly, verify
two retained integrators and one condition, then compare the complete operator
and nonzero boundary-induced load against analytic P1 mass coefficients.

Both component cases pass individually in optimized run1015, Debug run1018
and ASan run1020, with leak detection enabled and matching Eigen allocation
ABI. Support archives remain optimized and their counters are not captured.
An isolated header mutation drops the second integrator; both real and complex
cases reject it through the matrix and load oracles (run1022). Mutation,
discovery and sanitizer profiles are excluded. No production files changed.

Run1025 confirms both scalar overload instantiations execute four times each
across the normal configurations. Its raw reports and source/producer hashes
are retained under `build-coverage-linux/run1025-problem-body-inputs/`. The
first capture attempt used an incorrect test-source path and produced no
admitted report; the capture verifier was then corrected to check executable
statements rather than the closing brace, which remains a zero-count compiler
attribution even on successful returns.

Run1026 preserves the source resets and verifies every production file against
the unchanged snapshot. Five existing misses become hits, with no measured-line
expansion: 48,224 measured, 45,905 hit, approximately 95.1912 percent. Another
391 hits reach 96 percent at that denominator. This is diagnostic CI/local
evidence, not final exact-head CI confirmation.

### Instantiated scalar P0 gradient contracts

The source audit identified an exposed but unmeasured scalar P0 gradient
specialization. Instantiating real and complex point evaluation failed against
the existing implementation (run1027): the complex output was hardcoded real,
P0Element no longer provided getGradient, GridFunction no longer accepted a
coefficient index through getValue, and an Optional attribute was passed to a
set of concrete attributes. Preserve the failed compile, test source and exact
pre-fix header in `build-coverage-linux/run1027-p0-gradient-prefix-evidence/`;
no counters from that producer are admitted.

Use the operand scalar field for gradient storage, the existing basis
getDerivative interface, indexed coefficient access, and explicit optional
attribute handling. Preserve one-sided trace selection and Alert diagnostics.
The new real and complex tests assign distinct cell constants, including
nonzero imaginary coefficients, then check exact zero gradients on cells,
boundaries and each selected interface side across all seven positive-
dimensional geometries. Copies, moves and clones retain the original field;
missing and removed attributes reject the requested trace.

All 17 gradient cases pass individually in optimized run1030, Debug run1032
and ASan run1034. Leak detection is enabled with matching Eigen allocation ABI;
linked support archives remain optimized and their counters are not captured.
A private header mutation adds a spurious unit derivative, and both new tests
reject it through exact value checks (run1042). The first mutation harness
forced the specialization ahead of its required declarations and failed to
compile; that attempt is excluded. The successful harness substitutes the
private header after the public declarations without editing production files.

Refresh existing matrix tests sharing the header: 22 optimized matrix-range
cases, four optimized and four Debug matrix derivative lifecycle cases, and
12 Debug derivative-trace cases all pass individually. Fresh normal counters
replace every historical P0/Grad.h record; earlier BoundaryNormal and
LinearSystem resets remain in force. Source hashes, producer results and
executed variants are retained in the run1039, run1046 and run1053 input
folders. Tests, benchmarks, library test helpers, failed builds, discovery,
mutation and sanitizer counters are excluded.

The broad Debug matrix-quadrature refresh was narrowed after encountering an
unrelated potential quadrature case. The entire interrupted run1055 producer
is excluded, including its completed first case. Run1059 instead executes the
existing complex Jacobian and curved/embedded/point Jacobian cases individually
with fresh profiles; no tests or numerical tolerances are changed for this
refresh.

Two further published-head Debug CI artifacts add 29 previously missed hits
in run1043. Those reports reset all three locally changed headers before
normalization. The final local repair still requires exact-head CI after a
later authorized push.

Run1057 proves that all 117 previously hit matrix lines remain measured and
hit after their verified five-line source shift. No historical P0 header
counters are reused. The refreshed source audit still inventories 566 files,
6,189 definitions and 390 exception references, with 94 parser recoveries.

The final diagnostic union (run1056) has 48,294 measured and 45,998 hit lines,
approximately 95.2458 percent. It includes 19 current published-head CI
artifacts and 67 successful local normal reports, in addition to the older
verified inputs. Another 365 hits would reach 96 percent at this denominator;
this remains incomplete diagnostic evidence rather than exact-head CI proof.

### Invalid HDF5 current extents (run1067)

Audit existing dataset-helper misses before adding tests. Missing paths,
read-only files, duplicate datasets and incompatible numeric conversions
already have coverage. The focused guard inventory in run1066 distinguishes
backend-failure-only guards from a reachable invalid-current-extent path.

Exercise unlimited current row, column and depth extents with zero payload
cardinality. Native HDF5 rejects these current extents before dataset creation.
Require the Alert dataspace diagnostic, no new link or file-owned handles,
unchanged original matrix values and shape, and a successful subsequent write.
No production code or numerical tolerances change.

All 102 HDF5 component cases pass individually in optimized run1061, Debug
run1063 and ASan run1065. Leak detection is enabled with matching Eigen
allocation ABI. The first sanitizer-launch permission review expired before
process creation; the authorized retry completed normally. Support archives
remain optimized, and sanitizer/discovery counters are excluded.

Run1067 retains exact source/producer hashes and only the two successful
normal Test profiles. It verifies all four matrix-dataset scalar variants;
the new case checks three invalid extents in each normal configuration. The
diagnostic insertion line has six hits, while the Raise line is attributed
only by Debug (three hits). Do not interpret lexical line counters as a
one-to-one event count. Raw reports are retained in
`build-coverage-linux/run1067-hdf5-extents-inputs/`.

Run1068 preserves all source resets and adds two existing hits at HDF5.h
lines 1237 and 1238, with no measured-line expansion. The diagnostic union is
48,294 measured and 46,000 hit lines, approximately 95.2499 percent; another
363 hits reach 96 percent at this denominator. It is not final exact-head
CI confirmation. The paired delta artifact identifies both newly hit lines.

### Reject unsupported MMG dimensions before allocation (run1085)

The audit candidate in run1070 is confirmed by an isolated leak-enabled
probe: createMesh with a one-dimensional volume reaches the unsupported
backend exception after allocating a native mesh and leaks 568 bytes. The
exact new regression bodies reject three positive unsupported dimensions,
including the largest size_t, against the pre-fix implementation and expose
1,704 leaked bytes in three allocations (run1076). Both diagnostic producers
and the earlier failed compile are excluded from coverage. Their pre-fix
sources, bodies and leak reports are retained in the run1072/run1076 evidence
folders. A copied structured binding warning in the first regression compile
was corrected before validation.

Validate topology/embedding pairs before allocation or backend initialization.
Only planar 2D, surface 2D in 3D, and volume 3D are supported. Unsupported and
zero dimensions now receive the existing class/function-qualified MMG Alert
exception with supplied dimensions; the static factory has no instance for a
MemberFunctionException. Document native ownership and supported dimensions.
The tests also reject invalid embeddings and recover through all three
supported conversion paths. The accepted backend initialization and ownership
paths are unchanged.

All 122 MMG component cases pass individually in optimized run1082, Debug
run1083 and ASan run1084. Both the native MMG5.cpp and test translation units
are compiled freshly; their objects precede the static MMG archive at link.
ASan instruments both units with matching Eigen allocation ABI and leak
detection enabled. Vendor and other support archives remain optimized. The
exact regression bodies match the observed pre-fix failure. No numerical
tolerances change.

Discard every historical MMG5.cpp and MMG5.h record before the new capture;
the previous three source resets remain intact. Run1085 captures only the
current successful normal native and test objects. The first intermediate
union lacked reverse solution-copy template counters instantiated in Adapt.cpp.
Refresh that unchanged consumer, rerun all 122 cases individually (run1090),
and capture only its successful Debug object (run1091). Its reused native and
test objects have matching current source identities; their duplicate profiles
are not admitted. Do not use the incomplete intermediate percentage as a gain.

Run1093 verifies that every previously hit unchanged line in both files is
measured and hit by current-source profiles. MMG5.cpp now has 610 measured/562
hit lines and MMG5.h 44/40. The parser audit still inventories 566 files,
6,189 definitions and 390 exception references, with 94 recoveries; it is not
a whole-program instantiation proof. Tests, benchmarks, helpers, pre-fix
diagnostics, failed compile, discovery and sanitizer counters are excluded.

The final diagnostic union (run1092) has 48,293 measured and 46,008 hit lines,
approximately 95.2685 percent. Another 354 hits reach 96 percent at that
denominator. It remains incomplete diagnostic evidence, not final exact-head
CI confirmation.
