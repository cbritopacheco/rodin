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
  endpoint being discarded before the bracket check, returning no result for
  a linear function. Both new typed endpoint cases fail on the preserved
  pre-fix binary, while its fourteen other cases pass. Endpoint roots are now
  accepted before any interior iteration; the header also owns its required
  standard includes. Sixteen cases pass individually in Release, normal Debug
  and ASan/LSan, covering invalid initial guesses, zero/absent derivative
  safeguards and depleted work budgets. All 21 heart-model caller regressions
  pass individually in Release, including branch-flow derivatives and dynamic
  Jacobian finite differences. A doubled returned-root mutation fails all four
  new initial-guess/derivative cases; endpoint and zero-budget controls pass.
