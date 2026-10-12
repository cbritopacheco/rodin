# Benchmark module architecture

Library benchmarks mirror the owning source module under
`tests/benchmarks/Rodin/<Module>/`. MPI and PETSc specializations retain
their source-side hierarchy rather than being mixed into the local module.
Placement follows the operation being measured: spatial-value arithmetic is
in Math, while matrix-valued finite-element assembly is in Assembly.

```text
tests/benchmarks/
├── CMakeLists.txt
├── README.md
├── Workloads/                  Shared forms and numerical oracles
├── Rodin/
│   ├── Math/                   Spatial arithmetic and contractions
│   ├── QF/                     Reference quadrature
│   ├── Geometry/               Meshes, connectivity and transformations
│   ├── Variational/
│   │   ├── P1/                 P1-specific space operations
│   │   └── H1/                 High-order reference kernels
│   ├── Location/               AABB construction and queries
│   ├── IO/                     Mesh loading
│   ├── Assembly/               Local/global assembly and physical workloads
│   ├── Solver/                 Solves on fixed assembled systems
│   ├── MPI/
│   │   ├── Geometry/           Distributed construction
│   │   └── Location/           Distributed location
│   └── PETSc/
│       └── Assembly/           PETSc local/distributed assembly
└── convergence/                Verification-support profiling, not library kernels
```

Each module owns its CMake targets. The independently selectable
`RodinMathBenchmarks`, `RodinQFBenchmarks`, `RodinGeometryBenchmarks`,
`RodinVariationalBenchmarks`, `RodinLocationBenchmarks`,
`RodinIOBenchmarks` and `RodinAssemblyBenchmarks` partition the historical
registration set. Component profilers retain their existing target names.
For example, `RodinFieldOperationBenchmarks` belongs to Variational and
`RodinConstrainedSolveBenchmarks` belongs to Solver.

All module and component executables are opt-in (`EXCLUDE_FROM_ALL`), with
no CTest registration. The compatibility executable `RodinBenchmarks`
retains the original source ordering, benchmark names and publication
workflow; new standalone profilers are not added to that workflow.
Executables retain their default `build/tests/benchmarks/` output location,
independently of source placement. An explicitly configured runtime output
directory is respected.

A module executable has its own registration and cache-initialization history.
Timing comparisons must retain the executable, selection order and cold/warm
protocol; a module timing is not automatically comparable with a historical
aggregate timing.

For example:

```sh
cmake --build build --target RodinGeometryBenchmarks -j1
build/tests/benchmarks/RodinGeometryBenchmarks --benchmark_list_tests=true
cmake --build build --target RodinFieldOperationBenchmarks -j1
build/tests/benchmarks/RodinFieldOperationBenchmarks \
  --benchmark_filter='Field/Interpolation/H1/P2/Real/Scalar/Pyramid/9/'
```

Shared workloads contain mathematical forms and oracles, not assembly-backend
policy. References from `tests/convergence` remain test-only dependencies.
The pressure-spectrum oracle profiler is kept under `convergence/` because
its implementation belongs to verification support, not to `src/Rodin`.
Future benchmark directories should follow newly measured source components;
empty directories do not imply performance coverage.

Profiling proceeds from spatial algebra, reference elements and quadrature
through mapped geometry, fields, expressions and element integrators to
global assembly. Solving is measured separately on the assembled operators.
An isolated improvement is checked again through its consuming library
operation. Correctness regressions remain in unit/manufactured/convergence
suites; benchmark numerical guards and separately instrumented source coverage
do not replace convergence certification.

## Hierarchical profiling ledger

Benchmark availability, numerical qualification, repeatable measurement and
source profiling are separate deliverables. Registration counts are not a
percentage of this workplan. A row is complete only when its applicable cases
have qualified numerical guards, controlled measurements, attributed profiles,
and consumer-level confirmation of any repaired hotspot. The existing
pre-NAFEMS convergence certification is a correctness foundation, not evidence
that this performance matrix has been profiled.

| Layer / owning module | Available measurement | Remaining qualification or isolation |
| --- | --- | --- |
| Spatial values / Math | Real/complex vector, matrix and tensor arithmetic; spatial dot; fixed-capacity construction/copy/assignment and temporary/direct array fill | Remaining arithmetic consumer attribution and controlled scaling baselines beyond selected storage scopes |
| Reference metadata / Geometry | Topology, connectivity and reference transformation consumers | Isolated metadata access, embedding and shape-dependent work accounting |
| Reference elements / Variational | Scalar P0/P0g/P1 and H1 degrees 1–6; values, derivatives, prefetched descriptors | Vector/matrix component expansion, DOF functionals and first-use initialization |
| Reference quadrature / QF | Owned construction/copy, polynomial traversal, controlled eight-key hot and nine-key canonical-pool sweeps; table/fallback boundaries | Process-first allocation, concurrent pool contention and profiles beyond selected orders |
| Transformations / Geometry | Pointwise/bulk Jacobians; pointwise map consumer profiles | Independent geometry degree, curvature, embeddings and inverse-map stage matrix |
| Mapped points / Geometry | Isolated warm/rebound coordinates, Jacobian, inverse and measure; mapped-quadrature construction; field consumer profiles | General curved embeddings, cache-miss provenance, controlled timings and profiles beyond selected cases |
| Space setup / Variational, MPI | P1 construction and distributed construction | P0/P0g/H1 setup, DOF lookup/restriction, reconciliation and degree scaling |
| Fields / Variational | Interpolation and true mass projection; 1,008 native and OpenMP selections qualified; point/integration-point sweeps | P0/P0g, higher degrees, matrix evaluation, point cells, norms and boundary restriction |
| Expressions / Variational | Generic expressions, coefficients and named forms | Cold expression construction, binding, changes of state and complete range/backend matrix |
| Local integrators / Variational | Volume entry kernels and reference-loop comparisons | Boundary/interface kernels; mixed blocks and nonlinear residual/tangent stages |
| Operator setup / Assembly | Public repeated assembly and sequential stage scopes | Sparsity discovery, allocation and replacement separated from kernels |
| Scatter / Assembly, PETSc | Sequential triplets; complete native/PETSc insertion consumers | Isolated OpenMP buffers/merge, PETSc insertion and ownership-dependent MPI work |
| Constraints and completion / Assembly, PETSc | Constrained scalar assembly and fixed-system solves | Vector/mixed constraints, identification, finalization and communication attribution |
| Public physical assembly / Assembly | Diffusion, conductivity, LinearElasticity, mass/reaction and nonlinear Poisson | Helmholtz, coupled reaction–diffusion, Stokes, boundary variants and remaining degrees/backends |
| Solvers / Solver, PETSc | Real/complex fixed-system constrained CG solves; projection CG profile | Setup versus iteration/factorization, indefinite/mixed systems, backend scaling and memory |
| Complete workflows / owning modules | Consumer-level interpolation/projection and assembly profiles | Matched setup/assembly/solve/postprocessing accounting across physical contexts |

The required axes are scalar field (real/complex), supported scalar/vector/matrix
range, P0/P0g/P1 and H1 degrees 1–6, reference-cell family, field and geometry
degree, embedding dimension, and native/OpenMP/PETSc/MPI execution. Applicability
is mathematical: point cells do not have reference derivative sweeps, global
constants do not support diffusion solves, and Taylor–Hood requires a stable
mixed pair. Unsupported Cartesian products are not missing measurements.
At least three sizes are required for mesh-based scaling; reference-element
sweeps instead use three quadrature orders. Cold/warm cache state and changing
coefficient state are distinct experiments. Darcy and NAFEMS remain deferred.

For a sample count $N_q$, basis count $N_b$ and reference dimension $d$,
value and derivative sweeps write $N_qN_b$ and $dN_qN_b$ scalar entries,
respectively. Timings are reported alongside these counts, not interpreted as
dimension costs in isolation. Assembly additionally records cells, global DOFs,
local entry count and nonzeros; solver comparisons retain residuals and iteration
counts. A before/after repair must preserve the numerical observations and be
checked in its consuming operation, not only in an isolated microbenchmark.

### Scalar reference-element protocol

`RodinReferenceElementBenchmarks` measures all eight reference geometries for
real/complex P0, P0g, P1 and H1 degrees 1–6 at orders 2, 4 and 8. The two derivative
sweeps are omitted on points, leaving 2,484 registrations. P0/P0g use the common
reference-table utility explicitly; this does not imply that their production
integrators consume that cache. Direct and prefetched output entries are compared
exactly. An independent affine/constant reproduction check applies each family's
public DOF functionals instead of assuming that coefficients are nodal values.

`Values` and `Derivatives` include descriptor lookup in the sample loop;
`PrefetchedValues` and `PrefetchedDerivatives` exclude that lookup.
`WarmTable` measures lookup alone. `CapacityWarmMiss` cycles nine formula
identities through the eight-slot cache after allocation warmup: it measures
rebuilding, not process-cold initialization. The cache-hit case is warmed after
the eviction sequence. Output buffers, input formulas and reference guards are
outside timing. Derivative tables include reference derivatives, not mapped ones.
The untimed miss protocol verifies absence by logical formula identity before
each rebuild. Complex reconstruction uses a nonzero imaginary amplitude.
Warm lookup costs depend on cache-slot history; isolated process selections
must be used for matched comparisons rather than interpreting a full-suite
lookup ramp as geometry-dependent work.

```sh
cmake --build build --target RodinReferenceElementBenchmarks -j1
python3 dev/check_p1_descriptor_initialization.py \
  build/tests/benchmarks/RodinReferenceElementBenchmarks
build/tests/benchmarks/RodinReferenceElementBenchmarks \
  --benchmark_filter='ReferenceElement/P1/Real/(Values|PrefetchedValues)/Pyramid/8/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

The initialization diagnostic requires an unstripped binary and GNU/LLVM-style
demangled `nm` output. It checks all 66 real/complex P1 descriptor objects and
rejects runtime initialization guards; it is a logical binary-contract check,
not a timing threshold or a replacement for numerical tests.

### Field evaluation and pyramid attribution

`RodinFieldOperationBenchmarks` additionally measures `PointEvaluation` and
`IntegrationPointEvaluation`. The two scopes use P1 and H1 degrees one/two,
real/complex scalar/vector ranges, all seven positive-dimensional geometries
and grid sizes 3, 5 and 9: 504 additional selections, or 1,512 in total.
Both sweep the same preconstructed mapped samples, coefficients and output
buffers through the public `GridFunction::getValue` interface. Analytic field
values and physical-coordinate caches are prepared outside timing. Sample
construction, interpolation, assembly, solving and numerical checks are not
timed. Changing coefficients while preserving sample identities is checked
before timing, and a perturbed result must be rejected.

The integration-point sweep carries formula provenance but visits successive
points and cells; it is not a warm repeated-point cache-hit test. The point path
uses structured component contractions where the element provides them, whereas
the current integration-point cache constructs full range-valued basis entries.
`point_evaluations`, `scalar_basis_evaluations_per_point` and
`scalar_output_entries` make these different work counts visible. The scope is
local field evaluation even in an OpenMP-enabled build, not parallel evaluation
or distributed ghost synchronization.

`RodinH1PyramidBenchmarks`, owned by `Variational/H1`, isolates index decoding,
binomial coefficients, Bernstein values/derivatives and modal values/derivatives
at degrees one through six. The four sampled scopes use quadrature orders 2, 4
and 8; the two metadata scopes use their degree-dependent index sets, giving
84 selections. These are intrinsic real kernels also consumed by complex fields.
No nodal Vandermonde contraction, reference-table lookup, mesh or assembly is
timed. They complement, rather than replace, the nodal reference-element and
complete field-evaluation measurements. Mode indices and Bernstein arguments
are read from runtime arrays exposed through compiler barriers, preventing
compile-time loop bounds from folding metadata into precomputed outputs.
Reported sweeps include those barriers, argument reads and output writes;
they are not isolated instruction latencies.

The independent oracle uses Pascal addition for binomial coefficients and
polynomial powers for
$B_i^n(t)=\binom{n}{i}t^i(1-t)^{n-i}$. For $q=1-z$, $a=x/q$ and $b=y/q$,
the pyramid modal value is
$B_i^{K-k}(a)B_j^{K-k}(b)B_k^K(z)$; reference derivatives are obtained by
the chain rule. Sampled quadrature points are interior, so $q>0$. Mode triples
are enumerated independently in layer order. Metadata is checked exactly;
floating values use the stated dimensionless roundoff budget, with explicit
perturbation rejection. Raw times must be interpreted alongside output counts.

An optimized Apple LLVM 21 local measurement at H1 degree two and grid size
five illustrates the required accounting. At order eight, triangle, hexahedron
and pyramid sweeps visit 512, 8,000 and 18,048 points, with 6, 27 and 14 scalar
basis functions per point. Scalar point evaluation measured approximately
76, 199 and 1,190 ns per point, respectively, in five repetitions. Thus the
roughly 554-fold pyramid/triangle sweep-time ratio combines 35.25-fold sample
growth with approximately 15.7-fold per-point cost; neither ratio is a universal
dimension penalty. An Instruments profile of the pyramid consumer attributed
about 86% of selected CPU sample weight to modal evaluation, Bernstein factors,
binomial coefficients and index decoding. Each nodal basis currently evaluates
all $N_b$ modes independently, producing $N_b^2$ modal evaluations per point.
This identifies a source-level reuse investigation, not a completed repair or
an assembly/solver speedup. Warm reference-table consumers have a different
cost model and must be profiled separately.

The corresponding vector integration-point sweep measured about 3,773 ns per
pyramid point versus 1,217 ns for vector point evaluation. Its profile placed
about 79% of selected CPU sample weight in the same modal/Bernstein/index work.
The integration-point path evaluates 42 component bases where the structured
point path evaluates 14 scalar bases. This is a consumer-specific duplication,
not evidence for slow mapping, MPI sharding or a solver bottleneck. Profile
percentages describe selected source CPU samples, including untimed guards
executing those kernels, and must not be interpreted as end-to-end speedups.

```sh
cmake --build build --target RodinH1PyramidBenchmarks -j1
build/tests/benchmarks/RodinH1PyramidBenchmarks \
  --benchmark_filter='H1/Pyramid/(BernsteinDerivatives|ModalDerivatives)/P2/8/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
build/tests/benchmarks/RodinFieldOperationBenchmarks \
  --benchmark_filter='Field/(PointEvaluation|IntegrationPointEvaluation)/H1/P2/Real/Scalar/Pyramid/5/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

### Mapped-point protocol

`RodinMappedPointBenchmarks`, owned by Geometry, separates physical coordinates,
Jacobians, left inverses and measure factors. `Warm` sweeps reuse each selected
point cache. `Rebound` sweeps call `setPolytope` before each access, including
cache invalidation and lazy computation but excluding point allocation. Both
write preallocated outputs; their times are not pure accessor latencies. The
points in these scopes do not carry formula provenance, so rebound map and
Jacobian evaluation are pointwise rather than reference-table operations.

`Construction` measures fresh owned `PolytopeQuadrature` allocation, bulk
Jacobians and destruction/replacement with warmed reference tables.
`ConstructionAndCoordinates` additionally consumes all lazy physical
coordinates. These are separate scopes, not independent additive costs.
Mesh-owned cache lookup, first-use basis initialization, reference-table
eviction, assembly and solving are excluded. Zero-dimensional cases check
empty Jacobian dimensions and unit measure; their inverse sweeps are omitted.

Geometry degrees one through six, all eight reference families, every
admissible embedding through dimension three, and orders 2, 4 and 8 are
registered. Each workload has one cell; quadrature-order sweeps are sample-count
studies, not mesh refinement. Affine maps use $\alpha=0$, and quadratic maps use
$\alpha=1/8$ at degrees two through six. Quadratic maps are omitted for point
cells, leaving 4,062 registrations. The independent analytic map is

$$
x_j(\widehat x)=
\begin{cases}
\widehat x_j+\alpha\widehat x_j^2,&j<d,\\
(j+1)/4,&d\le j<s.
\end{cases}
$$

Its Jacobian is diagonal in the first $d$ rows, with entries
$1+2\alpha\widehat x_j>0$. The left inverse has reciprocal diagonal entries;
the measure is their original diagonal product. Geometry coefficients are
obtained through the finite-element DOF functionals. Every sampled result is
checked against these expressions, with a declared dimensionless roundoff
budget and deliberate output-perturbation rejection. Extra embedding coordinates
are constant: this does not certify general curved embedded surfaces or
ill-conditioned inverse maps. Those remain separate profiling workloads.
Both map policies use `ParametricTransformation` with H1 geometry; the affine
label describes the represented map, not the default affine transformation
implementation. Holding that map fixed while raising its representation degree
isolates degree-dependent evaluation overhead.

An optimized native study at quadratic map degree six and order eight measured
about 826 ns per pyramid sample and 2,008 ns per hexahedron sample for warmed
mapped-quadrature construction. Their basis counts are 140 and 343: the costs
normalize to 5.90 and 5.86 ns per basis/sample, respectively. An Instruments
profile attributed about 91% of selected source CPU sample weight to the bulk
Jacobian contraction and reference-gradient access. Initial table-building
samples were also captured outside benchmark timing and are not warm-path
modal recomputation. In contrast, rebound pointwise pyramid Jacobians measured
about 1.1 ms per sample at degree six. These are different operation scopes,
not interchangeable assembly timings or a newly achieved speedup. The broad
five-repetition study contained three noisy selections; isolated process repeats
reduced their coefficients of variation to approximately 0.8–1.6%. The retained
raw results distinguish elapsed time from CPU time rather than silently dropping
the original samples.

```sh
cmake --build build --target RodinMappedPointBenchmarks -j1
build/tests/benchmarks/RodinMappedPointBenchmarks \
  --benchmark_filter='MappedPoint/(Construction|Jacobian/Rebound)/Quadratic/P2/Pyramid/3/8/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

### Arithmetic-preservation gate for pyramid reuse

Pyramid nodal functions contract modal values with inverse Vandermonde columns.
For basis index $a$ and reference point $r$, the computed expression is

$$
\phi_a(r)=\sum_m (V^{-1})_{ma}\psi_m(r).
$$

Materializing each modal value once per point avoids repeated modal evaluation,
but preserving the written summation order alone does not guarantee identical
floating-point results. An isolated prototype passed polynomial reproduction
while failing exact direct-versus-table comparisons in eight of twelve
real/complex degree combinations under the optimized native build. Disabling
floating-point contraction in a diagnostic build removed those discrepancies.
Neither the prototype nor a floating-point flag change was retained; no
performance gain is claimed for that rejected implementation.

`RodinVariationalH1TabulationTest` now compares every cached value and first
derivative with direct evaluation for degrees one through six, both scalar
fields and all eight reference families. Samples include vertices, centroids
and order-eight interiors, including the existing pyramid-apex convention.
Separate tests reconstruct a degree-dependent polynomial through public DOF
functionals at orders two, four and eight, with a dimensionless
$10^{-8}\max(1,|q_{\mathrm{ref}}|)$ roundoff allowance. These are numerical
regressions, not timing gates. All 24 tests passed against the original source;
the retained exact-comparison test failed against the isolated prototype.
Future reuse changes must pass this gate without weakening arithmetic policy.

An integer-only lookup-table prototype was also evaluated and rejected. At
degrees two, four and six, isolated decoding sweeps improved from approximately
14.7, 70.7 and 201 ns to 11.8, 35.8 and 89.0 ns. Modal sweeps did not show a
consistent improvement. In alternating original/prototype process runs, real
scalar degree-two pyramid point-evaluation sweeps instead increased from
21.4–21.5 ms to 22.1–22.3 ms; integration-point sweeps increased from
22.1–22.2 ms to about 22.7 ms. Each sweep evaluated 18,048 points with 14 scalar
basis functions. Thus the isolated decoding result did not transfer to its
consumer and does not justify adding the table. These are native optimized
five-repetition measurements, not a backend-wide conclusion.

The original and prototype basis/derivative snapshots were byte-identical
across both scalar fields, degrees one through six and all eight reference
families. `RodinVariationalH1PyramidIndexTest` separately checks layer enumeration,
inverse indexing and invalid-index compatibility through degree ten using exact
logical comparisons. Its sensitivity was confirmed with an isolated perturbed
decoder. No lookup table or source arithmetic change was retained.

### Reference-quadrature protocol

`RodinReferenceQuadratureBenchmarks`, owned by `QF`, separates owned rule
construction/destruction, polymorphic copying/destruction, traversal, and
canonical access. All eight reference families are exercised at orders two,
eight and sixteen; triangle and wedge additionally use order 51, and pyramid
uses order 11. Thus the first conical-product fallback is explicit for each
tabulated factor family. Order is a degree of exactness, not a point count.

For a rule $Q_p(f)=\sum_q w_qf(x_q)$, the numerical guards compare all monomial
moments of total degree at most two with closed-form integrals from the existing
quadrature verification helper. The dimensionless acceptance budget is
$10^{-12}\max(1,|I(f)|)$. Positive finite weights, finite coordinates and
closed-cell membership are checked; deliberately perturbing a weight must fail.
Closed-cell membership is intentional: the published degree-six pyramid table
contains boundary nodes. No quadrature coefficients or moment tolerances are
changed. These low-degree guards do not replace the full exactness unit suite.
The timed traversal integrates

$$
f(x)=1+\sum_{j=1}^{d}(x_j+x_j^2).
$$

Owned construction is repeated after library initialization. It is not a
process-first allocation experiment. Traversal includes a compiler memory
barrier to prevent hoisting an invariant result out of the timing loop.

The current canonical dispatcher has an eight-entry round-robin thread-local
cache backed by a locked process-wide pool. Untimed disjoint point-rule keys
first remove prior hot entries. An eight-key measured working set then visits
all cache positions once per sweep, so insertion rotation cannot masquerade as
a geometry cost. A nine-key cyclic working set instead misses the hot cache on
every access while retaining canonical pool objects. Pointer identity is
checked after eviction; no private cache-reset API is added. Cache-capacity
changes require revisiting this explicitly source-dependent protocol.
Both scopes include size access and report requests and points in the working
set; their per-sweep timings are not directly comparable per-request costs.

All 135 selections passed the guards in the native optimized build and the
LLVM OpenMP-enabled build. These scopes execute locally; enabling OpenMP does
not make this evidence of concurrent lookup throughput or MPI scaling.

In a native five-repetition study, order-sixteen owned construction measured
approximately 0.187 microseconds on triangles (53 points), 6.58 microseconds on
pyramids (729 points) and 11.14 microseconds on hexahedra (729 points). Traversal
normalized to about 3.13, 4.31 and 4.33 ns per point, respectively. Construction
has different algorithms and cannot be explained by dimension or point count
alone. Instruments attributed about 81% of selected hexahedron source sample
weight to `TensorProduct::append` and repeated point-vector/size access. The
pyramid profile instead includes Gaussian root generation and the Jacobi
eigensolve; the order-51 triangle fallback is dominated by those factor stages.
These profiles motivate a separate factor-access investigation, not an already
achieved source speedup.

The hot-cache profile contains thread-local access and key comparisons; the
nine-key pool profile additionally contains flat-map lookup and mutex operations
(about 33% inclusive sample weight). This confirms the intended experimental
paths, not lock contention: only one caller executes. The focused study retained
one noisy copy selection with an elapsed-time outlier; an isolated repeat had
about 1% variation. Raw results are retained separately from the profile captures.

The resulting `TensorProduct::append` repair reads factor counts once, borrows
left coordinates/weight once per outer iteration, and borrows each right point
and its dimension once per product entry. Ownership, coordinate ordering and
weight-product association are unchanged. Empty products read no samples;
aliased factors remain valid. Deterministic query-count regressions failed on
the original implementation and pass after the repair, alongside exact output,
copy-lifetime and three-factor ordering checks.

Matched native five-repetition owned-construction measurements were:

| Rule | Points | Before, microseconds | After repeat, microseconds | Speedup |
| --- | ---: | ---: | ---: | ---: |
| Quadrilateral, order 16 | 81 | 1.659 | 1.231 | 1.35 |
| Hexahedron, order 16 | 729 | 11.037 | 6.298 | 1.75 |
| Wedge, order 51 | 17,576 | 254.260 | 159.430 | 1.59 |

Triangle and pyramid controls remained essentially unchanged. The original
and repaired rule snapshots were byte-identical for all eight families through
order 20 and the triangle/wedge order-51 fallbacks. Both builds passed all 135
benchmark guards and 14 quadrature exactness tests with unchanged counts and
moment diagnostics; the three new unit tests passed in both builds. Fourteen
native Poisson p-convergence tests passed over all seven mesh geometries.
No solver or whole-assembly speedup is inferred from these construction results.

A matched post-repair Instruments capture reduced `append` plus point-vector
access from approximately 68% to 26% of selected source sample weight; repeated
size access ceased to be a leading frame. Point movement/copying and storage
initialization are now exposed as follow-up source hotspots. Profile fractions
describe attribution within each capture, not absolute elapsed-time ratios.

A subsequent placement repair constructs each zero-initialized point directly
in its reserved owning slot rather than filling and moving a temporary. Relative
to the factor-access repair above, matched five-repetition measurements reduced
order-16 hexahedron construction from 6.305 to 3.361 microseconds (1.88 times)
and order-16 wedge construction from 4.274 to 2.514 microseconds (1.70 times).
The order-51 wedge fallback improved from 159.634 to 97.401 microseconds;
triangle and pyramid controls remained essentially unchanged. Repeat variation
coefficients were below 1.7%. These are incremental construction results, not
whole-workflow speedups.

The 1,423,992-byte rule snapshot remained identical. All 135 benchmark guards,
four tensor-product unit tests and 14 quadrature exactness tests passed in both
optimized builds; 14 native Poisson p-convergence tests also passed. The added
inactive-capacity regression passes before and after placement and rejects an
isolated padding corruption. Zero initialization, including unused capacity,
is retained; the change removes movement, not initialization or ownership.
The matched hexahedron Instruments capture no longer has memory-copy or
fixed-matrix move frames among its leading samples. Gaussian factor generation
now accounts for about 15% of selected source sample weight, with `append`
itself about 39%; these fractions identify remaining work, not elapsed gains.

```sh
cmake --build build --target RodinReferenceQuadratureBenchmarks -j1
build/tests/benchmarks/RodinReferenceQuadratureBenchmarks \
  --benchmark_filter='ReferenceQuadrature/(OwnedBuild|OwnedCopy|Traversal)/Pyramid/16/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

### Spatial-vector storage protocol

`RodinSpatialVectorStorageBenchmarks`, owned by `Math`, isolates zero
construction, default-construction followed by logical resize, copying into a
reserved array, copy assignment, temporary-fill-then-move and direct fill in the
owning slot. Allocation and input preparation are excluded. Append scopes clear
the preceding batch; assignment retains the destination array size. The operation
is fixed at compile time so per-vector runtime dispatch is not measured.
Input pointers and memory barriers prevent invariant-result hoisting.

Real and complex scalars, dimensions zero through three, and batch sizes 32,
256 and 2,048 give 144 selections. Every logical size and every fixed-capacity
component is checked exactly, including inactive zeros; deliberately corrupting
an active or inactive component must fail. The storage count is
$N\operatorname{sizeof}(\mathrm{SpatialVector})$, whereas the active entry count
is $Nd$ and capacity is $3N$. Thus a zero-dimensional vector still has storage
and initialization cost; it does not represent a missing physical component.

All selections passed in both optimized builds. In a native five-repetition
256-vector study, real zero construction remained approximately 0.74 ns per
vector across logical dimensions. Real three-dimensional temporary-fill/move
and direct-fill scopes measured about 3.66 and 2.33 ns per vector; complex
three-dimensional scopes measured about 4.61 and 2.53 ns. These are complete
reserved-array operation scopes, not isolated move-instruction latencies.
The largest variation coefficient in this focused study was about 2.7%.

Instruments profiles at 2,048 real three-dimensional vectors attributed about
73% of temporary-fill source sample weight to memory-copy frames. Direct fill
instead concentrated sampled work in the inlined fill/storage loop. This
supports investigating point placement in the tensor-product consumer, but
does not by itself establish a quadrature or assembly speedup. The zeroed
inactive-storage contract must remain intact in any subsequent change.

```sh
cmake --build build --target RodinSpatialVectorStorageBenchmarks -j1
build/tests/benchmarks/RodinSpatialVectorStorageBenchmarks \
  --benchmark_filter='SpatialVectorStorage/(TemporaryFillAppend|DirectFillAppend)/Real/3/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

### Fixed-system scalar solver protocol

`RodinConstrainedSolveBenchmarks`, owned by `Solver`, retains manufactured
Poisson, variable-conductivity and scalar reaction–diffusion systems. Real and
complex H1 degrees one through three, all seven mesh geometries, and grids with
3, 5 and 9 points per axis give 378 selectable measurements. Complex fields use
the same positive real coefficients; these are Hermitian positive-definite CG
problems, not indefinite Helmholtz or mixed Stokes systems.

Assembly and its original-kernel reference are untimed. Each measured public
solve includes preconditioner setup and CG iterations, but excludes zeroing the
initial guess, field writeback and verification. Every repetition starts from
zero; a warmup fixes the expected iteration count. The retained diagnostics are

$$
r=\frac{\lVert Ax-b\rVert_2}{\max(1,\lVert b\rVert_2)},
\qquad
e=\frac{\lVert x-x_\star\rVert_2}{\max(1,\lVert x_\star\rVert_2)}.
$$

For complex coefficients these norms include both components. The requested
CG tolerance is $10^{-13}$; acceptance budgets are $r<10^{-11}$ and
$e<10^{-9}$, with unchanged operator/load comparisons against the reference.
A deliberately wrong field must have a detectable residual. Timings carry
cells, DOFs, matrix rows, nonzeros and reported iterations, so a dimensional
comparison can distinguish larger systems from more expensive iterations.
These algebraic guards do not replace continuous convergence certification.
No CTest registration or CI timing gate is added.

All 378 selections passed in the native and two-thread OpenMP-enabled optimized
builds. Maximum residuals were below $10^{-13}$ and maximum coefficient error
below $2.0\times10^{-12}$. Iteration counts are required to remain fixed within
each measured case, not across compilers or assembly reductions: two of the
378 cross-build comparisons differed by one iteration while satisfying the
same numerical budgets. This qualification is not an OpenMP solver speedup
claim or coverage of PETSc/MPI, indefinite systems or higher degrees.

A native five-repetition Poisson degree-two comparison gave the following
median public-solve times. The maximum timing variation coefficient was 3.2%.

| Geometry | Points per axis | DOFs | Nonzeros | Reported CG iterations | Real, microseconds | Complex, microseconds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Triangle | 3 | 25 | 61 | 4 | 1.87 | 2.84 |
| Triangle | 5 | 81 | 437 | 22 | 9.82 | 23.52 |
| Triangle | 9 | 289 | 2,293 | 61 | 105.30 | 249.52 |
| Tetrahedron | 3 | 125 | 369 | 10 | 6.39 | 15.13 |
| Tetrahedron | 5 | 729 | 6,721 | 40 | 230.03 | 502.21 |
| Tetrahedron | 9 | 4,913 | 80,721 | 82 | 5,060.07 | 9,575.76 |

At the largest level, the real solve ratio is 48.1, while the sparse-work proxy
$\operatorname{nnz}(A)\,k$ gives a ratio of 47.3. The proxy is not an exact
operation count: reported iterations can omit the converging update, and setup
and vector operations are also timed. The complex ratio is 38.4, illustrating
that constants and storage effects differ. Equal points per axis do not mean
equal algebraic problem size; these samples do not establish a pathological
dimension-specific solver slowdown or MPI imbalance.

Four matched Instruments captures at the largest level place sparse-product
source frames on approximately 71%/83% of selected CG sample weight for real
triangle/tetrahedron solves and 69%/81% for complex solves. Complex leaf samples
often resolve to standard-library multiplication/conversion helpers nested in
those products; leaf names alone must not be interpreted as a separate stage.
These inclusive fractions are source attribution, not stage timings or a proof
that setup has zero cost. Sparse products, reductions and preconditioner setup
remain candidates for isolated Math/Solver measurements. No solver source or
floating-point arithmetic policy was changed in this batch.

# Spatial value benchmarks

The spatial arithmetic benchmarks are part of the existing `RodinBenchmarks`
executable and benchmark publication workflow. Build an optimized configuration
with `RODIN_BUILD_BENCHMARKS=ON`, then run:

```sh
build/tests/benchmarks/RodinBenchmarks \
  --benchmark_filter='^Spatial(Vector|Matrix|Tensor[34])/' \
  --benchmark_min_time=0.1s --benchmark_repetitions=3 \
  --benchmark_out=spatial-algebra.json --benchmark_out_format=json
```

Each result measures one construction, copy, addition, scalar multiplication,
inner product, or contraction. Operand initialization, counters and independent
reference checks run outside the timed loop. Compiler barriers expose input
objects and materialize outputs on each iteration. Timings include these
barriers, so very small operations should be compared under the same compiler,
optimization flags and hardware. `entries` reports the active input entries;
`object_bytes` reports the complete object size, including unused capacity and
metadata. Throughput counts active entries (multiply-accumulate terms for
`MatMat`), not floating-point operations.

Both real and complex values are registered. Vectors cover sizes zero through
three. Matrices cover the empty shape and all nine nonempty shapes, with all
27 compatible matrix-product dimension combinations. Rank-three and rank-four
tensors cover the empty shape, uniform axis extents one through three, and two
rectangular shapes each. These tensor cases are representative, not an exhaustive
axis-extent sweep or coverage of every possible rank.

`Contract` measures matrix-vector, rank-three tensor-vector, or rank-four
tensor-matrix multiplication. Products use ordinary linear contractions, checked
against independent sums before timing. `Dot` measures each type's existing
method: complex vector `dot` is bilinear, while complex matrix and tensor `dot`
conjugate the second operand. Their complex timings therefore measure different
arithmetic contracts.

`SpatialTensor` owns a fixed-capacity `std::array` of entries. Its capacity is
three to the power of its compile-time rank under the current configuration;
changing runtime extents does not allocate storage. A local tensor stores that
array on the stack; a tensor embedded in another object shares that object's
storage duration. Rank-three and rank-four real tensors reserve 27 and 81 scalar
entries, respectively. The object size counter also includes their extents and
active-entry count.

The separate `RodinSpatialMatrixBenchmarks` executable measures finite-element
assembly against nine scalar assemblies. It uses total process CPU time, including
OpenMP workers, and separates complete assembly, triplet generation, element
binding, cached entry scans and sparse conversion. Those assembly costs should
not be inferred from these value-operation timings.

# Assembly performance methodology

The assembly benchmarks measure computational cost separately from the
multi-resolution accuracy assertions in `tests/convergence`. The performance
track complements convergence certification: timings do not replace accuracy,
rate, or backend-equivalence checks. The initial suite introduces no assembly
optimization; later comparisons must retain the same mathematical workload.

These standalone assembly executables are opt-in local profiling tools, not
CTest tests or CI gates. Build them explicitly with
`RODIN_BUILD_BENCHMARKS=ON` and the desired target. Their internal checks guard
local measurements; they do not contribute to convergence certification.
The historical `RodinBenchmarks` publication workflow remains separate and
does not build these standalone assembly targets.

## Implemented workload

`RodinPhysicsAssemblyBenchmarks` measures warmed real `H1<1>` and `H1<2>`
bilinear-form assembly on unit-box affine meshes. Each context covers segment,
triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge with
`n=3,5,9` grid points per axis, $h=1/(n-1)$: 126 registered cases.

| Context | Timed form $a(u,v)$ | Affine oracle field | Expected $a(u,u)$ |
| --- | --- | --- | --- |
| Poisson | $\int_\Omega\nabla u\cdot\nabla v\,\mathrm{d}x$ | $u=\sum_jx_j$ | $d$ |
| Conductivity | $\int_\Omega\gamma\nabla u\cdot\nabla v\,\mathrm{d}x$, $\gamma=1+\sum_jx_j$ | $u=\sum_jx_j$ | $d(1+d/2)$ |
| Linear elasticity | $\int_\Omega\lambda\operatorname{div}u\operatorname{div}v+2\mu\varepsilon(u):\varepsilon(v)\,\mathrm{d}x$ | $u_i=(i+1)\sum_jx_j$ | $\lambda C^2+\mu(dS+C^2)$ |

Here $d=\dim\Omega$, $\lambda=1.5$, $\mu=0.5$,
$\varepsilon(u)=(Du+Du^T)/2$, $C=d(d+1)/2$, and
$S=d(d+1)(2d+1)/6$. Component indices are zero-based. Quadrature order
6 is held fixed and recorded; this is not the order-12 convergence workload.
Affine fields are interpolated through the space's DOF functionals rather
than assuming nodal coefficient layouts.

Before timing, the matrix $A$ and an original-loop reference matrix
$A_{\mathrm{baseline}}$ are assembled through the same backend. The relative energy error
$|z^TAz/a(u,u)-1|$ must be below $10^{-9}$, where $z$ represents the oracle
field. Operator comparison and subsequent repeat assembly must satisfy
$\|A-A_{\mathrm{baseline}}\|_F<10^{-12}\max(1,\|A_{\mathrm{baseline}}\|_F)$.
Both checks are repeated after the timed loop. A failed check marks the
benchmark as erroneous and makes the executable return nonzero, including in
optimized builds. These checks detect incorrect energy and accumulation;
they are not full operator equivalence proofs for future optimizations.

The timed public `BilinearForm::assemble` operation includes traversal,
kernel evaluation, scatter, allocation/replacement as performed by that
backend, and sparse-matrix completion. Mesh/space/form construction,
interpolation, warmup, matrix comparisons, and energy evaluation are excluded.
No boundary constraints, loads, solve, or error integration are included.
This is repeated complete operator assembly, not isolated kernel timing or
cold-cache construction. Existing integrator benchmarks remain separate.

Eigen sequential or OpenMP assembly is selected by `RODIN_USE_OPENMP` and
reported in the output context. Wall time is measured with `UseRealTime`;
cell count, global DOF count, nonzeros, components, degree, quadrature order,
and cell throughput are emitted as counters. Timing comparisons must retain
the same geometry, degree, quadrature, build, and thread configuration.

Configure with `RODIN_BUILD_BENCHMARKS=ON` and an optimized build. For example:

```sh
cmake --build build --target RodinPhysicsAssemblyBenchmarks -j4
OMP_NUM_THREADS=1 build/tests/benchmarks/RodinPhysicsAssemblyBenchmarks \
  --benchmark_filter='/3/real_time$' --benchmark_min_time=0.01s
OMP_NUM_THREADS=4 build/tests/benchmarks/RodinPhysicsAssemblyBenchmarks \
  --benchmark_filter='LinearElasticity/P2/Pyramid' \
  --benchmark_repetitions=5 --benchmark_out=/tmp/assembly.json \
  --benchmark_out_format=json
```

The first command is a small-mesh numerical smoke run, not a performance
baseline. Run performance studies in isolation, without concurrent CTest or
build jobs. Record CPU, compiler, optimization flags, dependencies, sanitizer
status, affinity, thread counts, and repetitions. Report dispersion; vary
registration order to detect cache-order artefacts. Sanitized timings are not
interchangeable with optimized unsanitized measurements. No timing threshold
is imposed on general-purpose CI runners.

## Sequential stage measurements

`RodinPhysicsStageBenchmarks` registers four overlapping scopes for the same
three contexts, P1/P2 spaces, seven geometries, and `n=3,5,9`: 504 cases.
The stage policy is explicitly sequential even in an OpenMP-enabled build.
The canonical full operator supplies a comparison baseline; independently
contracting local entries with the affine coefficients must recover the
analytic energy before and after timing. The matrix generated from triplets
must satisfy the same energy and Frobenius-norm checks as full assembly.

| Scope | Timed operations | Excluded operations |
| --- | --- | --- |
| Binding | Cell traversal and each local integrator's `setPolytope` | Entry extraction, global scatter, sparse completion |
| Kernel | Binding plus extraction and checksum of every local entry | Global scatter and sparse completion |
| Triplets | Library sequential triplet assembly: binding, entry extraction, nonzero filtering, scatter and vector growth | Sparse completion |
| Finalize | `setFromTriplets` into a fixed-size matrix, including duplicate reduction and replacement | Triplet generation and local kernels |

Binding is not merely cache preparation: the current quadrature rules assemble
their local matrix inside `setPolytope`. Consequently, a costly Binding scope
can identify local integration work, while `integrate(m,l)` reads an already
computed entry. With $T_B,T_K,T_T,T_F$ denoting the measured scopes, no identity
$T_A=T_B+T_K+T_T+T_F$ is asserted. In particular, $T_B$ is contained in both
$T_K$ and $T_T$. Differences between independently timed scopes are diagnostic
estimates, not exact stage accounting. Kernel checksum work is not performed
by the production triplet assembler.

For the generic quadrature path with $N_Q$ points and $n_r,n_t$ trial/test
basis functions, the original test-side expression evaluation occurred
$N_Q n_r n_t$ times per cell. Point-local owning tabulation now reduces this
count to $N_Q n_t$; the $N_Q n_r n_t$ dot products remain. Values are refreshed
after every `setIntegrationPoint`, not reused across points, cells, or changed
coefficients. Scratch capacity is retained per integrator; OpenMP integrator
copies own independent scratch. Existing H1 tabulations/mapped-gradient caches
and specialized quadrature dispatch are unchanged. Scalar–spatial-matrix
products now promote real/complex scalar types in both operand orders.

Whole-process sampling traces also include framework startup, CPU-frequency
estimation, and untimed numerical checks. Attribute samples to the assembly
call path before reporting stage percentages. In particular, Google
Benchmark's frequency-estimation PRNG is not an assembly workload. Sampling
weights and unprofiled benchmark wall times are different measurements.

Counters report cells, DOFs, nonzeros, local-entry count summed across terms,
stored triplets, components, degree, quadrature order, and points per cell.
All stages are warmed; repeated Triplets retains vector capacity. Initial allocation,
mesh/space construction, constraints, PETSc insertion, and MPI communication
remain outside this executable's coverage. The 168 smallest-mesh stage cases
are available for opt-in local runs, without timing thresholds.

```sh
cmake --build build --target RodinPhysicsStageBenchmarks -j4
build/tests/benchmarks/RodinPhysicsStageBenchmarks \
  --benchmark_filter='LinearElasticity/P2/Hexahedron/5/real_time$' \
  --benchmark_repetitions=5 --benchmark_out=/tmp/stages.json \
  --benchmark_out_format=json
```

## PETSc local and MPI workload

`RodinPETScPhysicsAssemblyBenchmarks` uses the same `PhysicsForm` definitions,
affine energy oracles, degrees, geometries, sizes, and quadrature order as the
Eigen executable. The target requires `RODIN_USE_PETSC=ON` and
`RODIN_USE_MPI=ON`, including for the local-context mode. Local assembly uses
PETSc sequential/OpenMP according to `RODIN_USE_OPENMP`; `--assembly_mpi`
selects distributed assembly. This batch measures real fields with real-scalar
PETSc; complex-scalar PETSc builds and complex fields are not claimed as
verified coverage.

Each case uses three fixed iterations to keep all ranks on the same collective
sequence. With rank elapsed times $t_r$, the manual wall-time sample is
$T_A=\max_r t_r$. A barrier precedes each sample; neither that barrier nor the
timing reduction is included in $T_A$. Communication and matrix finalization
performed by `assemble` are included. Google Benchmark reports CPU time from
its own timer, which is not this maximum-rank assembly metric. Only manual
wall time and counters derived from it are used for distributed comparisons.

Workload metadata includes global owned-cell count $N_K$, global DOFs and
nonzeros, rank count $R$, minimum/maximum owned cells, and cell imbalance
$I_K=R\max_r N_{K,r}/N_K$. Ghost cells are excluded. The matrix difference
is measured in a global Frobenius norm and the energy by the distributed form
action; neither compares rank-local coefficient layouts. Empty owned-cell
partitions are included in the smallest segment cases at three/four ranks.

Only rank zero writes benchmark reports/JSON; other ranks still execute every
case. Failures may emit rank-local diagnostics. Random
benchmark interleaving and user-specified adaptive warmup are rejected to
prevent rank-dependent collective order. Local mode requires one process.
An empty benchmark selection returns nonzero rather than passing a smoke job.

```sh
build/tests/benchmarks/RodinPETScPhysicsAssemblyBenchmarks \
  --benchmark_filter='/3/iterations:3/manual_time$'
mpiexec -n 4 build/tests/benchmarks/RodinPETScPhysicsAssemblyBenchmarks \
  --assembly_mpi --benchmark_filter='LinearElasticity/P2' \
  --benchmark_repetitions=5 --benchmark_out=/tmp/assembly-mpi.json \
  --benchmark_out_format=json
```

Opt-in local runs support local and MPI (1–4 ranks) smallest-mesh cases in
sequential/OpenMP builds and JSON output without timing thresholds. Full
three-size registrations support fixed-global-work comparisons across ranks;
no weak-scaling sequence or speedup claim follows merely from registration.

Full Eigen/PETSc operators are compared against a test-only `ReferenceIntegral`
which retains the original quadrature/trial/test evaluation loop and passes
through the same assembly backend. Thus the comparison includes backend
scatter, constraints present in the workload, and MPI ownership/reduction;
the current workload still has no boundary constraints. The original loop is
not a second production dispatch path.

## Paired generic-expression benchmarks and regression coverage

`RodinQuadratureExpressionBenchmarks` measures local binding only, with paired
`Original` and `Reuse` registrations. Both retain expression and matrix storage
across binds, own mapped quadrature only for the bound cell, and exclude setup and
entry comparisons from timing. The original loop lives in the shared test-only
`tests/QuadratureReference.h`. Every local entry must agree exactly before and
after timing; a mismatch returns nonzero. Comparisons are separate measurements,
not simultaneous timers. Registration order, frequency variation and profiling
overhead must still be controlled before claiming speedup.

The matrix includes P0, P0g, P1 and H1 orders 1–3, real/complex scalar/vector
fields, all seven geometries, and `n=2,3,5`: 1008 registrations. Vector P1/H1
uses the compound symmetric-gradient form; scalar and constant spaces use
compound mass forms. Scalar constant cases have $n_r=1$ and are explicit
overhead controls rather than expected evaluation-count improvements. These
are sequential kernel measurements even in an OpenMP build. The 336
smallest-mesh registrations are available for opt-in local measurements.

The unit regression suite compares local entries against the original loop,
assembled operators through logical DOF maps, and exact integer callback
counts. It covers changing points/cells/orders/coefficients, copies and moves,
mixed trial/test sizes, complex phases/conjugation, a curved triangle and facet
integrators. The scalar-promotion regression checks both operand orders and
all spatial matrix dimensions from 0 through 3, including inactive-storage
zeros. Explicit `Conjugate` nodes are tested on scalar expressions; existing
explicit vector `Conjugate` support is outside this change. Complex vector
forms test the implicit test-side conjugation in `Math::dot` with non-real
phases. Real-scalar PETSc does not certify complex PETSc assembly.

An isolated Apple M1 Pro measurement with Apple Clang 21, Release `-O3`,
no sanitizers, sequential Eigen, order-6 quadrature and 64 hexahedra (`n=5`)
gave the following paired binding medians over five repetitions:

| H1 P2 expression | Original loop (ms) | Owning reuse (ms) | Ratio | CV, original/reuse |
| --- | ---: | ---: | ---: | ---: |
| Real symmetric gradient | 846.9 | 190.3 | 4.45 | 3.54% / 0.41% |
| Complex symmetric gradient | 1135.9 | 260.9 | 4.35 | 1.47% / 0.74% |

The complex reference uses the repaired scalar promotion; it is not a timing
of the formerly uncompilable expression. Scalar H1 P2 compound mass decreased
from 0.9673 to 0.7929 ms. Scalar P0 and P0g controls increased by approximately
1.6% and 1.3%, respectively (0.06568 to 0.06675 ms and 0.01615 to 0.01637 ms).
These small costs are reported rather than treated as evaluation-count gains:
both constant scalar cases have only one local trial function. The measurements
are workload-specific observations, not CI thresholds or universal speedups.

These historical paired measurements predate the reference oracle's bounded
mapped-point ownership. Fresh comparisons must use the current ownership
contract on both sides; the historical numbers are not a baseline for that
lifecycle change.

With matched original-header and repaired-header builds of the full physics
benchmark, real P2 linear-elasticity assembly decreased from 886.8 to 208.7 ms
(five-repeat medians; CV 4.54% and 0.45%). The specialized Poisson control
measured 5.08–5.16 ms in the original-header build and 5.24–5.34 ms in the
repaired build across both execution orders. The observed 1.5–5.1% shift and
varying dispersion do not establish performance neutrality for that control;
its specialized quadrature dispatch and numerical results remain unchanged.

## Full mass/reaction operators and projection loads

`RodinMassAssemblyBenchmarks` extends local-kernel coverage to full Eigen
sequential/OpenMP operator and load assembly. The matrix contains P0, P0g, P1,
H1 P1/P2/P3, real/complex scalar/vector fields, seven affine geometries, and
`n=3,5,9` grid points per axis: 2016 registrations. Each context has separate
`Operator` and `Load` registrations; neither includes a projection solve.

For a scalar or $r$-component field, define
$a_\rho(u,v)=\int_\Omega\rho\sum_{j=1}^r u_j\overline{v_j}\,\mathrm{d}x$
and $L(v)=a_\rho(c,v)$. Mass uses $\rho=1$; reaction uses
$\rho=s\alpha(1+\sum_{k=1}^d x_k)$, with $s=1$ during timing and
$\alpha=1$ over the reals or $\alpha=1+i/2$ over the complex numbers.
The constant field has components $c_j=j\eta$, using one-based component
indices and $\eta=1$ or $1+2i$, respectively. Scalars use $r=1$ and vectors
use $r=d$. These fields belong to every tested space, including P0g.

If $z$ is obtained by interpolation through the space's DOF functionals,
the assembled operator and load must satisfy $Mz=b$. On the unit box the
independent analytic action is
$z^*Mz=|\eta|^2\sum_{j=1}^r j^2$ for mass and
$z^*Mz=s\alpha(1+d/2)|\eta|^2\sum_{j=1}^r j^2$ for reaction.
Complex reaction is not claimed to be Hermitian. The unweighted mass equation
is the ordinary $L^2$ projection equation; its right-hand-side assembly is
measured here, not its solution time or approximation error.

Before and after timing, the analytic action has relative tolerance $10^{-10}$.
With $M_0,b_0$ the initial operator and load, the checks also require
$\|M-sM_0\|_F<10^{-11}\max(1,s\|M_0\|_F)$,
$\|Mz-b\|_2<10^{-11}\max(1,\|b\|_2)$, and
$\|b-sb_0\|_2<10^{-11}\max(1,\|b\|_2)$.
Here $M_0$ is assembled through the original quadrature loop, while $b_0$
is a snapshot of the initial load; $s=1$ throughout the mass workload.
Reaction coefficients are doubled after timing, and both outputs must equal
twice their initial values. This checks replacement and changing coefficient
state, not merely a second assembly with identical data. Checks, interpolation,
reference assembly, and coefficient updates are untimed. A failed check or a
filter matching no case returns nonzero.

Quadrature order 8 is fixed and recorded, together with the actual number of
points. Degree-three polynomial products with an affine coefficient have
degree at most seven on simplices; pyramidal bases and mapped tensor-product
spaces require their own exactness analysis. No blanket claim of exact mass
entries on every geometry follows from this order. The constant-field
identities provide an analytic check independent of the basis and geometry
enumeration, while the original loop checks the implemented quadrature matrix.

`RodinPETScMassAssemblyBenchmarks` uses the same mathematical workload with
PETSc local/MPI storage and sequential/OpenMP-enabled builds. It registers
1008 cases for the PETSc installation's scalar type, uses three aligned manual
iterations, and reports maximum-rank time with owned-cell imbalance metadata.
The installed real-scalar PETSc certifies only real workloads; a complex-scalar
PETSc build must be compiled and tested separately. P0g here uses per-cell
integrals assembled into shared global constant DOFs; it does not cover the
separate global-integrator lists or mixed multiplier couplings.

Smallest-mesh filters support opt-in local Eigen and native-scalar PETSc runs
locally and at MPI ranks one through four, in both OpenMP configurations.
These runs are not CI gates; registrations do not establish passing results
or three-size performance qualification. Larger-size runs should use
`--benchmark_min_time=1x` for numerical checks, or controlled repeated timings
for performance studies; H1 P3 vector reference matrices can be costly.

The smallest segment mesh has two cells, so three- and four-rank runs
explicitly exercise empty shards. A P0g shard can be empty while retaining
access to shared global constant DOFs whose unique owner is rank zero.
Distributed interpolation selects an eligible entity using logical distributed
indices and evaluates its finite-element functionals before updating storage.
Only DOF owners commit the resulting coefficients; ghost copies are refreshed
from their owners. For $P_{0g}$, the selected source can be on another rank,
including when the DOF owner's shard is empty. Empty selections preserve the
existing coefficients. Source selection and evaluation are implemented in the
MPI layer independently of PETSc storage. The mesh dimension is the maximum
over shard dimensions, while an empty shard retains local dimension zero.
Numerical smoke runs establish correctness checks, not isolated performance
baselines.

## Full constrained scalar problems

`RodinConstrainedAssemblyBenchmarks` measures complete Eigen problem assembly,
including volume loads, nonhomogeneous Dirichlet evaluation/elimination, and
matrix/vector completion. `RodinPETScConstrainedAssemblyBenchmarks` measures
the same forms in PETSc local/MPI contexts. Mesh/space/problem construction,
interpolation, reference assembly, checks, and solves are excluded. There is
no solve in this workload. The PETSc target requires MPI support even for
single-rank local mode. Eigen uses wall time; PETSc uses three aligned
iterations and maximum-rank wall time, excluding the preceding barrier and
timing reduction as in the existing distributed benchmarks.

The unit-box problem is

$$
-\nabla\cdot(\gamma\nabla u)+\rho u=f\quad\text{in }\Omega=(0,1)^d,
\qquad u=g\quad\text{on }\partial\Omega.
$$

Poisson uses $\gamma=s$, $\rho=0$; conductivity uses
$\gamma=s(1+\sum_jx_j)$, $\rho=0$; reaction--diffusion uses the latter
diffusion coefficient and $\rho=s$. The manufactured field is

$$
g=u=\eta t\left[1+\sum_{j=1}^{d}(x_j+\chi x_j^2)\right],
\qquad f=-\nabla\cdot(\gamma\nabla u)+\rho u,
$$

where $\chi=0$ at degree one and $\chi=1$ at degrees two/three;
$\eta=1$ for real fields and $\eta=1+2i$ for complex fields. Thus degree-one
Poisson has zero volume forcing, while conductivity and reaction--diffusion
exercise nonzero loads at every degree. The other Poisson cases have nonzero
forcing. Fields are interpolated through DOF functionals, not coefficient
layout assumptions. All seven positive-dimensional geometries, H1 degrees
one through three, and `n=3,5,9` are registered: 378 Eigen cases across both
scalar fields, and 189 PETSc cases for the installed scalar field. P0/P0g
are excluded because these conforming diffusion problems require H1 spaces;
point/0D physics and curved maps remain separate workplan items.

The production constrained operator is compared against the original-loop
quadrature reference, assembled through the same backend and constraints.
The load reference shares the production load integrator; its independent
check is the manufactured residual, not an independent entrywise load oracle.
With $z=I_hu$, each initial, repeated, timed-final and changed-state assembly
must satisfy

$$
\frac{\|Az-b\|_2}{\max(1,\|b\|_2)}<10^{-10},\qquad
\|A-A_{\mathrm{ref}}\|_F<10^{-11}\max(1,\|A_{\mathrm{ref}}\|_F),
$$

and $\|b-b_{\mathrm{ref}}\|_2<10^{-11}\max(1,\|b\|_2)$.
An intentionally incorrect coefficient vector $z+\mathbf{1}$ must give a
normalized residual greater than $10^{-4}$; this checks sensitivity of the
oracle, not a physically constant field in a nonnodal basis. After timing,
$s$ changes from $1$ to $1.5$ and $t$ from $1$ to $2$, exercising both
coefficient re-evaluation and changing essential data. Quadrature order is
eight at degree one and twelve otherwise; these choices are part of the
workload and are reported, not a blanket exactness claim on rational bases.
Opt-in local runs cover the smallest-mesh registrations, including PETSc local
and MPI at one through four ranks in both OpenMP configurations. Small segment
partitions include empty shards. An empty benchmark filter or failed oracle
returns nonzero. Performance claims require separate isolated repeated runs.

## Prescribed-state nonlinear Poisson assembly

`RodinNonlinearPoissonAssemblyBenchmarks` isolates warmed native residual and
tangent assembly at a prescribed field. The active native backend is sequential
or OpenMP according to the build configuration. Degrees one through three,
all seven positive-dimensional geometries, grid-point counts $n=3,5,9$, and
quadrature orders eight and sixteen give 252 registrations per configuration.
An additional 216 `CurvedQ2` registrations use degree-two geometry on all six
two- and three-dimensional geometries, with the same field degrees, sizes and
quadrature orders. The existing `CurvedGeometry` utility installs
$\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1}$ on cells and traces before spaces are built.
For $d\ge2$, $\det D\Phi=1$ and $x_0=\xi_0$. Consequently the physical affine
state, its gradient, and both analytic actions below are unchanged on
$\Omega=\Phi((0,1)^d)$. This comparison changes geometry evaluation while
holding the physical integrand and domain moments fixed. Segment shear is
excluded from this comparison because neither identity holds in dimension one.
The order-sixteen cases match the volume quadrature order of the nonlinear
natural-boundary convergence studies; these are not complete boundary-problem
or Newton-solve timings.

For $\Omega=(0,1)^d$, the interpolated state is $q_s(x)=s(1+x_0)$. The timed
forms are

$$
R(q_s;v)=\int_\Omega\left[\nabla q_s\cdot\nabla v+(q_s+q_s^3)v\right]\,dx,
\qquad
J(q_s;w,v)=\int_\Omega\left[\nabla w\cdot\nabla v+(1+3q_s^2)wv\right]\,dx.
$$

The affine field is interpolated through DOF functionals, without assumptions
about coefficient numbering. Before and after timing, the analytic actions

$$
R(q_s;q_s)=\frac{10}{3}s^2+\frac{31}{5}s^4,
\qquad
J(q_s;q_s,q_s)=\frac{10}{3}s^2+\frac{93}{5}s^4
$$

must agree to relative tolerance $10^{-9}$. The tangent is also compared with
the original test-only quadrature loop to relative Frobenius tolerance
$10^{-11}$, using denominator $\max(1,\|J_{\mathrm{ref}}\|_F)$.
An assembled central difference in the direction $w=1+x_0$, with step
$\varepsilon=10^{-5}$, must agree with $J(q_1)w$ to tolerance $10^{-7}$,
normalized by $\max(1,\|J(q_1)w\|_2)$. These are numerical oracle policies,
not performance thresholds. After timing, $s$ changes from one to two and
both forms are reassembled, checking state-update visibility and replacement
rather than accumulation. Omitting the cubic contribution would change the
analytic residual action by more than ten percent.

Setup, interpolation, independent checks and state changes are untimed. The
residual and tangent scopes exclude source terms, natural-boundary integrals,
constraints, SNES/KSP, and error norms. PETSc/MPI parity and separate timings for those excluded operations
remain required before attributing a full distributed solve's cost to this
benchmark. Counters report cells, global native DOFs, nonzeros, field and geometry degrees,
quadrature order and points per cell. Smallest-mesh cases are available for
opt-in local numerical checks, without timing thresholds.

### Physical-coordinate-dependent source

Separate `Source` registrations measure replacement assembly of

$$
L_s(v)=s\int_\Omega (1+x_0+x_{d-1}^2)v\,dx.
$$

These loads exercise physical-coordinate evaluation, which the prescribed-state
residual and tangent do not require. Field degrees one through three, all seven
positive-dimensional geometries, affine and exact quadratic geometry,
$n=3,5,9$, and quadrature orders eight and sixteen give 252 source registrations
per backend configuration. Setup and reference integration remain untimed.

With the quadratic amplitude $a=0.1$ (or $a=0$ for affine geometry), the
constant-test action in dimensions two and three is

$$
L_s(1)=s\left(\frac{11}{6}+\frac{a}{3}+\frac{a^2}{5}\right).
$$

The curved segment is included in this scope: its image is $(0,\ell)$ with
$\ell=1+a$, so

$$
L_s(1)=s\left(\ell+\frac{\ell^2}{2}+\frac{\ell^3}{3}\right).
$$

Every coefficient is compared with an independent pointwise mapping, Jacobian
and basis integration loop, using relative vector tolerance $10^{-11}$ with
denominator $\max(1,\|b_{\mathrm{ref}}\|_2)$. The analytic action uses relative
tolerance $10^{-9}$. Gates run before and after timing and after changing
$s=1$ to $s=2$. Omitting the quadratic source term changes the analytic action
by more than ten percent. These are accuracy gates, not timing thresholds.
The source workload is not a nonlinear solve or a distributed assembly benchmark.

The independent quadrature oracle owns mapped points only for its currently
bound cell, just as production assembly does. It does not request or evict
mesh-owned borrowed quadratures. This prevents untimed validation from retaining
all cells' quadrature points during timed assembly. Integer request-count
regressions enforce that lifecycle separately from numerical agreement.

Mapped-quadrature construction prepares the current cell's Jacobians in one
bulk transformation call. H1 geometry reuses its existing reference derivative
table, but mapped matrices remain owned by each geometric point. The same
assembly benchmark therefore measures this change without a separate assembly
path. Pointwise-versus-bulk exact comparisons, logical formula-lifetime checks,
and point copy/move/rebinding tests gate the geometry operation independently
of the assembled numerical oracles. Assembly introduces no MPI collective for
this operation and retains mapped points only for its currently bound cells;
the explicit mesh-owned quadrature cache remains a separate public facility.

```sh
build/tests/benchmarks/RodinNonlinearPoissonAssemblyBenchmarks \
  --benchmark_filter='^NonlinearPoisson/(Residual|Tangent)/P1/Tetrahedron/' \
  --benchmark_min_time=0.1s --benchmark_repetitions=5 \
  --benchmark_out=nonlinear-volume.json --benchmark_out_format=json
```

Paired optimization measurements must hold the compiler, assembly backend,
thread count, mesh, degree, quadrature and benchmark scope fixed. Independent
repetitions are required; a smallest-mesh numerical smoke run is not evidence
of a speedup or of complete convergence certification.

For the curved pyramid geometry-evaluation workload, use the same executable:

```sh
build/tests/benchmarks/RodinNonlinearPoissonAssemblyBenchmarks \
  --benchmark_filter='^NonlinearPoisson/CurvedQ2/(Residual|Tangent)/P1/Pyramid/3/16/real_time$' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5 \
  --benchmark_out=nonlinear-curved-pyramid.json --benchmark_out_format=json
```

The pyramid modal derivative evaluates only the Bernstein factors required by
the selected coordinate direction. This removes unused evaluations without
changing the returned expression, apex policy, or nodal transformation.
`PyramidDirectionalDerivativesPreserveArithmetic` compares every modal
derivative exactly with the unpruned expression for orders one through six,
at quadrature points, vertices and near-apex points. Independent analytic
gradient tests and the assembly oracles above remain separate checks.
Compare baseline and candidate binaries in separate processes, reversing their
execution order, and retain unchanged-geometry controls. A single selected
case per process fixes initial cache-population order; random interleaving of
multiple cases is a different experiment. Neither check supplies a portable
wall-time threshold or certifies the complete nonlinear convergence matrix.
If an unchanged-geometry control shifts, inspect both generated instructions
and helper placement before attributing the shift to its numerical kernel.
A diagnostic relink with matched helper placement can distinguish binary-layout
effects from added arithmetic; it is an experimental control, not a production
linker policy. Repeat optimization comparisons after changes to shared caches.

## Remaining performance workplan

| Extension | Required measurement or evidence |
| --- | --- |
| Extended workload/constraint parity | Full constrained scalar Poisson/conductivity/reaction--diffusion assembly is implemented; extend to vector physics, mixed boundary conditions and identification constraints |
| Complex Helmholtz, coupled reaction–diffusion, Taylor–Hood Stokes, nonlinear Poisson | Loads/block forms and residual/tangent assembly; nonlinear Poisson has prescribed-state native volume assembly through degree three, with PETSc/MPI and boundary parity remaining; stable mixed spaces |
| Projection solves, global couplings, and higher-order H1 physics | Mass/reaction and constrained scalar diffusion through degree three are implemented; extend to projection solves, mixed/global-integrator couplings, degree-four/vector physics, and meaningful point/0D forms |
| Remaining setup and stage isolation | Cold allocation, mesh/space setup, PETSc/MPI insertion and completion, constraints, solve, and norm timings; warmed sequential binding/kernel/triplet/finalization scopes are implemented |
| Curved geometry and boundary variants | Map degree/regularity, quadrature-point counts, constraints and boundary workload metadata |
| Scaling and regression baselines | At least three sizes; fixed-global-work strong scaling and fixed-per-rank-work weak scaling; controlled runner and established variance before thresholds |

All meaningful combinations of existing contexts, supported spaces, the seven
positive-dimensional geometries, and assembly backends remain the intended
coverage matrix. Unsupported combinations are stated explicitly, not counted
as completed measurements. Darcy remains deferred. Any subsequent optimization
requires independently checked operators/loads/constraints and baseline-
equivalent solution behavior, residuals, and solver iterations.
