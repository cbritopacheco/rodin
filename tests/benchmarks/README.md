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

Exported Instruments `time-profile` XML can be summarized with
`dev/summarize_xctrace.py`. Repeated stack selectors are combined by logical
AND, so an owning workload frame and the relevant source kernel can both be
required. A process selector distinguishes MPI participants or unrelated
processes. For example:

```sh
python3 dev/summarize_xctrace.py profile.xml \
  --process 'RodinHelmholtzSolveBenchmarks' \
  --require 'HelmholtzSolve<.*::run\(' \
  --require 'Eigen::internal::gmres<' \
  --category 'sparse-product=sparse_time_dense_product' \
  --category 'orthogonalization=Householder' \
  --limit 12
```

The report retains its process, stack and category selectors, weighted sample
counts, leaf attribution and inclusive categories. Inclusive categories may
overlap; their percentages are
not additive. Sampled CPU weights are not elapsed timings, and stack selection
is not a timed-region boundary: untimed warmup calls to the same kernel can
remain selected. Fixture kernels without the required solver frame are excluded.
For assembly captures, `--exclude 'ReferenceIntegral|QuadratureReference'` removes
the original-reference kernel stacks; this still is not timed-region bracketing.
The parser retains referenced weights even when first defined in a row without
a stack, rejects unresolved references and discards processed XML rows. Its
offline parser checks are not benchmark tests, timing gates or CI registrations.
Category patterns must identify executing methods, not merely a type appearing
in a demangled template argument. For example, `DiagonalPreconditioner` also
appears in GMRES expression signatures; its method pattern
`Eigen::DiagonalPreconditioner<(?:[^<>]|<[^<>]*>)*>::(?:factorize|_solve_impl)<`
distinguishes setup
and application from this incidental type occurrence.

### Verification-support attribution

The finite pressure-spectrum oracle is not a production Stokes solver. Its
global collection and result broadcast intentionally surround one rank-zero
spectral calculation. A complete native three-rank $P_4/P_3$ pyramid hierarchy
passed all affine/curved and order-16/order-20 comparisons in 1,276.57 seconds.
The observed algebraic dimensions and global-spectrum diagnostic intervals were:

| Points per axis | Free velocity DOFs | Pressure DOFs | Global-spectrum wall interval |
| --- | --- | --- | --- |
| 3 | 2,565 | 559 | 1.50–1.64 s |
| 4 | 9,177 | 1,729 | 27.85–28.54 s |
| 5 | 22,413 | 3,925 | 276.65–290.13 s |

One physical mean constraint reduces the pressure spectral dimensions to
558, 1,728 and 3,924. Fine-grid Schur-column construction took 81–91 seconds
and incremental QR took 95–99 seconds. These are stage wall observations from
an instrumented test, not controlled benchmark baselines. Two short Instruments
captures were attached sequentially to the existing run: the root capture
selected generalized eigenspectrum work, while the peer capture selected the
result broadcast. Neither capture covers the complete hierarchy.
The observed rank CPU imbalance therefore does not establish a partitioning
defect. These support profiles do not count as production assembly/solve
profiling; the opt-in `RodinPressureOracleBenchmarks` retains that distinction.

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
| Spatial values / Math | Real/complex vector, matrix and tensor arithmetic; single/batched spatial dot; fixed-capacity construction/copy/assignment and temporary/direct array fill | Remaining arithmetic consumer attribution and controlled scaling baselines beyond selected storage scopes |
| Sparse algebra / Math | Direct/transpose sparse products, Hermitian dots and squared norms; 1,512 selections qualified in both builds with selected source profiles | P0/P0g mass operators, mixed operators, remaining degrees and distributed representations |
| Reference metadata / Geometry | Topology, connectivity and reference transformation consumers | Isolated metadata access, embedding and shape-dependent work accounting |
| Reference elements / Variational | Scalar P0/P0g/P1 and H1 degrees 1–6; values, derivatives, prefetched descriptors | Vector/matrix component expansion, DOF functionals and first-use initialization |
| Reference quadrature / QF | Owned construction/copy, polynomial traversal, controlled eight-key hot and nine-key canonical-pool sweeps; table/fallback boundaries | Process-first allocation, concurrent pool contention and profiles beyond selected orders |
| Transformations / Geometry | Pointwise/bulk Jacobians; pointwise map consumer profiles | Independent geometry degree, curvature, embeddings and inverse-map stage matrix |
| Mapped points / Geometry | Isolated warm/rebound coordinates, Jacobian, inverse and measure; mapped-quadrature construction; field consumer profiles | General curved embeddings, cache-miss provenance, controlled timings and profiles beyond selected cases |
| Space setup / Variational, MPI | P1 construction and distributed construction | P0/P0g/H1 setup, DOF lookup/restriction, reconciliation and degree scaling |
| Fields / Variational | Interpolation and true mass projection; 1,008 native and OpenMP selections qualified; point/integration-point sweeps | P0/P0g, higher degrees, matrix evaluation, point cells, norms and boundary restriction |
| Expressions / Variational | Generic expressions, coefficients and named forms; native real/complex sparse/dense coupled Problem lifetimes at degree two | Cold expression construction, remaining degrees, changes of state and complete range/backend matrix |
| Local integrators / Variational | Volume entry kernels and reference-loop comparisons; mixed Stokes vector diffusion and both pressure/divergence orientations, with 42 selections qualified per build; separate weighted-viscosity binding scopes | Remaining mixed degrees/geometries/sizes, boundary/interface kernels and nonlinear residual/tangent stages |
| Operator setup / Assembly | Public repeated assembly and sequential stage scopes | Sparsity discovery, allocation and replacement separated from kernels |
| Scatter / Assembly, PETSc | Sequential triplets; complete native/PETSc insertion consumers | Isolated OpenMP buffers/merge, PETSc insertion and ownership-dependent MPI work |
| Constraints and completion / Assembly, PETSc | Constrained scalar assembly and fixed-system solves | Vector/mixed constraints, identification, finalization and communication attribution |
| Public physical assembly / Assembly | Diffusion, conductivity, LinearElasticity, mass/reaction and nonlinear Poisson; complex Helmholtz and real/complex coupled reaction–diffusion, each with 126 cases qualified in each build and selected native profiles; mixed Stokes with 20 selections qualified in native and two-thread OpenMP assembly, plus explicit weighted-viscosity scopes | Stokes mesh/backend expansion, boundary variants and remaining degrees/backends; constructor/lifecycle isolation |
| Solvers / Solver, PETSc | Real/complex fixed-system CG; isolated diagonal compute/apply; projection CG profile; complex Helmholtz GMRES and coupled reaction–diffusion CG, each with 126 cases qualified in each build and selected native profiles; native fixed-system Stokes SparseLU with 20 selections qualified | Direct factorization stages, remaining setup/lifetimes, mixed-system mesh/backend expansion, backend scaling and memory; repeatable 3D GMRES baselines |
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

### Mixed Stokes measurement protocol

`RodinStokesAssemblyBenchmarks` and `RodinStokesSolveBenchmarks` share
`Workloads/StokesWorkload.h`. For velocity degree $k\in\{2,3,4\}$, pressure
degree $k-1$ and one global $P_0^g$ multiplier, they retain the certified form

$$
(\nabla u_h,\nabla v_h)-(p_h,\nabla\cdot v_h)
+(\nabla\cdot u_h,q_h)+\lambda_h(1,q_h)+\mu_h(p_h,1)
=(f,v_h).
$$

On the unit box, $u=x_1e_0$, $p=x_0-1/2$ and $f=e_0$ supply independent
analytic data, with prescribed velocity traces and zero pressure mean.
Public field interpolation and Problem's block gathering provide the exact
coefficient vector; no interleaved DOF numbering is assumed. Qualification
checks the patch residual, solved coefficients and true residual, and rejects
both a wrong coefficient vector and the operator missing the momentum pressure
term. Those guards are outside timing and do not replace mixed-space convergence
or stability certification. The pressure-spectrum oracle is not timed here.

Each executable registers 108 opt-in selections: six two-/three-dimensional
cell families, three degree pairs, $n\in\{3,5,9\}$ points per axis and two
expression dispatches. Names containing `WeightedViscosity` explicitly multiply
the trial Jacobian by unit viscosity; the original names retain the unweighted
expression. Although the forms are mathematically identical, their specialized
local kernels differ and their costs must be measured independently.
Assembly includes public repeated assembly, constraints and completion. Solving
uses the fixed assembled operator and public SparseLU, including a fresh
factorization and one residual correction on every iteration. It excludes
assembly, field scatter and norm integration; no triangular-solve-only timing
is claimed. Cells, separate velocity/pressure/mean DOFs, nonzeros and quadrature
points accompany the measurements. Larger registered cases are not implicitly
qualified by a smaller run.

The initial native qualification covers all degree pairs and geometries at
$n=3$, plus the $P_2/P_1$ triangle at $n=5,9$: 40 distinct assembly/solve
selections. Five-repeat native $P_2/P_1$ medians at $n=3$ are:

| Geometry | Assembly | SparseLU, solve and correction |
| --- | --- | --- |
| Triangle | 0.179 ms | 0.0191 ms |
| Quadrilateral | 0.194 ms | 0.0164 ms |
| Tetrahedron | 27.54 ms | 0.201 ms |
| Pyramid | 39.78 ms | 1.654 ms |
| Hexahedron | 16.52 ms | 0.155 ms |
| Wedge | 12.07 ms | 0.191 ms |

These Apple LLVM 21 RelWithDebInfo baselines are not portable timing thresholds.
Assembly coefficients of variation range from 0.24% to 5.06%, and solving
from 0.27% to 5.19%. Matching triangular levels $n=3,5,9$ have assembly medians
0.179, 0.720 and 2.861 ms and solve medians 0.0191, 0.238 and 1.496 ms.
The tetrahedron/triangle assembly ratio must retain the work counts: 48 versus
8 cells and 729 versus 53 quadrature points per cell at order 16, as well as
different vector and pressure basis sizes. It is not an isolated dimension
penalty. Remaining qualifications include larger three-dimensional meshes,
OpenMP and distributed/PETSc scopes, direct-solver stage isolation and curved
or boundary variants. No production repair or speedup is claimed by this batch.

Matching Instruments captures select the workload and executing assembly or
SparseLU method together. Of selected assembly CPU samples, `setPolytope`
accounts for 88.39% on triangles and 97.73% on tetrahedra; this method includes
eager element-matrix construction, not merely cursor binding. Bulk geometric
Jacobians account for 17.14% and 17.00%, respectively, inside that work.
`setFromTriplets` accounts for 1.82% and 0.16%. Tetrahedral leaf samples identify
vector-Jacobian local construction (28.66%), bulk geometric Jacobians (13.68%)
and divergence/pressure local construction (12.43%) as the next isolation targets.
The corresponding solve captures attribute 54.68%/62.72% to numeric factorization,
23.73%/24.11% to symbolic analysis and 11.44%/7.72% to substitution, for
triangle/tetrahedron. These are inclusive CPU-sample categories, can overlap,
include matching untimed qualification calls and are not elapsed-time fractions.
The retained raw captures and JSON summaries are local certification artifacts;
only files prefixed `final-` represent the canonical Stokes expression. Earlier
unit-weighted pressure-expression experiments are not its baseline.

### Mixed Stokes local-kernel isolation

`RodinStokesBindingBenchmarks` belongs to Variational and isolates the three
volume operators consumed by Stokes assembly: vector diffusion
$(\nabla u_h,\nabla v_h)$, divergence $(\nabla\cdot u_h,q_h)$ and pressure
$(p_h,\nabla\cdot v_h)$. Both mixed orientations are measured independently.
An explicit unit-viscosity diffusion expression selects a fourth kernel under
`WeightedVectorDiffusion`. The 216 registrations span four terms, three degree
pairs, six geometries and $n\in\{3,5,9\}$. Formula order 16, spaces and geometry
match the complete
assembly consumer. Each timed sweep includes cell traversal, owned mapped
quadrature construction, physical basis preparation, eager local contraction
and one entry read, but no insertion, constraints or solver work.

Before timing, every matrix entry on every selected cell is compared with
`QuadratureReference`, the original generic expression-entry loop. This oracle
also accepts a directly supplied expression: the vector-Jacobian specialization
does not expose its owned integrand, and no production accessor is added merely
for instrumentation. Perturbing an entry must fail qualification. A final sweep
must reproduce all initially retained entries exactly. Reference-loop work and
retained qualification matrices are outside timing; they remain part of process
memory and must not be interpreted as production assembly storage.

The first qualification covers $P_2/P_1$ on all geometries at $n=3$ and all
three degree pairs on triangles at $n=3,5,9$: 42 selections in each of the
native and OpenMP-enabled builds. The matching 20 complete assembly selections
also qualify with two OpenMP workers. Local sweeps remain serial in both builds;
OpenMP-build qualification is not a parallel speedup measurement. Five-repeat
native $P_2/P_1$, $n=3$ medians are:

| Geometry | Vector diffusion | Divergence | Pressure |
| --- | --- | --- | --- |
| Triangle | 0.0504 ms | 0.0306 ms | 0.0285 ms |
| Quadrilateral | 0.0596 ms | 0.0346 ms | 0.0320 ms |
| Tetrahedron | 11.090 ms | 4.853 ms | 5.085 ms |
| Pyramid | 17.317 ms | 7.177 ms | 7.279 ms |
| Hexahedron | 7.963 ms | 3.270 ms | 3.210 ms |
| Wedge | 5.529 ms | 2.271 ms | 2.198 ms |

These are isolated baselines, not additive stage fractions of a separately
measured assembly. Local rows, columns, total dense entries, cells and quadrature
points are reported with the timings. Dense-entry counts include structural
component zeros and are not a count of executed multiply-adds. The triangle and
tetrahedron have respectively $12\times12$ and $30\times30$ vector matrices,
in addition to the cell/quadrature differences stated above. Vector diffusion
is the largest of these three terms. Its unweighted kernel now shares physical
gradient scratch only for identical borrowed reference tables, as described below.

Embedded-cell qualification is a separate prerequisite for this reuse experiment.
It exposed ambient-column truncation in H1 vector Jacobian evaluation and P1
shape Jacobians, as well as incorrectly shaped global-constant zero derivatives.
The generic expression oracle must retain the same ambient derivatives as the
specialized integrator before it can qualify an embedded performance change.
The physical dimension contract is stated in
[the numerical contracts](../../doc/agents/numerical-contracts.md#physical-differential-dimensions).
Analytic tangent-projector regressions cover real/complex P1, H1 degrees 1–6 and
global constants with component counts 1–3, on embedded segments, triangles,
quadrilaterals and points. A separate all-entry quadrature regression includes
all volume families, both mixed-degree directions and curved embedded triangles.
These checks were observed failing before repair. They establish correctness
outside the previous volume benchmark matrix; they are not a speedup or a new
performance-coverage percentage. That correctness repair did not change
full-dimensional contraction arithmetic or specialized assembly kernels.

### Shared vector-Jacobian preparation

Unweighted H1 vector-Jacobian binding compares the two borrowed reference-table
object identities during each bind. Identical tables imply identical gradients
under the current geometric inverse and share one existing scratch buffer;
distinct tables retain independent buffers. Physical gradients are recomputed
on every cell and quadrature point, with unchanged mapping expressions and
entry-summation order. No reference-table pointer survives the bind and no new
cross-cell cache is introduced. The original entry-loop oracle checks every
matrix entry. Mixed-degree regressions include both $P_1/P_2$ and $P_2/P_3$
directions, real/complex fields, changed orders and curved embedded cells.
Forcing unconditional sharing was observed to fail the regression before the
identity guard was restored.

Fresh before/after Apple LLVM 21 RelWithDebInfo measurements use five repetitions,
$P_2/P_1$, $n=3$ and quadrature order 16. Both executables retain their
independent numerical qualifications outside timing:

| Geometry | Local before (ms) | Local after (ms) | Assembly before (ms) | Assembly after (ms) |
| --- | --- | --- | --- | --- |
| Triangle | 0.0504 | 0.0437 | 0.1790 | 0.1688 |
| Quadrilateral | 0.0597 | 0.0530 | 0.1886 | 0.1835 |
| Tetrahedron | 11.1589 | 8.6975 | 26.8680 | 25.1601 |
| Pyramid | 16.9699 | 14.4913 | 38.4778 | 36.1106 |
| Hexahedron | 7.9732 | 6.9986 | 16.3942 | 15.2367 |
| Wedge | 5.5219 | 4.6350 | 12.0195 | 11.1148 |

The selected local-binding speedup is 1.13–1.28 times; complete native assembly
speedup is 1.03–1.08 times, with the smallest changes requiring cautious
interpretation. Matching Instruments captures perform 800 tetrahedral local
sweeps: selected kernel CPU sample weight decreases from 8.851 to 6.901 seconds.
These weights are not elapsed timings and include matching untimed binds.

The timed complete consumer is the public unit-viscosity Stokes expression with
an **unweighted** Jacobian integral. The configurable-viscosity convergence helper
uses an explicitly weighted Jacobian expression, even at unit viscosity; its
separate specialization was not modified by this gradient-sharing change and no
speedup is claimed for it here. Its separate measurement protocol follows below.
Serial local
binding in the LLVM 19 OpenMP-enabled build improves by 1.12–1.21 times. The
two-worker whole-assembly measurements do not establish a parallel gain or
regression: substantial elapsed variability coincides with competing SWIFT and
Docker CPU workloads. An unchanged-binary control also changes substantially
under that contention. Parallel scaling needs a quieter-host remeasurement;
local OpenMP-build qualification is not evidence of parallel speedup.

### Weighted vector-Jacobian contract and measurement

The weighted scopes retain the same spaces, meshes, formula order, analytic
patch data and qualification guards as their unweighted controls. Their local
operator is

$$
A_{(b,i),(a,j)}
=\sum_q w_q\,\delta_K(\widehat x_q)\,
C_{ij}(x_q)\,\nabla_K\phi_a(x_q)\cdot\overline{\nabla_K\psi_b(x_q)}.
$$

Here $C=cI$ for a scalar viscosity, and $\delta_K$ is the mapped measure
factor. For a reference dimension $d$ and ambient dimension $s$, the physical
gradients retain all $s$ components of $(J_K^\dagger)^T\widehat\nabla\phi$.
Point cells have zero derivatives. A general component matrix $C$ does not
justify triangular mirroring; only equal-basis scalar-coefficient operators
use that reduction. These logical contracts are independently regressed for
real/complex fields, equal and mixed degrees, volume and embedded geometries,
curved triangles and point cells. Before repair the entry oracle rejects the
embedded and nonsymmetric-matrix operators, while assertion-enabled point
binding aborts. No performance improvement is inferred from correcting an
incorrect operator. Full-dimensional scalar controls retain their original
contraction arithmetic.

The initial weighted qualification selects 14 local bindings in each build:
all six geometries at $P_2/P_1$, $n=3$, and all three degree pairs on triangles
at $n=3,5,9$, counting overlapping selections once. Complete native assembly
and SparseLU each qualify 20 selections: all degree pairs and geometries at
$n=3$, plus the $P_2/P_1$ triangle at $n=5,9$. The same 20 assembly selections
qualify with two OpenMP workers. These are 88 selected scope/build combinations,
not qualification of every registered size or backend. The correctness checks
also pass 11 new geometry cases in native, assertion-enabled and OpenMP builds,
seven existing weighted-complex regressions and seven native Stokes
exactness/convergence tests. The assertion-enabled check recompiles the test
translation unit with `-UNDEBUG`, not the complete library in Debug mode.

Corrected native local-binding medians for $P_2/P_1$, $n=3$ are shown below.
Each measurement has three repetitions of three sweeps; these short controls
are not a claimed optimization or parallel-scaling result.

| Geometry | Cells | Points per cell | Vector local rows | Binding sweep |
| --- | ---: | ---: | ---: | ---: |
| Triangle | 8 | 53 | 12 | 0.0382 ms |
| Quadrilateral | 4 | 81 | 18 | 0.0480 ms |
| Tetrahedron | 48 | 729 | 30 | 7.2922 ms |
| Pyramid | 48 | 729 | 42 | 12.1257 ms |
| Hexahedron | 8 | 729 | 81 | 6.2197 ms |
| Wedge | 16 | 477 | 54 | 4.0486 ms |

The tetrahedral sweep performs $48\times729$ mapped-point visits, compared
with $8\times53$ on triangles: an 82.5-fold count ratio before accounting
for larger bases and component counts. The scalar gradient-dot work proxy
$N_KN_qN_b(N_b+1)d/2$ has a 324.2-fold ratio, whereas the observed sweep-time
ratio is 190.9. This proxy excludes mapping and scatter and is not an exact
instruction count. A matching 800-sweep tetrahedral Instruments capture selects
5.635 seconds of sampled CPU weight in the weighted binding call tree: 57.1%
is attributed directly to `QuadratureRule::setPolytope`, 21.0% inclusively to
mapped-quadrature/Jacobian operations, and 16.2% inclusively to spatial dot
products. Inlined gradient preparation and contraction still require finer
source attribution. These categories can overlap; samples include matching
untimed binds and are not a partition of benchmark elapsed time. The raw
before/after controls, regressions and capture are retained in
`build-convergence/certification/weighted-jacobian-20261013/`.

### Weighted physical-gradient reuse

The weighted Jacobian kernel now applies the same reference-table identity
criterion as the unweighted kernel. If both borrowed tables are the same
object, trial and test physical-gradient scratch aliases within the current
quadrature point, and the duplicate test-gradient map is omitted. Different
tables retain independent storage. No new cache, retained table pointer,
coordinate comparison or collective is introduced. Coefficient evaluation,
matrix contraction and accumulation order are unchanged.

The all-entry regression additionally covers both $P_2/P_3$ orientations for
real and complex fields. Forcing reuse across different tables was observed
failing before restoring the identity predicate. Native, assertion-enabled
and OpenMP regressions pass, as do the selected weighted-complex and native
Stokes convergence controls: 47 test executions in total. Before/after
qualification covers 14 bindings and 20 complete assemblies per build, with
two workers for OpenMP assembly. Selected numerical and work counters match
exactly; this is 68 scope/build combinations, not full backend certification.

Native $P_2/P_1$, $n=3$ three-dimensional medians are:

| Geometry | Binding before / after | Assembly before / after |
| --- | ---: | ---: |
| Tetrahedron | 7.177 / 6.856 ms | 23.572 / 22.862 ms |
| Pyramid | 11.976 / 11.502 ms | 34.287 / 32.871 ms |
| Hexahedron | 6.129 / 6.052 ms | 14.787 / 14.211 ms |
| Wedge | 3.996 / 3.880 ms | 10.719 / 10.406 ms |

These measurements use five repetitions of ten sweeps. Repeating the unchanged
baseline after the candidate yields assembly median drift of 0.8–2.0% on these
four geometries. The observed gains are modest and do not resolve the larger
three-dimensional cost. Longer low-degree 2D binding windows, with 10,000
sweeps per repetition, give 0.03615 / 0.03297 ms on triangles and
0.04398 / 0.04069 ms on quadrilaterals. No material 2D whole-assembly gain is
established. Two initial higher-degree 2D control windows, only 1–2.5 ms long,
contained large wall/CPU timing outliers; they are not used as stable controls.
Repeating those cases with 1,000 sweeps reduces elapsed variability below 1.4%
and shows about 1–2% differences. Short-window noise must not be interpreted as
a higher-degree optimization. Serial local binding in the LLVM 19 OpenMP build
shows 1.02–1.09-fold changes across the selected cases; this is not a parallel
speedup claim.

The broad sweep initially measured a 3.4% slowdown for $P_4/P_3$ wedge
assembly. A focused five-repetition comparison gives 161.03 / 162.01 ms,
with elapsed coefficients of variation of 0.54% / 0.96%. The unchanged
baseline also shifts relative to the broad sweep. No material regression or
gain is established for that case; its original and focused observations
remain in the artifacts rather than being replaced by a favorable selection.

Matching 800-sweep native tetrahedral captures select 5.635 seconds before
and 5.468 seconds after in the weighted binding call tree. These are sampled
CPU weights, not elapsed timings. Matching 250-assembly captures select
1.813 / 1.698 seconds specifically in the consuming weighted kernel. Profiles
include matching untimed binds, and changing inlining can change operation
labels; absent named frames do not imply absent arithmetic. Exported kernel
source-line metadata collapses
the inlined body to one line; it cannot independently split the gradient loops.
The isolated source change and matched consumer measurements test the reuse
hypothesis without inferring those missing source-line fractions. Profiles and
controls are retained under
`build-convergence/certification/weighted-gradient-sharing-20261013/`.

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

### Native P1 bulk Jacobians and their assembly consumer

`RodinBulkP1JacobianBenchmarks`, under `Rodin/Geometry`, isolates the native
`RealP1Element` transformation used by mesh geometry. This is not an H1 element
of degree one. For control coefficients $X_{ja}$, reference basis functions
$\phi_a$, reference dimension $d$ and embedding dimension $s$, it computes

$$
J_{ji}(\widehat x_q)=\sum_{a=0}^{n_b-1}
X_{ja}\,\partial_i\phi_a(\widehat x_q),
\qquad 0\le j<s,\quad 0\le i<d.
$$

`BulkJacobian` overwrites retained output matrices through the public bulk
transformation method, including the existing table lookup, vector resize,
zeroing and contraction. `OwnedQuadrature` constructs/replaces an owned
`PolytopeQuadrature`, including allocation, bulk Jacobians and point-cache
initialization. Coordinates, inverses and measures remain lazy and are not
included in this construction scope. Neither scope includes setup, first-use
reference tabulation, guards, assembly or solving.

All eight reference families, all admissible embeddings through dimension
three, orders 2, 8 and 16, and affine/perturbed control clouds are registered:
168 selections. Each workload contains one cell. Controls are obtained through
the element DOF functionals. Affine cases represent $x=b+A\widehat x$ with mixed
components; perturbed simplex P1 maps remain affine, whereas no constant-Jacobian
assumption is made for tensor-product, wedge or pyramid maps. Every active
entry must agree exactly with the unchanged pointwise contraction. Affine cases
also reproduce $A$ within a dimensionless absolute roundoff budget of
$10^{-12}$. Deliberately wrong entries or dimensions must be rejected.

An Apple LLVM 21 optimized native study used five repetitions per selection.
The following order-16 CPU medians exclude setup and guards. The last column
normalizes the bulk time by $N_q n_b d s$ scalar products; it is a work proxy,
not an instruction count.

| Geometry / embedding | $N_q$ | $n_b$ | Bulk, ns | Owned quadrature, ns | Bulk ns/product |
| --- | ---: | ---: | ---: | ---: | ---: |
| Triangle / 2 | 53 | 3 | 537 | 1,084 | 0.845 |
| Tetrahedron / 3 | 729 | 4 | 14,947 | 23,662 | 0.570 |
| Pyramid / 3 | 729 | 5 | 18,427 | 26,847 | 0.562 |
| Hexahedron / 3 | 729 | 8 | 29,776 | 39,183 | 0.567 |
| Wedge / 3 | 477 | 6 | 14,400 | 20,087 | 0.559 |

The tetrahedral raw bulk cost is about 27.8 times the triangular cost, while
the scalar-product proxy is about 41.3 times larger. This comparison does not
establish dimension-dependent inefficiency: the selected three-dimensional
families have similar normalized costs. An Instruments capture selected by
both benchmark owner and bulk-kernel frames attributed 88.2% of selected CPU
sample weight to the contraction method itself and 9.6% to reference-gradient
access. The warmed table lookup contributed less than 0.1%; a second cache is
not motivated by this profile.

The existing Assembly-owned coupled reaction-diffusion benchmark supplies the
consumer comparison: two real H1 degree-two fields, two diffusion terms, four
reaction terms, native P1 geometry, order 16, and grid sizes 3, 5 and 9. Native
assembly CPU medians were 0.220, 0.889 and 3.551 ms on triangles, versus 30.202,
254.350 and 1,959.505 ms on tetrahedra. At size 9, the latter contains 3,072
cells, 729 points per cell and 10 local field functions, versus 128 cells,
53 points and 6 functions on triangles. The approximately 552-fold raw ratio
therefore combines different cell, quadrature and local-matrix counts, rather
than isolating spatial dimension. A size-5 tetrahedral Instruments capture
selected within actual assembly attributed 20.2% inclusively to bulk P1
Jacobians, including 16.5% leaf weight in the contraction method. This confirms
an assembly consumer of the measured kernel; inclusive weights are not additive
or predictions of end-to-end speedups.

Two candidate changes were rejected. Copying packed controls once per local
basis preserved exact numerics but increased triangular order-16 bulk CPU time
from 537 to 595 ns, with essentially unchanged tetrahedral time. An Eigen
column update failed exact bulk/pointwise equality in four geometry regression
tests and was not timed as an eligible repair. The canonical production kernel
is unchanged; no source or assembly speedup is claimed. The retained regression
poisons packed input padding and output storage, varies dimensions and formula
sizes, and enforces exact active entries and zero inactive output storage.
All 168 selections and all 14 geometry unit tests passed in both native and
OpenMP-enabled optimized builds. The 42 real/complex degree-two assembly
consumer selections retained exactly the native baseline numerical/work counters.

```sh
cmake --build build --target RodinBulkP1JacobianBenchmarks -j1
build/tests/benchmarks/RodinBulkP1JacobianBenchmarks \
  --benchmark_filter='BulkP1Geometry/(BulkJacobian|OwnedQuadrature)/Affine/Tetrahedron/3/16/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

### Local quadrature binding and physical assembly

`RodinQuadratureBindingBenchmarks`, under `Rodin/Variational`, measures public
`Integral::setPolytope` sweeps for the four distinct operator expressions in
the coupled reaction-diffusion consumer:

$$
m(u,v)=\int_\Omega u\,\overline v\,\mathrm{d}x,\qquad
m_\alpha(u,v)=\int_\Omega(\alpha u)\overline v\,\mathrm{d}x,\qquad
a(u,v)=\int_\Omega\nabla u\cdot\overline{\nabla v}\,\mathrm{d}x,\qquad
a_2(u,v)=\int_\Omega(2\nabla u)\cdot\overline{\nabla v}\,\mathrm{d}x,
\quad \alpha=0.2.
$$

Complex cases follow Rodin's trial-first, test-second sesquilinear convention;
the overbars reduce to the identity over real fields. Each scope includes
mesh traversal, formula selection,
owned mapped-quadrature construction, basis preparation and local contraction.
Reading one prepared entry prevents dead-code elimination. Global scatter,
boundary constraints and solving are excluded; complete physical assembly
remains in its Assembly-owned benchmark. Geometry remains native P1, while
field H1 degrees one through three, real/complex scalars, seven volume families
and sizes 3, 5 and 9 give 504 registered selections. This registration count is
not a claim that every selection has been profiled or numerically qualified.

Before timing, every cell matrix is checked against the retained original
entry-loop algorithm with the existing relative budget $10^{-11}$ and scale
$\max(1,\|M_{\mathrm{ref}}\|_F)$. Deliberately perturbing one entry by that
scale must fail. An untimed sweep after timing must reproduce every initially
retained entry exactly. The entry-loop oracle shares mapped geometry; it is
not an independent mapping certification. Reference-table initialization,
oracle evaluation and retained-matrix allocation are excluded from timing.

The initial qualification covered all 56 degree-two real/complex cases at
size 3 in both native and OpenMP-enabled optimized builds. Local binding is
sequential in both builds; this is not a parallel scaling study.
Selected size-5 real studies used five repetitions. Initial CPU
medians on triangle/tetrahedron meshes were respectively 0.100/29.113 ms for
mass, 0.122/38.871 ms for weighted mass, 0.108/33.540 ms for diffusion and
0.116/43.219 ms for weighted diffusion. An elapsed-time variability anomaly in
the initial mass sweep prompted isolated fixed-work repeats, rather than a
source-performance conclusion. Those repeats measured 0.090/27.011 ms for
mass and 38.539 ms for tetrahedral weighted mass, with CPU coefficients of
variation 0.56%, 1.18% and 0.65%, respectively; elapsed coefficients were
0.62%, 1.57% and 0.75%. The original samples are retained. No speedup is claimed
between these unchanged-source runs.

At size 5, triangles contain 32 cells with 53 quadrature points and six field
basis functions per cell; tetrahedra contain 384 cells with 729 points and ten
functions. The local matrix accumulation proxy
$N_{\mathrm{cells}}N_qn_b^2$ is about 458.5 times larger on tetrahedra. The
observed real CPU ratios, approximately 300--372 depending on the operator
and process repeat, do not demonstrate an additional dimension-only penalty.
This proxy excludes mapped geometry, gradient transforms, traversal and
allocation, so equal normalized values are not expected.

Separate Instruments captures of size-5 tetrahedral weighted mass and
diffusion selected stacks containing both the benchmark owner and
`QuadratureRule::setPolytope`. Weighted mass placed 72.0% of selected CPU
sample weight in the local preparation/contraction method itself and 16.7%
inclusively in bulk P1 Jacobians. Diffusion placed 52.5% in its own method,
19.4% inclusively in bulk Jacobians, and 5.9% in spatial dot products.
Allocation-related frames contributed less than 0.6% in either capture.
Inlined operations remain attributed to their enclosing method; these
profiles do not yet separate every source expression inside that method.
Untimed guards also invoke production kernels, so selected CPU samples are
source attribution, not timed-stage fractions. The within-cell gradient reuse
described below removes one identified duplication without changing contraction
arithmetic; remaining scalar contraction and gradient-expression costs still
require source attribution and complete-assembly measurements.

```sh
cmake --build build --target RodinQuadratureBindingBenchmarks -j1
build/tests/benchmarks/RodinQuadratureBindingBenchmarks \
  --benchmark_filter='CoupledReactionDiffusion/Binding/.*/H1/P2/Real/Tetrahedron/5/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

### Physical-gradient preparation and within-cell reuse

For reference dimension $d$, embedding dimension $s$, cached inverse/left-inverse
$B_q\in\mathbb{R}^{d\times s}$ and reference basis gradient
$\widehat g_a(q)$, physical basis preparation evaluates

$$
g_a(q)=B_q^{\mathsf T}\widehat g_a(q).
$$

`RodinPhysicalGradientBenchmarks` measures this preparation through the public
H1 gradient expression, including integration-point rebinding and basis
extraction. One affine cell and its mapped quadrature are retained; inverse
Jacobians and reference tables are warmed. Every sample changes the quadrature
index, so this is not a repeated-point cache-hit measurement. Mapped-quadrature
construction, matrix contraction, scatter and solving are excluded. Specialized
diffusion rules own their own transformation loops; expression timings are not
an additive decomposition of those rules.

Degrees one through three, real/complex scalar fields, seven volume families
and orders 2, 8 and 16 give 126 cases. Each case reconstructs the analytic
gradient of $f(x)=\eta(1+\sum_j(j+1)x_j)$ from DOF-functional coefficients,
with $\eta=1$ or $1+i/2$, within a dimensionless absolute budget $10^{-10}$.
A contributing basis perturbation must be rejected. Every prepared basis
component is retained and must repeat exactly after timing. All 126 cases
passed in both optimized native and OpenMP-enabled builds. Curved maps,
embedded cells and matrix-valued gradient-expression timings are not certified
by this benchmark matrix.

Native degree-two, order-16 real CPU medians were approximately 3.20 microseconds
per triangular sweep and 76.23 microseconds per tetrahedral sweep. These contain
53/729 points and 6/10 basis functions, respectively. Normalizing by
$N_qn_bd^2$ gives approximately 2.52/1.16 ns per scalar-product proxy, not a
dimension-only slowdown. Isolated five-repetition fixed-work repeats had CPU
coefficients of variation near 1.1%; the initial wider sweep and its noisier
selections remain in the artifacts. In the selected tetrahedral preparation
profile, matrix-vector multiplication carried 16.0% of CPU sample weight;
`bzero` and explicit spatial-vector zero-storage frames together carried
10.7%. Inlined work remains attributed to the enclosing gradient method.
These observations motivate the next Math-owned matrix-vector and storage
benchmarks, not unconditional removal of padding initialization.

The unweighted H1 diffusion rule previously transformed both trial and test
reference gradients even when both operands borrowed the same table object.
It now shares one physical-gradient buffer only under that exact object
identity during the current bind. Distinct tables retain independent buffers.
Every cell and quadrature point still recomputes physical gradients with its
own inverse Jacobian; physical results are not reused across cells or meshes,
scratch storage remains thread-local, and no collective is introduced.
Arithmetic expressions, matrix
summation order and quadrature requirements are unchanged.

Five-repetition native degree-two real binding measurements at grid size 5
gave 32.577 to 30.041 ms on tetrahedra and 19.528 to 18.826 ms on hexahedra.
Additional tetrahedral size-3 studies gave speedup ratios 1.06/1.27 for degree
one real/complex and 1.09/1.19 for degree three real/complex. The LLVM 19
OpenMP-enabled build gave 33.312 to 29.929 ms for degree-two real tetrahedral
binding; this operation is sequential, not an OpenMP scaling measurement.
All matched numerical/work counters remained identical. Native complete
two-field degree-two real tetrahedral assembly measured 250.118 to 250.495 ms,
within CPU variation of 1.3--1.7%; no end-to-end assembly speedup is established.

The exact same-/mixed-degree entry regression failed when table sharing was
forced unconditionally. The scalar/vector/matrix component regression failed
when mass components were omitted. Restored implementations passed these tests,
the matrix-range suite and all fourteen Poisson polynomial-convergence cases
in both builds. Diffusion binding additionally qualified degrees one through
three, real/complex fields and seven geometries at size 3 in both builds; the
fourteen native degree-two full-assembly cases at that size retained their
original numerical/work counters. These are selected recertification results,
not full performance coverage of the convergence workplan.

A separate mass-kernel prototype replaced runtime component-count division
with the existing space accessor. Its repeated tetrahedral weighted-mass CPU
time increased from 38.714 to 44.114 ms, so it was rejected. The production mass
kernel is unchanged. The benchmark and rejected-prototype measurements remain
available; no compiler flags or numerical budgets were relaxed.

```sh
cmake --build build --target RodinPhysicalGradientBenchmarks -j1
build/tests/benchmarks/RodinPhysicalGradientBenchmarks \
  --benchmark_filter='PhysicalGradient/H1/P2/Real/Tetrahedron/16/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

### Spatial dot batch protocol

`RodinSpatialDotBenchmarks` retains its 32 single-contraction scopes and adds
192 paired batch scopes. Both scalar orders over $\mathbb R$ and $\mathbb C$,
logical dimensions $d\in\{0,1,2,3\}$ and vector counts $N\in\{4,16,64\}$
are represented. Each batch stores

$$
D_{ba}=\sum_{j=0}^{d-1}u_{bj}\overline{v_{aj}}.
$$

`Full` evaluates $N^2$ pairs and `Lower` evaluates $N(N+1)/2$ pairs. The latter
is a loop shape, not a symmetry assumption: left and right arrays differ, and
no omitted entries are reconstructed by mirroring. Barriers surround the whole
batch, permitting normal caller optimization within the pair loop. Vector
construction, allocation and reference checks are outside timing; pair-loop
control and output stores remain included. This is a Math operation benchmark,
not a replacement for mapped-gradient, weighted-kernel or assembly measurements.

The historical dynamic reduction and its reference use the same const Eigen
view types as the original public implementation, and obtain dimension from
the left operand rather than a separately supplied loop constant. Mutable views
produce different Eigen expression types even when their numerical values agree.
Compile-time guards reject that comparator drift; a deliberately mutable view
was observed failing compilation before restoration. Inactive vector entries
are poisoned with NaNs. All output entries must equal the dynamic reference
before and after timing, and changing one entry must be rejected, including
empty contractions. Existing scalar-product unit regressions additionally
retain bitwise, cancellation, signed-zero and infinity checks.

All 224 scopes qualify in both native and LLVM 19 OpenMP-enabled builds: 448
scope/build combinations, with seven scalar-product regression tests passing
in each build. The Math batches themselves remain serial and opt-in, without
new CTest or publication-workflow registrations. Counters distinguish vector
count, active pair count, field-scalar product terms, retained input/output
bytes and bytes written; field-scalar terms are not a real-flop count.

For $N=64$, the lower real/real batch has 2,080 pairs. Selected medians, in
nanoseconds per pair, compare the existing public operation with the historical
dynamic reduction; they are not a speedup from a new library repair:

| Build | Dimension | Public | Dynamic reference |
| --- | ---: | ---: | ---: |
| Native Apple LLVM 21 | 2 | 0.615 | 1.130 |
| Native Apple LLVM 21 | 3 | 0.797 | 1.519 |
| LLVM 19, OpenMP enabled | 2 | 0.613 | 1.154 |
| LLVM 19, OpenMP enabled | 3 | 0.827 | 1.489 |

These use five repetitions with a 20 ms minimum timing window. Isolated
per-pair timings are not additive fractions of separately compiled assembly.
In particular, the weighted Stokes profile's Eigen `coeff` label resolves,
after checking the image UUID and rebasing its address, to an inline
third-component load inside `QuadratureRule::setPolytope`, not a separate
function call. A method-shaped profile label alone does not justify forced
inlining or an interposition repair. Raw traces, generated instructions,
compiler optimization records and the explicitly noncanonical initial
mutable-view measurements are retained in
`build-convergence/certification/spatial-dot-batches-20261014/`.
### Complex one-component dispatch findings

The one-component complex difference is a calling-context dispatch cost, not
an additional complex-algebra operation or an out-of-line Eigen getter. In
the native lower-pair loop, sampled instructions in the public path compare
logical dimension with three, two and one before entering the one-component
arithmetic. The historical dynamic reduction reaches its first complex packet
after the empty-range check, followed by its reduction-loop control. Both
paths use one-complex-value packets; their Eigen reductions differ in complete
versus runtime unrolling. Changing only that unrolling instantiation leaves
the surrounding dimension checks and does not remove the observed difference.

Controlled source experiments retained the existing reduction formulas:

| Experiment | Observation | Disposition |
| --- | --- | --- |
| Dynamic reduction inside the existing one-component case | No stable removal of the one-component cost | Reverted |
| Wider dimension metadata | No consistent removal; the dimension load is already outside the inner loop | Reverted; original representation and serialization retained |
| Separate zero/one complex range before fixed-extent dispatch | Removes the one-component penalty but slows the two-component and empty cases | Reverted |
| Separate exact one-component condition | Optimizes back to effectively the original dispatch; no stable gain | Reverted |

No shared Math repair is retained merely to improve one dimension at the
expense of another. The benchmark scopes and all experimental observations
remain available. Producer sites with invariant dimension require their own
consumer-level evidence; a global branch-order or representation change is
not justified by this microbenchmark alone.

Source regression coverage now includes full complex pair batches of 4, 16
and 64 operands after shrinking three-entry scratch to one active entry.
Inactive entries are poisoned, both result components are compared bit-for-bit
with the original const-view reduction, and an unconjugated result is rejected.
The wrong-conjugation mutation was observed failing. The older randomized
reference now also retains the original const expression types rather than
using mutable Eigen views. All eight source tests pass under native Apple LLVM
21, LLVM 19 and GCC 15, and all 224 benchmark scopes qualify in each Clang build.
The GCC check is a local optimized standalone build of the same test source,
not evidence for GCC 10/14 CI completion.

The physical baseline additionally qualifies 24 local-binding and 12 coupled
reaction--diffusion assembly selections: complex segment degrees 1–3 at
$n\in\{3,5,9\}$, with real-segment and complex triangle/tetrahedron controls.
No end-to-end improvement is claimed because no production repair was accepted.
Profiles, dispatch diagnostics, rejected experiments and final regressions are
retained under `build-convergence/certification/complex-one-dispatch-20261014/`.

### Spatial matrix-vector initialization and its assembly consumer

`RodinSpatialMatrixVectorBenchmarks`, owned by Math, measures the canonical
runtime-shape product $y=Ax$ including result construction. All shapes
$0\le r,c\le3$ and four real/complex operand pairs give 64 registrations.
Mixed scalar promotion and genuinely complex products are included. Input
padding is poisoned, while every inactive result entry must be zero.
Independent row sums use bounded dyadic operands and at most three products,
so equality is exact, without a numerical tolerance. Deliberately wrong values
or dimensions must be rejected. Operand/result addresses escape on every
iteration, preventing constant folding or discarded padding stores.

The source repair removes only the second active-entry clear in the product.
`SpatialVector` already initializes its complete storage in its constructor;
that initialization remains mandatory and unchanged. Empty contractions retain
zero active values, and empty/shrinking results retain zero inactive capacity.
The new unit regression covers all 64 shape/scalar combinations, including
poisoned destination storage. A padding-corruption mutant failed; both the
original implementation and the restored cleanup pass. The retained
`RodinSpatialVectorStorageBenchmarks` additionally passed all 144 native scopes;
no constructor/storage-throughput improvement is claimed from this repair.

Five-repetition native CPU medians for real $2\times2$/$3\times3$ products were
3.926/4.290 ns before and 1.352/1.511 ns after removal. Selected $3\times3$
mixed/complex products showed ratios approximately 1.45--1.46. The empty real
control was unchanged within variation. Empty complex products also improved,
although their runtime path never entered the removed loop: disassembly shows
the common function prologue now saves six fewer registers and uses a
144-byte rather than 176-byte stack frame. Their improvement therefore includes
code-generation effects, not arithmetic or runtime clearing on an empty result.
LLVM19 controls compiled with the original Git header show the same qualitative
pattern; these are separate compiler/build measurements, not parallel scaling.

Degree-two physical-gradient sweeps on triangles/tetrahedra retained all
numerical/work counters and showed native CPU speedup ratios 1.24--1.50 across
real/complex fields. The Assembly-owned coupled executable now also registers
126 explicitly labeled `GenericGradientAssembly` cases, alongside its original
126 canonical cases. Only the second diffusion block places the real factor
two on the test operand:

$$
\int_\Omega\nabla w\cdot\overline{2\nabla z}\,\mathrm{d}x
=\int_\Omega(2\nabla w)\cdot\overline{\nabla z}\,\mathrm{d}x.
$$

The physical problem, sources, coupling, constraints, degrees, geometries and
three mesh sizes are unchanged. This benchmark-only expression selects the
existing generic gradient path; it does not modify convergence formulations or
add a production dispatch path. Canonical specialized diffusion transforms
gradients directly and is retained as a negative timing control.

Selected native real tetrahedral size-3 generic assembly measured 36.463 to
33.686 ms CPU, with coefficients of variation 0.50%/0.25%. Canonical real and
complex assembly and generic complex assembly remained within variation.
The two-thread LLVM19/OpenMP generic assembly elapsed ratio was approximately
1.03, within elapsed variation of 3.6%; no parallel assembly gain is established.
Its default CPU clock measures the calling thread, not aggregate worker CPU.
Before/after selected numerical/work counters were identical. Full-assembly
profiles selected within the owning assembler show gradient preparation at
19.2% before and 12.4% after of sampled CPU weight. The matrix-vector method
was separately visible before and inlined afterward; missing separate frames
must not be interpreted as zero remaining arithmetic cost.

All 64 primitive cases and 126 gradient cases passed in each optimized build.
Generic degree-two real/complex assembly qualified all seven geometries at
size 3 natively; the OpenMP build additionally recertified those fourteen
canonical cases. Both builds passed the exact unit regression and fourteen
Poisson polynomial-convergence cases; fourteen existing native tensor/scalar
product tests also passed. Wider timing, curved/matrix-field consumers and
distributed assembly profiling remain separate workplan items. No overall
performance completion is inferred from these selected measurements.

```sh
cmake --build build --target RodinSpatialMatrixVectorBenchmarks -j1
build/tests/benchmarks/RodinSpatialMatrixVectorBenchmarks \
  --benchmark_filter='SpatialMatrixVector/.*/3/3/' \
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

### Sparse algebra protocol

`RodinSparseAlgebraBenchmarks`, owned by `Math`, isolates $y\leftarrow Ax$,
$y\leftarrow A^Tx$, $x^*z$ and $\lVert x\rVert_2^2$ on the same manufactured
constrained operators
used by the assembly and CG workloads. The matrix is the actual compressed
column-major `Math::SparseMatrix`; no alternative storage or synthetic stencil
is substituted. Input vectors vary by entry, with non-collinear real/imaginary
parts so complex reductions exercise conjugation. Products overwrite an already
allocated output. Mesh/space setup, assembly, input/output allocation and
verification are excluded; these are warmed arithmetic scopes, not cold
allocation or complete solver measurements.

The installed Eigen full-matrix CG uses a transpose view for real column-major
operators and a direct view for complex operators. Both product paths therefore
have distinct scopes. `TransposeProduct` is a plain transpose, not an adjoint;
its independent oracle swaps row/column accumulation indices. No owning layout
conversion is performed, and no equivalence is assumed for arbitrary matrices.

Independent scalar traversals check all products and reductions, using
`long double` accumulators. The precision is reported as
`reference_mantissa_bits`: on the native Apple ARM build it is 53, so extended
precision must not be assumed. An original-kernel matrix/load comparison and
the manufactured equation separately guard the input operator. Wrong products,
zero dot products and zero norms must be rejected. Acceptance uses a relative
budget of $10^{-11}$ with a unit scale floor; the manufactured patch uses
$10^{-10}$. These checks are untimed and are not timing gates.

Four operations, three scalar equations, real/complex H1 degrees one through
three, seven mesh geometries and three grid levels give 1,512 selections.
Reports include cells, matrix rows/nonzeros, scalar size and compressed storage
bytes. The latter counts values and indices, not allocator overhead, resident
memory or memory traffic. The `nnz` counter is Eigen's stored-entry count,
including any explicit zeros retained by assembly. Source profiles and numerical
qualification must be distinguished from this availability count. This does not cover P0/P0g mass
operators, mixed systems, degrees four through six or PETSc/MPI sparse algebra.

All 1,512 selections passed in the native and two-thread OpenMP-enabled
optimized builds, with identical case lists and three fixed iterations per
case. Product checks matched the independent traversal exactly; maximum
relative dot/norm differences were below $3.2\times10^{-14}$ and
$5.4\times10^{-14}$, respectively. Fixed-iteration qualification samples are
not performance baselines. They do not establish an OpenMP speedup or complete
the remaining space/operator/backend axes listed above.

A native five-repetition study covers all four operations over three Poisson
degree-two grid levels in triangle and tetrahedron meshes. Largest-level
medians, in microseconds, are:

| Geometry / field | Product | Transpose product | Hermitian dot | Squared norm |
| --- | ---: | ---: | ---: | ---: |
| Triangle / real | 1.455 | 1.199 | 0.051 | 0.045 |
| Tetrahedron / real | 43.360 | 48.881 | 1.184 | 1.154 |
| Triangle / complex | 2.751 | 3.237 | 0.171 | 0.214 |
| Tetrahedron / complex | 90.029 | 101.685 | 2.899 | 4.746 |

The short complex triangle dot initially varied by about 7%; an isolated
seven-repetition, longer-duration repeat retained the median with 0.7%
variation. All other initial selections varied by at most 2.9%. Both raw runs
are retained. Direct and transpose paths are not interchangeable performance
models: even their relative costs differ between the selected 2D/3D systems.

Six Instruments captures distinguish the real-CG transpose products, complex
direct products, and complex reductions. Attribution selects stacks containing
the operation's kernel, excluding untimed FE fixture construction visible in
the raw captures. About 86%/88% of selected real product sample weight lies in
Eigen's row traversal for triangle/tetrahedron cases; about 83%/89% of complex
product sample weight lies in standard-library multiply/conversion helpers
nested in the sparse product. Dot and squared-norm profiles instead expose
conjugating packet reductions and real/imaginary component access. These are
source hotspots to investigate, not evidence that replacing standard complex
arithmetic or changing floating-point policy would preserve correctness.
The native binary has no standalone definitions of those complex
multiply/conversion helpers: sampled inline source frames must not be mistaken
for out-of-line function-call overhead. Nor does a transpose view uniformly
improve throughput; the measured tetrahedral direct products were faster.

```sh
cmake --build build --target RodinSparseAlgebraBenchmarks -j1
build/tests/benchmarks/RodinSparseAlgebraBenchmarks \
  --benchmark_filter='^SparseAlgebra/(Product|TransposeProduct)/Poisson/H1/P2/(Real|Complex)/Tetrahedron/' \
  --benchmark_min_time=0.2s --benchmark_repetitions=5
```

### Diagonal-preconditioner protocol

`RodinDiagonalPreconditionerBenchmarks`, owned by `Solver`, isolates the
diagonal preconditioner selected by the current native CG backend. It uses the
same constrained manufactured operators as the sparse-algebra and CG scopes.
This is backend-stage attribution, not an additional Rodin solver or inspection
API. With $D=\operatorname{diag}(A)$, application computes

$$
y=D^{-1}x,\qquad y_j=\frac{x_j}{A_{jj}}.
$$

| Scope | Timed work | Excluded |
| --- | --- | --- |
| `OwnedCompute` | New instance, diagonal compute and destruction | FE assembly, inputs, verification |
| `RetainedCompute` | Diagonal recomputation with retained, already-sized storage | Initial allocation, application |
| `Apply` | Preconditioner application into a preallocated vector | Compute, output allocation |

Owned construction is not a process-first or allocator-cold guarantee. Setup
scopes do not include CG matrix binding, stopping logic or iterations. Independent
componentwise division verifies the result with a $10^{-11}$ relative budget and
unit scale floor. Fresh/retained instances are checked outside timing; the last
timed application is checked before its output is overwritten. Wrong outputs
are rejected. Recomputing on $1.5A$ must yield $y/1.5$, then recomputing on $A$
must restore $y$. This tests changed numerical state, not only repeated setup
on an immutable matrix.

Three lifetimes across real/complex degrees one through three, seven geometries,
three scalar equations and three grid levels give 1,134 selections. All 378
smallest-grid cases passed in both optimized builds; maximum relative error
was below $6.8\times10^{-17}$. Reports distinguish matrix stored entries from
`diagonal_search_entries`, the entries visited through each column's diagonal,
and inverse-diagonal storage bytes. These are work/storage counts, not traffic
or instruction counts. Zero/missing-diagonal fallback, mixed/block operators,
higher degrees and PETSc/MPI preconditioners remain separate workloads.

A pre-integration native baseline contains 36 five-repetition Poisson
degree-two measurements across triangle/tetrahedron meshes and three levels.
At the largest level, retained real setup measured approximately 21–23
microseconds in tetrahedra; complex setup was approximately 62 microseconds.
Application measured approximately 0.84 microseconds real and 2.6 microseconds
complex. Compared with the corresponding fixed-system CG measurements, setup
is below 1% of the solve cost. These are separately measured scopes, not an
exact additive time decomposition.

Initial variation reached 9%; longer repeats of six noisy largest-level cases
still reached about 9% and showed drift. Raw runs are retained, and no precise
owned-versus-retained speedup or source optimization is claimed. Six Instruments
captures isolate real/complex retained setup in both dimensions and application
in tetrahedra. Setup resolves to compressed-column diagonal lookup and division;
the complex division runtime helper accounts for about 61%/48% of selected
setup samples in triangle/tetrahedron cases. Application is predominantly the
dense packet assignment loop. Kernel-stack selection excludes untimed FE setup.
Owned construction/destruction attribution remains pending; the baseline
predates the subsequent `develop` integration. After that integration, all
1,134 selections passed in each optimized build, with maximum relative error
below $1.1\times10^{-16}$. All 164
affected numerical tests passed in each build, including reference-element,
quadrature, pyramid-apex and Poisson degree-refinement regressions.

### Matched Helmholtz assembly and GMRES protocol

`RodinHelmholtzAssemblyBenchmarks` and `RodinHelmholtzSolveBenchmarks` belong
to Assembly and Solver, respectively. `Workloads/HelmholtzWorkload.h` defines
their common equation on the unit box:

$$
-\Delta u+\alpha u=f,\qquad u|_{\partial\Omega}=u_\star.
$$

The affine and quadratic fields from convergence certification are used for
degree one and higher degrees, respectively. The source is shifted analytically
from the existing Helmholtz source, not obtained by multiplying a fabricated
coefficient vector by the assembled matrix. The two mass coefficients are
$\alpha=-1/4$ and $\alpha=-64+i$. The first matches the low-wave-number
convergence context; the continuous real part of the second is indefinite in
dimensions one through three. This does not assert that every coarse discrete
interior block has eigenvalues of both signs.

Assembly measures warmed full `Problem::assemble`, including the matrix, load
and Dirichlet processing but excluding setup and solving. GMRES measures its
public fixed-system solve, including diagonal-preconditioner setup and restarted
Householder iterations. Zero initial guesses are restored outside timing. The
restart is fixed at 50, the requested tolerance at $10^{-13}$ and the iteration
cap at 20,000; no geometry-specific tuning is performed. The true residual is
checked independently of the backend's reported error estimate:

$$
r=\frac{\|Ax-b\|_2}{\max(1,\|b\|_2)},\qquad
e=\frac{\|x-x_\star\|_2}{\max(1,\|x_\star\|_2)}.
$$

The solver budgets are $r<10^{-11}$ and $e<10^{-9}$. Assembly compares full
operators and loads with the original local kernels at a $10^{-11}$ relative
budget and checks the manufactured patch at $10^{-10}$. Wrong coefficients
must be rejected. The solver additionally checks an omitted mass action and
records a smooth probe's real Rayleigh quotient. Matrix rows/nonzeros, local
entry counts, quadrature points, solver iterations and restart are reported
together; nonzeros alone do not account for GMRES's dense Householder sweeps.

For the mass-omission control, let $I$ be the free DOF index set obtained from
the assembled boundary constraints and let $\widetilde M$ retain the $I\times I$
mass block with zeros elsewhere. The deliberately wrong operator is
$A-\alpha\widetilde M$. With the nonzero action
$q=\alpha\widetilde M x_\star$, its sensitivity is

$$
s_M=\frac{\|Ax_\star-b-q\|_2}{\|q\|_2}>10^{-4}.
$$

This denominator measures the term deliberately removed, rather than the
boundary-dominated load. The earlier load-normalized guard falsely rejected
five fine-grid low-wave-number cases despite residuals below
$3.7\times10^{-14}$. The corrected control retains the rejection threshold
and uses logical DOF indices, not coordinate comparisons or matrix-entry
thresholds. On those five cases, $s_M$ differs from one by less than
$8\times10^{-12}$; the solver residual and coefficient budgets are unchanged.

Each executable registers complex H1 degrees one through three on all seven
mesh geometries at 3, 5 and 9 points per axis: 126 selections per operation.
All 126 selections passed in each native and two-thread OpenMP build for both
operations: 504 qualified operation/build selections. The largest true solver
residual was $3.3\times10^{-13}$ and the largest coefficient error was
$6.1\times10^{-11}$. The largest assembly patch residual was
$3.1\times10^{-15}$. Fixed-iteration qualification is not a timing baseline.
Four unit cases independently compare public GMRES diagnostics with Eigen for
sparse/dense real/complex systems, including failed solves, recovery and exact
initial guesses; all pass in each build. An isolated header overlay returning
zero for both diagnostics makes all four cases fail; the production getters
pass without changing the solve algorithm.

Matched native degree-two triangle/tetrahedron measurements use all three
registered sizes. Eight Instruments captures isolate both operations, both
coefficients and both geometries at nine points per axis. In selected assembly
stacks, quadrature accounts for 88–96% of sampled CPU weight. Complex
multiplication alone accounts for 33–48% of leaf weight. Geometry and basis
evaluation are nested contributors, not additional disjoint percentages.
GMRES selected stacks instead spend 79–86% in Householder operations and
11–19% in sparse products; diagonal setup/application accounts for less than
0.8%. The source profiles therefore do not justify treating either operation
as uniformly sparse-product dominated. These are sampled stack weights,
including matching warmup calls, not elapsed-time decompositions or speedups.

The native assembly repeats at nine points per axis used seven repetitions
with at least 0.8 seconds per repetition. Observed median CPU times were:

| Coefficient | Triangle | Tetrahedron | CPU coefficient of variation |
| --- | --- | --- | --- |
| $-1/4$ | 1.352 ms | 265.163 ms | 0.65%, 1.13% |
| $-64+i$ | 1.341 ms | 267.333 ms | 1.18%, 0.88% |

At this size, triangles have 128 cells, 32 quadrature points per cell and six
local basis functions; tetrahedra have 3,072 cells, 120 points and ten basis
functions. Thus the entry-quadrature work proxy
$N_{\mathrm{cell}}N_qN_b^2$ increases by 250, whereas the measured CPU time
increases by approximately 196–199. This proxy excludes gradient dimension,
load terms, constraints and insertion; it is not an exact operation count.
The global row counts are 289 and 4,913 and the constrained matrix nonzero
counts are 2,293 and 80,721. The large raw assembly gap is accompanied by a
larger increase in actual local integration work, rather than an unexplained
equal-work dimensional penalty.

The GMRES comparison cannot use that assembly proxy. At the same size the
low-wave case takes 112/110 iterations on triangles/tetrahedra; the shifted
case takes 234/288, with restart 50 throughout. Dense orthogonalization depends
on rows, restart and position within each cycle, not only matrix nonzeros.
Longer repeats retain isolated 3D CPU-time spikes: coefficients of variation
are approximately 20% and 11%, despite unchanged iteration counts, residuals
and coefficient errors. The cause of those spikes remains unresolved; precise
3D solver timing baselines and optimization claims are withheld. Raw initial
measurements, longer repeats and traces are retained separately under
`build-convergence/certification/helmholtz-gmres-20261013/`.

P0/P0g, higher degrees, curved geometry, boundary variants, direct solves and
PETSc/MPI policies are not covered by these scopes.
The benchmarks remain outside CTest and the publication workflow; only the
ordinary GMRES correctness tests are registered with CTest.

### Coupled reaction--diffusion protocol

`RodinCoupledReactionDiffusionAssemblyBenchmarks` and
`RodinCoupledReactionDiffusionSolveBenchmarks` are opt-in Assembly and Solver
targets. Their common workload retains two scalar trial fields and two scalar
test fields, rather than replacing the coupled problem by a single scalar
reaction term. On the unit box the physical equations are

$$
-\kappa_i\Delta u_i+u_i+\alpha u_{1-i}=f_i,
\qquad \kappa=(1,2),\quad \alpha=0.2,\quad i\in\{0,1\}.
$$

Full manufactured Dirichlet traces are prescribed. Affine fields are selected
at degree one and quadratic fields otherwise from `ReactionDiffusionData`.
The real reaction matrix has eigenvalues 0.8 and 1.2; the complex-valued
extension retains the same Hermitian operator, not a non-Hermitian shift.
For real certified fields $u_i^0$ and sources $f_i^0$, the complex data are
defined analytically by

$$
u_i=\eta_i u_i^0,\qquad
f_i=\eta_i f_i^0+\alpha(\eta_{1-i}-\eta_i)u_{1-i}^0,
\qquad (\eta_0,\eta_1)=(1+i,2-i).
$$

The distinct complex multipliers exercise both components without assuming a
common phase. Original local kernels qualify all six operator terms. The two
loads retain their production integration path; comparisons between formulations
are supplemented by the analytic manufactured patch, not by an independent
load algorithm. The complete exact coefficient vector is gathered by public block
assembly; the component ordering is not fabricated in the benchmark.

The negative control drops both coupling terms but retains analytic sources
and prescribed traces. Its constrained system is denoted $(A_0,b_0)$; even
with unchanged physical sources, essential elimination can change $b_0$.
Let $c$ gather the independent scalar mass actions
$\alpha Mx_{\star,1-i}$ at the public test offsets, with essential output
rows masked using boundary DOF indices. Prescribed input coefficients remain
in this action. The control requires a finite nonzero $\|c\|_2$ and rejects
the omitted reaction through

$$
s_C=\frac{\|A_0x_\star-b_0\|_2}{\|c\|_2}>10^{-4}.
$$

Thus the sensitivity denominator measures the omitted physical term, not a
boundary-dominated load. The correct assembled patch and full operator/load
comparisons retain budgets $10^{-10}$ and $10^{-11}$, respectively. Public
CG solves use tolerance $10^{-13}$, a 20,000-iteration cap, independent true
residual budget $10^{-11}$ and coefficient budget $10^{-9}$. Zero initial
guesses are restored outside timing; iteration counts are checked after every
measured solve. Assembly, interpolation, field scatter and validation are
excluded from the solve scope. Full assembly includes block bookkeeping,
constraints and gathering the initial guess, but no solving.

Each target registers real/complex H1 degrees one through three on all seven
mesh geometries at 3, 5 and 9 points per axis: 126 operation selections.
Quadrature order 16 is held fixed. Cell, quadrature, per-field basis, local
entry, global row and nonzero counts distinguish the two diffusion and four
reaction terms. All 126 selections per operation passed in each optimized build:
504 assembly/solve selections across native and two-thread OpenMP builds. The
126 native real-valued assembly/solve selections match the original-header
controls exactly in all recorded numerical and work counters. Across both
builds the largest patch residual was $3.6\times10^{-15}$, the largest solve
true residual was $1.1\times10^{-13}$ and the largest joint coefficient error
was $1.1\times10^{-12}$. The omitted-coupling normalized residual stayed near
one, well above its rejection floor. Native/OpenMP iteration counts can differ;
cross-build bitwise equality is not claimed. These targets are not CTest tests
or publication-workflow additions.

### Coupled source attribution and dimensional scaling

Matched native real-valued degree-two measurements used five repetitions at
each of the three sizes, fixed iteration counts selected by a pilot, and
alternating original/repaired order between geometries. The constructor repair
is outside both timed operations, so these measurements do not establish a
speedup from that repair. At nine points per axis the repaired CPU medians and
work counts were:

| Quantity | Triangle | Tetrahedron |
| --- | --- | --- |
| Cells | 128 | 3,072 |
| Quadrature points per cell | 53 | 729 |
| Scalar basis functions per field | 6 | 10 |
| Global rows | 578 | 9,826 |
| Matrix nonzeros | 9,044 | 319,808 |
| CG iterations | 94 | 103 |
| Full assembly CPU median | 3.610 ms | 1,965.046 ms |
| Fixed-system CG CPU median | 0.581 ms | 26.170 ms |

For cell count $N_T$, quadrature count $N_q$ and scalar local basis size $N_b$,
the local-entry work proxy $N_TN_qN_b^2$ grows by approximately 917, compared
with an observed assembly CPU ratio of 544. This proxy does not account for
different expression costs, load terms, geometry evaluation or insertion;
it is not an operation count or a predicted timing. The solve CPU ratio is 45,
while nonzeros grow by 35 and iterations by $103/94$. Dimension alone therefore
does not isolate either cost difference.

Eight sequential Instruments captures covered assembly and CG, real and complex
fields, and triangles and tetrahedra at this degree and size. Required executing
assembly/solver frames select production consumers; original-reference kernels
are excluded from assembly selection. Matching untimed production warmups and
controls can remain selected. Inclusive sampled CPU categories were:

| Operation / field | Triangle | Tetrahedron |
| --- | --- | --- |
| Assembly / real: quadrature | 87.9% | 97.7% |
| Assembly / complex: quadrature | 91.1% | 98.3% |
| CG / real: sparse products | 82.2% | 92.0% |
| CG / complex: sparse products | 75.8% | 89.9% |

Geometry contributed 11.5–23.2% of selected assembly samples, overlapping the
quadrature category. Complex multiplication contributed 35.9–47.1% of complex
assembly samples and 63.4–78.2% of complex CG samples, overlapping their owning
kernels. These observations locate source work; they neither prove an
optimization nor justify changing scalar semantics. Sample percentages are not
elapsed-time fractions and must not be added.

CPU coefficients of variation in the matched measurements were at most 2.9%,
but selected elapsed-time coefficients reached 22.2%. The elapsed/CPU discrepancy
has not been causally resolved; elapsed baselines are withheld rather than
attributed to scheduling, allocation or memory waits without further evidence.
Remaining work includes construction/lifecycle isolation, the other physical
contexts, higher degrees, curved geometry and PETSc/MPI counterparts.

Complex block-problem construction initially failed because `ProblemUsBase`
created private sequential tuple assemblers even though execution already
uses the canonical `Assembly::Default` problem assembler. The unused tuple
assemblers restricted scalar/storage types independently of that active path.
The repair removes their aliases, allocations, copy/move handling and members;
the active full/targeted assembly and solve code is unchanged. An original-header
overlay demonstrates the compile regression for complex sparse and real dense
systems. Mixed P0/H1 unit cases cover real/complex, sparse/dense storage, three
H1 degrees and seven geometries, including logical block offsets, independent
mass blocks, repeated assembly, true residuals and component solution scatter.
No alternate complex assembler, numerical tolerance change or speedup claim
is introduced by this construction repair.

### Problem construction and assignment lifetimes

`RodinProblemLifecycleBenchmarks` belongs to Variational and is opt-in, outside
CTest and benchmark publication. It uses the same degree-two real/complex
coupled reaction–diffusion formulation as the assembly and CG targets, with
sparse and dense native storage. Seven geometries, three sizes (2, 3 and 4
points per axis), four field/storage combinations and six lifetimes provide
504 selections. All 504 passed in each optimized native/OpenMP build (1,008
operation/build selections). This is degree-two native lifetime coverage, not qualification
of the remaining degree, vector-range or PETSc/MPI matrix.

Empty construction includes destruction; Defined additionally clones the
formulation body. Copy copies an already assembled body and system, including
their destruction. Move transfers a fresh untimed staging copy and times result
destruction. CopyAssignment and MoveAssignment retain an existing destination;
copy assignment can reuse fixed-size system buffers, whereas move assignment
replaces the previous value. Staging copies and moved-from cleanup are excluded.
The allocator is warmed; none of these scopes is process-first allocation.
Untimed guards require exact operator/load reproduction after reassembly and
reject a wrong body through logical integrator counts. Assignment also checks
the actual final timed destination before and after reassembly. Source spaces
and mesh setup, interpolation and PDE assembly remain outside timing.

The original four-lifetime executable is preserved locally. All its 336
selections passed. The applicable repaired numerical/work counters match those
controls exactly; object size is intentionally different. In the native Apple
LLVM build, removing two unused
linear forms and four unused bilinear forms reduces the two-field object from
1,888 to 720 bytes for sparse storage and from 1,648 to 672 bytes for dense
storage. Those forms contained initially empty matrices/vectors: this is not
removal of dense per-block matrix allocations or a total-memory measurement.

Matched five-repetition native real/sparse measurements at all three sizes on
triangles and tetrahedra observed empty-construction CPU medians of
0.527–0.546 microseconds before and 0.232–0.235 microseconds after repair.
The ratio was 2.25–2.35. Definition and assembled-copy changes were smaller;
no complete-workflow or solver speedup follows from the construction result.
Four Instruments captures compare empty construction and assembled copy at
fixed iteration counts. Residual empty-construction samples locate allocation,
deallocation and identity work; samples do not count allocations. Copy samples
remain dominated by memory copying (47–48%) and owned cloning (23–24%), with
overlapping categories. Their similar sampled CPU totals do not establish a
whole-copy gain. Matching untimed preparation can remain selected; sampled
CPU weights are not elapsed-time fractions.

Move timings include a substantial timer-transition floor. An offline control
of the same two pause/resume pairs measured a CPU median of 1.258 microseconds,
with 1.6% variation. It is not a Rodin benchmark or a CI test. Raw move numbers
are retained, but are not reported as isolated move latency; subtracting this
control would not establish a precise kernel time. Batched instrumentation
remains a separate measurement improvement before fine-grained move claims.
The installed Boost 1.88 bimap also exposes copy construction/assignment, not
move overloads: the logical UUID maps are copied even when passed through
`std::move`. Thus a public problem move is not a constant-time pointer transfer.
This source finding is separate from timer overhead; batched move attribution
and map-lifetime measurements are the next isolation step, not an adopted
container replacement.

The assignment repair retains destination identity, transfers the owned body,
system, borrowed function bindings and execution policy, and initializes mixed
offsets from logical space sizes at construction. Only complete-problem
assembler bases support identity-preserving assignment; FormLanguage root
assignment remains deleted. The original assignment fails to compile, so there
is no valid original assignment timing or assignment speedup claim. Logical
regressions reject stale bodies and zeroed cumulative offsets, alongside exact
reassembly and initial-guess gathering across P0/H1 degrees one through three,
all seven geometries and real/complex sparse/dense storage. Remaining work
includes other degrees/ranges, backend-specific handle lifetimes and the
unmeasured physical-context matrix.

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
