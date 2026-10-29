# Library performance benchmarks

This module measures the cost of library operations used by finite-element
assembly and solution. Its hierarchy follows the production modules, from
spatial arithmetic and reference data to mapped geometry, element kernels,
global assembly and solvers. Interpolation and projection are measured as
library operations in their own right.

The purpose is to locate work, distinguish its mathematical and implementation
dependencies, and evaluate proposed repairs without changing the discrete
problem. Local measurements are observations for a specified configuration,
not architecture-independent performance guarantees. This README therefore
defines measurement scopes and coverage gaps rather than recording timings,
speedup rankings or the chronology of individual investigations.

## Module architecture

```text
tests/benchmarks/
├── Workloads/                 Shared forms, inputs and numerical oracles
├── Rodin/
│   ├── Math/                  Spatial arithmetic and sparse algebra
│   ├── QF/                    Reference quadrature
│   ├── Geometry/              Connectivity, transformations and mapped points
│   ├── Variational/
│   │   ├── P1/                Low-order space-specific operations
│   │   └── H1/                High-order basis and local quadrature operations
│   ├── Solid/                 LinearElasticity element kernels
│   ├── Assembly/              Assembly stages and complete physical operators
│   ├── Solver/                Fixed-system solves and preconditioning
│   ├── Location/              Spatial search
│   ├── IO/                    Mesh input
│   ├── MPI/                   Distributed module specializations
│   └── PETSc/                 PETSc module specializations
└── Stokes/
    └── InfSupStability/        Finite-mesh pressure-spectrum profiling
```

Placement follows the operation being measured, not the physical context of
its caller. A spatial-matrix product belongs to Math; an element integral
belongs to Variational or Solid; global insertion belongs to Assembly or its
backend specialization. Shared workloads contain forms and reference checks,
not a second assembly policy. Qualified `Rodin/...` includes resolve to
production headers before mirrored benchmark fixtures.

Module and component executables are opt-in (`EXCLUDE_FROM_ALL`) and are not
registered with CTest. Every enabled executable must have an explicit benchmark
registration; configuration fails for an omitted registration. The old
`RodinBenchmarks` aggregate and its historical dashboard are retired. Its source
cases remain in their owning component executables. Benchmark sources and
measurements do not enter the production-coverage denominator or test lanes.

The pressure-spectrum profiler under `Stokes/InfSupStability/` measures verification
support, not a production Stokes solver. A root-owned spectral calculation and
its result communication must not be interpreted as distributed assembly
imbalance. Such measurements do not count as production solver coverage.

## Registration and display

An executable is declared in its owning `CMakeLists.txt`. The adjacent
`registrations.benchmarks.cmake` records its public path, authored description,
report format, timing policy, execution profile, ranks and case-batch limit:

```cmake
rodin_register_benchmark(RodinMathSpatialMatrixBenchmarks
  PATH Rodin/Math/SpatialMatrix
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/SpatialMatrix.html"
  FORMAT googlecpp PROFILE slow METRIC cpu_time RANKS 1
  CASE_BATCH 256)
```

`PATH` follows the owning module/class or physical context. `DESCRIPTION` is an
HTML fragment explaining the operation, mathematics, algorithm, correctness
checks and timing boundary. Inline mathematics uses `\(...\)`; descriptions
become more specific down the hierarchy. Parent descriptions can be supplied in
`tests/benchmarks/<publication path>/index.html`. Do not commit generated timing
reports, investigation files or performance rankings.

`PROFILE` is `regular` or `slow`. Slow is a scheduling variant under the same
owner, not a separate module. `METRIC cpu_time` and `METRIC real_time` select a
Google Benchmark field explicitly; `auto` follows the driver's `UseRealTime` or
`UseManualTime` declaration, otherwise CPU time. Manually supplied values remain
labeled `manual_time`: their scientific meaning comes from the driver and is
not assumed to be wall time. `RANKS` preserves the required launcher counts.
`ARGS` applies to every invocation; `MPI_ARGS` additionally applies above one
rank.

The Stokes assembly, SparseLU and GMRES/block-Schur registrations use `CASE_BATCH 1`: individual
large cases can take many minutes outside their reported timed operations.
Separate lanes retain every case while limiting work lost to a runner shutdown.
The [runner resource budget](../README.md#github-runner-resource-budget) applies
to the complete invocation, including untimed setup and numerical qualification.

Stokes assembly qualification checks the exact affine patch residual and both
wrong-field and missing-pressure sensitivity controls without constructing a
solver. Every repeated assembly retains its exact-field residual check. The
solver workloads additionally require solver success, coefficient reproduction
and the true solved residual; all existing thresholds are unchanged.

`RodinStokesGMRESBenchmarks` uses the production native
`Solver::GMRES<System>` API, constructed with a configured
`StokesBlockPreconditioner`. The preconditioner type is retained inside the
solver's private Eigen backend. Its benchmark-owned
preconditioner applies incomplete Cholesky to the velocity block and an LDLT
pressure mass solve for the Schur approximation, including the explicit mean
constraint. Pressure mass storage has shared immutable ownership across solver
configuration copies; factors remain independent. Restart is 60, the iteration
limit is 1200 and the backend tolerance is 1e-13. Each solve starts from zero,
and timing includes preconditioner setup and the solve. Pressure mass assembly
is untimed. This variant has no residual-correction solve, unlike SparseLU.
Its 108 cases are separately labelled; SparseLU cases are not replaced.

The largest pyramid cases currently fail the existing missing-pressure
sensitivity floor despite passing the affine patch and linear-solver checks.
Their failure remains explicit; changing the solver does not waive that check.

Non-root display reporters must redirect their output stream as well as suppress
measurement callbacks: Google Benchmark writes `--benchmark_list_tests` output
directly to that stream. Only the root catalogue should reach launcher stdout;
every rank still participates in measurement and its numerical checks.

`CASE_BATCH` is the maximum number of discovered cases per execution lane,
not a sample reduction. The single-argument selector-size budget may split a
batch further. Expensive spectral cases use one case per Slow lane.

The configured `RodinBenchmarkTargets-Release.tsv` contains executable and
library targets. Documentation-only configuration reads the same registrations
into `doc/RodinBenchmarkDescriptions.tsv` without building executables. The
rendered `rodin-module-catalogue.json` supplies the exact API module prose,
representative formulas/code and module groups; no separate catalogue is
maintained. Lists are alphabetical. Physical contexts are siblings of `Rodin/`;
the public label for that source tree is **Modules**. There is no convergence
benchmark category: the pressure spectrum is a finite-mesh stability diagnostic.

Navigation places each measured operation directly under its owning module or
class. There is no Pipeline category. For example:

```text
Modules / Assembly / Mass
Modules / Math / SpatialMatrix / MatrixProduct
Modules / MPI / Variational / Construction
Modules / PETSc / Assembly / Mass
Modules / Solver / Stokes
Stokes / InfSupStability
```

Class-level children remain where they distinguish operations, such as matrix
products and contractions. Complete assembly and solve workloads sit directly
under Assembly and Solver. InfSupStability contains its spectral algorithms as
sections on one benchmark page.

## Shared builds, execution and pooling

The main-branch `Benchmarks` workflow runs this graph. A manual dispatch can
validate a branch. Numerical tests, sanitizer and production-coverage workflows
remain separate.

```text
Benchmarks
├── Infrastructure checks
├── Documentation ── Doxygen XML + m.css + benchmark descriptions
│   └── Publish documentation immediately
│       └── Retain the latest complete benchmark histories
├── Local / Release / GCC 14
│   └── Shared libraries (Rodin + Google Benchmark, compiled once)
│       ├── Modules / Assembly ── shared component build ── Regular/Slow lanes
│       ├── Modules / Geometry ── shared component build ── Regular/Slow lanes
│       ├── Modules / IO       ── shared component build ── Regular/Slow lanes
│       ├── Modules / Location ── shared component build ── Regular/Slow lanes
│       ├── Modules / Math     ── shared component build ── Regular/Slow lanes
│       ├── Modules / QF       ── shared component build ── Regular/Slow lanes
│       ├── Modules / Solid    ── shared component build ── Regular/Slow lanes
│       ├── Modules / Solver   ── shared component build ── Regular/Slow lanes
│       ├── Modules / Variational ─ shared component build ─ Regular/Slow lanes
│       └── Stokes ── shared component build ── InfSupStability / Slow
└── MPI and PETSc / Release / GCC 14
    └── Shared libraries (compiled once for this backend configuration)
        ├── Modules / MPI / Variational ── shared build ── NP2, NP4 lanes
        ├── Modules / MPI / Location ── shared build ── NP2, NP4 lanes
        └── Modules / PETSc / Assembly ─ shared build ── NP1, NP2, NP4 lanes
                           │
                           ▼
                 Validate complete campaign
              exact discovery + every raw report
                           │
                           ▼
                 Publish complete benchmark histories
                  (after documentation publication)
```

Documentation deploys as soon as its own build and navigation checks pass;
benchmark failures cannot block that deployment. Fresh pages reconnect to the
last published complete histories, retaining their measured commit and date.
Once every benchmark contribution is validated, a second deployment updates
the histories without replacing the documentation. Both publishers and Static
Analysis use the same Pages publication lock and verify their deployed commit.

CMake's configured dependency graph selects only libraries needed by that
configuration's benchmark targets, including linked third-party libraries.
Unrelated configured libraries are not built. Each group starts after its
configuration's common build. Its execution lanes
start immediately after **that group's** component build; there is no barrier
waiting for sibling groups. The same compiled binary serves all its case
partitions and rank counts. Runtime jobs restore an immutable archive and do no
compilation. Archives validate commit, full configuration, compiler version,
architecture, checksums and runtime dependencies. Stable source mtimes prevent a
fresh checkout from rebuilding the restored library objects.

The two Release configurations distinguish native and MPI/PETSc compilation.
OpenMP and numerical-library thread counts are one in this campaign. These are
not replacements for the test workflow's configuration matrix. Within a lane,
executables run sequentially. Separate hosted runners can execute independent
lanes concurrently. Google Benchmark's existing iteration/repetition policies
and all discovered cases are retained; the runner requests JSON output. The
AABB workload retains its CSV driver and exposes its exact configuration list
with `--list`.

Every lane uploads its raw report and a checksummed receipt. Pooling requires
both configuration plans, every build group and every discovered case at every
registered rank. Missing, duplicate, skipped, failed, stale or altered
contributions fail the campaign. Raw Google Benchmark repetitions, counters, context and metrics remain
available; chart values use the median of iteration records, excluding aggregate
statistics such as the coefficient of variation. The CSV driver already emits
its paired-block medians; its internal per-block samples are not exported. Different configurations,
metrics and ranks are never averaged together. The final complete package also
checksums every component result file.

## Serving the pages

Doxygen/m.css generates the benchmark section with the same themes and sidebar
as the API documentation. `RodinDoxygen` requires no numerical build for these
pages. For example:

```sh
cmake -S . -B build-doc -DRODIN_BUILD_SRC=OFF -DRODIN_BUILD_DOC=ON \
  -DRODIN_USE_MCSS=ON -DRODIN_INSTALL_RESOURCES=OFF
cmake --build build-doc --target RodinDoxygen
```

The publisher waits for both documentation and the complete benchmark package,
then makes one GitHub Pages commit. For `develop`, navigation is served at
`docs/refs/heads/develop/benchmarks/index.html`; a leaf is, for example,
`docs/refs/heads/develop/benchmarks/Rodin/Math/SpatialMatrix/index.html`.
Per-component, per-configuration histories live separately under
`benchmark-data/refs/heads/develop/`, so regenerating documentation cannot erase
them. A rerun replaces the same commit's sample rather than duplicating it.
The generated legacy `benchmarks/` dashboard is removed only during successful
publication of the complete replacement. Earlier aggregate history is not
imported into the new component series.

A leaf offers configuration and case selectors, then a Chart.js timing history.
It loads only the selected case’s series and retains the latest 20 commits per
case. This bounded window limits Pages storage; complete campaign reports remain
CI artifacts. Histories are sharded by case to avoid monolithic class-sized JSON
files, and configurations have separate indexes.
Tooltips identify date, run and sample count; points link to source commits.
Unmeasured pages say so explicitly. A local documentation build contains the
navigation and scientific descriptions without fabricated timings.

## Adding a benchmark or lane

1. Add the driver under its owning module. Keep independent class benchmarks
   separate from whole-pipeline workloads. Retain numerical guards outside the
   timed region and describe all included work.
2. Add an opt-in executable and link its actual dependencies. Include the
   adjacent `registrations.benchmarks.cmake` from that directory's CMake file.
3. Register a unique path and description, format, metric, profile, ranks and
   batch size. Optional MPI/PETSc executables stay in their backend CMake guards;
   their descriptions remain visible in documentation-only builds.
4. Configure and build the target, list its cases, then run a bounded selection.
   Run one numerical benchmark process at a time. Test all required MPI ranks
   sequentially.
5. Run the infrastructure checks below and build documentation. A new target
   automatically joins its owning group's shared build and discovered execution
   lanes. To split a slow suite further, reduce `CASE_BATCH`; no workflow matrix
   duplication is required. MPI rank changes require an explicit registration
   change, not a workflow-only selector.

To create a new top-level module or context, register its path and add parent
prose. `group_path()` in `dev/benchmark_ci.py` controls the common hierarchy cutoff:
`Rodin/<Module>`, `Rodin/MPI/<Module>`, `Rodin/PETSc/<Module>`, or a physical context.
Only split this cutoff when a group exceeds the hosted matrix limit; do not
repeat common compilation to create execution parallelism.

The workflow implementation is in `BenchmarkConfiguration.yml` and
`BenchmarkGroup.yml`; pooling, publication and site generation are in
`dev/benchmark_ci.py`, `dev/benchmark_publish.py` and `dev/benchmark_site.py`.

Intermediate benchmark pages present their children as compact linked blocks,
with two columns when the content area is wide and one on narrow screens.
Each block reuses the opening paragraph of its child's authored description;
keep that paragraph a concise explanation of the measured operation. The full
mathematics and timing contract remain on the destination page.

Benchmark pages use the same content grid, heading sizes and sidebar placement
as the API documentation. Its top bar reuses the rendered documentation navigation,
including dropdowns, responsive controls and theme selection. Only search changes
to a combined benchmark-and-search icon control. Both search icons
open the same search overlay and share hover and keyboard focus behavior. The sidebar supplies
section links and scroll-following highlights. Its controls retain the chosen
side across pages without shifting the content column. Visibility follows the
same viewport rules as the API sidebar. Clicking a subsection scrolls to that
heading. Navigation updates do not steal focus.

Benchmark search reuses the pinned m.css engine and overlay. Its generated index
contains module, collection and benchmark names, paths, CMake target names and
authored descriptions. Exact names rank first; description queries may contain
several words, and slash or C++-style scope paths are accepted. Each result shows
its owning path and opening description, with no separately maintained catalogue.
The index works in a local file preview as well as on Pages.

The search icon, Tab or T opens the overlay; arrow keys select a result, Enter
opens it and Escape closes it. A visible close button supports touch screens.
On published measurement pages, **Filter cases** narrows the case selector by
operation, geometry, degree or other case-name fields. The selected configuration,
case and active filter are reflected in the URL and restored when reopened.

Initialize the search dependency before running site-generation checks:

```sh
git submodule update --init third-party/m.css
```

Run infrastructure checks sequentially:

```sh
python3 -m unittest discover -s dev -p 'test_benchmark_*.py' -v
node --check dev/benchmarks/charts.js
node dev/test_benchmark_navigation.cjs build-doc/doc/html/benchmarks
node dev/test_benchmark_charts.cjs build-doc/doc/html/benchmarks
node dev/test_benchmark_search.cjs build-doc/doc/html/benchmarks
```

The browser check requires the `playwright` Node package and a Chromium install.
`RODIN_BROWSER_EXECUTABLE` can select an installed Chromium-compatible browser.
It checks destination icons, sidebar positioning, direct subsection links,
manual-scroll highlights, theme selection and narrow-screen layout. Run one
numerical benchmark at a time, retaining the registered MPI ranks for a
distributed benchmark.

## Spatial vector and matrix arithmetic

`RodinMathSpatialVectorBenchmarks` and `RodinMathSpatialMatrixBenchmarks` compare the
public spatial APIs (including `Math::sum`, `Math::minus`, `Math::mult`,
`Math::division`, `Math::conj` and `Math::dot`) with bounded-capacity Eigen values, dynamic Eigen vector/matrix representations, and selected compile-time fixed Eigen values.
`SpatialArithmetic.h` owns the common input, execution and numerical guards.
Each opt-in executable has its own component page and case histories.

The vector suite covers zero construction, copying, arithmetic and compound
assignments, scalar multiplication in both orders, division, conjugation,
transpose, dot conventions, ordinary/stable/Blue/p-norms, normalization and
cross products. The matrix suite covers all supported rectangular shapes,
arithmetic, norms, contractions, conjugation, transpose/adjoint, row/column
extraction, matrix/vector products, trace, determinant, inverse, solve and
full-rank pseudoinverse. Empty shapes are included where the operation is
defined; solve/inverse require nonempty square matrices, row/column extraction
requires nonempty shapes, and normalization requires nonempty vectors.
The `ScalarDot` comparator uses a transparent scalar loop and the same right
conjugation as `MathDot`. It supplies an algorithm comparison without changing
production arithmetic. Fixed comparisons cover selected operations on dimensions one through three
(square matrices). Rectangular and empty shapes retain runtime extents in the
other three representations. Each case uses either one or 32 independent inputs.
The matched `MatrixProduct` cases multiply $r\times k$ and $k\times r$
operands with independently populated entries. The existing
`RodinSpatialMatrixProductBenchmarks` supplies the complete $r,k,c$ shape
sweep, including mixed scalar fields; these scopes complement one another.

The names identify different mathematical contracts explicitly:

| Operation | Definition | Implementation distinction |
| --- | --- | --- |
| `ScalarDot` | $\sum_i u_i\overline{v_i}$ | Explicit scalar-loop comparator |
| `MathDot` | $\sum_i u_i\overline{v_i}$ | `Math::dot`, including Eigen-backed values |
| `RightDot` | $\sum_i u_i\overline{v_i}$ | Spatial vector/matrix member dot, or Eigen vector `v.dot(u)` |
| `BilinearProduct` | $\sum_i u_i v_i$ | Explicit entrywise reduction without conjugation |
| `Cross` | Ordinary bilinear cross product | Eigen's complex cross is conjugated back to match this contract |
| `MatrixProduct` | $C=AB$ | Ordinary product without conjugation |
| `Adjoint` | $A^*=\overline{A}^{\mathsf T}$ | Distinct from transpose |

Real and complex inputs have nonuniform bounded coefficients. Independent
scalar-loop oracles check contractions, products and conjugation. Full-pivot LU
checks solves/inverses, and complete orthogonal decomposition checks full-rank
pseudoinverses outside timing. NaN padding detects inactive spatial reads;
perturbation checks verify that the guard rejects changed results. All outputs
are checked before and after timing; numerical failures return a nonzero status.

Timing includes an owned result and assignment into retained output. Compound
assignments include copying their left value; dynamic result allocations are
included. Setup, reference calculations and array allocation are excluded.
The `BoundedDestination`, `DynamicDestination` and `FixedDestination` variants
instead evaluate safe arithmetic/product expressions directly into prepared,
disjoint output via Eigen's `noalias()` assignment. They include expression
construction and evaluation but avoid an explicitly materialized intermediate.
They are a different output-lifecycle contract from the owned-result cases;
compare them separately. Spatial operations always produce their public owned
result. This distinguishes allocation/temporary costs from arithmetic costs.

Pseudoinverse comparators retain the same normal-equation algorithm, rather
than comparing it with a more general rank-revealing algorithm inside timing.
Results report calling-thread CPU time per sweep: divide by the `operations`
counter for time per value. Counts one and 32 expose barrier/loop overhead.
Fixed shapes give the compiler additional information and must be identified
when comparing them with runtime dimensions. Run optimized, uninstrumented
builds for performance conclusions and retain repetition/dispersion data outside
the repository. No timing or investigation reports belong in version control.

```sh
cmake --build build --target RodinMathSpatialVectorBenchmarks RodinMathSpatialMatrixBenchmarks
build/tests/benchmarks/RodinMathSpatialVectorBenchmarks \
  --benchmark_filter='/rows:3/count:32$' \
  --benchmark_min_time=0.1s --benchmark_repetitions=5
build/tests/benchmarks/RodinMathSpatialMatrixBenchmarks \
  --benchmark_filter='/rows:3/cols:3/count:32$' \
  --benchmark_min_time=0.1s --benchmark_repetitions=5
```

Run these commands sequentially. The existing scalar-dot, spatial-dot,
matrix-contraction, matrix/vector-product, matrix-product and vector-storage
executables retain mixed scalar orders, basis-sized batches and storage
lifecycles. They complement the matched API suite; a primitive speedup still
requires confirmation in its assembly consumer.

## Sparse elimination

`RodinSparseEliminationBenchmarks` measures the public native sparse
`LinearSystem::eliminate` operation on P1 and degree-three stiffness-plus-mass
systems across Segment, Triangle, Quadrilateral, Tetrahedron, Hexahedron,
Pyramid and Wedge meshes. Sizes range from roughly one thousand to one hundred
thousand DOFs (the segment sweep also includes smaller systems). Random and
contiguous subsets prescribe 1%, 2%, 5%, 10%, 40%, 50%, 80%, or 90% of the
coefficients. The `placement` argument distinguishes random (`0`), contiguous
(`1`), and the complete assembled boundary (`2`, with `percent:0`). Boundary
indices come from `DirichletBC`; their values are coefficients of the exact
algebraic vector, not samples of an analytic boundary function.

The `Direct` family constructs directed upper-band CSC matrices of 1,000,
10,000 and 100,000 rows, with widths 8 and 64 and 1–10% scattered constraints.
The `missing:1` cases also omit columns and diagonals. They exercise the
absence of transpose storage and include diagonal insertion in the timed
operation. These are algebraic systems, not manufactured PDE convergence tests.

Assembly, subset construction, independent projector guards and matrix/load
restoration are excluded. FEM calls retain their compressed allocation. Direct
calls receive an exact-capacity copy so structural growth cannot disappear
through allocation reuse; copying is outside timing. All algorithm-owned
scratch and storage growth are inside timing. Manual timings report calling-
thread CPU seconds, including clock-read overhead; they are not wall-time or
whole-process throughput measurements. Each result is checked against the
matrix/load projector formula and the exact algebraic residual. Numerical guard
failures make the executable return failure.

The target is opt-in, serial, and not registered with CTest. For example:

```sh
cmake --build build --target RodinSparseEliminationBenchmarks
build/tests/benchmarks/Rodin/Math/RodinSparseEliminationBenchmarks \
  --benchmark_filter='SparseElimination/Direct/dofs:10000/width:64/percent:5/missing:1' \
  --benchmark_repetitions=7
```

Select a bounded family before running under a sampling profiler. Run one
benchmark process at a time. Use the same binary configuration, constraints,
coefficient restoration and timing scope for comparisons; retain repetition
and dispersion information on shared machines. The mathematical elimination
contract and algorithm are documented in the boundary-conditions guide and
in the sparse `LinearSystem::eliminate` Doxygen block.

## Mathematical measurement model

Let $\mathcal{T}_h$ be a mesh of $\Omega\subset\mathbb{R}^{D}$, with
reference cells $\widehat K\subset\mathbb{R}^{d}$ and maps
$F_K:\widehat K\to K$. The embedding dimension $D$ and reference dimension
$d$ are distinct. Let $V_h,W_h$ be finite-element trial and test spaces over
$\mathbb{F}\in\{\mathbb{R},\mathbb{C}\}$, with scalar, vector or matrix
values as supported by the operation. Their local basis counts are $N_u,N_v$.

For a quadrature rule $\{(\widehat x_q,w_q)\}_{q=1}^{N_q}$, a local
bilinear operator is represented by

$$
A^K_{ba}=\sum_{q=1}^{N_q}w_q\mu_K(\widehat x_q)
\mathcal{B}\bigl(\phi_a,\psi_b;c\bigr)(F_K(\widehat x_q)),
\qquad
\mu_K=\sqrt{\det(J_K^{\mathsf T}J_K)},\quad J_K=DF_K.
$$

Here $c$ denotes prescribed coefficients; $\mathcal{B}$ includes the required
pullbacks, derivatives and contractions. For full-dimensional cells the
measure factor is $|\det J_K|$. Point cells use the discrete point-measure
convention rather than this positive-dimensional derivative formula.

If $R_K^u,R_K^v$ restrict global coefficients to a cell, unconstrained assembly
has the algebraic structure

$$
A_h=\sum_{K\in\mathcal{T}_h}(R_K^v)^*A^K R_K^u.
$$

The adjoint $*$ is transpose-conjugation over $\mathbb{C}$. Constraint
expansion, insertion and completion are additional operations, not part of a
reference-basis microbenchmark. A component benchmark isolates one dependency
of this construction; an assembly consumer checks whether its cost matters to
the complete operation.

Comparisons record actual work, including $N_q,N_u,N_v$, component counts,
cells, global DOFs and stored nonzeros where applicable. A basis-value sweep
produces $N_qN_u$ scalar entries; a complete reference-derivative sweep produces
$dN_qN_u$. A general local matrix contains $N_vN_u$ entries. Equal mesh
parameters or polynomial degrees do not imply equal work across geometries.
Large dimensional differences require attribution to these counts and to
measured stages before they are explained as an architectural limitation.

## Physical contexts

| Context | Quantity measured and purpose | Important exclusions or gaps |
| --- | --- | --- |
| Mass, reaction and projection | Operator/load assembly and selected projection solves; distinguishes basis preparation, weighting and global insertion | Complete projection-solve/backend matrix, mixed/global couplings and component-storage attribution |
| Poisson, conductivity and reaction–diffusion | Diffusion/reaction operators and constrained manufactured problems; distinguishes geometry, coefficient and boundary-data evaluation | Vector/mixed constraints, broader boundary variants and distributed stage profiles |
| Complex Helmholtz | Complex assembly and fixed-system GMRES; separates operator construction from indefinite-system solution | Broader mapped/backend matrix and repeatable size-dependent solver attribution |
| Coupled reaction–diffusion | Block assembly and fixed-system solution; exposes coupling and scalar-promotion costs | Broader mixed spaces, backend parity and block-preconditioner stages |
| Stokes | Stable mixed velocity/pressure operators, pressure-mean coupling, local bindings and selected fixed-system solves | Wider degree/mesh/backend profiles; production solver timings remain distinct from the verification spectrum |
| LinearElasticity | Element and complete assembly kernels; measures symmetric-gradient/constitutive contractions | Higher degrees, complex and mapped/backend parity where supported |
| Nonlinear Poisson | Residual, tangent and source assembly at prescribed fields; isolates state-dependent kernels | Full Newton workflows, natural-boundary terms, PETSc/MPI parity and solver stages |

Darcy and NAFEMS remain deferred. Available mass benchmarks do not establish
completion of the remaining performance matrix.

### Coverage axes and applicability

The intended matrix distinguishes real/complex scalar fields; supported
scalar/vector/matrix ranges; P0, P0g, P1 and high-order H1 spaces; field degree;
reference-cell family; geometry degree and embedding; and native sequential,
OpenMP, PETSc-local and MPI execution. Reference-element drivers include H1
degree sweeps through six, while complete physical workloads generally cover
a smaller degree range. These are not interchangeable coverage claims.

The eight reference families are Point, Segment, Triangle, Quadrilateral,
Tetrahedron, Hexahedron, Pyramid and Wedge. Point algebra and basis checks do
not constitute a positive-dimensional PDE experiment. Constant spaces do not
provide conforming diffusion approximations; mixed problems require applicable
stable pairs. Unsupported combinations are excluded with a mathematical reason,
not silently counted as measurements or described as unimplemented physics.
A real-scalar PETSc installation does not certify complex PETSc behavior.

Mesh-based scaling requires at least three sizes; reference-element scaling
uses at least three quadrature orders. Field degree, geometry degree and
quadrature order remain independent inputs. Fixed quadrature defines a discrete
operator, not an exactness or convergence guarantee. In particular, rational
pyramidal bases and nonaffine maps require their own integration analysis.

## Measurement boundaries: what, how and why

Each benchmark defines its input state, timed operation and excluded work in
the owning fixture documentation. The following distinctions are essential.

| Scope | Measurement boundary | Interpretation |
| --- | --- | --- |
| Construction/lifecycle | Creation, copy, move, assignment or destruction of a named object | Separates setup/ownership costs from evaluation; construction is not implicitly included in a warm kernel |
| Reference preparation | Basis/rule construction or production at specified reference samples | Isolates degree, shape and quadrature dependence without global scatter |
| Mapped preparation | Coordinates, Jacobians, pullbacks, measure and cell binding | Exposes geometry-dependent work and distinguishes reference reuse from physical-state reuse |
| Local kernel | Preparation and/or entry evaluation for one cell, as explicitly named | Measures contractions before backend insertion; complete binding and repeated-entry reads are distinct scopes |
| Global assembly | Operator/load replacement, with setup/constraints/completion included only when stated | Determines whether an isolated change benefits its actual assembly consumer |
| Fixed-system solve | Solve on an already assembled operator, with preconditioner setup/reuse explicitly specified | Prevents assembly cost from being attributed to the solver or reuse from being confused with first setup |
| Interpolation/projection | Public DOF-functional interpolation or mass-system projection | Distinguishes coefficient construction from solving an approximation problem |

Cold/first-use, retained-reference and repeated-use states are separate
experiments. Changes of cell, rule, coefficient and geometry must remain visible
to the operation. Prepared-reference measurements do not represent ordinary
unprepared public calls. Setup, reference integration and numerical checks are
excluded from timing unless their inclusion is the stated measurement scope.

For distributed scopes, wall time is reported using the maximum participating
rank when the driver implements that protocol. Alignment and timing reduction
are excluded from the operation interval and documented separately. Owned and
halo work, empty shards and imbalance metadata are distinguished. No collective
is introduced into production code merely to simplify a local benchmark.
Strong scaling fixes global work; weak scaling fixes work per rank. Neither is
established by changing rank count and problem size without declaring which
quantity is held fixed.

## Supporting documents

[Measurement procedures](methodology.md) contain qualification, comparison and
execution instructions. The [scope inventory](coverage-ledger.md) tracks
available measurements and remaining work separately from the module design.
