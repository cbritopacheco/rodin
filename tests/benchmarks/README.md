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
remain outside this executable's coverage. CI exercises all 168 smallest-mesh
stage cases without timing thresholds.

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

Dedicated PR smoke jobs run local and MPI (1–4 ranks) smallest-mesh cases in
sequential/OpenMP builds and archive JSON without timing thresholds. Full
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
are sequential kernel measurements even in an OpenMP build. CI checks all
336 smallest-mesh registrations without timing thresholds.

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

CI runs all 672 smallest-mesh Eigen cases and all 336 native-scalar PETSc
cases locally and at 1–4 MPI ranks, in both OpenMP configurations. These are
numerical smoke gates with no timing thresholds. Larger-size runs should use
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
CI checks every smallest-mesh registration, including PETSc local and MPI
at one through four ranks in both OpenMP configurations. Small segment
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

Setup, interpolation, independent checks and state changes are untimed. Source
terms, natural-boundary integrals, constraints, SNES/KSP, and error norms are
excluded. PETSc/MPI parity and separate timings for those excluded operations
remain required before attributing a full distributed solve's cost to this
benchmark. Counters report cells, global native DOFs, nonzeros, field and geometry degrees,
quadrature order and points per cell. CI executes all smallest-mesh cases as
numerical checks, without timing thresholds.

The independent quadrature oracle owns mapped points only for its currently
bound cell, just as production assembly does. It does not request or evict
mesh-owned borrowed quadratures. This prevents untimed validation from retaining
all cells' quadrature points during timed assembly. Integer request-count
regressions enforce that lifecycle separately from numerical agreement.

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
