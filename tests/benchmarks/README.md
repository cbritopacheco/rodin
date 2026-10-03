# Assembly performance methodology

The assembly benchmarks measure computational cost separately from the
multi-resolution accuracy assertions in `tests/convergence`. The performance
track is developed independently of the convergence-certification PR; both
target `develop`. No assembly optimization is introduced by the initial suite.

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
across binds, use the same mesh-owned mapped quadrature, and exclude setup and
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

## Remaining performance workplan

| Extension | Required measurement or evidence |
| --- | --- |
| Extended workload/constraint parity | Full constrained scalar Poisson/conductivity/reaction--diffusion assembly is implemented; extend to vector physics, mixed boundary conditions and identification constraints |
| Complex Helmholtz, coupled reaction–diffusion, Taylor–Hood Stokes, nonlinear Poisson | Loads/block forms and residual/tangent assembly at a prescribed state; supported scalar/backend paths; stable mixed spaces |
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
