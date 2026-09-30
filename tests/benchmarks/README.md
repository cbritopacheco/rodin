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

Before timing, the matrix $A$ is assembled twice. The relative energy error
$|z^TAz/a(u,u)-1|$ must be below $10^{-9}$, where $z$ represents the oracle
field. Repeat assembly must satisfy
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

## Remaining performance workplan

| Extension | Required measurement or evidence |
| --- | --- |
| Extended PETSc workload/constraint parity | Loads, boundary constraints, and independent operator-action comparisons beyond affine energy; retain maximum-rank timing and ownership metadata |
| Complex Helmholtz, coupled reaction–diffusion, Taylor–Hood Stokes, nonlinear Poisson | Loads/block forms and residual/tangent assembly at a prescribed state; supported scalar/backend paths; stable mixed spaces |
| P0/P0g real/complex scalar/vector mass and projection forms, higher H1 orders | Global and local integrator coverage; meaningful point/0D forms where supported; explicit mathematical exclusions |
| Cold/setup and stage isolation | Mesh, space, sparsity/allocation, kernels, insertion, constraints, finalization, solve, and norm timings reported separately |
| Curved geometry and boundary variants | Map degree/regularity, quadrature-point counts, constraints and boundary workload metadata |
| Scaling and regression baselines | At least three sizes; fixed-global-work strong scaling and fixed-per-rank-work weak scaling; controlled runner and established variance before thresholds |

All meaningful combinations of existing contexts, supported spaces, the seven
positive-dimensional geometries, and assembly backends remain the intended
coverage matrix. Unsupported combinations are stated explicitly, not counted
as completed measurements. Darcy remains deferred. Any subsequent optimization
requires independently checked operators/loads/constraints and baseline-
equivalent solution behavior, residuals, and solver iterations.
