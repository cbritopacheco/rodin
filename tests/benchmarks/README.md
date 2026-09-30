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

## Remaining performance workplan

| Extension | Required measurement or evidence |
| --- | --- |
| PETSc sequential/OpenMP/MPI for the implemented physics | Same forms and numerical oracles; support checks; MPI maximum-rank wall time; explicit ownership and partition imbalance |
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
