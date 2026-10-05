# PETSc conductivity h-convergence

On the unit box $\Omega=(0,1)^d$, $d\in\{1,2,3\}$, the problem is

$$
-\nabla\cdot(\gamma\nabla u)=f\quad\text{in }\Omega,
\qquad u=g\quad\text{on }\partial\Omega,
\qquad \gamma(x)=1+\sum_{j=1}^d x_j.
$$

The coefficient satisfies $1\le\gamma\le d+1$. The weak form is
$a(u,v)=\int_\Omega\gamma\nabla u\cdot\nabla v\,\mathrm{d}x$
for $v\in H_0^1(\Omega)$, with the exact trace prescribed by Dirichlet
elimination. Uniform positivity makes the constrained real scalar problem
coercive and permits CG. These tests use Rodin's PETSc assembly and solver
with either a local mesh or an MPI mesh. The local and distributed test
fixtures share the same manufactured data, forms, and assertions.

## Exact fields and references

The affine P1 patch is $u=1+\sum_j x_j$, with $f=-d$ and
$\nabla u=(1,\ldots,1)$. The quadratic P2 patch is
$u=1+\sum_j x_j^2$, with

$$
f=-2d\gamma-2\sum_j x_j,\qquad \partial_j u=2x_j.
$$

Both patches are exactly representable on the affine cell maps used here.
Their nonzero traces exercise lifting and boundary elimination as well as
the variable coefficient. They use `n=5` grid points per axis and require
both L2 and H1-seminorm errors below $10^{-9}$. This bound accounts for
floating-point assembly, quadrature, and iterative solution; it is not a
topological-index tolerance.

The smooth rate study uses $u=1+w$, where
$w=\prod_j\sin(\pi x_j)$, so $g=1$. Its source is

$$
f=d\pi^2\gamma w-\pi\sum_{j=1}^d\cos(\pi x_j)
\prod_{k\ne j}\sin(\pi x_k).
$$

The shared `ConductivityData` derives $f=-\gamma\Delta u-\sum_j\partial_j u$.
The affine and quadratic formulas above supply separate polynomial checks
of that realization.

## Levels, norms, and acceptance

| Field degree | Grid points per axis | Expected L2 / H1 orders | Accepted adjacent orders |
| --- | --- | --- | --- |
| P1 (`H1<1>`) | `5→9→17` | $2/1$ | $1.6<r_0<2.4$, $0.75<r_1<1.4$ |
| P2 (`H1<2>`) | `3→5→9` | $3/2$ | $2.6<r_0<3.4$, $1.75<r_1<2.4$ |

Here $h=1/(n-1)$. Every adjacent interval checks finite positive errors,
strict reduction, and the stated rate bounds. The standard elliptic orders
assume the smooth manufactured field, conforming spaces, regular meshes,
consistent assembly, and sufficient dual regularity. Two observed intervals
provide case-specific numerical evidence, not an arbitrary-mesh theorem.

Forms and independent error norms use quadrature order 12. CG uses relative
tolerance $10^{-13}$, absolute tolerance $10^{-14}$, divergence threshold
$10^5$, and at most 20,000 iterations. Each solve checks a positive PETSc
convergence reason and a finite reported residual below $10^{-9}$. This is
the solver's residual norm, which can depend on its preconditioner; it is not
reported as an independently recomputed unpreconditioned residual.

L2 and H1-seminorm errors are evaluated through `ErrorNorm`, separately from
the assembled residual. For MPI, only owned cells contribute and squared
norms are globally summed before taking square roots. A constant unit error
on the unit box checks the exact integral $E_0^2=|\Omega|=1$, detecting
duplicate halo integration. A sensitivity test on `n=5`, P2 raises quadrature
to 14 and tightens relative solver tolerance to $10^{-14}$; both errors must
change by less than $10^{-6}$ relative to the baseline.

The negative control replaces $\gamma$ in the stiffness by 1 while keeping
the correct affine forcing and trace. Its L2 error must exceed $10^{-3}$ and
its H1-seminorm error must exceed $10^{-2}$, demonstrating that the patch
oracle detects an omitted variable coefficient.

## Geometry and backend scope

All cases cover segment, triangle, quadrilateral, tetrahedron, pyramid,
hexahedron, and wedge. Distributed executions use 1, 2, 3, and 4 ranks.
The coarsest P2 segment mesh has two cells, so the four-rank case also
exercises ranks without owned cells. Point entities have no PDE h-rate here.
The suite covers real scalar conductivity, affine geometry, and h refinement;
complex/vector fields, curved geometry and p/hp paths are not claimed by
this Dirichlet target. The additional boundary target is specified below.

`DistributedUniformGrid` in `tests/convergence/MPIConvergence.h` constructs
the root mesh, computes the explicit incidence requirements, partitions it,
and distributes it through the existing Sharder. Poisson reuses this factory
and the same global norm evaluator. No library assembly or ownership code is
changed by this suite.

The `PETScConvergence` CI job builds this target and selects its CTest entries
along with the existing PETSc Poisson suites. Configure with
`RODIN_BUILD_CONVERGENCE_TESTS=ON`, `RODIN_USE_PETSC=ON`, and, for distributed
cases, `RODIN_USE_MPI=ON`. Run
`ctest --test-dir build/tests -R 'RodinConvergenceHPETSc(MPI)?Conductivity' --output-on-failure`.

## Natural and Robin boundary target

`RodinConvergenceHPETScConductivityBoundary` instantiates the shared scalar
diffusion boundary fixture with the variable coefficient above. Its
manufactured field is $u=\exp(\sum_jx_j)$, so

$$
\nabla u=u\boldsymbol{1},\qquad
f=-d(\gamma+1)u,\qquad
g_N=\gamma\nabla u\cdot n.
$$

The essential partition is the face $x_0=0$; the remaining faces carry
either $g_N$ or the Robin datum $g_R=g_N+2u$. Boundary attributes are
assigned before partitioning. The mixed weak form, for test functions
vanishing on the essential partition, is

$$
\int_\Omega\gamma\nabla u\cdot\nabla v\,\mathrm{d}x
+\alpha\int_{\Gamma_N}uv\,\mathrm{d}s
=\int_\Omega fv\,\mathrm{d}x
+\int_{\Gamma_N}(g_N+\alpha u)v\,\mathrm{d}s,
\qquad \alpha\in\{0,2\}.
$$

With MUMPS available, a pure-Neumann variant uses all boundary faces and
$u=\exp(\sum_jx_j)-(e-1)^d$, enforcing the zero integral by a constant
Lagrange multiplier. The gradient and forcing are unchanged. Independently
integrated mean and compatibility multiplier must satisfy absolute budgets
$10^{-10}$; omitted flux must produce a compatibility defect above $10^{-2}$.

Each boundary condition also checks the affine P1 and quadratic P2 fields
specified above on `n=3`, requiring both independent field errors below
$10^{-9}$. For pure Neumann data their unit-box means are respectively
$1+d/2$ and $1+d/3$; subtracting these constants retains the exact gradient,
forcing and flux while enforcing the same zero-mean gauge. Thus polynomial
reproduction, smooth-field rates and the compatibility constraint have
distinct oracles.

Each condition uses P1 on `5→9→17` and P2/P3 on `3→5→9` grid points per
axis. Every adjacent interval must reduce both errors and exceed the
L2/H1-seminorm floors $1.65/0.75$, $2.45/1.55$, and $3.45/2.35$,
respectively. These are case-specific acceptance floors, not a claim of
universal mixed-boundary dual regularity. Assembly order 16 and independent
norm order 18 are checked separately against orders 18 and 20. Tightening
the solver tolerance from $10^{-13}$ to $10^{-14}$ must also change each
error by less than $10^{-6}$ relatively. The independently recomputed
coefficient residual satisfies
$\|A u_h-b\|_2/\max(1,\|b\|_2)<10^{-11}$.

Removing the natural flux while retaining the forcing and essential/Robin
data must increase each field error by a factor greater than 5. This control
distinguishes a genuinely exercised natural boundary from an irrelevant
zero-flux example. All seven positive-dimensional geometries are registered
for local meshes and MPI ranks 1–4, in both sequential and OpenMP CI jobs.
Registration describes intended coverage; passing numerical evidence is
required before certification. The target has its own CI runtime budget.
