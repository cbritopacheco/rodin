# PETSc linear-elasticity h-convergence

On $\Omega=(0,1)^d$, $d\in\{1,2,3\}$, the real displacement
$u:\Omega\to\mathbb{R}^d$ satisfies

$$
-\nabla\cdot\sigma(u)=f,\qquad
\sigma(u)=\lambda(\nabla\cdot u)I+2\mu\varepsilon(u),\qquad
\varepsilon(u)=\tfrac12(Du+Du^T),
$$

with the exact displacement prescribed on the full boundary. Here
$\lambda=1.5$ and $\mu=0.5$. The weak form is

$$
a(u,v)=\int_\Omega\lambda(\nabla\cdot u)(\nabla\cdot v)
+2\mu\varepsilon(u):\varepsilon(v)\,\mathrm{d}x
=\int_\Omega f\cdot v\,\mathrm{d}x.
$$

The homogeneous test space is $H_0^1(\Omega;\mathbb{R}^d)$. Korn's
inequality and positive elastic coefficients give coercivity and remove the
rigid-motion nullspace under full Dirichlet conditions, permitting CG.
The segment case is the one-component reduction with stiffness
$\lambda+2\mu$.

## Manufactured fields

Let $s=\sum_j x_j$, $c_i=i+1$ for zero-based component indices, and
$C=\sum_i c_i=d(d+1)/2$. The shared `ManufacturedSolution` supplies:

| Case | Displacement | Jacobian entries | Body force |
| --- | --- | --- | --- |
| Affine P1 patch | $u_i=1+c_i s$ | $D_j u_i=c_i$ | $f_i=0$ |
| Quadratic P2 patch | $u_i=1+c_i s^2$ | $D_j u_i=2c_i s$ | $f_i=-2[\mu d c_i+(\lambda+\mu)C]$ |
| Smooth rate study | $u_i=c_i e^s$ | $D_j u_i=c_i e^s$ | $f_i=-e^s[\mu d c_i+(\lambda+\mu)C]$ |

These sources follow from
$f=-\mu\Delta u-(\lambda+\mu)\nabla(\nabla\cdot u)$.
Both diagonal and off-diagonal derivatives are nonzero. The quadratic and
smooth fields have a nonconstant divergence, so the volumetric coupling is
tested along with the shear term. The trace is nonzero in all cases.

Affine and quadratic patches use `n=5` grid points per axis. They require
L2 and H1-seminorm errors below $10^{-9}$, accounting for floating-point
assembly and iterative solution of exactly representable fields. This is a
numerical patch bound, not a tolerance for mesh or DOF indices.

## Refinement and acceptance

| Space | Levels | Expected L2 / H1 orders | Accepted adjacent orders |
| --- | --- | --- | --- |
| Vector `H1<1>` | `n=5→9→17` | $2/1$ | $1.65<r_0<2.35$, $0.75<r_1<1.25$ |
| Vector `H1<2>` | `n=3→5→9` | $3/2$ | $2.45<r_0<3.55$, $1.55<r_1<2.45$ |

With $h=1/(n-1)$, every adjacent interval checks finite positive errors,
strict reduction, and both rate bounds. The expected estimates assume a
smooth solution, conforming spaces, regular affine meshes, consistent
integration, and the required dual regularity. The suite uses the same
Dirichlet rate bounds as the existing Eigen elasticity h studies.

Independent quadrature measures

$$
E_0^2=\int_\Omega\|u-u_h\|_2^2\,\mathrm{d}x,\qquad
E_1^2=\int_\Omega\|Du-Du_h\|_F^2\,\mathrm{d}x.
$$

For MPI, `ErrorNorm::computeVector` integrates owned cells only and sums
squared contributions globally before taking square roots. A zero discrete
field compared with $(1,\ldots,1)$ on the unit box checks $E_0=\sqrt d$
and $E_1=0$ through both the vector and value-only norm interfaces.

All forms and norms use quadrature order 12. CG uses relative tolerance
$10^{-13}$, absolute tolerance $10^{-14}$, divergence threshold $10^5$,
and at most 50,000 iterations. Each solve requires a positive PETSc
convergence reason and a finite reported residual below $10^{-8}$. That
reported residual can depend on preconditioning; it is not presented as an
independently recomputed unpreconditioned residual.

At `n=5`, the P2 sensitivity check raises quadrature to 14 and tightens the
relative solve tolerance to $10^{-14}$. Each measured error must change by
less than $10^{-6}$ relative to the baseline. The negative control omits the
volumetric bilinear term while retaining the correct quadratic forcing and
trace. Its L2 and H1-seminorm errors must exceed $10^{-3}$ and $10^{-2}$,
respectively, demonstrating that the patch oracle rejects this operator error.

## Geometry and backend scope

Local PETSc and distributed PETSc cases cover segment, triangle,
quadrilateral, tetrahedron, pyramid, hexahedron, and wedge. MPI executions
use 1, 2, 3, and 4 ranks; the coarsest P2 segment case has two cells and thus
also exercises ranks without owned cells. Mesh construction uses the shared
`DistributedUniformGrid` factory. Local and MPI fixtures use the same data,
formulation, degree families, and assertions.

This suite covers real vector displacement, full Dirichlet conditions,
moderate compressibility, affine geometry, and h refinement. PETSc mixed
traction, near-incompressibility, p/hp, and curved geometry remain separate
coverage tasks. Point entities have no PDE h-rate here.

The executable is built by `PETScConvergence` CI. Configure with
`RODIN_BUILD_CONVERGENCE_TESTS=ON`, `RODIN_USE_PETSC=ON`, and, for distributed
cases, `RODIN_USE_MPI=ON`. Run
`ctest --test-dir build/tests -R 'RodinConvergenceHPETSc(MPI)?LinearElasticity' --output-on-failure`.
CTest entries are split by geometry and, for MPI, rank count. Each carries
the `slow` label and a 600-second timeout; MPI entries also declare their
process count for CTest scheduling. The dedicated
PETSc job selects these entries by name, so that label does not remove their
CI execution.
