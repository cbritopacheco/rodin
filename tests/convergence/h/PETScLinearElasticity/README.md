# PETSc linear-elasticity h-convergence

On $\Omega=(0,1)^d$, $d\in\lbrace 1,2,3\rbrace$, the real displacement
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
+2\mu\varepsilon(u):\varepsilon(v)\thinspace \mathrm{d}x
=\int_\Omega f\cdot v\thinspace \mathrm{d}x.
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
E_0^2=\int_\Omega\Vert u-u_h\Vert_2^2\thinspace \mathrm{d}x,\qquad
E_1^2=\int_\Omega\Vert Du-Du_h\Vert_F^2\thinspace \mathrm{d}x.
$$

For MPI, `ErrorNorm::computeVector` integrates owned cells only and sums
squared contributions globally before taking square roots. A zero discrete
field compared with $(1,\ldots,1)$ on the unit box checks $E_0=\sqrt d$
and $E_1=0$ through both the vector and value-only norm interfaces.

All forms and norms use quadrature order 12. CG uses relative tolerance
$10^{-13}$, absolute tolerance $10^{-14}$, divergence threshold $10^5$,
and at most 50,000 iterations. Each solve requires a positive PETSc
convergence reason and a finite reported residual below $10^{-8}$. That
reported residual can depend on preconditioning. Separately, the shared
workload recomputes the unpreconditioned residual and requires
$\lVert Ax-b\rVert_2/\max(1,\lVert b\rVert_2)<10^{-11}$.

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

The original executable covers real vector displacement, full Dirichlet
conditions, moderate compressibility, affine geometry, and h refinement.
The boundary executable below supplies mixed traction and nearly
incompressible counterparts. PETSc [p](../../p/PETScLinearElasticity/README.md),
[hp](../../hp/PETScLinearElasticity/README.md) and
[curved](../../isoparametric/LinearElasticity/README.md) studies have separate
specifications. Point entities have no PDE h-rate here.

The executable is built by `PETScConvergence` CI. Configure with
`RODIN_BUILD_CONVERGENCE_TESTS=ON`, `RODIN_USE_PETSC=ON`, and, for distributed
cases, `RODIN_USE_MPI=ON`. Run
`ctest --test-dir build/tests -R 'RodinConvergenceHPETSc(MPI)?LinearElasticity' --output-on-failure`.
CTest entries are split by geometry and, for MPI, rank count. Each carries
the `slow` label and a 600-second timeout; MPI entries also declare their
process count for CTest scheduling. The dedicated
PETSc job selects these entries by name, so that label does not remove their
CI execution.

## Mixed traction and nearly incompressible counterparts

`RodinConvergenceHPETScLinearElasticityBoundary` reuses the common PETSc
workload, manufactured fields and vector norms. For mixed conditions,
$\Gamma_D=\lbrace x_0=0\rbrace$ and $\Gamma_T=\partial\Omega\setminus\Gamma_D$.
The natural datum is the analytic traction $t=\sigma(u_\ast)n$, so that

$$
a(u_h,v_h)=\int_\Omega f\cdot v_h\thinspace \mathrm{d}x
+\int_{\Gamma_T}t\cdot v_h\thinspace \mathrm{d}s,
\qquad v_h\rvert_{\Gamma_D}=0.
$$

The nonempty displacement partition removes rigid motions; no penalty
replaces the essential condition. Boundary attributes are prescribed on the
complete mesh before MPI partitioning. Physical normals use `BoundaryNormal`
in dimensions two and three and the explicit endpoint signs in one dimension.
The traction uses the closed-form manufactured stress, independently of the
discrete Jacobian used to measure the solution error.

Let $V_D=\lbrace v\in H^1(\Omega;\mathbb R^d):v\rvert_{\Gamma_D}=0\rbrace$ and
$e_h=u_\ast-u_h$. For conforming approximation of degree $p$ to a smooth solution,
with compatible approximation of the essential data, the energy estimate is
$\lVert e_h\rVert_{H^1}\le C h^p\lVert u_\ast\rVert_{H^{p+1}}$.
The L2 estimate additionally depends on the adjoint problem: for
$g\in L^2(\Omega;\mathbb R^d)$, let $z_g\in V_D$ satisfy
$a(v,z_g)=\langle g,v\rangle_{L^2}$ for all $v\in V_D$. If
$\lVert z_g\rVert_{H^{1+s}}\le C_{\mathrm{reg}}\lVert g\rVert_{L^2}$,
$0<s\le1$, and the boundary approximation is dual-consistent, the expected
L2 order is $p+s$. The familiar order $p+1$ corresponds to $s=1$; that
regularity is not established here for the mixed boundary partition.
The constants may depend on the Lamé ratio and boundary partition.

With $\lambda=1.5$ and $\mu=0.5$, mixed P1/P2 studies retain the native
exponential field and three-level sequences. P1 uses $n=5,9,17$, except
tetrahedra, which use $n=9,17,33$: the coarsest interval in the native study
did not meet the stated L2 floor. P2 uses $n=3,5,9$. Every interval requires positive
finite errors, strict decrease and L2/H1 rate floors $(1.65,0.75)$ for P1
and $(2.45,1.55)$ for P2. Expected orders $(p+1,p)$ require the stated
regularity and stable conforming approximation; the measured floors are
finite-workload observations, not unconditional mixed-boundary estimates.

At $n=3$, an asymmetric affine vector field verifies P1 traction
reproduction, including nonsymmetric Jacobian entries; the quadratic
manufactured field verifies P2 reproduction. Both field errors must remain
below $10^{-9}$. Removing the traction from the P2 patch while retaining
the volume load and displacement trace must give L2 error above $10^{-3}$
and H1-seminorm error above $10^{-2}$. This control isolates the natural
boundary contribution rather than merely changing the constitutive law.

The nearly incompressible case retains full Dirichlet data and uses
$\lambda=10^4$, $\mu=1$, with the exact divergence-free shear field

$$
u_\ast(x)=(\sin(\pi x_1),0,\ldots,0),\qquad
\nabla\cdot u_\ast=0,\qquad f=\mu\pi^2u_\ast.
$$

It is defined in dimensions two and three. Its exact field, source and
stress do not depend on $\lambda$; a shared-data regression checks that
identity and the known component derivatives on all six applicable geometries.
P2 uses $n=3,5,9$, with L2/H1 floors $(2.25,1.35)$, matching the native
case. Segment rate/sensitivity cases are explicitly skipped because a
nonconstant divergence-free displacement does not exist in one dimension.
This is verification of one resolved divergence-free workload at a fixed
large Lamé ratio. It does not prove uniform locking-free behavior as
$\lambda/\mu\to\infty$, or stability for arbitrary loads and boundary data.

Both variants use assembly order 16, independent norm order 18 and CG/Jacobi
relative/absolute tolerances $10^{-13}/10^{-14}$, with at most 50,000
iterations. Solver status is checked separately from
$\lVert Ax-b\rVert_2/\max(1,\lVert b\rVert_2)<10^{-11}$.
At P2 and $n=3$, assembly order 18, norm order 20 and relative solve
tolerance $10^{-14}$ are varied independently; relative changes in each
field error must remain below $10^{-6}$.

Jacobi is explicitly selected in this executable because its positive
diagonal scaling respects CG's definiteness requirement. The default local
incomplete-factorization policy failed the nearly incompressible
hexahedral case with `KSP_DIVERGED_INDEFINITE_PC`; selecting Jacobi resolves
that failure without modifying the operator, loads or error bounds.
This is a verification-workload solver choice, not a new library solver
default or a change to the existing p/hp studies.

All seven geometry families have mixed-condition registrations for local
and MPI ranks 1–4; the nearly incompressible cases have the dimension
restriction above. Sequential and OpenMP executions are separate evidence.
CTest uses the shared slow label, 1,800-second per-geometry timeout and
pyramid resource lock. A separate boundary CI runtime partition builds and
executes this target without removing refinement levels.
