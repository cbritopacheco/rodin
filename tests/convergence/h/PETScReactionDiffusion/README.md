# PETSc coupled reaction–diffusion h-convergence

On $\Omega=(0,1)^d$, $d\in\lbrace 1,2,3\rbrace $, two real fields satisfy

$$
-\kappa_i\Delta u_i+\sum_{j=1}^{2}R_{ij}u_j=f_i,
\qquad \kappa=(1,2),\qquad
R=\begin{pmatrix}1&0.2\\0.2&1\end{pmatrix},
$$

with full manufactured Dirichlet traces. The test space is
$H_0^1(\Omega)^2$ and the weak form is

$$
a(u,v)=\sum_{i=1}^{2}\int_\Omega\kappa_i\nabla u_i\cdot\nabla v_i
+\sum_{j=1}^{2}R_{ij}u_jv_i\,\mathrm{d}x
=\sum_{i=1}^{2}\int_\Omega f_iv_i\,\mathrm{d}x.
$$

Positive diffusion and the symmetric reaction matrix, with eigenvalues
$0.8$ and $1.2$, give a coercive symmetric form. Two scalar H1 trial/test
pairs are assembled into a coupled PETSc problem and solved by CG. This
exercises block offsets, both off-diagonal terms, and both Dirichlet traces;
it is not two independently solved scalar equations.

## Manufactured data and oracles

Let $s=\sum_jx_j$. The shared `ReactionDiffusionData` supplies the following
fields, with zero-based component index $i\in\lbrace 0,1\rbrace $ and $c_i=i+1$:

| Case | Fields | Gradient entries | Laplacian |
| --- | --- | --- | --- |
| Affine P1 patch | $u_i=c_i(1+s)$ | $\partial_j u_i=c_i$ | $0$ |
| Quadratic P2 patch | $u_i=c_i(1+s^2)$ | $\partial_j u_i=2c_i s$ | $2dc_i$ |
| Smooth rates | $u_0=e^s$, $u_1=2e^{-s}$ | $\partial_j u_0=u_0$, $\partial_j u_1=-u_1$ | $du_i$ |

Sources are $f_i=-\kappa_i\Delta u_i+u_i+0.2u_{1-i}$. Each field's
L2 and H1-seminorm errors are integrated separately:

$$
E_{0,i}^2=\int_\Omega|u_i-u_{h,i}|^2\,\mathrm{d}x,\qquad
E_{1,i}^2=\int_\Omega\Vert\nabla u_i-\nabla u_{h,i}\Vert_2^2\,\mathrm{d}x.
$$

Patches use `n=5` grid points per axis and require both errors below
$10^{-9}$ for both fields. The negative control sets the two off-diagonal
reaction coefficients to zero but retains the correct affine sources and
traces. Both fields must then have $E_{0,i}>10^{-3}$ and
$E_{1,i}>10^{-2}$, establishing that the oracle rejects missing coupling.

## Refinement and numerical budget

| Degree | Levels | Expected L2/H1 orders | Accepted adjacent orders |
| --- | --- | --- | --- |
| P1 | `n=5→9→17` | $2/1$ | $1.6<r_0<2.4$, $0.75<r_1<1.4$ |
| P2 | `n=3→5→9` | $3/2$ | $2.6<r_0<3.4$, $1.75<r_1<2.4$ |

With $h=1/(n-1)$, every interval checks finite positive errors, strict
reduction, and both rate bounds for each component. Expected rates assume
smooth solutions, conforming regular affine meshes, consistent integration,
and the required dual regularity. Observed rates are numerical evidence
under these hypotheses, not an unconditional theorem about every geometry.

All forms and error integrals use order 12. CG uses relative tolerance
$10^{-13}$, absolute tolerance $10^{-14}$, divergence threshold $10^5$,
and at most 50,000 iterations. Each solve requires a positive PETSc
convergence reason and a finite reported residual below $10^{-8}$; this
reported residual is distinct from the independently recomputed
$\Vert Ax-b\Vert_2/\max(1,\Vert b\Vert_2)<10^{-11}$ also required by the shared
`PETScReactionDiffusionProblem` workload. The workload constructs a fresh
fixed-layout two-field problem for each measurement, with independently
selectable norm quadrature and a solve-scoped const observer. It is shared
with the [p](../../p/PETScReactionDiffusion/README.md) and
[hp](../../hp/PETScReactionDiffusion/README.md) counterparts.
At `n=5`, P2 sensitivity raises quadrature to 14 and tightens relative
tolerance to $10^{-14}$; every component error must change by less than
$10^{-6}$ relative to its baseline.

## Natural-boundary extension

`RodinConvergenceHPETScReactionDiffusionBoundary` uses the same coupled
operator and manufactured fields. Let $\Gamma_D=\lbrace x_0=0\rbrace $ and
$\Gamma_N=\partial\Omega\setminus\Gamma_D$. Mixed Neumann and Robin
cases prescribe the exact trace on $\Gamma_D$ and

$$
\kappa_i\partial_nu_i+\beta u_i=g_i,
\qquad g_i=\kappa_i\nabla u_i\cdot n+\beta u_i,
\qquad \beta\in\lbrace 0,1\rbrace ,
$$

on $\Gamma_N$. Pure Neumann cases instead use $\Gamma_D=\varnothing$,
$\Gamma_N=\partial\Omega$, and $\beta=0$. The two-field test space is
$V=\lbrace v\in H^1(\Omega)^2:v|_{\Gamma_D}=0\rbrace $, and the weak form is

$$
a(u,v)+\beta\sum_{i=1}^{2}\int_{\Gamma_N}u_iv_i\,\mathrm{d}s
=\sum_{i=1}^{2}\left(\int_\Omega f_iv_i\,\mathrm{d}x
+\int_{\Gamma_N}g_iv_i\,\mathrm{d}s\right).
$$

In contrast to pure diffusion, constant fields do not form a nullspace:

$$
a(v,v)\ge\sum_{i=1}^{2}\Vert\nabla v_i\Vert_{L^2(\Omega)}^2
+0.8\sum_{i=1}^{2}\Vert v_i\Vert_{L^2(\Omega)}^2.
$$

Therefore pure Neumann cases require neither a compatibility projection nor
a mean multiplier. The flux of the second field includes its diffusion
coefficient $\kappa_2=2$. Segment endpoints use the signed outward normals
$-1$ and $+1$; higher-dimensional faces use `BoundaryNormal`.

For each boundary variant, P1 uses `n=5→9→17`; P2/P3 use `n=3→5→9`.
Both fields and both adjacent intervals require positive finite, strictly
decreasing errors. The L2/H1-seminorm rate floors are respectively
$1.65/0.75$, $2.45/1.55$, and $3.45/2.55$. These are acceptance policies
under the regularity hypotheses above, not substitutes for those hypotheses.
Affine P1 and quadratic P2 patches at `n=3` require every field error below
$10^{-9}$. Removing coupling while retaining affine manufactured loads
must give both errors above $10^{-3}$ for each field. Independently removing
normal flux, while retaining the Robin reaction term and its exact load,
must increase each smooth P2 error by more than a factor of five.

At `n=3`, independent sensitivity checks change assembly order $16\to18$,
norm order $18\to20$, or relative solver tolerance $10^{-13}\to10^{-14}$,
one at a time. Each field error must change by less than $10^{-6}$ relative
to its baseline. Residual budgets are unchanged. All seven geometries have
local and MPI rank 1–4 registrations: 24 cases per geometry/context, split
into 35 CTest entries with `slow` labels and 1800-second limits. Run
`ctest --test-dir build/tests -R '^RodinConvergenceHPETScReactionDiffusionBoundary_' --output-on-failure`.

Distributed construction, assembly, solving, residuals, and error norms
require the mesh communicator's participation. Manufactured data, normal
contractions, and local metadata queries do not introduce collectives.

## Geometry and backend scope

Local PETSc and MPI PETSc cover segment, triangle, quadrilateral,
tetrahedron, pyramid, hexahedron, and wedge, with MPI ranks 1–4.
The coarsest P2 segment mesh has two cells, so some ranks have no owned
cells. `DistributedUniformGrid` supplies partitioned meshes. Error norms
integrate owned cells only and globally reduce squared contributions before
taking square roots. For each discrete field, a zero function compared with
the unit constant checks $E_0=1$ and $E_1=0$ through both norm interfaces.

Configure with `RODIN_BUILD_CONVERGENCE_TESTS=ON`, `RODIN_USE_PETSC=ON`,
and `RODIN_USE_MPI=ON` for distributed cases. Run
`ctest --test-dir build/tests -R 'RodinConvergenceHPETSc(MPI)?ReactionDiffusion' --output-on-failure`.
Entries are split by geometry/rank count, labelled `slow`, and have a
600-second timeout; MPI entries declare their process count. The dedicated
PETSc CI job selects them explicitly by name.

The original target covers real two-field full-Dirichlet h studies on affine
meshes; the boundary target extends it as specified above. Complex fields
and nonsymmetric reaction are outside these suites. The p/hp and curved studies have separate
specifications. Point/0D has no PDE h-rate here.
