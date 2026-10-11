# PETSc-backed Poisson h-convergence

This suite runs the same smooth sine-product Poisson problem as the local
Eigen-backed suite, but assembles the weak form into PETSc matrices and vectors
and solves with PETSc CG. In dimension $d\in\lbrace 1,2,3\rbrace$, the exact field is
$u_\ast(x)=\prod_{j=1}^d\sin(\pi x_j)$, the load is
$f=d\pi^2u_\ast$, and the boundary trace is prescribed from $u_\ast$. It checks
independently integrated $L^2$ and $H^1$-seminorm errors on all seven
`UniformGrid` cell geometries.

P1 uses $n=5,9,17$ points per axis and verifies respective orders two and
one; P2 uses $n=3,5,9$ and verifies orders three and two. Both sequences
give two adjacent rate checks. Assembly and error norms use quadrature order
12. The PETSc CG relative and absolute tolerances are $10^{-12}$ and
$10^{-14}$, with at most 20,000 iterations. This is local-context backend
evidence; distributed PETSc convergence remains a separate global-norm suite.

As a negative control, multiplying the manufactured load by $1.25$ while
keeping the exact field and boundary data fixed made the P2 rate checks fail
on all seven geometries. The correct load was restored before the passing
run.

## Natural-boundary counterparts

`RodinConvergenceHPETScPoissonBoundary` adds real-PETSc local and MPI
counterparts of the native boundary studies. On $\Omega=(0,1)^d$, let
$U(x)=\exp(\sum_jx_j)$ and $f=-dU$. For the mixed cases,
$\Gamma_D=\lbrace x_0=0\rbrace$ and $\Gamma_N=\partial\Omega\setminus\Gamma_D$.
The prescribed trace is $U$ and the outward flux is
$g=\nabla U\cdot n$. The Neumann weak form is

$$
\int_\Omega\nabla u_h\cdot\nabla v_h
=\int_\Omega fv_h+\int_{\Gamma_N}gv_h.
$$

The Robin case uses $\partial_nu+2u=r$, with $r=g+2U$, so that
$2\int_{\Gamma_N}u_hv_h$ is added on the left and
$\int_{\Gamma_N}rv_h$ replaces the flux term on the right.
Boundary attributes are assigned on the complete grid before partitioning;
the distributed assembly then selects physical boundary entities and owned
DOFs through the existing topology and index contracts. In one dimension,
the lower and upper endpoint normals are explicitly $-1$ and $+1$.

For pure Neumann data, the exact field is
$u_\ast=U-(e-1)^d$ and the entire boundary receives $g$. A global P0g
multiplier $\lambda_h$ removes the constant nullspace:

$$
\int_\Omega\nabla u_h\cdot\nabla v_h+\lambda_h\int_\Omega v_h
=\int_\Omega fv_h+\int_{\partial\Omega}gv_h,
\qquad \int_\Omega u_h=0.
$$

The exact compatibility identity is
$\int_\Omega f=-d(e-1)^d=-\int_{\partial\Omega}g$.
The computed mean and multiplier must each have absolute magnitude below
$10^{-10}$ on the unit-volume box. The multiplier is a discrete compatibility
diagnostic, not an independent field-error norm. Pure Neumann tests require
MUMPS and are explicitly excluded when PETSc lacks that solver; mixed cases
remain available.

All three boundary conditions use degree $p=1,2,3$, with three levels:
$n=5,9,17$ for P1 and $n=3,5,9$ for P2/P3. Every interval requires
finite positive errors, strict decrease, and respective L2/H1 floors
$(1.65,0.75)$, $(2.45,1.55)$ and $(3.45,2.35)$.
The usual target orders are $(p+1,p)$ under smooth-solution, mesh-regularity
and dual-regularity hypotheses. These finite-sequence floors provide
case-specific verification, not a theorem for arbitrary mixed boundaries.

Assembly order is 16; independent norm order is 18. Mixed cases use CG
with relative tolerance $10^{-13}$, absolute tolerance $10^{-14}$ and at
most 50,000 iterations. Pure Neumann uses PREONLY/LU/MUMPS on the augmented
system. Every solve requires a positive PETSc convergence reason and
$\lVert Ax-b\rVert_2/\max(1,\lVert b\rVert_2)<10^{-11}$.
At P2 and $n=3$, assembly order 18, norm order 20 and CG relative tolerance
$10^{-14}$ are varied separately; relative changes in both field errors
must remain below $10^{-6}$. Tightening CG does not change a direct
factorization in the pure-Neumann case; its two quadrature checks remain
independent of the assembled residual.

Removing the manufactured flux while retaining the source, trace and Robin
coefficient is the boundary negative control. Both field errors must exceed
five times their correct counterparts; for pure Neumann, the compatibility
multiplier must also exceed $10^{-2}$. All seven cell families are registered
for local execution and MPI ranks 1–4, separately in sequential/OpenMP builds.
These registrations are tagged slow and retain the shared 1,800-second
per-geometry timeout. Remote CI evidence remains separate from local runs.

The direct solver retains the Stokes policy: MUMPS `ICNTL(20)=0` selects a
centralized right-hand side, avoiding the distributed-RHS scatter failure
observed with the local installation, and `ICNTL(14)=100` reserves extra
factor workspace. Explicit PETSc options take precedence. Neither setting
modifies the assembled operator or right-hand side.
