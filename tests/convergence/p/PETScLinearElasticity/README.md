# PETSc vector linear-elasticity degree refinement

On $\Omega=(0,1)^d$, the real displacement $u:\Omega\to\mathbb R^d$
satisfies the isotropic small-strain equation

$$
-\nabla\cdot\sigma(u)=f,\qquad
\sigma(u)=\lambda(\nabla\cdot u)I+2\mu\varepsilon(u),\qquad
\varepsilon(u)=\tfrac12(\nabla u+\nabla u^T),\quad
\lambda=1.5,\quad\mu=0.5.
$$

Full manufactured Dirichlet traces eliminate rigid motions; positivity
and Korn's inequality give a coercive symmetric form. With
$s=\sum_jx_j$ and zero-based component index $i$, the analytic field is
$u_i=(i+1)e^s$, and

$$
f_i=-e^s\left[\mu d(i+1)+(\lambda+\mu)\sum_{j=0}^{d-1}(j+1)\right].
$$

The fixed grid has $n=2$ points per coordinate; degrees are
$p=1\to2\to3\to4$. These fields, Lamé parameters and levels match the
native degree study. Independently integrated vector $L^2$ and
$H^1$-seminorm errors include every displacement component and every
physical derivative. For either error $E_p$, every adjacent interval
requires finite positive errors, strict reduction and
$\log(E_{p-1}/E_p)>0.25$. This is a finite-path analytic decay policy,
not a universal exponential constant or an h order.

## Architecture and controls

`LinearElasticity::ManufacturedSolution` supplies physical displacement,
Jacobian, forcing, strain and stress. `PETScLinearElasticityProblem` creates
a fresh fixed-layout vector space/system per solve and exposes only
solve-scoped const views. `PETScLinearElasticityRefinement` reuses
`FieldConvergence` for every-interval acceptance. The p/hp suites share
their entry point and geometry/rank registration, and all norms use the
common owned-cell integration and global MPI reduction.

On $n=2$, the quadratic field $u_i=1+(i+1)s^2$ has nonzero degree-one
errors above $10^{-3}$ and $10^{-2}$, while degree two reproduces it
within $10^{-9}$. Degree four also reproduces it on $n=3$. A resolved
control then removes $\lambda(\nabla\cdot u)(\nabla\cdot v)$ from the
operator, retaining the correct forcing and trace; both errors must exceed
the same separation floors, after the correct operator passes its patch.

A separate degree-four affine patch uses $u=\mathbf1+Ax$, with
$A_{ij}=(i+1)(j+1)+\delta_{i0}\delta_{j,d-1}$. In dimensions two and
three its Jacobian is not symmetric. Displacement and its full Jacobian,
plus independently integrated Frobenius strain and stress errors, must
be below $10^{-9}$. This checks component ordering and constitutive
postprocessing beyond symmetric-gradient data alone.

On $n=3$ at degree four, assembly order $16\to18$, norm order $18\to20$,
and CG relative tolerance $10^{-13}\to10^{-14}$ vary separately.
Each displacement error changes by less than $10^{-6}$ relatively. CG
requires a positive reason, a finite reported residual below $10^{-8}$,
and an independently recomputed residual
$\Vert Ax-b\Vert _2/\max(1,\Vert b\Vert _2)<10^{-11}$.

All seven positive-dimensional cell families run locally and on MPI ranks
one through four, with sequential and OpenMP assembly. Complex PETSc is
not registered for this real-displacement suite. Point/0D has no spatial
elasticity degree rate. Tests are slow-labelled, use 30-minute safety
timeouts, and serialize pyramid workloads through the shared resource lock.
