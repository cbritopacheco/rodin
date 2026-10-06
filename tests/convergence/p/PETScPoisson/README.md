# PETSc Poisson degree refinement

On $\Omega=(0,1)^d$, the real field satisfies

$$
-\Delta u=f,\qquad u|_{\partial\Omega}=u_*,\qquad
u_*(x)=e^{s(x)},\quad s(x)=\sum_{j=0}^{d-1}x_j,\quad f=-d e^s.
$$

Conforming $H^1$ spaces of degrees $p=1\to2\to3\to4$ are solved
on the fixed grid with $n=2$ points per coordinate. These are the analytic
field and levels of the native degree suite. For the independently
integrated $L^2$ and $H^1$-seminorm errors $E_p$, every adjacent interval
requires finite positive errors, strict decrease, and
$\log(E_{p-1}/E_p)>0.25$. This is a finite-resolution decay policy under
analytic regularity, conformity and consistent integration, not a universal
exponential constant.

## Architecture and controls

`ConductivityData` derives physical fields, gradients and sources for the
selected coefficient. `PETScDiffusionProblem` creates a fresh fixed-layout
system per solve; `PETScDiffusionRefinement` owns the refinement policy and
uses `FieldConvergence` for every-interval acceptance. The
[conductivity counterpart](../PETScConductivity/README.md) and both hp
suites share these classes and one test entry point. A solve-scoped const
observer permits mapped-domain comparisons without new norm implementations.

The polynomial field $u_*=1+\sum_jx_j^2$ has $f=-2d$.
On $n=2$, degree one must have errors above $10^{-3}$ and $10^{-2}$,
while degree two reproduces the field and gradient within $10^{-9}$.
Degree four independently reproduces this patch on $n=3$. A resolved
wrong-operator control doubles the stiffness, retaining the correct source
and trace: both errors must then exceed the degree-one separation bounds,
after the correct operator has passed its exact-patch budget.

On $n=3$ at degree four, assembly order $16\to18$, norm order $18\to20$,
and CG relative tolerance $10^{-13}\to10^{-14}$ vary separately.
Each norm must change by less than $10^{-6}$ relatively. CG must return a
positive convergence reason, a finite reported residual below $10^{-8}$,
and an independent residual
$\Vert Ax-b\Vert _2/\max(1,\Vert b\Vert _2)<10^{-11}$.

All seven positive-dimensional cell families run with local meshes and
MPI ranks one through four, under sequential and OpenMP assembly. MPI
norms count owned cells once before globally reducing squared errors.
A point has no nonconstant spatial degree rate. Complex PETSc is not
registered for this real-field suite. Shared slow labels, 30-minute safety
timeouts and pyramid resource locks affect scheduling, not error budgets.
