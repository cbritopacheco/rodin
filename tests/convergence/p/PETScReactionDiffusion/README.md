# PETSc coupled reaction–diffusion degree refinement

On $\Omega=(0,1)^d$, the two real fields satisfy

$$
-\kappa_i\Delta u_i+u_i+\alpha u_{1-i}=f_i,\qquad
\kappa=(1,2),\quad\alpha=0.2,\quad i\in\lbrace 0,1\rbrace .
$$

Both fields have their manufactured Dirichlet trace. Positive diffusion and
reaction eigenvalues $0.8,1.2$ give a symmetric coercive weak problem.
With $s=\sum_jx_j$, the analytic fields are $u_0=e^s$, $u_1=2e^{-s}$,
and $f_i=(1-d\kappa_i)u_i+\alpha u_{1-i}$. This is the coefficient/data
variant shared with the h and hp suites; the native p suite also retains
its original equal-diffusion variant.

The grid is fixed at $n=2$ points per coordinate; degrees are
$p=1\to2\to3\to4$. For each field, both independently integrated
$L^2$ and $H^1$-seminorm errors must be finite, positive and strictly
decreasing on every adjacent interval. Each error $E_{i,p}$ must obey

$$
\log(E_{i,p-1}/E_{i,p})>0.25.
$$

This is a finite-resolution decay policy for analytic fields, not a
universal exponential constant. A point has no nonconstant spatial rate.

## Shared workload and independent controls

`PETScReactionDiffusionProblem` owns the manufactured data and constructs a
fresh fixed-layout two-field system per measurement. CG must terminate
with a positive reason, a finite reported residual below $10^{-8}$, and
an independently recomputed residual
$\Vert Ax-b\Vert _2/\max(1,\Vert b\Vert _2)<10^{-11}$. `FieldConvergence` checks
each field and every interval without averaging errors across fields.
The p/hp entry point and geometry/rank registration are shared.

On $n=3$, degree four reproduces the quadratic patch
$u_i=(i+1)(1+s^2)$ within $10^{-9}$ in both norms. A separate affine
patch $u_i=(i+1)(1+s)$ first satisfies the same absolute budget; then
omitting both coupling terms, while retaining the original source and
trace, must give each field $E_{L^2}>10^{-3}$ and $E_{H^1}>10^{-2}$.
This distinguishes incorrect coupling from unresolved approximation.

At degree four on $n=3$, assembly quadrature $16\to18$, norm quadrature
$18\to20$, and CG relative tolerance $10^{-13}\to10^{-14}$ vary
separately. Both norms of both fields must change by less than $10^{-6}$
relatively. These controls isolate integration and algebraic contamination.

All seven positive-dimensional geometry families run with local meshes
and MPI ranks one through four, under sequential and OpenMP assembly.
MPI norms integrate only owned cells and globally reduce squared errors.
Complex PETSc is not registered for this real-field suite. Tests are
slow-labelled, have 30-minute safety timeouts, and lock pyramid workloads.
