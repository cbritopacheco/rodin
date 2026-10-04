# PETSc nonlinear Poisson degree refinement

This suite verifies real-PETSc/SNES solutions of

$$
-\Delta u+u+u^3=f\quad\text{in }\Omega=(0,1)^d,\qquad
u|_{\partial\Omega}=0,\qquad
u(x)=A\prod_{j=0}^{d-1}\sin(\pi x_j),\qquad
f=(d\pi^2+1)u+u^3.
$$

All seven positive-dimensional cell families are tested with sequential and
OpenMP assembly, local meshes and MPI ranks one through four. A point has
no nonconstant spatial approximation rate and is not included here.

The fixed grid has three points per axis; field degrees are
$p=1\to2\to3\to4$. Both independently integrated L2 and H1-seminorm
errors must be finite, positive and decrease on every adjacent interval.
For either error $E_p$, the finite-resolution policy is

$$
\log(E_{p-1}/E_p)>0.1.
$$

Analytic data support degree-refinement decay under the approximation and
solver hypotheses. The test does not identify a universal exponential
constant or equate a degree increment with an h-refinement order.

The shared `PETScNonlinearRefinement` class creates fresh fixed-layout
`PETScNonlinearPoissonProblem` instances for each measurement. Residual,
Jacobian, state synchronization and SNES solve logic remain in the common
problem class. MPI norms count owned cells once before global reduction.
The p and hp registrations compile one common test entry point, rather
than duplicating assertions or PDE assembly.

For every degree, the centered residual-difference tangent defect must be
below $10^{-6}$. Replacing $3u^2$ by $u^2$ in the highest-degree tangent
must produce a defect above $10^{-3}$. An amplitude-$4$ control retains
the original source but omits the cubic reaction from the solved operator:
the correct L2/H1 errors must be below $0.05/0.2$, while the incorrect
errors must exceed these same bounds. This separates unresolved fields
from incorrect physics.

Assembly quadrature order $16\to18$, norm order $18\to20$, and SNES
tolerance $10^{-11}\to10^{-12}$ are varied separately at degree four.
Each norm must change by less than $10^{-6}$ relatively. Residual and
convergence-reason checks remain enforced inside the problem class.
Tests are slow-labelled with 30-minute safety timeouts; pyramid
registrations share a resource lock.
