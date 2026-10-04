# Combined hp Poisson refinement

The exact solution is the analytic sine product
$u=\prod_i\sin(\pi x_i)$, with $-\Delta u=d\pi^2u$ and homogeneous
Dirichlet data. Three runs change the grid and representing-space degree
together: $(h,p)=(1,1),(1/2,2),(1/4,3)$. They cover every seven
UniformGrid cell geometries.

Analyticity gives rapid polynomial approximation and mesh refinement reduces
the length scale simultaneously. Since both $h$ and $p$ change, this suite
does not call the observed log-ratio a fixed-order h rate. Instead it checks
strict reduction in independently integrated L2 and H1-seminorm errors and
requires the reduction to beat the minimum P1 h-only powers: two in L2 and
one in H1. This is a conservative combined-refinement certification, not a
claim of a universal hp exponent.
