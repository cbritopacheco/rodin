# PETSc-backed Poisson h-convergence

This suite runs the same smooth sine-product Poisson problem as the local
Eigen-backed suite, but assembles the weak form into PETSc matrices and vectors
and solves with PETSc CG. In dimension $d\in\{1,2,3\}$, the exact field is
$u_*(x)=\prod_{j=1}^d\sin(\pi x_j)$, the load is
$f=d\pi^2u_*$, and the boundary trace is prescribed from $u_*$. It checks
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
