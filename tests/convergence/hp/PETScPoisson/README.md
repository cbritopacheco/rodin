# PETSc Poisson combined refinement

The analytic equation, shared workload, degree-threshold patch, residual
budgets and wrong-stiffness control are specified by the
[degree suite](../../p/PETScPoisson/README.md).
The three-level mesh/degree path is

$$
(n,p)=(2,1)\to(3,2)\to(5,3),\qquad h=1/(n-1).
$$

For each independently integrated error $E_i^{(0)}$ in $L^2$ and
$E_i^{(1)}$ in the $H^1$ seminorm, every adjacent interval requires
finite positive errors, strict reduction and

$$
r_i^{(s)}=\frac{\log(E_{i-1}^{(s)}/E_i^{(s)})}{\log(h_{i-1}/h_i)},
\qquad r_i^{(0)}>1.9,\quad r_i^{(1)}>0.9.
$$

These are effective slopes along the stated combined path, not
fixed-degree h orders or a universal hp theorem. Degree three on $n=3$
supplies the highest-degree patch, wrong-operator control and three
separate quadrature/solver checks; their budgets are those of the p suite.
All seven positive-dimensional geometries run locally and with MPI
ranks one through four, in sequential and OpenMP assembly configurations.
