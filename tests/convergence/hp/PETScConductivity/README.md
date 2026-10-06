# PETSc variable-conductivity combined refinement

The equation, analytic data and wrong-coefficient control are those of
the [degree suite](../../p/PETScConductivity/README.md), with
$\gamma=1+\sum_jx_j$ and $u_\ast=e^{\sum_jx_j}$.
The three-level path is

$$
(n,p)=(2,1)\to(3,2)\to(5,3),\qquad h=1/(n-1).
$$

Every adjacent interval checks finite positive and strictly decreasing
$L^2$ and $H^1$-seminorm errors. The effective slopes are

$$
r_i^{(s)}=\frac{\log(E_{i-1}^{(s)}/E_i^{(s)})}{\log(h_{i-1}/h_i)},
\qquad r_i^{(0)}>1.9,\quad r_i^{(1)}>0.9.
$$

This is a combined-path policy, not a fixed-degree h order. It uses the
same levels and analytic problem as the native hp suite. The shared
`PETScDiffusionRefinement` checks every interval and reuses the degree-two
reproduction threshold. Degree three on $n=3$ supplies the quadratic
patch, resolved affine wrong-coefficient control and independent assembly,
norm-quadrature and solver-tolerance checks, with the p suite's budgets.
All seven positive-dimensional geometries have local and MPI rank-one-
through-four registrations under sequential and OpenMP assembly.
