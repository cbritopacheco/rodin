# Complex PETSc Helmholtz combined refinement

The coercive wave-number regime, analytic complex plane wave, shared
workload, residual budgets and wrong-mass control are specified in the
[degree suite](../../p/PETScHelmholtz/README.md).
The three-level path matches the native hp study:

$$
(n,p)=(2,1)\to(3,2)\to(5,3),\qquad h=1/(n-1).
$$

For both independently integrated $L^2$ and $H^1$-seminorm errors,
every adjacent interval requires finite positive errors, strict reduction
and the effective slopes

$$
r_i^{(s)}=\frac{\log(E_{i-1}^{(s)}/E_i^{(s)})}{\log(h_{i-1}/h_i)},
\qquad r_i^{(0)}>1.9,\quad r_i^{(1)}>0.9.
$$

These are combined-path policies, not fixed-degree h orders or a
universal hp theorem. Degree three on $n=3$ supplies the highest-degree
quadratic patch, resolved affine omitted-mass control and separate
assembly, norm and solver sensitivity checks, with the p suite's budgets.
The degree-one/two reproduction threshold is independently retained.
All seven positive-dimensional geometries have local and MPI rank-one-
through-four cases in sequential and OpenMP native-complex PETSc builds.
