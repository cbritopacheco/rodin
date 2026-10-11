# PETSc vector linear-elasticity combined refinement

The isotropic equation, vector manufactured field, Lamé parameters, shared
workload and independent numerical controls are specified by the
[degree suite](../../p/PETScLinearElasticity/README.md).
The three-level path matches the native combined study:

$$
(n,p)=(2,1)\to(3,2)\to(5,3),\qquad h=1/(n-1).
$$

Each interval checks finite positive and strictly decreasing vector
$L^2$ and $H^1$-seminorm errors. Their effective slopes obey

$$
r_i^{(s)}=\frac{\log(E_{i-1}^{(s)}/E_i^{(s)})}{\log(h_{i-1}/h_i)},
\qquad r_i^{(0)}>1.9,\quad r_i^{(1)}>0.9.
$$

These are combined-path policies, not fixed-degree h orders or a
universal hp theorem. Degree three on $n=3$ supplies the quadratic patch,
resolved omitted-volumetric-term control, asymmetric affine displacement/
strain/stress patch and three separate quadrature/solver checks, with the
p suite's budgets. The degree-one/two reproduction threshold is also
retained. All seven positive-dimensional geometries have local and MPI
rank-one-through-four cases under sequential and OpenMP real-PETSc assembly.
