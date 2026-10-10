# PETSc nonlinear Poisson combined refinement

This suite uses the same real-PETSc/SNES manufactured equation, shared
problem class, tangent checks and incorrect-reaction oracle as the
[degree-refinement suite](../../p/PETScNonlinearPoisson/README.md).
The shared study and test entry point keep those contracts identical
across refinement axes and local/MPI assembly.

The three-level path is

$$
(n,p)=(2,1)\to(3,2)\to(5,3),\qquad h=1/(n-1).
$$

For L2 error $E_i^{(0)}$ and H1-seminorm error $E_i^{(1)}$, every adjacent
interval has a finite, positive, strictly decreasing error and obeys

$$
r_i^{(s)}=\frac{\log(E_{i-1}^{(s)}/E_i^{(s)})}
 {\log(h_{i-1}/h_i)},\qquad
r_i^{(0)}>1.9,\qquad r_i^{(1)}>0.9.
$$

These are effective slopes along the stated combined path, not
fixed-degree h-orders or a universal exponential hp theorem.
On the coarsest P1 mesh, six families have no unconstrained DOFs and
legitimately require zero Newton iterations; cube-centred pyramids retain
one free vertex. This is decided by the global logical free-dimension
count in the shared problem, not by a residual-size threshold. Tangent
controls use $n=3$, where every family has free DOFs.
All seven positive-dimensional geometries have local and MPI rank
one-through-four registrations under sequential and OpenMP assembly.

Tangent consistency is checked at degrees one through three, with the
wrong derivative rejected at degree three. The amplitude $4$ omitted-cubic
control uses degree three on the fixed $n=5$ grid, so the correct field
is resolved below the same rejection bounds. Independent assembly,
norm-quadrature and SNES-tolerance controls use degree three on $n=3$,
with the orders, relative budget and solver checks of the p suite.
Slow labels, 30-minute safety timeouts and the pyramid resource lock
bound execution without relaxing mathematical acceptance.
