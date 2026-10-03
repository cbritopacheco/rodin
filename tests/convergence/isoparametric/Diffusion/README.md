# Poisson and conductivity on quadratic geometry

## Continuous problems and reference fields

The fixed physical domain is $\Omega=\Phi((0,1)^d)$, where
$\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1}$ and $d\in\{1,2,3\}$.
The map is regular: its determinant is one for $d=2,3$ and
$1+0.2\xi_0$ for $d=1$. Consequently, $|\Omega|=1$ in dimensions two
and three, and $|\Omega|=1.1$ in dimension one. Exact quadratic geometry
is installed on cells, traces and MPI halos before constructing field spaces.
This is field refinement on a fixed exact curved domain, not a geometry-error
study. A point geometry has no spatial refinement rate.

Two real scalar problems are considered:

$$
-\nabla\cdot(\gamma\nabla u)=f\quad\text{in }\Omega,
\qquad u=g\quad\text{on }\partial\Omega,
\qquad
\gamma(x)=\begin{cases}1&\text{Poisson},\\1+\sum_jx_j&\text{conductivity}.
\end{cases}
$$

Full manufactured Dirichlet data are imposed; both coefficients are positive
on the stated domains. For $V_{h,0}=V_h\cap H^1_0(\Omega)$, the discrete
problem is to find $u_h\in g_h+V_{h,0}$ such that

$$
\int_\Omega\gamma\nabla u_h\cdot\nabla v_h\,dx
=\int_\Omega f v_h\,dx\qquad(v_h\in V_{h,0}).
$$

The common smooth reference, its gradient and its Laplacian are

$$
u(x)=1+\prod_{j=0}^{d-1}\sin(\pi x_j),\qquad
\partial_j u=\pi\cos(\pi x_j)\prod_{k\ne j}\sin(\pi x_k),
\qquad\Delta u=-d\pi^2(u-1).
$$

Sources are derived from the selected continuous coefficient:

$$
f_{\mathrm P}=d\pi^2(u-1),\qquad
f_{\mathrm C}=\gamma d\pi^2(u-1)-\sum_j\partial_j u,
\qquad g=u|_{\partial\Omega}.
$$

The constant patch $u=1$ has $f=0$ for both problems. The physical affine
patch $u=1+\sum_jx_j$ has $f_{\mathrm P}=0$ and $f_{\mathrm C}=-d$.
On quadratic geometry, a physical affine function has a quadratic pullback:
it is tested in P2, not asserted to be exactly representable in P1.
The constant patch is tested in P1.

## Refinement and independent acceptance

Native P1 and H1 degree-two spaces are exercised locally; PETSc uses H1
degrees one and two. Each problem has three refinement levels:
$n=5,9,17$ for P1 and $n=3,5,9$ for P2, with $n$ grid points per
coordinate and nominal spacing $h=(n-1)^{-1}$.
All seven positive-dimensional geometries are included: Segment, Triangle,
Quadrilateral, Tetrahedron, Pyramid, Hexahedron and Wedge.

Physical errors are integrated independently of the assembled equations:

$$
E_0(h)=\|u_h-u\|_{L^2(\Omega)},\qquad
E_1(h)=\|\nabla u_h-\nabla u\|_{L^2(\Omega)}.
$$

Under smoothness, stable conforming approximation, shape regularity and
the dual regularity needed for the L2 estimate, degree $p$ gives
$E_0=O(h^{p+1})$ and $E_1=O(h^p)$. Both adjacent intervals require
finite positive decreasing errors and rates

$$
p+1-0.55<r_0<p+1+0.55,\qquad
p-0.45<r_1<p+0.45,\qquad
r_i=\frac{\log(E_i(h_{\ell-1})/E_i(h_\ell))}
{\log(h_{\ell-1}/h_\ell)}.
$$

These intervals are finite-resolution policies, not proofs of arbitrary-level
asymptotics. Representable patches require $E_0,E_1<10^{-9}$.
Independent map sampling on every local entity, including halos, checks
quadratic transformation order, map error below $10^{-11}$, positive finite
metric factors, and the analytic volume within $10^{-12}$.

## Controls, budgets and architecture

The wrong Poisson operator doubles stiffness without changing source or trace.
For the smooth P2 case at $n=5$, its errors must satisfy
$E_0>0.02$ and $E_1>0.05$. The wrong conductivity operator replaces
$\gamma$ by one while retaining the affine patch source and trace;
it must give $E_0>10^{-3}$ and $E_1>10^{-2}$.
Each wrong error must also exceed twice its correct counterpart.
The absolute bounds are dimensionless policies on the prescribed unit-scale
domains, chosen to separate a model error from discretization and roundoff.

CG uses relative tolerance $10^{-13}$ and at most $50000$ iterations.
Native solver success or a positive PETSc convergence reason is required.
An independently recomputed residual must satisfy

$$
\frac{\|A_hU_h-b_h\|_2}{\max(1,\|b_h\|_2)}<10^{-11}.
$$

PETSc absolute tolerance is $10^{-14}$ and divergence tolerance is $10^5$.
Assembly order eleven, norm order thirteen and solver tolerance $10^{-13}$
are varied separately to sixteen, eighteen and $10^{-14}$, respectively.
At $n=5$, both positive error norms must change relatively by less than
$10^{-6}$ under each variation, for P1 and P2 and for both equations.
No polynomial-exactness claim is made for the mapped trigonometric integrands.

ConductivityData supplies shared physical fields, gradients and sources;
its existing variable-coefficient default is unchanged. Workload owns the
mapped mesh and creates fresh spaces, fields, matrices and solvers per solve.
Both equations share immutable geometry setup within each study, not algebraic
state or a new cache. CurvedGeometry, ErrorNorm and ErrorHistory supply the
existing mapping, physical integration and rate machinery.

Native local and real-PETSc local/MPI configurations are registered.
MPI runs use one to four ranks, including sparse and empty shards on coarse
segment meshes. The continuous dimension is taken from the declared grid
family, independently of the dimension of an empty local shard.
Owned cells alone contribute to globally reduced squared
errors; an independent constant-norm test checks $E_0=\sqrt{|\Omega|}$
and $E_1=0$. Sequential/OpenMP assembly is selected by the build.
Complex-PETSc builds do not register this real-field suite; complex scalar
or vector diffusion coverage is not claimed here.
Tests are labelled slow, have 1800-second timeouts and share a pyramid lock.
