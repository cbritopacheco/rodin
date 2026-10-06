# Poisson and conductivity on exact and approximated curved geometry

## Continuous problems and reference fields

The fixed physical domain is $\Omega=\Phi((0,1)^d)$, where
$\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1}$ and $d\in\lbrace 1,2,3\rbrace $.
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
E_0(h)=\Vert u_h-u\Vert_{L^2(\Omega)},\qquad
E_1(h)=\Vert\nabla u_h-\nabla u\Vert_{L^2(\Omega)}.
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
\frac{\Vert A_hU_h-b_h\Vert_2}{\max(1,\Vert b_h\Vert_2)}<10^{-11}.
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

## Field refinement on approximated sine geometry

The opt-in `Sine` cases use the same equations, analytic physical fields,
sources, solvers and independent physical norms, but change the geometry data
to

$$
\Phi(\xi)=\xi+0.1\sin(\pi\xi_0)e_{d-1},\qquad
\Omega_{h,2}=\Phi_{h,2}((0,1)^d),
$$

where $\Phi_{h,2}$ is the nodal degree-two geometry approximation installed
on cells, traces and halos. The exact determinant is one for $d=2,3$ and
$1+0.1\pi\cos(\pi\xi_0)>0$ for $d=1$. The
[geometry specification](../GeometryApproximation/README.md) separately
defines map and derivative errors against this analytic map at degrees
$q=1,2,3$; it does not reuse a field error as a geometry oracle.

For every represented domain, the manufactured extension
$u(x)=1+\prod_j\sin(\pi x_j)$ defines the continuous comparison problem

$$
-\nabla\cdot(\gamma\nabla u)=f\quad\text{in }\Omega_{h,2},
\qquad g=u|_{\partial\Omega_{h,2}},
$$

with the same independently derived $f_{\mathrm P}$ and $f_{\mathrm C}$
as above. Thus the measured field errors are

$$
E_{0,h}=\Vert u_h-u\Vert_{L^2(\Omega_{h,2})},\qquad
E_{1,h}=\Vert\nabla u_h-\nabla u\Vert_{L^2(\Omega_{h,2})}.
$$

The domain changes with refinement. These quantities do not compare solutions
on different domains and do not include the error of replacing the exact
domain by $\Omega_{h,2}$. The separate lifted affine studies below prescribe
that comparison explicitly. In one dimension, the endpoint
images are fixed at zero and one: only the interior parametrization changes.

Field degrees $p=1,2$ are tested separately from the fixed geometry degree
$q=2$. Three levels remain $n=5,9,17$ for $p=1$ and $n=3,5,9$ for
$p=2$, with $h=(n-1)^{-1}$. Both adjacent intervals require decreasing
positive errors and the same degree-dependent windows above. The sine-map
segment with $p=2$ instead uses $n=5,9,17,33$, testing three intervals:
the two-cell $n=3$ parametrization is pre-asymptotic for this smooth field.
The finer levels retain the original acceptance windows, not relaxed bounds.
Interpretation as $E_{0,h}=O(h^{p+1})$, $E_{1,h}=O(h^p)$ requires uniform approximation,
shape regularity, and the necessary dual regularity across this changing-domain
family; the finite-resolution tests do not prove these hypotheses.

Each physics retains the constant P1 and physical-affine P2 patches, incorrect
operator controls, and separate assembly/norm/solver sensitivity checks with
unchanged bounds. All seven geometries, native and real-PETSc local execution,
sequential/OpenMP assembly, and MPI ranks 1–4 are included. Complex-PETSc and
point geometry retain the exclusions above.

An analytic sine-node oracle checks the degree-two map values independently
of the installation helper on every local cell, trace and halo entity. The
node discrepancy must be below $10^{-11}$ on the unit-scale domains; no
whole-cell exactness is asserted for the sine map. This check must reject
silently using the quadratic default: a field patch or rate alone does not
validate the domain map.

An independent constant-field norm also requires $E_0=1$, $E_1=0$ on the
represented sine domains. Their volume is one: the 1D interval endpoints are
unchanged; in 2D/3D, matching upper and lower graph translations leave the
unit-box cross-section heights unchanged. This checks physical metric weights
and MPI owned-cell reduction, including sparse/empty segment partitions.
The shared workload owns the opt-in map choice; the quadratic default and its
existing exact-map assertions are unchanged. No additional solver, quadrature
cache or parallel mesh factory is introduced.
The approximated-domain cases have distinct CTest registrations, retaining
the slow label, 1800-second per-geometry/rank budget and shared pyramid lock.

## Lifted field, geometry and total errors on the exact domain

Let $\widehat\Omega=(0,1)^d$, $\Omega=\Phi(\widehat\Omega)$ and
$\Omega_{h,q}=\Phi_{h,q}(\widehat\Omega)$ for the sine map above. Define
the lift by correspondence of the reference coordinates, not by inverse
point location or coordinate matching:

$$
u_h^\ell(\Phi(\xi))=u_h(\Phi_{h,q}(\xi)),\qquad
B(\xi)=D\Phi(\xi)^{-T}D\Phi_{h,q}(\xi)^T.
$$

For the same manufactured physical extension $u$ on both domains, define
three functions on $\Omega$, with $x=\Phi(\xi)$ and
$x_h=\Phi_{h,q}(\xi)$:

$$
e_F(x)=u_h(x_h)-u(x_h),\qquad
e_G(x)=u(x_h)-u(x),\qquad
e_T(x)=u_h(x_h)-u(x)=e_F(x)+e_G(x).
$$

Their physical gradients are

$$
\nabla e_F=B(\nabla u_h(x_h)-\nabla u(x_h)),\qquad
\nabla e_G=B\nabla u(x_h)-\nabla u(x),\qquad
\nabla e_T=B\nabla u_h(x_h)-\nabla u(x).
$$

For $X\in\lbrace F,G,T\rbrace $, report $E_{0,X}=\Vert e_X\Vert_{L^2(\Omega)}$ and
$E_{1,X}=\Vert\nabla e_X\Vert_{L^2(\Omega)}$. These norms do not add. Both
the triangle and reverse-triangle inequalities are checked:

$$
|E_{i,F}-E_{i,G}|\le E_{i,T}\le E_{i,F}+E_{i,G},\qquad i\in\lbrace 0,1\rbrace .
$$

The implemented geometry-limited studies solve both equations with the
physical affine field $u(x)=1+\sum_jx_j$, field degree $p=2$ and geometry
degrees $q=1,2$. Its pullback is representable, so $E_{i,F}<10^{-9}$ and
$|E_{i,T}-E_{i,G}|<10^{-9}$. The geometry defects must be positive,
decrease at both adjacent intervals and satisfy the same rate windows with
$q+1$ for the L2 norm and $q$ for the gradient seminorm. Levels are
$n=3,5,9$ except Segment, which uses $n=5,9,17$. In 1D this measures
parametrization/lift error even though both physical domains are $(0,1)$;
it is not a domain-boundary error. These affine checks isolate geometry error;
the smooth-field studies below exercise interacting field and geometry defects.

An independent metric oracle deliberately sets $\Phi_{h,q}(\xi)=\xi$
while keeping the exact sine map and affine reference field. With amplitude
$a=0.1$, its closed-form geometry and total norms are

$$
E_{0,G}=E_{0,T}=\frac{a}{\sqrt{2}},\qquad
E_{1,G}=E_{1,T}=
\begin{cases}
\left((1-(a\pi)^2)^{-1/2}-1\right)^{1/2},&d=1,\\
a\pi/\sqrt{2},&d=2,3.
\end{cases}
$$

This nonvanishing oracle uses $n=3$, norm order eighteen and absolute
dimensionless tolerance $10^{-9}$, independently of the rate checks.
It tests the inverse-transpose gradient chain rule, the exact-domain metric
weight, and MPI reduction. Assembly order, norm order and solver tolerance
are varied independently at $n=5$ for each geometry degree and both physics;
total and geometry norms must change relatively by less than $10^{-6}$,
while field errors remain below the absolute patch budget.

`LiftedErrorNorm` is shared under `tests/convergence`. The opt-in workload
keeps a copy of the original mesh before installing the geometry; ordinary
tests do not incur that copy. Cell correspondence, ordered vertex indices,
and ownership are checked by exact logical indices. Both mesh chart
Jacobians come from cached `Geometry::Point` objects; `SineMap` supplies the
independent analytic exact Jacobian. Integration uses the original cell
distortion multiplied by $\det D\Phi$, not the represented-domain metric.
MPI includes owned original cells, sums six squared contributions globally,
and only then takes square roots. This is scalar real-field coverage on all
seven positive-dimensional geometries, native/real-PETSc local and MPI
ranks 1–4, with sequential/OpenMP assembly. The new `LiftedDiffusion`
registrations retain the slow label, timeout and pyramid lock above.

### Smooth-field lifted convergence and independent controls

The same lift is applied to $u(x)=1+\prod_j\sin(\pi x_j)$ for Poisson
and variable conductivity, using field degrees $p=1,2$ and fixed geometry
degree $q=2$. Unlike the affine patches, all three error components are
nonzero. Under uniform map regularity, approximation stability, solution
smoothness and the dual regularity needed for the L2 estimate, the field
and geometry contributions satisfy the formal estimates

$$
E_{0,F}=O(h^{p+1}),\qquad E_{1,F}=O(h^p),\qquad
E_{0,G}=O(h^{q+1}),\qquad E_{1,G}=O(h^q).
$$

The triangle inequality gives the total-error upper bounds

$$
E_{0,T}=O(h^{\min(p,q)+1}),\qquad
E_{1,T}=O(h^{\min(p,q)}).
$$

These are upper bounds, not lower bounds excluding cancellation. The tests
add a finite-resolution requirement: each of the six measured norms is
finite, positive and decreasing, and every adjacent rate lies in the
two-sided rate windows stated above. The target orders are $p+1,p$ for
the field, $3,2$ for geometry and $\min(p,2)+1,\min(p,2)$ for the total.
All components are measured independently, and both norm inequalities remain
checked at each level. No deduction of total-error rates from field-error
rates alone is made.

Levels are $n=5,9,17$ for P1 and $n=3,5,9$ for P2, except P2 Segment,
which retains $n=5,9,17,33$ and three intervals. The spacing is
$h=(n-1)^{-1}$. These are the same sequences as the represented-domain
sine-field studies, so the new comparison changes the metric and reference
domain, not the discrete PDE or its refinement path.

At $n=5$, assembly order eleven, norm order thirteen and solver tolerance
$10^{-13}$ are varied separately to sixteen, eighteen and $10^{-14}$.
Every field, geometry and total norm must change relatively by less than
$10^{-6}$ for each physics and field degree. The existing residual budget
and convergence-reason checks are retained. The shared sensitivity method
also retains the affine absolute field-error checks, avoiding relative
normalization by a near-zero patch error.

Independent incorrect-operator controls are solved with smooth P2 fields at
$n=5$: Poisson stiffness is doubled, or variable conductivity is replaced
by one, without changing the manufactured source or boundary trace. Their
field and total errors must exceed twice the correct errors and the
absolute wrong-operator bounds stated above. The geometry norms must remain
exactly unchanged, since their definitions contain no numerical solution
or discrete operator. Thus changing the PDE is distinguished from changing
the geometry contribution; neither is hidden inside a total norm.

`SmoothLiftedDiffusion` has separate CTest registrations to preserve the
affine batch's runtime budget. The same seven geometries, real scalar
spaces, native and real-PETSc local/MPI backends, sequential/OpenMP assembly,
rank counts, slow label, timeout and pyramid lock apply. No new solver,
mesh factory, integration cache or library-source path is introduced.
