# Constant-space projection on quadratic geometry

## Spaces and mapped mass problem

The physical domain is $\Omega=\Phi((0,1)^d)$, with
$\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1}$. Exact quadratic maps are installed on
cells, traces and MPI halos. The geometry determinant is one in dimensions
two and three; in dimension one it is $1+0.2\xi_0>0$ and $|\Omega|=1.1$.
This is fixed curved-geometry field refinement, not a geometry-error study
or a strictly isoparametric degree-zero discretization.

For $\mathbb K\in\{\mathbb R,\mathbb C\}$ and scalar or two-component vector
fields, the discontinuous space is

$$
V_h=\{v\in L^2(\Omega;\mathbb K^m):v|_K
\text{ is constant for each cell }K\},\qquad m\in\{1,2\}.
$$

The assembled mass problem computes the orthogonal projection:

$$
\int_\Omega u_h\cdot\overline{v_h}\,dx
=\int_\Omega f\cdot\overline{v_h}\,dx,\qquad v_h\in V_h.
$$

Consequently, $u_h|_K=|K|^{-1}\int_K f\,dx$.
Assigning a function directly to a P0 grid function instead evaluates its
DOF functional at the mapped reference centroid; that is interpolation,
not the volume-weighted projection on a curved cell.

The Workload class owns the mapped mesh and creates a fresh space and system
per projection. Native and PETSc versions share physical data, norm integration,
moment checks and acceptance policies. Backend-specific solver checks remain
explicit. CurvedGeometry, ErrorNorm and NormHistory supply the common machinery;
the scalar/vector/complex squared magnitude is reused from ErrorNorm.
Within each study, value types share the immutable mapped mesh and its existing
geometry quadrature cache. Spaces, coefficients, matrices and solvers are fresh
for every solve; no new cache or production fast path is introduced.

## Manufactured data and independent moments

For real cases, define $c(a,b)=a$; for complex cases, $c(a,b)=a+ib$.
The scalar reference is

$$
f(x)=c(1.25,0.5)+c(1,-0.75)x_0+c(0.3,0.2)x_{d-1},
$$

and the vector reference is

$$
f(x)=\begin{pmatrix}
c(1.25,0.5)+c(1,-0.75)x_0\\
c(-0.75,0.25)+c(0.8,0.2)x_{d-1}
\end{pmatrix}.
$$

Both vector components vary; complex references have nonzero imaginary parts.
Physical error integration and independently evaluated cell moments measure

$$
E_h=\|u_h-f\|_{L^2(\Omega)},\qquad
M_h=\left(\sum_K\frac{1}{|K|}
\left\|\int_K(u_h-f)\,dx\right\|_{\mathbb K^m}^2\right)^{1/2}.
$$

Projection must have $M_h<10^{-10}$ on every refinement level, while
$E_h$ decreases. Moments are evaluated through physical field values and
a higher-order rule, independently of assembled coefficients.
Vector magnitudes sum component moduli; complex errors are not measured
with a nonconjugated square.

All seven positive-dimensional geometries are tested: Segment, Triangle,
Quadrilateral, Tetrahedron, Pyramid, Hexahedron and Wedge.
Three levels use $n=5,9,17$ grid points, or $4,8,16$ subdivisions per
coordinate, with $h=(n-1)^{-1}$. Under map regularity, shape-regular refinement
and smooth data, the expected estimate is $E_h=O(h)$.
Every adjacent interval requires finite positive decreasing errors and

$$
0.8<\frac{\log(E_{h_{\ell-1}}/E_{h_\ell})}
{\log(h_{\ell-1}/h_\ell)}<1.2.
$$

The interval is a finite-resolution acceptance policy, not an asymptotic proof.
No H1-seminorm rate is claimed for discontinuous P0 fields.
A spatial refinement rate is not defined for a point geometry.

## Global constants, controls and numerical budgets

The P0 and P0g spaces must reproduce the constant parts of the references
within $10^{-10}$ in physical L2 norm. P0g is the single global constant
space, so no refinement-driven rate is claimed for nonconstant data.
Its projected value must equal the analytic volume mean of the references.
The coordinate means are

$$
\overline{x_0}=\overline{x_{d-1}}=0.55\quad(d=1),\qquad
\overline{x_0}=0.5,\quad
\overline{x_{d-1}}=0.5+\frac{0.1}{3}\quad(d=2,3).
$$

Substituting these values into the references gives a constant physical
oracle. Its L2 distance from the computed P0g field must be below $10^{-10}$.

Two deliberate wrong computations establish the acceptance gates:
centroid interpolation must give $M_h>10^{-7}$, and doubling the mass
operator without changing the constant load must give $E_h>0.1$.
The bounds are absolute, dimensionless policies on the stated unit-scale
domain. The interpolation control matters because an incorrect centroid
method can still exhibit the expected first-order rate.

Assembly order eight is compared separately with fourteen; norms and moments
use order twelve. CG relative tolerance $10^{-13}$ is compared separately
with $10^{-14}$. Each positive error changes relatively by less than $10^{-6}$,
and moments retain their absolute bound. Every solve also checks

$$
\frac{\|A_hU_h-b_h\|_2}{\max(1,\|b_h\|_2)}<10^{-11}.
$$

Native solver success and positive PETSc convergence reason are required.
CG permits at most $20000$ iterations. PETSc absolute tolerance is $10^{-14}$
and divergence tolerance $10^5$, remote from the tested residual budget.

Native local configurations exercise all four value types. Real-PETSc builds
exercise real scalar/vector fields; native-complex PETSc builds exercise complex
scalar/vector fields. Each PETSc configuration registers local and MPI contexts,
with one to four ranks. Sequential/OpenMP assembly is selected by the build.
Owned cells alone contribute to MPI squared norms and cell moments before
global reduction; halo cells are not counted twice. P0g global DOF ownership
is handled by the existing distributed space, not inferred from cell ownership.

Registrations are labelled slow, use 1800-second timeouts, and share a pyramid
resource lock. No production finite-element, geometry or assembly code is changed.
