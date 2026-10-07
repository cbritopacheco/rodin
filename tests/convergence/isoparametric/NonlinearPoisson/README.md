# Semilinear Poisson on exact and approximated geometry

## Problem, boundary lifting and architecture

On $\Omega=\Phi((0,1)^d)$ with the exact quadratic map
$\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1}$, the physical problem is

$$
-\Delta u+u+u^3=f,\qquad
u=g\ \text{on }\partial\Omega.
$$

The common physical reference is $u(x)=\prod_{j=0}^{d-1}\sin(\pi x_j)$,
with analytic gradient and $f=(d\pi^2+1)u+u^3$. Its trace is generally
nonzero on the mapped boundary. The discrete lifting $g_h=I_hu$ is the
finite-element interpolant, not a mass projection. The unknown is a
homogeneous correction $w_h$ and $u_h=g_h+w_h$.

For $v_h$ with homogeneous trace, the residual and tangent are

$$
F_h(w_h)[v_h]=\int_\Omega
\nabla(g_h+w_h)\cdot\nabla v_h+
(g_h+w_h+(g_h+w_h)^3-f)v_h\thinspace dx,
$$

$$
J_h(w_h)[z_h,v_h]=\int_\Omega
\nabla z_h\cdot\nabla v_h+
(1+3(g_h+w_h)^2)z_hv_h\thinspace dx.
$$

The reaction derivative is positive, so the tangent is coercive.
Native Newton begins from $g_h$ and uses homogeneous increments.
PETSc SNES begins from zero correction; each callback reconstructs
$g_h+w_h$, including synchronized MPI halos. Neither path changes the
prescribed trace during iteration. The shared test workloads are extended
with opt-in lifting and representable data; existing flat-mesh defaults
remain homogeneous sine data with no lifting. Production solvers are unchanged.

Workload owns the mesh and installs exact P2 maps before constructing
spaces. Native and PETSc use scalar H1 elements of degree $k\in\lbrace 1,2\rbrace$.
P1 fields on quadratic geometry are superparametric; P2 is strictly
isoparametric. Maps are installed on cells, boundary traces and MPI halos.

## Refinement and acceptance

All seven positive-dimensional geometries are covered: Segment, Triangle,
Quadrilateral, Tetrahedron, Pyramid, Hexahedron and Wedge. A spatial
diffusion refinement study is not defined on a point geometry.

| Field degree | Grid points per coordinate | Subdivisions |
| --- | --- | --- |
| $P_1$ | $5,9,17$ | $4,8,16$ |
| $P_2$ | $3,5,9$ | $2,4,8$ |

With $h=(n-1)^{-1}$, independent physical integration measures

$$
E_0(h)=\Vert u_h-u\Vert_{L^2(\Omega)},\qquad
E_1(h)=|u_h-u\rvert_{H^1(\Omega)}.
$$

Under smoothness, uniform map regularity, shape-regular refinement and
the requisite dual regularity, the expected orders are $k+1$ and $k$.
Every adjacent interval must have finite, positive, decreasing errors and

$$
r_{m,\ell}=
\frac{\log(E_m(h_{\ell-1})/E_m(h_\ell))}
{\log(h_{\ell-1}/h_\ell)},\qquad
|r_{0,\ell}-(k+1)|<0.55,\quad |r_{1,\ell}-k|<0.45.
$$

These are finite-resolution acceptance windows, not a proof of the
asymptotic theorem. MPI norms sum squared contributions on owned cells
only before taking the global square root.

## Patches, derivatives and independent controls

The constant $u=1$ is a P1 patch. Physical affine
$u=1+\sum_jx_j$ has quadratic pullback and is a P2 patch.
For either, $f=u+u^3$, and both physical error norms must be below $10^{-9}$.
These patches may converge without a Newton correction because their
interpolants already solve the discrete problem.

The tangent check uses the physical state $I_hu$ and the homogeneous
direction $z_h=\tfrac14 I_h(\prod_j\sin(\pi\xi_j))$. The original mesh chart
supplies $\xi$ independently of the physical field data. With
$\epsilon=10^{-5}$, the central-difference defect is

$$
D=\frac{\Vert J_hz_h-
(F_h(w_h+\epsilon z_h)-F_h(w_h-\epsilon z_h))/(2\epsilon)\Vert_2}
{\Vert(F_h(w_h+\epsilon z_h)-F_h(w_h-\epsilon z_h))/(2\epsilon)\Vert_2}.
$$

Both P1 and P2 require $D<10^{-6}$. Replacing $3u_h^2$ by $u_h^2$
in the P2 tangent must yield $D>10^{-3}$.

Two separate affine-patch controls retain the exact physical source:
omitting the cubic reaction, or omitting the nonzero boundary lifting.
Each wrong problem must have $E_0>10^{-3}$ and $E_1>10^{-2}$.
The quantities are dimensionless on the stated unit-scale domains.
These controls separate physical correctness from a self-consistent tangent.

Quadrature order twelve is compared separately with sixteen, and nonlinear
tolerance $10^{-11}$ with $10^{-12}$ at both degrees. Relative changes in
each positive norm must be below $10^{-6}$. Mapped/nonpolynomial integrands
are not claimed polynomial-exact. Native Newton uses SparseLU; PETSc SNES
uses Newton line search with CG/Jacobi tangents. Both require nonlinear
convergence and independently reassemble the final residual.
Native residual normalized by $\max(1,\Vert F_h(w_{h,0})\Vert_2)$ is below $10^{-9}$;
PETSc requires $10^{-10}$. Linear solver success is checked whenever a
correction was computed; exact initial patches need no linear solve.

CMake registers native local and real-PETSc local/MPI suites with one to four
ranks. Sequential/OpenMP is selected by the build. Complex-PETSc builds do
not register this real suite. Registrations are labelled slow, have
1800-second timeouts, and serialize pyramid cases through a resource lock.

## Approximated geometry and exact-domain errors

Separate cases interpolate $\Phi(\xi)=\xi+0.1\sin(\pi\xi_0)e_{d-1}$
with degree-two geometry $\Phi_h$, defining $\Omega_h=\Phi_h((0,1)^d)$.
The same physical sine field, source and essential trace are prescribed on
$\Omega_h$. The exact-domain comparison uses $x=\Phi(\xi)$,
$x_h=\Phi_h(\xi)$ and the defects

$$
e_F(x)=u_h(x_h)-u(x_h),\qquad
e_G(x)=u(x_h)-u(x),\qquad e_T=e_F+e_G.
$$

The shared lift integrates on the exact domain, transforms gradients by
$D\Phi^{-T}D\Phi_h^T$, and counts owned reference cells only in MPI.
Both triangle inequalities are checked with an absolute allowance of
$10^{-11}$; norms are not assumed to add. Logical cell indices and ordered
vertices preserve chart correspondence without coordinate-tolerance matching.
The regularity argument for these exact and quadratic-interpolated sine maps
is given in the [coupled diffusion specification](../ReactionDiffusion/README.md).

Under the stated regularity hypotheses, represented and lifted-field errors
have orders $(p+1,p)$ in L2/H1 seminorm; geometry defects have orders $(3,2)$;
total errors have orders $(\min(p,2)+1,\min(p,2))$. Every adjacent interval
checks positive, finite, decreasing errors within the original $0.55$ and $0.45$
rate windows. P1 uses $n=5,9,17$ and P2 uses $n=3,5,9$, except
Segment uses $n=5,9,17,33$ to resolve its coarse pre-asymptotic regime.

Native Newton and PETSc SNES share a post-solve measurement callback. It
receives the converged physical state and analytic data while their lifetimes
are valid; SNES synchronizes the state before invoking it. This supplies
lifted norms without reimplementing nonlinear assembly or retaining a field
beyond its space lifetime. Existing solve defaults are unchanged.

At $n=5$, assembly order $12\to16$, norm order $14\to18$ and
nonlinear tolerance $10^{-11}\to10^{-12}$ are varied independently;
all represented/lifted components must change by less than $10^{-6}$
relatively. A physical-affine P2 patch requires represented and lifted-field
errors below $10^{-9}$. Omitting the cubic term while retaining the correct
source and trace must exceed the dimensionless field-error floors above and
increase both total norms by more than a factor of two; geometry errors must
remain identical. At $n=3$, the homogeneous reference-coordinate direction
also checks residual/tangent consistency on the sine-map mesh for both field
degrees, retaining the deliberately incorrect P2 tangent control.

All seven geometries use separate `Approximated` CTest registrations in
native and real-PETSc local/MPI configurations, including ranks one through
four and sequential/OpenMP builds. Slow labels, timeouts and pyramid locking
are retained. These checks concern the stated finite hierarchies, not a
uniform asymptotic result for arbitrary curved meshes.

## Linear and cubic geometry: representable-field studies

The `LocalQ1Test`, `LocalQ3Test`, `MPIQ1Test` and `MPIQ3Test` suites
interpolate the same sine map with geometry degree $q\in\lbrace1,3\rbrace$.
For field degree $k=\max(2,q)$, the physical affine solution

$$
u_\ast(x)=1+\sum_{j=0}^{d-1}x_j,\qquad
f(x)=u_\ast(x)+u_\ast(x)^3,\qquad g=u_\ast\rvert_{\partial\Omega_h}
$$

is represented on the discrete geometry: $u_\ast\rvert_{\Omega_h}\in V_h$.
Equivalently, $u_\ast\circ\Phi_h\in\widehat V_{h,k}$, where
$\widehat V_{h,k}$ denotes the field space on the original unit-box mesh.
Consequently, the represented-domain and lifted-field defects vanish in
exact arithmetic, whereas

$$
e_G(\Phi(\xi))=\sum_{j=0}^{d-1}
\left(\Phi_{h,j}(\xi)-\Phi_j(\xi)\right),\qquad e_T=e_F+e_G
$$

isolates the map interpolation defect. The same physical problem is used
at every level. The shared `LiftedConvergence` representable-field path
requires both field norms below the dimensionless absolute budget $10^{-9}$,
checks the triangle inequalities, and compares total and geometry norms
within this budget. Field errors are not assigned a logarithmic rate.

| Geometry degree $q$ | Field degree $k$ | Segment grid points | Other geometry grid points |
| --- | --- | --- | --- |
| $1$ | $2$ | $5,9,17$ | $3,5,9$ |
| $3$ | $3$ | $5,9,17$ | $3,5,9$ |

For $h=(n-1)^{-1}$, the geometry and total defects have expected orders
$q+1$ in $L^2$ and $q$ in the $H^1$ seminorm under the stated map
regularity and nonvanishing-leading-defect hypotheses. Both adjacent
intervals require finite, positive, decreasing errors and the same
$0.55$ and $0.45$ rate windows. The linear/cubic map regularity argument
and the distinction between changing domains and one-dimensional
parametrization comparisons are given in the
[Helmholtz specification](../Helmholtz/README.md).

At $n=5$, assembly order $12\to16$, norm order $14\to18$ and
nonlinear tolerance $10^{-11}\to10^{-12}$ are varied independently.
Every solve retains the field-error budgets; each positive geometry and
total norm changes by less than $10^{-6}$ relatively.

An exactly represented affine field can satisfy the discrete residual at
the initial iterate, without a Newton correction. A separate case at
$n=3$ therefore exercises the actual residual and tangent at field degree
$k$, with a homogeneous reference-coordinate bubble direction. The
central-difference defect must be below $10^{-6}$, while replacing
$3u_h^2$ by $u_h^2$ must produce a defect above $10^{-3}$. This is a
derivative consistency check, not a claim about Newton iteration counts.

All seven positive-dimensional geometries are registered for native local
and real-PETSc local/MPI configurations, with ranks one through four.
Sequential/OpenMP remains a build choice. These are finite hierarchy
certificates; arbitrary map degrees and arbitrary curved meshes are not
certified by these cases.
