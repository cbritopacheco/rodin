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

## Natural boundaries on the exact quadratic domain

`RodinConvergenceIsoparametricPETScNonlinearPoissonBoundary` reuses the
natural-boundary workload on the exact map

$$
\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1},\qquad
\Omega=\Phi((0,1)^d),\qquad d\in\lbrace 1,2,3\rbrace.
$$

Boundary attributes are assigned in reference coordinates before installing
the map. The mixed cases prescribe the exact trace on
$\Gamma_D=\Phi(\lbrace \xi_0=0\rbrace)$; the remaining boundary is $\Gamma_N$.
The pure-Neumann case has $\Gamma_D=\varnothing$. Physical coordinates,
physical outward normals, and mapped measures define all manufactured data.
For $\beta=0$ (Neumann) or $\beta=1$ (Robin), the problem is

$$
-\Delta u+u+u^3=f\quad\text{in }\Omega,\qquad
u=u_\star\quad\text{on }\Gamma_D,\qquad
\partial_nu+\beta u=g\quad\text{on }\Gamma_N,
$$

where $f=-\Delta u_\star+u_\star+u_\star^3$ and
$g=\nabla u_\star\cdot n+\beta u_\star$. The residual and its derivative
for test and increment functions $v,z$ with zero trace on $\Gamma_D$ are

$$
F(u;v)=\int_\Omega[\nabla u\cdot\nabla v+(u+u^3-f)v]\thinspace\mathrm{d}x
+\int_{\Gamma_N}(\beta u-g)v\thinspace\mathrm{d}s,
$$

$$
D_uF(u)[z;v]=\int_\Omega[\nabla z\cdot\nabla v+(1+3u^2)zv]\thinspace\mathrm{d}x
+\beta\int_{\Gamma_N}zv\thinspace\mathrm{d}s.
$$

Since $1+3u^2\ge1$, the reaction controls constants even for pure Neumann
data: no mean constraint or nullspace gauge is introduced.

The smooth field is $u_\star(x)=\prod_{j=0}^{d-1}\sin(\pi x_j)$.
P1, P2, and P3 check both adjacent refinement intervals in physical $L^2$
and $H^1$ seminorms. P1 uses $n=5,9,17$, except curved tetrahedra,
which use $n=9,17,33$ for all three boundary conditions; P2 and P3 use
$n=3,5,9$. With $h=(n-1)^{-1}$,
the expected orders are $h^{K+1}$ and $h^K$. The observed-rate floors are
$(1.65,0.75)$ for P1 and $(K+0.45,K-0.45)$ for P2/P3. These are tests
of the stated finite hierarchies, not a proof of regularity on arbitrary domains.

The finer tetrahedral P1 path retains the original rate floors. On the coarse
interval $n=5\to9$, mixed Neumann and Robin data gave $L^2$ rates about
$1.611$ and $1.637$, respectively, below the $1.65$ policy. Independent
coarse-level assembly-order, norm-order and solver-tolerance variations changed
the errors by less than $1.5\times10^{-14}$ relatively. Separate local
three-level studies on $n=9,17,33$ recovered $L^2$ rates approximately
$1.877,1.965$ for mixed Neumann and $1.887,1.968$ for Robin; their
$H^1$-seminorm rates were approximately $0.937,0.981$ for both conditions.
These observations support selecting a resolved finite hierarchy rather than
lowering its acceptance bounds. They do not establish the mechanism of every
coarse-grid transient or certify the corresponding MPI and OpenMP runs;
those remain separate verification gates. The flat mixed-Neumann and Robin
hierarchies are unchanged.

Exact physical patches use $u_\star=a$, $a(1+\sum_jx_j)$, and
$a(1+\sum_jx_j^2)$ with $a=1/4$. Their pullbacks are represented by
P1, P2, and P4 respectively; both norm errors must be below $10^{-9}$.
Separate negative controls omit the cubic source response or the natural
flux. Finite-difference residual/tangent comparisons exercise P1/P2 and
reject missing cubic and Robin derivatives; the direction $z(x)=x_0$
vanishes on the mapped Dirichlet face but has nonzero natural-boundary trace.
Assembly order $16\to18$, norm order $18\to20$, and nonlinear tolerance
$10^{-12}\to10^{-13}$ are varied independently, requiring relative error
changes below $10^{-6}$.

There are 37 case selections per geometry and context, covering all seven
positive-dimensional geometries. CTest separates local and MPI ranks one
through four; real-PETSc sequential and OpenMP configurations use the same
cases. Geometry installation is rank-local after explicit mesh distribution;
assembly, solver operations, and global norm reductions retain their global
contracts. The workload does not add collectives to local geometry queries.
Each of the nine boundary/order rate hierarchies has a separate process;
the remaining 28 fixed-mesh controls are grouped separately. All three
levels and both adjacent intervals remain in each rate case. The tetrahedral
P1 entries have a three-hour watchdog: the resolved local pure-Neumann
hierarchy took approximately 6836 seconds. Other entries retain 1800-second
watchdogs. These are scheduling allowances, not numerical budgets or a
performance certificate. MPI processor counts and the common pyramid
resource lock are retained. Full runtime certification of the revised
hierarchy and process partition remains separate from its registration.
The CI workflow additionally assigns each tetrahedral P1 boundary hierarchy
and local/MPI rank count to a separate job. The six P2/P3 hierarchies and
28 fixed-mesh controls form another partition for each context. Both thread
configurations retain every selection; a four-hour job budget accommodates
the three-hour P1 test watchdog and dependency/build time. This scheduling
structure does not assert that the hosted jobs have passed.

## Cubic fields on quadratic approximated geometry

The `ApproximatedP3Q2` extension retains the physical sine solution of
$-\Delta u+u+u^3=f$, its full essential trace and the represented sine-map
domain. Field degree $p=3$ and geometry degree $q=2$ are independent.
Three levels are $n=3,5,9$, except Segment ($n=5,9,17$).
Represented-domain errors and the lifted field, geometry and total defects
are integrated independently. Under the preceding regularity, coercivity
and map assumptions, field errors have $L^2/H^1$ orders $4/3$ and geometry
errors have orders $3/2$.

For norm index $j\in\lbrace0,1\rbrace$, let $F_j,G_j,T_j$ denote field,
geometry and total errors and let $\rho=h_f/h_c<1$. The shared
`LiftedConvergence` assertion checks component-rate windows, both norm
triangle inequalities and the adjacent-level envelope

$$
T_{j,f}<T_{j,c},\qquad
T_{j,f}\le F_{j,c}\rho^{4-j-\delta_j}
             +G_{j,c}\rho^{3-j-\delta_j}+10^{-11},
\qquad\delta_0=0.55,\quad\delta_1=0.45.
$$

These margins, monotonicity and the dimensionless floor are finite-hierarchy
acceptance policies. No two-sided total-error rate or established geometry
dominance is inferred from the component estimates.
At $n=5$, assembly order $12\to16$, norm order $14\to18$, and nonlinear
tolerance $10^{-11}\to10^{-12}$ are varied separately, retaining the
$10^{-6}$ relative sensitivity budget for every positive error component.
A physical affine field is representable at these degrees; represented and
lifted-field errors must remain below $10^{-9}$. Omitting the cubic reaction
while retaining the original source and trace must violate the existing
field-error floors $10^{-3}/10^{-2}$ and double the total errors, without
changing geometry defects.

A fourth case tests the cubic-field residual derivative on the represented
sine-map mesh at $n=3$. The existing homogeneous reference-bubble direction
and central-difference oracle with step $10^{-5}$ are reused. The correct
tangent must have relative defect below $10^{-6}$; changing the cubic
derivative coefficient from $3u^2$ to $u^2$ must produce a defect above
$10^{-3}$. This derivative consistency check is
separate from convergence of the solved field.

Separate registrations cover all seven geometries, native and real-PETSc
local execution, MPI ranks one through four and both thread configurations.
Slow labels, 1800-second watchdogs, processor counts and pyramid resource
locks are retained. The complete finite matrix is locally verified: 84
registrations select 336 configurations and produce 672 rank reports across
the seven geometries, native and real-PETSc local/MPI execution, and both
thread modes. Dependency freshness, case-selection uniqueness and every
participant report are checked independently. This is not hosted-CI
certification or a statement about arbitrary meshes and degrees.
