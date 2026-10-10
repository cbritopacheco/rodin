# Coupled reaction–diffusion on exact and approximated geometry

## Mathematical problem and representation

Let $\widehat\Omega=(0,1)^d$ and $\Omega=\Phi(\widehat\Omega)$, where
$\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1}$. The exact quadratic map is installed
on cells, traces and MPI halos. Its determinant is one in dimensions two
and three; in dimension one it is $1+0.2\xi_0$, giving $|\Omega|=1.1$.
Geometry remains fixed while the field mesh is refined.

The two physical fields satisfy

$$
-\kappa_i\Delta u_i+\sum_{j=0}^1R_{ij}u_j=f_i
\quad\text{in }\Omega,\qquad u_i=g_i\quad\text{on }\partial\Omega,
\qquad i\in\lbrace 0,1\rbrace,
$$

with $\kappa=(1,2)$ and

$$
R=\begin{pmatrix}1&0.2\cr 0.2&1\end{pmatrix}.
$$

The eigenvalues of $R$ are $0.8$ and $1.2$. The coercive weak form is

$$
a(u,v)=\sum_{i=0}^1\int_\Omega\kappa_i\nabla u_i\cdot\nabla v_i\thinspace dx
+\sum_{i,j=0}^1\int_\Omega R_{ij}u_jv_i\thinspace dx
=\sum_{i=0}^1\int_\Omega f_iv_i\thinspace dx.
$$

Boundary data and sources are prescribed in physical coordinates independently
of discrete coefficients. Shared ReactionDiffusionData supplies
$s(x)=\sum_{\ell=0}^{d-1}x_\ell$, $u_0=e^s$, $u_1=2e^{-s}$ and

$$
f_i=-\kappa_i d u_i+u_i+0.2u_{1-i}.
$$

Both fields use degree $k\in\lbrace 1,2\rbrace$ with exact quadratic geometry.
The first case is superparametric; the second is strictly isoparametric.
Native degree one uses P1; degree two and PETSc use H1.
Workload owns the mapped mesh and creates a fresh two-field system per solve.
Data and observables are shared; solver/storage policy remains explicit.
No production assembly or solver implementation is changed.

## Refinement and error oracles

All seven positive-dimensional geometries are exercised: Segment, Triangle,
Quadrilateral, Tetrahedron, Pyramid, Hexahedron and Wedge. A spatial diffusion
rate has no counterpart on a point geometry.

| Field degree | Grid points per coordinate | Subdivisions per coordinate |
| --- | --- | --- |
| $P_1$ | $5,9,17$ | $4,8,16$ |
| $P_2$ | $3,5,9$ | $2,4,8$ |

With $h=(n-1)^{-1}$, each component has independent physical error measurements

$$
E_{i,0}(h)=\left(\int_\Omega|u_{i,h}-u_i|^2\thinspace dx\right)^{1/2},
\qquad
E_{i,1}(h)=\left(\int_\Omega|\nabla u_{i,h}-\nabla u_i|^2\thinspace dx\right)^{1/2}.
$$

For smooth solutions, regular maps, shape-regular refinement and the requisite
dual regularity, the expected orders are $E_{i,0}=O(h^{k+1})$ and
$E_{i,1}=O(h^k)$. These hypotheses motivate the numerical gates; finite runs
do not prove an asymptotic theorem. Every adjacent interval and both components
must have finite, positive, strictly decreasing errors. The observed rate is

$$
r_{i,m,\ell}=
\frac{\log(E_{i,m}(h_{\ell-1})/E_{i,m}(h_\ell))}
{\log(h_{\ell-1}/h_\ell)}.
$$

The finite-resolution policy windows are
$|r_{i,0,\ell}-(k+1)|<0.55$ and $|r_{i,1,\ell}-k|<0.45$.
They are test acceptance margins, not mathematical bounds.

## Independent controls and backends

Constant fields $u_i=i+1$ form representable $P_1$ patches.
Physical affine fields $u_i=(i+1)(1+s)$ have quadratic pullbacks and
form representable $P_2$ patches. Both norms must be below $10^{-9}$.
A physical quadratic field is not claimed to be a mapped quadratic patch.

On the affine patch, both off-diagonal reaction terms are omitted while
sources and boundary data remain unchanged. Each component must then have
$E_{i,0}>10^{-3}$ and $E_{i,1}>10^{-2}$, whereas the correct coupled solve
remains below the patch tolerance. These absolute thresholds concern
dimensionless fields and gradients on the stated unit-scale domain.

Assembly order eleven and norm order thirteen are compared separately with
orders sixteen and eighteen. CG relative tolerance $10^{-13}$ is compared
separately with $10^{-14}$; each positive error changes relatively by less
than $10^{-6}$. Mapped/nonpolynomial integrands are not claimed to be
polynomial-exact. Every solve also checks

$$
\frac{\Vert A_hU_h-b_h\Vert_2}{\max(1,\Vert b_h\Vert_2)}<10^{-11}.
$$

Native solver success and PETSc positive convergence reason are required.
The iteration cap is $50000$; PETSc absolute tolerance is $10^{-14}$ and
divergence tolerance is $10^5$, remote from the tested residual budget.

CMake registrations cover native local and real-PETSc local/MPI contexts,
with one to four MPI ranks. Sequential/OpenMP assembly is selected by the
build configuration. Complex-PETSc builds do not register this real suite.
In MPI, squared norms integrate owned cells only and are globally summed
before taking square roots. An analytic-volume control measures the error
of the zero field against each constant $i+1$, requiring
$E_{i,0}=(i+1)\sqrt{|\Omega|}$ within $10^{-12}$ and $E_{i,1}=0$.
This detects double-counted halo cells independently of PDE convergence.

The suite is labelled slow; each geometry registration has a 1800-second
timeout, and pyramid registrations share a resource lock.

## Represented and lifted sine-map studies

The approximated-domain cases instead use the exact map
$\Phi(\xi)=\xi+0.1\sin(\pi\xi_0)e_{d-1}$ and its quadratic
interpolant $\Phi_h$. Physical sources and traces retain the same exact
fields on $\Omega_h=\Phi_h((0,1)^d)$. This separates field approximation
from geometric displacement of the domain.
In dimensions two and three, the first coordinate is unchanged and the
exact-map determinant equals one. The interpolant's determinant is not
assumed to be identically one on every cell family; the shared geometry
study independently checks positive represented Jacobians at its integration
points. Uniform regularity remains a hypothesis of the rate interpretation.
In one dimension,
$\Phi'=1+0.1\pi\cos(\pi\xi)>0$. The quadratic interpolant satisfies
$\Phi_h'\ge1-0.3\pi>0$: its endpoint derivatives are twice a
half-interval secant minus the full-interval secant, bounded by $3\pi$
for the sine contribution. The piecewise-linear derivative of the quadratic
interpolant lies between those endpoints. Thus the specified hierarchy
retains regular maps; no assumption about interpolation overshoot is needed.

For $x=\Phi(\xi)$ and $x_h=\Phi_h(\xi)$, each component has defects

$$
e_{i,F}(x)=u_{i,h}(x_h)-u_i(x_h),\qquad
e_{i,G}(x)=u_i(x_h)-u_i(x),\qquad
e_{i,T}=e_{i,F}+e_{i,G}.
$$

The shared `LiftedErrorNorm` integrates these quantities on the exact domain
and applies $D\Phi^{-T}D\Phi_h^T$ to physical gradients before forming
derivative defects. It pairs cells through retained logical indices and
ordered vertices, without inverse point location. Only owned reference cells
contribute to MPI norms. `LiftedConvergence` applies the same decomposition,
adjacent-rate and sensitivity policies to each scalar field; the two-field
solve and manufactured data remain specific to this formulation.

For each field, the expected displacement/gradient orders are $(p+1,p)$
for represented-domain and lifted-field errors, $(3,2)$ for geometry defects,
and $(\min(p,2)+1,\min(p,2))$ for total errors, under the same regularity
and geometry assumptions as above. Both triangle inequalities are checked
with absolute roundoff allowance $10^{-11}$; errors are not assumed to add
in norm. Every adjacent interval retains the margins $0.55$ and $0.45$.
P1 uses $n=5,9,17$; P2 uses $n=3,5,9$, except Segment uses
$n=5,9,17,33$ to resolve the coarse pre-asymptotic regime.

At $n=5$, assembly order $11\to16$, norm order $13\to18$ and solve
tolerance $10^{-13}\to10^{-14}$ are changed independently. All four
error components for both fields must change by less than $10^{-6}$
relatively. The physical-affine P2 patch retains represented and lifted-field
errors below $10^{-9}$. Omitting both off-diagonal reaction terms, while
retaining the correct sources and traces, must violate the absolute field
error floors above and increase each total norm by a factor greater than two.
The geometry defect must remain identical to the correct solve.

All seven geometries have separate `Approximated` registrations for native,
real-PETSc local and MPI ranks one through four. Sequential/OpenMP builds
use identical acceptance logic. These cases retain slow-test labels,
1800-second timeouts and the shared pyramid resource lock.

## Linear and cubic geometry with representable coupled fields

Additional sine-map studies use geometry degree
$q\in\lbrace 1,3\rbrace$ and field degree $p=\max(2,q)$.
The affine physical fields are

$$
u_i(x)=(i+1)\left(1+\sum_jx_j\right),\qquad
f_i(x)=u_i(x)+0.2u_{1-i}(x),\qquad i\in\lbrace 0,1\rbrace.
$$

Both pullbacks belong to the represented field family, whereas the exact
sine map is not represented by any of these finite geometry degrees.
The coupled operator, diffusion coefficients and reaction matrix remain
unchanged. Each field must independently satisfy represented and lifted
field errors below $10^{-9}$; its total and geometry norms must agree
within the same absolute budget. Agreement of a combined two-field norm
is not used as a substitute for either component check.

| Geometry degree $q$ | Field degree $p$ | Grid points $n$ | Nominal geometry L2 / H1-seminorm orders |
| --- | --- | --- | --- |
| 1 | 2 | $3\to5\to9$; segment $5\to9\to17$ | $2/1$ |
| 3 | 3 | $3\to5\to9$; segment $5\to9\to17$ | $4/3$ |

For both geometry and total errors, each adjacent interval requires
finite, positive, decreasing norms, with L2 rates within $0.55$ of
$q+1$ and H1-seminorm rates within $0.45$ of $q$. The gradients
of the chosen affine fields are nonzero in the perturbed coordinate,
so this experiment measures a geometry defect rather than a constant-field
null case. The finite rate windows are numerical policies, not general
two-sided interpolation estimates.

At $n=5$, assembly order $11\to16$, norm order $13\to18$ and
CG tolerance $10^{-13}\to10^{-14}$ are varied separately.
Field reproduction remains required at every setting. Each component's
geometry and total norms must change by less than $10^{-6}$ relatively;
near-zero field errors are not used as relative denominators.
The shared representable-field path in `LiftedConvergence` owns these
decomposition, rate and sensitivity checks; this fixture owns the coupled
solve and analytic data. The independent omitted-coupling controls remain
separate.

The [segment regularity argument](../Helmholtz/README.md#linear-and-cubic-geometry-with-a-representable-complex-field)
applies to these linear/cubic sine interpolants. In one dimension the
physical interval is unchanged, and the lift compares parametrizations.
In higher dimensions, regular represented maps remain a hypothesis of
the geometry-rate interpretation.

Degree-specific entries cover all seven positive-dimensional geometries
in native and real-PETSc local/MPI ranks one through four, separately
under sequential and OpenMP assembly. Slow labels, 1800-second watchdogs,
MPI processor counts and the pyramid resource lock are retained.
Registration alone does not certify a numerical run.

## Natural boundaries on the exact quadratic domain

The real-PETSc natural-boundary target reuses the flat boundary fixture and
physical reaction–diffusion data, with the exact quadratic map installed after
reference-boundary classification. Attributes and MPI ownership retain their
logical identities; no physical-coordinate matching or new collective query
is used. Global assembly, solution and error norms retain collective semantics.

Mixed Neumann and Robin cases prescribe essential data on the mapped image
of $\lbrace \xi_0=0\rbrace$ and use the complementary boundary $\Gamma_N$. Pure
Neumann uses $\Gamma_N=\partial\Omega$. With $\beta=0$ for Neumann and
$\beta=1$ for Robin, the physical data and weak form are

$$
g_i=\kappa_i\nabla u_i\cdot n+\beta u_i,\qquad
a(u,v)+\beta\sum_{i=0}^1\int_{\Gamma_N}u_iv_i\thinspace ds
=\sum_{i=0}^1\left(\int_\Omega f_iv_i\thinspace dx+
\int_{\Gamma_N}g_iv_i\thinspace ds\right).
$$

The positive reaction bound $\lambda_{\min}(R)=0.8$ controls constants,
including pure Neumann cases; no artificial mean constraint is introduced.
Physical unit normals and surface measures come from the mapped geometry.
Each component is checked independently in both error norms.

Degrees $p=1,2,3$ use respectively $n=5,9,17$ and $n=3,5,9$ for each
higher degree. Both adjacent intervals require strict error reduction and
L2/H1-seminorm rate floors $1.65/0.75$, $2.45/1.55$, and $3.45/2.55$.
Physical affine and quadratic patches use field degrees two and four:
their pullbacks under the degree-two map have degrees at most two and four,
respectively. Correct patches require $E_{i,0},E_{i,1}<10^{-9}$.
Omitting either cross-coupling or normal flux on the affine patch must give
$E_{i,0},E_{i,1}>10^{-3}$; the flux control also retains the fivefold
separation from the correct error, without coarse smooth-field contamination.
The Robin exact-value load is retained in the omitted-normal-flux control.

Assembly order $16\to18$, norm order $18\to20$, and CG relative tolerance
$10^{-13}\to10^{-14}$ are varied independently on the smooth P2 field,
with relative changes below $10^{-6}$ in both component norms. Registrations
cover all seven geometries locally and at MPI ranks one through four,
separately under sequential and OpenMP assembly. These registrations specify
the validation matrix; their presence alone does not certify a passing run.

The natural-boundary targets retain slow labels and the common pyramid
resource lock. Pyramid contexts have a 14400-second CTest budget and are
separated into local and individual MPI rank-count jobs, each with a
300-minute CI budget including configuration and compilation. Other
natural-boundary contexts retain 1800 seconds and their existing CI
partitions. This scheduling policy preserves all 27 cases per context,
refinement levels, quadrature and numerical acceptance. The completed
sequential local pyramid group took approximately 139 minutes with other
validation work active; MPI-1 exceeded 157 minutes before its final boundary
family. These are scheduling observations, not isolated benchmarks. The
complete local matrix has passed under both assembly modes, with Local and
MPI ranks one through four on all seven geometries. Hosted-CI success and
runner-specific budget headroom remain separate checks.

An additional P1 quadrature-adequacy study uses $n=5$ for each natural
boundary variant and geometry. It compares assembly order 12 and norm
order 12 against the existing 16/18 reference, varying assembly and norms
separately before testing their combined change. For each component and
each nonzero norm $E$, it requires

$$
\left|\frac{E_{12,18}}{E_{16,18}}-1\right|<10^{-6},\qquad
\left|\frac{E_{16,12}}{E_{16,18}}-1\right|<10^{-6},\qquad
\left|\frac{E_{12,12}}{E_{16,18}}-1\right|<10^{-6}.
$$

CG tolerance remains $10^{-13}$. This independent candidate study does
not change the quadrature settings or acceptance bounds of the rate tests,
and a coarse-mesh comparison alone does not establish adequacy at every
refinement level. It adds three cases, giving 27 case selections per
geometry/context. Any subsequent adoption of cheaper settings must retain
the three-level convergence checks and justify their integration budgets.

Polynomial moment exactness is not a sufficient criterion for these
integrands. For example, the reference pyramid approximation space contains
the rational mode

$$
\psi(r,s,t)=\frac{rs}{1-t},\qquad
0\le r,s\le1-t,\quad0\le t<1.
$$

Its gradients and their mapped products are not polynomial moments in
$(r,s,t)$. Quadrature adequacy is therefore checked on the actual field
errors, independently of the nominal polynomial order of the rule.

## Cubic fields on quadratic approximated geometry

The `ApproximatedP3Q2` extension retains the coupled smooth fields, unequal
diffusion coefficients, reaction matrix, sources and essential traces, but
uses field degree $p=3$ on quadratic sine-map geometry ($q=2$). Three
grid-point levels are $n=3,5,9$, except Segment ($n=5,9,17$). Existing
P1/P2 sequences and acceptance rules remain unchanged.

Each field $i\in\lbrace 0,1\rbrace$ has independent represented-domain,
lifted field, geometry and total error measurements and a separate
`LiftedConvergence` history. Under the preceding coercivity, regularity and
approximation hypotheses, the expected $L^2/H^1$ orders are $4/3$ for
field errors and $3/2$ for geometry errors. The shared mixed-order rule in
the [general methodology](../../README.md#acceptance-and-reproducibility)
checks every adjacent interval, both component-rate windows and both norm
triangle inequalities. Its total-error envelope uses the independently
measured field and geometry errors; geometry dominance is not presumed.
The rate margins remain $0.55/0.45$ and the dimensionless absolute floor
remains $10^{-11}$. Strict total-error decrease is an additional policy,
not a consequence of the triangle inequality. A combined two-field norm
cannot replace either field's checks.

At $n=5$, assembly order $11\to16$, norm order $13\to18$, and solver
tolerance $10^{-13}\to10^{-14}$ are varied separately. Every positive
error component retains the $10^{-6}$ relative sensitivity budget. Physical
affine cubic-field patches remain representable on quadratic geometry and
retain the $10^{-9}$ absolute reproduction budget. Omitting the off-diagonal
coupling with unchanged sources and traces must violate the existing field
error floors and increase each total norm by a factor greater than two;
geometry errors must remain unchanged.

Separate registrations cover seven geometries in native and real-PETSc
local contexts, MPI ranks one through four, and both thread configurations.
Slow labels, 1800-second watchdogs and pyramid locks are retained. The older
approximated groups exclude these additions. The complete finite matrix is
locally verified: 84 registrations select 252 configurations, with 504
successful rank-level reports. Both thread configurations have been built
with Clang and syntax-checked with GCC. Registration selection, source
identity, runtime reports and build dependency freshness are independently
checked. Sampled RSS guards are not continuous peak-memory measurements;
hosted CI and behavior outside these finite hierarchies remain separate.
The natural-boundary matrix uses a different target and is unaffected.
