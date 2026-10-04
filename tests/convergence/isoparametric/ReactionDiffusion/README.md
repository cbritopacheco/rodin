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
\qquad i\in\{0,1\},
$$

with $\kappa=(1,2)$ and $R=\begin{pmatrix}1&0.2\\0.2&1\end{pmatrix}$.
The eigenvalues of $R$ are $0.8$ and $1.2$. The coercive weak form is

$$
a(u,v)=\sum_{i=0}^1\int_\Omega\kappa_i\nabla u_i\cdot\nabla v_i\,dx
+\sum_{i,j=0}^1\int_\Omega R_{ij}u_jv_i\,dx
=\sum_{i=0}^1\int_\Omega f_iv_i\,dx.
$$

Boundary data and sources are prescribed in physical coordinates independently
of discrete coefficients. Shared ReactionDiffusionData supplies
$s(x)=\sum_{\ell=0}^{d-1}x_\ell$, $u_0=e^s$, $u_1=2e^{-s}$ and

$$
f_i=-\kappa_i d u_i+u_i+0.2u_{1-i}.
$$

Both fields use degree $k\in\{1,2\}$ with exact quadratic geometry.
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
E_{i,0}(h)=\left(\int_\Omega|u_{i,h}-u_i|^2\,dx\right)^{1/2},
\qquad
E_{i,1}(h)=\left(\int_\Omega|\nabla u_{i,h}-\nabla u_i|^2\,dx\right)^{1/2}.
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
\frac{\|A_hU_h-b_h\|_2}{\max(1,\|b_h\|_2)}<10^{-11}.
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
