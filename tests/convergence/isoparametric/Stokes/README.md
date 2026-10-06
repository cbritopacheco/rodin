# Curved Taylor–Hood Stokes verification

This suite tests the velocity/pressure pair $P_2/P_1$, together with a
global $P_{0g}$ pressure-mean multiplier, on exact quadratic geometry and
quadratic approximations of a nonpolynomial geometry map. In the exact-map
study, the physical domain is fixed across all refinement levels:

$$
\Omega=\Phi((0,1)^d),\qquad
\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1},\qquad d\in\lbrace 2,3\rbrace .
$$

All six applicable cell families are covered: triangle, quadrilateral,
tetrahedron, pyramid, hexahedron and wedge. The incompressible spatial
velocity/pressure problem is not asserted in zero or one dimension.
Native Eigen and real-PETSc local/MPI workloads reuse the shared Stokes
problem classes; sequential and OpenMP assembly consume the same forms.
MPI registrations cover one through four ranks.

## Physical formulation and gauge

With unit viscosity, the manufactured problem is

$$
-\Delta u+\nabla p=f,\qquad \operatorname{div}u=0,
\qquad u|_{\partial\Omega}=g,\qquad \int_\Omega p\,\mathrm{d}x=0.
$$

Since $\det D\Phi=1$ and $x_0=\xi_0$, the pressure fields retain their
analytic mean under this map. In particular, the rate data are

$$
u(x)=x_1^3e_0,\qquad p(x)=x_0^2-\tfrac13,
\qquad f(x)=(-6x_1+2x_0)e_0.
$$

The pressure gauge follows by change of variables:

$$
\int_\Omega p\,\mathrm{d}x
=\int_{(0,1)^d}(\xi_0^2-\tfrac13)\,\mathrm{d}\xi=0.
$$

The exact volume and pressure mean are independently integrated on the
mapped cells. A fresh global mean-multiplier space enforces the computed
pressure gauge; no pressure point is pinned. Every solve checks the
algebraic residual and computed pressure mean before field-error integration.

## Refinement and acceptance

| Dimension | Grid points per axis | Nominal spacing | Fixed field/geometry degrees |
| --- | --- | --- | --- |
| 2 | $3\to5\to9$ | $h=1/(n-1)$ | Velocity 2, pressure 1, geometry 2 |
| 3 | $3\to4\to5$ | $h=1/(n-1)$ | Velocity 2, pressure 1, geometry 2 |

Each sequence has three levels, and every adjacent interval is checked.
The smaller three-dimensional grids bound the cost of mixed direct solves;
the rate denominator uses the actual spacing ratio, not an assumed factor two.

Subject to mapped-pair stability and the requisite solution/dual regularity,
the expected velocity L2/H1-seminorm orders are three/two, and the pressure
L2/H1-seminorm orders are two/one. Every field error must be finite,
positive and strictly decreasing. The L2 acceptance window is within 0.55
of its expected order, and the derivative window is within 0.45.
These are finite-resolution acceptance policies, not uniform inf-sup proofs
for arbitrary curved meshes or degrees.

The independently integrated divergence is controlled through

$$
\Vert \operatorname{div}u_h\Vert _{L^2(\Omega)}
\le\sqrt d\,|u-u_h|_{H^1(\Omega)}.
$$

This follows from $\operatorname{div}u=0$ and the Frobenius trace bound.
It ties the divergence error to the decreasing velocity derivative error.
A strict divergence rate is not imposed: cancellation can make this
quantity vanish or reach roundoff before the field errors do.
An MPI interpolation control with $u_h=x_0e_0$ independently checks
$\Vert \operatorname{div}u_h\Vert _{L^2(\Omega)}=1$, counting owned cells once.

## Patch, negative and quadrature controls

The physical-affine patch uses $u=x_1e_0$ and $p=x_0-1/2$.
In two dimensions its velocity pulls back to
$(\xi_1+0.1\xi_0^2)e_0$, which belongs to the mapped P2 space.
The pressure pulls back to P1 in both dimensions. All velocity, pressure
and divergence errors must be below $10^{-9}$. A physical quadratic
velocity is not used as an exact P2 patch, because its pullback can have
degree four in two dimensions.

A wrong-viscosity control solves with viscosity two but retains the
unit-viscosity forcing and trace for $u=x_1^2e_0$, $p=x_0-1/2$.
Its pressure L2 error must exceed 0.1 and its pressure H1-seminorm error
must exceed one. In three dimensions the exact wrong pressure has gradient
$3e_0$ instead of $e_0$; in two dimensions an additional field-approximation
error is present. The control establishes rejection of an incorrect operator,
not a rate claim for that operator.

Assembly and field-norm integration use quadrature order 12. A separate
order-16 solve/integration must change each nonzero velocity and pressure
norm by less than $10^{-6}$ relatively. Mapped derivative integrands are
not assumed polynomial-exact. The coefficient residual obeys
$\Vert Az-b\Vert _2/\max(1,\Vert b\Vert _2)<10^{-11}$, where the coefficient vector
collects velocity, pressure and mean-multiplier unknowns. The computed pressure integral
has absolute magnitude below $10^{-10}$, in this dimensionless setting.
Independent analytic-volume/gauge checks use an absolute $10^{-12}$ bound.

Tests are labelled `convergence;slow`, with `petsc` and `distributed`
where applicable. Geometry-specific registrations have 30-minute safety
timeouts; pyramid registrations share a resource lock. Performance benchmarks
remain separate from these convergence assertions.

## Approximated domains and exact-domain lifts

The second study uses the exact reference map and its nodal quadratic
interpolant,

$$
\Phi(\xi)=\xi+0.1\sin(\pi\xi_0)e_{d-1},\qquad
\Phi_h=I_2\Phi,\qquad \Omega_h=\Phi_h((0,1)^d).
$$

Manufactured data are evaluated in physical coordinates on $\Omega_h$:

$$
u(x)=x_{d-1}^3e_0,\qquad p(x)=x_0^2-\tfrac13,\qquad
f(x)=(-6x_{d-1}+2x_0)e_0.
$$

Using the deformed coordinate in the shear profile gives a nonzero velocity
geometry defect in both dimensions. The reference domain has unit volume
and zero analytic pressure mean. Independent mapped-cell integrals check
these quantities on each represented domain too; no cellwise identity
$\det D\Phi_h=1$ is assumed. Map regularity is checked independently by
the shared geometry study and remains a hypothesis of the rate estimates.

Let $\ell_h=\Phi_h\circ\Phi^{-1}$ and $w_h^\ell=w_h\circ\ell_h$.
For either field $w\in\lbrace u,p\rbrace $, the measured defects on $\Omega$ are

$$
F_w=w\circ\ell_h-w_h^\ell,\qquad
G_w=w-w\circ\ell_h,\qquad T_w=w-w_h^\ell=F_w+G_w.
$$

Their L2 norms and H1 seminorms use the exact-domain Jacobian and the
chain rule $D(w_h^\ell)=(Dw_h\circ\ell_h)D\ell_h$.
Both triangle and reverse-triangle inequalities are checked numerically.
Velocity represented-domain and lifted field errors have orders three/two;
velocity geometry and total errors have orders three/two for this quadratic
geometry. Pressure represented-domain, lifted field and total errors have
orders two/one. Since $\ell_h$ preserves $x_0$, $G_p=0$ analytically;
its two measured norms must be below $10^{-10}$, rather than being assigned
a rate or used in a relative comparison.

| Cell family | Grid points per axis in the sine-map rate study |
| --- | --- |
| Triangle | $9\to17\to33$ |
| Quadrilateral | $3\to5\to9$ |
| Tetrahedron | $9\to11\to13$ |
| Wedge | $5\to7\to9$ |
| Pyramid, hexahedron | $3\to4\to5$ |

Every adjacent interval uses the actual spacing ratio. The finer triangle,
tetrahedron and wedge levels avoid coarse pressure transients; the same
fixed acceptance windows stated above apply without relaxation. These
finite-hierarchy observations do not certify a uniform inf-sup constant.

For each lifted velocity defect $E\in\lbrace F_u,G_u,T_u\rbrace $, independently
integrated divergence obeys

$$
\Vert \operatorname{tr}DE\Vert _{L^2(\Omega)}\le\sqrt d\,|E|_{H^1(\Omega)}.
$$

A lift of a physically divergence-free field on $\Omega_h$ need not remain
divergence-free on $\Omega$. The affine patch checks zero represented and
lifted field errors, but explicitly retains a nonzero velocity geometry
derivative and divergence defect. The wrong-viscosity control uses pressure
as its rejection oracle: for a shear flow, pressure can absorb an incorrect
viscosity without necessarily increasing velocity error.

Assembly order $12\to16$ and norm order $14\to18$ are varied separately
on the fixed $n=5$ mesh. Each nonzero field norm must change by less than
$10^{-6}$ relatively. Existing direct-solver residual and pressure-gauge
checks remain in force. Native and PETSc solvers expose const solution
views only during a solve-scoped observer callback; the common lifted-norm
helper evaluates owned quadrature points and reduces MPI contributions.
It neither duplicates the mixed solve nor retains solution references.
