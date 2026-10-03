# Curved Taylor–Hood Stokes verification

This suite tests the velocity/pressure pair $P_2/P_1$, together with a
global $P_{0g}$ pressure-mean multiplier, on exact quadratic geometry.
The physical domain is fixed across all refinement levels:

$$
\Omega=\Phi((0,1)^d),\qquad
\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1},\qquad d\in\{2,3\}.
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
\|\operatorname{div}u_h\|_{L^2(\Omega)}
\le\sqrt d\,|u-u_h|_{H^1(\Omega)}.
$$

This follows from $\operatorname{div}u=0$ and the Frobenius trace bound.
It ties the divergence error to the decreasing velocity derivative error.
A strict divergence rate is not imposed: cancellation can make this
quantity vanish or reach roundoff before the field errors do.
An MPI interpolation control with $u_h=x_0e_0$ independently checks
$\|\operatorname{div}u_h\|_{L^2(\Omega)}=1$, counting owned cells once.

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
$\|Az-b\|_2/\max(1,\|b\|_2)<10^{-11}$, where the coefficient vector
collects velocity, pressure and mean-multiplier unknowns. The computed pressure integral
has absolute magnitude below $10^{-10}$, in this dimensionless setting.
Independent analytic-volume/gauge checks use an absolute $10^{-12}$ bound.

Tests are labelled `convergence;slow`, with `petsc` and `distributed`
where applicable. Geometry-specific registrations have 30-minute safety
timeouts; pyramid registrations share a resource lock. Performance benchmarks
remain separate from these convergence assertions.
