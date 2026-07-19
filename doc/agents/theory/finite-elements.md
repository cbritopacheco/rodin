# Finite elements — theory ↔ code

The element-level mathematics behind `Variational/FiniteElement.h`, the
FES families, and the quadrature module.

## The Ciarlet triple

A finite element is a triple $(K, P, \Sigma)$: a reference geometry $K$, a
polynomial space $P$ on it, and a set of degrees of freedom $\Sigma$
(linear functionals: point values, moments, ...) that is *unisolvent* —
$\dim P = |\Sigma|$ and the interpolation problem is uniquely solvable.
The dual basis $\{\phi_i\} \subset P$ with $\sigma_j(\phi_i) =
\delta_{ij}$ is the shape-function basis.

In code: `Polytope::Type` is $K$; an element class (e.g.
`H1Element<K, Range>`, `P1Element`, `P0Element`, `P0gElement`) bundles $P$
and $\Sigma$, exposing `getCount()` ($=\dim P$), `getBasis(i)` ($\phi_i$
and its derivatives on reference coordinates), and `getNode(i)` (the point
functional's location, when $\Sigma$ is nodal). Shape-function evaluation
returns a `TensorBasis` — the whole $\{\phi_i(rc)\}$ at once — because
assembly always needs all of them.

## Reference → physical: the pullback

Physical basis functions are defined by composition with the geometric
map $x_\tau : K \to \tau$ (`PolytopeTransformation`):
$\phi_i^\tau = \hat\phi_i \circ x_\tau^{-1}$. Derivatives pull back through
the Jacobian: $\nabla_x \phi = J^{-T} \nabla_{\hat x}\hat\phi$ with
$J = \partial x_\tau/\partial \hat x$ = `Geometry::Point::getJacobian()`.
Integrals change variables with the metric factor:
volume terms pick up $|\det J|$; facet terms pick up the Gram factor
$\sqrt{\det(J^T J)}$ (length/area of the pushed-forward reference facet).

This is why the sanctioned pattern is *evaluate geometry through
`Geometry::Point`* — it is the chart. Reconstructing $J$ from vertex
coordinates silently assumes an affine chart and breaks the moment a
`ParametricTransformation` (curved cell) appears.

**Isoparametric principle:** representing the geometry itself in the same
element space as the fields ($x_h \in [V_h]^d$) is what "curved P2 mesh"
means. Then $A = \nabla_{\hat x} x_h + \nabla_{\hat x} u_h$ for a deformed
configuration — geometry Jacobian plus field Jacobian, each obtained from
its own object, never mixed by hand.

## DOF layout over the incidence complex

Global conformity is achieved by *sharing* DOFs on shared sub-polytopes.
For degree-$K$ H1 on simplices: 1 DOF per vertex, $K-1$ per edge,
$\binom{K-1}{2}$ per triangular face, interior bubbles per cell. This is
literally how `H1<K>` builds its local→global map over
`Geometry::Connectivity` — which is why constructing a space *requires*
the incidence chain to be computed, and why FES code never invents DOF
orderings: `getGlobalIndex({d, i}, local)` is the only contract.

Vector-valued spaces are $[V_h]^m$ with `vdim = m`; the local ordering of
(node, component) pairs is an implementation detail that has changed
before — discover it through the element (`getBasis(L)` component
structure) if you ever need it, never hard-code `node*vdim + c`.

## Nodal vs modal bases (the K ≥ 2 trap)

Two ways to realize $P$:

- **Nodal (Lagrange)** at good point sets. Equispaced nodes are
  ill-conditioned for growing $K$ (Runge phenomenon / exploding Lebesgue
  constant); good sets are **Fekete points** (maximize the Vandermonde
  determinant) and **Gauss–Lobatto–Legendre (GLL)**; `WarpBlend.h` is the
  Warburton warp-&-blend construction of near-optimal simplex nodes.
- **Modal (hierarchical/orthogonal)**: the **Dubiner basis** — orthogonal
  polynomials on the simplex built from Jacobi polynomials through a
  collapsed-coordinate map (`Dubiner.h`, `JacobiPolynomial.h`,
  `LegendrePolynomial.h`). Well-conditioned mass/stiffness matrices at
  high order.

Rodin's `H1<K>` machinery (`H1/LagrangeBasis.h`, `Fekete.h`, `GLL.h`)
combines these: the basis is *constructed through* the modal/orthogonal
apparatus at Fekete-type nodes. Consequence that has caused real bugs:
**at $K\ge2$ a coefficient vector is not a list of nodal values.**
$u_h(x_j) = \sum_i c_i \phi_i(x_j)$ requires evaluating the basis. Always
go through `GridFunction::getValue` / the element basis; "coefficient =
value at node i" is only true for P1/P0.

## Projection and interpolation

`GridFunction::operator=(expr)` realizes the interpolant
$\Pi_h f = \sum_i \sigma_i(f)\,\phi_i$ (apply the DOF functionals to the
expression). This is *interpolation*, not $L^2$ projection — they differ
at $K\ge1$ for non-polynomial data; if you need the $L^2$-orthogonal
projection, write the mass-matrix problem
`Integral(u, v) - Integral(f, v)` and solve it.

## Quadrature

Element integrals are approximated by rules on $K$:
$\int_K f \approx \sum_q w_q f(\hat x_q)$, chosen by (geometry, order):

- `QF::GaussLegendre` — exact to order $2n-1$ on tensor domains;
- `QF::GaussLobatto` — includes endpoints (useful for lumping/collocation);
- `QF::GrundmannMoller` — arbitrary-order simplex rules;
- `QF::Centroid` — the 1-point rule (exact for affine × constant).

Correctness rule: the rule must integrate the *integrand as mapped* —
polynomial degree of trial × test × coefficient **plus** the degree of
$|\det J|$ (non-constant on curved cells). Under-integration shows up as
lost convergence order (manufactured tests) or spurious zero-energy modes;
Rodin defaults are chosen per term (`Variational/QuadratureRule.h` and the
per-FES fast paths), and order-4 is the house sampling for curved-validity
checks. When adding a term, match the existing choice for like terms and
let a manufactured test confirm the rate.

Two structural quadrature facts:

- Rules live on *reference* polytopes; `PolytopeQuadratureFormula`
  dispatches per `Polytope::Type`; `Geometry::PolytopeQuadrature` caches
  per-polytope data.
- The per-cell hot path is `bind(polytope)` (setPolytope) then
  `integrate()`; anything that depends only on (cell, quadrature) —
  weights × det, reference-basis values/Jacobians — belongs in bind. This
  is the caching pattern with the measured ~7× precedent
  (conventions.md; solvers-assembly.md).

## Element families in one line each

- `P0` — $(K, \mathbb P_0, \{\text{cell mean}\})$: $L^2$, cellwise data.
- `P0g` — one global constant DOF: the Lagrange-multiplier space.
- `P1` — $(K, \mathbb P_1, \{\text{vertex values}\})$: the $H^1$ workhorse;
  nodal in the honest sense.
- `H1<K>` — degree-$K$ continuous Lagrange with Fekete/Dubiner internals;
  nodal *interface*, non-nodal coefficients.
