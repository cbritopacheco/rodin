# Quadrature formulas

The `QF` module supplies quadrature formulas on Rodin reference polytopes.
Its runtime architecture has two levels: published positive-interior tables
cover the degrees for which compact rules are available, and positive
Gaussian product rules provide deterministic arbitrary-order fallbacks. The
runtime library contains no nonlinear rule generator.

Read [numerical-contracts.md](numerical-contracts.md) before changing default
orders, fast paths, or mapped-integrand assumptions; quadrature exactness is a
cross-module numerical contract, not only a `QF` implementation detail.

## Contract

For a requested order $p$, `PolytopeQuadratureFormula` returns a rule that is
exact for every polynomial of total degree at most $p$ on the selected Rodin
reference polytope. Rules used by the default dispatcher have strictly
positive weights and points inside the reference polytope.

The order is a polynomial degree, not a number of points. A one-dimensional
Gaussian factor therefore uses

$$
n_p=\left\lceil\frac{p+1}{2}\right\rceil
$$

points. In integer arithmetic this is `(p + 2) / 2`.

## Default dispatch

| Reference polytope | Tabulated range | Arbitrary-order rule |
|---|---:|---|
| Point | all orders | `Centroid` |
| Segment | all orders | Gauss--Legendre |
| Triangle | Xiao--Gimbutas, $1\le p\le50$ | positive conical product |
| Quadrilateral | all orders | Gauss--Legendre tensor product |
| Tetrahedron | Xiao--Gimbutas, $1\le p\le15$ | positive conical product |
| Hexahedron | all orders | Gauss--Legendre tensor product |
| Wedge | XG triangle $\times$ GL, $1\le p\le50$ | conical triangle $\times$ GL |
| Pyramid | Witherden--Vincent, $1\le p\le10$ | positive conical product |

`XiaoGimbutas` and `WitherdenVincent` remain directly constructible. Their
complete table ranges are:

| Family | Triangle | Quadrilateral | Tetrahedron | Hexahedron | Wedge | Pyramid |
|---|---:|---:|---:|---:|---:|---:|
| Xiao--Gimbutas | 1--50 | -- | 1--15 | -- | -- | -- |
| Witherden--Vincent | 1--20 | 1--20 | 1--10 | 1--10 | 1--10 | 1--10 |

The dispatcher uses WV only on pyramids. Tensor-product Gaussian rules are
more direct on quadrilaterals and hexahedra, while the wedge is composed from
the selected triangle rule and a segment rule.

## Conical-product fallbacks

The simplex fallbacks use the Duffy transformation and carry its Jacobian in
Gauss--Jacobi weights. For the triangle,

$$
(r,s)=\bigl(u,(1-u)v\bigr),
\qquad dr\,ds=(1-u)\,du\,dv,
$$

which gives `GJ(alpha=1) x GL`. For the tetrahedron,

$$
(r,s,t)=\bigl(u,(1-u)v,(1-u)(1-v)w\bigr),
$$

with Jacobian $(1-u)^2(1-v)$, which gives
`GJ(alpha=2) x GJ(alpha=1) x GL`. For the pyramid,

$$
(x,y,z)=\bigl((1-z)u,(1-z)v,z\bigr),
$$

with Jacobian $(1-z)^2$, which gives
`GL x GL x GJ(alpha=2)`. These constructions are positive, interior,
arbitrary-order rules. They are not fully symmetric and generally use more
points than optimized tables.

## Runtime structure

- `QuadratureFormulaBase` is the value-access interface: size, weight, point,
  and polymorphic copy.
- `XiaoGimbutas` and `WitherdenVincent` expose immutable tabulated rules.
- `GaussLegendre` constructs one-dimensional Gaussian rules and the weighted
  conical products.
- `TensorProduct` composes two or three formulas and multiplies their weights.
- `GrundmannMoller` remains an explicitly selectable simplex family. It is not
  a default because its weights are signed above degree one.
- `PolytopeQuadratureFormula::build` constructs an owned formula.
- `PolytopeQuadratureFormula::get` returns a canonical formula from a
  process-wide pool. An eight-entry thread-local cache serves repeated hot-path
  lookups without locking; cache misses enter the locked canonical pool.

The Xiao--Gimbutas coefficients are taken from the authors' `triasymq`
distribution. The Witherden--Vincent coefficients are taken from PyFR's
published quadrature tables. The transformed coefficients, exact source
revisions, and complete attribution licenses live in `XiaoGimbutasData.h` and
`WitherdenVincentData.h`.
Numerical search, fitting, orbit enumeration, and table-emission code do not
belong to the runtime module.

## Verification

The tests under `tests/unit/Rodin/QF/` enforce:

- every tabulated monomial moment through the claimed degree;
- exact dispatch through degree 20 on every reference polytope;
- the first untabulated triangle and tetrahedron fallback degrees;
- positive weights, reference measure, and interior points;
- agreement between independent published families;
- affine covariance, subdivision additivity, random-polynomial exactness,
  deterministic construction, and rejection of perturbed coefficients;
- published point counts and external published-rule samples.

Tests use closed-form moments on Rodin's reference polytopes. Runtime tables
must never be validated by regenerating the coefficients they contain.
