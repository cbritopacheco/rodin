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

For a requested order $p$, `QuadratureFormula` returns a rule that is
exact for every polynomial of total degree at most $p$ on the selected Rodin
reference polytope. Rules used by the default dispatcher have strictly
positive weights and points inside the reference polytope.

The order is a polynomial degree, not a number of points. A one-dimensional
Gaussian factor therefore uses

$$
n_p=\left\lceil\frac{p+1}{2}\right\rceil
$$

points. In integer arithmetic this is `p / 2 + 1`, avoiding overflow in the addition.

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
  without polymorphic copying.
- `XiaoGimbutas` and `WitherdenVincent` expose immutable tabulated rules.
- `GaussLegendre` constructs one-dimensional Gaussian rules and the weighted
  conical products.
- `GaussLobatto` includes element boundaries with positive weights, using
  Jacobi--Lobatto conical products on collapsed elements and merging coincident
  nodes. Its constructor takes a point count, with exactness `2*n-3`.
- `GaussLobatto::get(geometry, count)` caches canonical rules directly by
  geometry and point count. The generic dispatcher remains unchanged.
  Positive reference weights require a nonsingular physical mapping at the nodes
  to remain positive after mapping. P1 pyramid gradients use the symmetric
  centreline trace at the apex, preserving affine reproduction. General rational
  fields can have direction-dependent apex limits; this trace does not certify
  all directional quality limits.
- SWIFT uses its own cached `QualityCovering` for hinges and actual-quality
  checks: closed-form shifted/contracted reference lattices, with a homothetic
  tetrahedral construction. `sampling.subdivision` selects the generating
  lattice resolution; point counts are unchanged, but vertices are not retained.
  Equal positive base weights sum to reference volume; adaptive hinge weights
  preserve mapped discrete mass. This sampling policy makes no polynomial
  exactness or continuous quality-certification claim. It does not replace
  ordinary surface/volume integration rules.
- `TensorProduct` composes two or three formulas and multiplies their weights.
- `GrundmannMoller` remains an explicitly selectable simplex family. It is not
  a default because its weights are signed above degree one.
- `QuadratureFormula::get(order, geometry)` selects a concrete family out of
  line and obtains it from the persistent factory. Tensor products receive
  persistent component formulas. No dispatcher forwarding object or shared
  ownership is retained. Default and explicit-family requests with identical
  normalized constructor arguments return the same object; different requested
  orders may share a formula when they select the same Gaussian point counts.
- Default requests first use an eight-entry non-owning thread-local lookup of
  order and geometry, preserving the fast path before family selection. Each
  concrete-family/key-type pool also has an eight-entry thread-local cache.
  Misses enter its locked canonical FlatMap. A published one-point rule needs
  one storage allocation when map capacity is available. Other rules may also
  allocate sample arrays; tensor products can publish missing component rules.

Reference-formula reuse and mapped-point ownership are separate lifetimes.
Variational quadrature rules own a `Geometry::PolytopeQuadrature` for their
currently bound polytope; rebinding releases the previous mapped points.
This bounds mapped-point storage by the active rules and quadrature sizes,
not by the number of traversed cells. Explicit `mesh.getQuadrature(...)`
calls still return mesh-owned, stable borrowed objects until geometry is
flushed. Assembly does not evict these objects. Both lifetimes are
backend-independent and introduce no MPI collective.

Persistent APIs accept only `const QuadratureFormula&`, a final noncopyable,
nonmovable value-access implementation. `QuadratureFormula::get<Family>(args...)`
constructs the concrete family and its read-only interface together in one
pool allocation outside its mutex, without copying reference samples, and
publishes that record under a tuple key. Optional `getSamples()` views expose
contiguous points and weights with a scalar stride (zero for uniform weights).
Publication retains these views to avoid pointwise virtual dispatch; families
without views retain ordinary pointwise access. Grundmann–Möller, Xiao–Gimbutas and Witherden–Vincent tables are
retained throughout static destruction because published rules borrow them. Each family/key-type instantiation
has its own shared FlatMap pool retained throughout static destruction. Integer
and enum arguments are values; published-formula arguments use a private
`Reference` containing a typed `pointer` and total pointer ordering. There are
no per-family key specializations or identity methods. Argument value types
matter; public family factories normalize orders/counts to `size_t`.

The lookup capacity is `RODIN_QF_QUADRATUREFORMULA_LOOKUP_CAPACITY` (default eight,
positive). Thread-local entries borrow shared published objects and require no
lock. Shared map lookup/insertion require the pool mutex. Concurrent misses may
build duplicate candidates but return the same retained object. Failed generation
publishes nothing. Constructors must depend on their arguments and synchronize
any shared state of their own; cyclic construction is unsupported.

Mapped and basis caches use persistent formula addresses. A formula cannot be
assigned, moved, destroyed with its producing thread, or reconstructed at a cached
address. Local concrete generators implement `QuadratureFormulaBase` and remain
usable by immediate `element.tabulate(table, generator)` operations, which retain
no formula reference. IntegrationPoint construction from mapped quadrature
establishes reference-chart association; raw point/formula/index construction
records provenance without authorizing reference-table reuse.

Mapped-quadrature construction prepares geometric Jacobians through the
transformation's bulk `jacobian` overload. Parametric elements reuse their
existing reference tables, while generic transformations
retain the pointwise implementation. The resulting matrices belong to the
mapped points; no borrowed reference-table storage escapes the bulk call.
Temporary matrix storage is proportional to the current formula size, not
the number of mesh cells. Copying or moving a point preserves its owned
geometry, and rebinding it invalidates that geometry normally.

`ReferenceTabulationCache<Element, Scalar>` owns the common bounded thread-local
store used by scalar P1 and H1 elements. Each element type retains its own basis
formulas and ordering, including the rational pyramid basis. The stored arrays
contain only reference values and first derivatives; physical Jacobians are
recomputed from each cell's control coordinates. The eight-entry eviction policy
and persistent formula addresses apply equally to both element families. Warm
hits compare the formula address and geometry without querying size or nodes;
published immutability makes those queries redundant. Cold
misses size and fill numerical storage before publishing an entry. The unit regressions
count formula queries and verify that published samples remain stable across generator and producer-thread lifetimes, including real and complex H1 degrees one through six.

Lookup and construction share one already obtained thread-local cache reference.
The private formula search serves both public access operations; table filling
and publication remain separate stages of this same owner. Replacement retires
the selected entry's previous formula association before changing its storage and publishes
the new formula association only after filling completes. Once replacement begins, the
round-robin successor is selected before filling and is not rewound on failure.
An interrupted fill therefore exposes neither formula association, while other entries remain
unchanged. The regression checks these formula keys and exact surviving arrays.

Physical-coordinate evaluation reuses the basis values from the same reference
table without making coordinates eager. A mapped point owns its reference
coordinates and stores only a persistent formula pointer and quadrature index. The
parametric transformation borrows a matching table during evaluation and
accumulates

$$
x_{K,j}(\widehat x_q)=\sum_a X_{K,ja}\phi_a(\widehat x_q)
$$

in the same local-basis order as direct evaluation. The table is never retained
by the point. If the entry has been evicted or evaluation occurs on another
thread, the direct map is evaluated at the owned reference coordinates instead.
The formula pointer remains valid until process termination. A missing table
still selects direct evaluation, preserving lazy coordinates for integrands that
need Jacobians but not physical coordinates.

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
