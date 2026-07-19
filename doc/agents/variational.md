# Variational — the form language in depth

The largest module (~127 headers). Read philosophy.md first for the
underlying algebra (typed expression trees, grading); this doc is the
concrete inventory and the contracts.

## FormLanguage substrate

- `FormLanguage::Base` — root of every node: `Copyable` (deep `copy()`),
  `Identifiable` (UUID), documented not thread-safe.
- `FormLanguage::Traits<T>` — open type-computation layer (ScalarType,
  RangeType, FESType, SpaceType...); scalar promotion via
  `FormLanguage::Mult<S1,S2>::Type` etc.
- `FormLanguage::List<T>` — polymorphic owning container of nodes (used by
  ProblemBody for integrator lists).

## Operator taxonomy (one header each)

- Arithmetic: `Sum`, `Minus`, `Mult`, `Division`, `UnaryMinus`, `Pow`,
  `Sqrt`, `Abs`, `Min`, `Max`.
- Linear algebra: `Dot` (vector dot / Frobenius A:B), `Trace`,
  `Transpose` (also on shape functions), `Frobenius` (norm),
  `IdentityMatrix`, `Component` (u(i), u(i,j), .x()/.y()/.z()).
- Calculus: `Grad`, `Div`, `Jacobian`, `Derivative` (directional). Each FES
  family provides its own specializations (see H1/, P1/, P0/, P0g/
  subdirectories) so `Grad(u)` dispatches to the right basis derivative.
- Transcendental: `Sine`, `Cosine`, `Tangent`, `Sinh`, `Cosh`, `Exp`.
- Comparison → BooleanFunction: `EQ`, `NEQ`, `LT`, `GT`, `LEQ`, `GEQ`;
  logic: `AND`, `OR`. Used for region predicates and conditional
  coefficients.
- Complex: `ComplexFunction`, `Re`, `Im`, `Conjugate`.
- Geometry probes: `F.h` (coordinate projections x,y,z as functions),
  `BoundaryNormal`, `FaceNormal`.
- DG/interface: `Jump`, `Average` (facet traces), used with
  `FaceIntegral`/`InterfaceIntegral`.
- Misc: `Zero`, `RelativeError`, `Restriction` (subdomain), `TraceOperator`
  (boundary restriction), `Flow` (flow maps), `Potential` (kernel
  potentials — the integral-equations entry point; examples in
  examples/IntegralEquations).

Function value categories: `RealFunction`, `ComplexFunction`,
`VectorFunction`, `MatrixFunction`, `BooleanFunction` — each a CRTP base
(`RealFunctionBase<Derived>` etc.) over `FunctionBase<Derived>`. Plain
values and lambdas lift implicitly (`RealFunction f = 1;`,
`VectorFunction v = {f, g};`).

## Shape functions and grading

`ShapeFunctionBase<Derived, FES, SpaceType>` with
`SpaceType ∈ {TrialSpace, TestSpace}` (compile-time enum).
`TrialFunction`/`TestFunction` are the leaves; `getLeaf()` recovers the
underlying unknown from any wrapped expression (needed to identify which
unknown an integrator contributes to). Combination rules (enforced by
overloads — ill-graded expressions do not compile):

- Function ∘ Function → Function
- Function ∘ ShapeFunction → ShapeFunction (same space)
- Trial ∘ Test (via Dot) → bilinear integrand

`getDOFs(polytope)` gives local size; evaluation produces a `TensorBasis`
(one value per basis function) rather than a single value.

## Finite element spaces

`FiniteElementSpace` = mesh association + element type + local↔global DOF
map (`getGlobalIndex({d, cellIdx}, local)`), global size, vector dimension.
Families:

- **P1** (`P1.h`, `P1/`) — piecewise linear Lagrange, scalar or vector
  (vdim from constructor). The workhorse. Note internal vector layout is
  interleaved but NEVER rely on it — use the local→global map.
- **P0** (`P0.h`, `P0/`) — piecewise constants (cellwise data,
  classifications, DG0 internal variables).
- **P0g** (`P0g/`) — *global constant* space: a single DOF shared by the
  whole mesh. This is the Lagrange-multiplier / global-constraint space
  (its integrators are the "global" integrator lists in Assembly inputs).
- **H1<K>** (`H1.h`, `H1/`) — degree-K continuous Lagrange, any K.
  Constructed `H1 vh(std::integral_constant<size_t, K>{}, mesh [, vdim])`.
  DOFs distribute over the incidence complex: vertices, then K−1 per edge,
  then face/cell interiors. **High-order bases are Fekete/Dubiner (with
  GLL/WarpBlend node machinery, Jacobi/Legendre polynomials under H1/) —
  NOT plain nodal Lagrange: coefficients ≠ nodal values at K ≥ 2.** Never
  reconstruct coordinates/values from coefficients; evaluate the
  GridFunction or the element basis.
- FES-specific operator specializations live in the family subdirectory
  (H1/Grad.h, P1/Div.h, ...) plus family QuadratureRule.h with optimized
  integration paths.

`FiniteElement.h` defines the element side: `getBasis(local)`,
`getNode(local)` (reference nodes), counts per `Polytope::Type`.

## GridFunction

`GridFunction(fes)` = coefficient vector bound to a space.
- Data access: `getData()`/`setData()` (Math::Vector), `getValue(Point)` /
  `getValue(IntegrationPoint)` for FES-agnostic evaluation (this is the
  sanctioned way — works at any order/basis).
- `operator=` from any FunctionBase expression **projects** onto the space
  (interpolation at DOF functionals); `projectOnBoundary` variants exist.
- Components `.x()/.y()/.z()` return component views.
- I/O: `load`/`save` in MFEM/MEDIT formats; XDMF/HDF5 via IO module.
- `GridFunctionBaseReference` / LazyEvaluator wraps a GridFunction for use
  inside expression trees without copying its data (expressions otherwise
  deep-copy operands).

## Integrals and integrators

`Integrator` hierarchy: `Type::Linear` (vector) vs `Type::Bilinear`
(matrix); under each, *local* (per-polytope, quadrature) and *global*
integrators. The user-facing constructors decide the class:

- `Integral(f, v)` / `Integral(v)` → LinearFormIntegrator over cells.
- `Integral(A(u), B(v))` (or a grade-2 integrand) → BilinearFormIntegrator.
- `BoundaryIntegral` — codim-1 exterior facets; `FaceIntegral` — all
  facets; `InterfaceIntegral` — interior facets (DG). All accept
  `.over(attr...)` restriction; functions accept `.traceOf(attr)` to pick
  the one-sided trace on interfaces.
- `QuadratureRule<...>` specializations are the default evaluation engines
  behind these (bind(polytope) → integrate()); FES-specific fast paths
  live in the family subdirectories.
- Assembly-facing contract: an integrator exposes the polytopes it ranges
  over (region: cells/boundary/faces) and fills a local matrix/vector for
  a bound polytope. Hot integrators should hoist per-cell work into the
  bind/setPolytope step (see conventions on caching).

## Boundary conditions

- `DirichletBC(u, g)` — value-prescribing essential BC (g any
  FunctionBase; scalar/vector). `.on(attr...)` selects facets; default =
  whole exterior boundary. Works on tagged *interior* facets too.
- `DirichletBC(u, A(v))` — **identification BC**: algebraically identifies
  u's boundary DOFs with a linear expression in another trial function v
  (component picks, rotations f·v, sums...). Both sides stay unknowns;
  assembly eliminates via the ConstraintMap. This is the mechanism for
  multi-field coupling constraints (e.g. mortar-like gluing) — do not
  emulate it with penalties.
- `PeriodicBC` — DOF identification across periodic pairs.

## Problem

`Problem(u1, ..., v1, ...)` takes the unknowns and test functions (multi-
field supported). Assignment of a **ProblemBody** — any `+`/`-` chain of
integrals and BCs — defines the system; `assemble()` builds the
`Math::LinearSystem` (`SparseProblem`/`DenseProblem` choose the matrix
type); `solve(solver)` assembles if needed, solves, and writes the
solution back into each `TrialFunction` (`u.getSolution()`).
`getLinearSystem()` exposes the assembled A, b for direct inspection.
Signs follow the math: `Integral(Grad(u), Grad(v)) - Integral(f, v)`
states a(u,v) − L(v) = 0.
