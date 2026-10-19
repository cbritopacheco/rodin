# Variational — the form language in depth

The largest module (~127 headers). Read philosophy.md first for the
underlying algebra (typed expression trees, grading); this doc is the
concrete inventory and the contracts. For extension work, keep
[workflows.md](workflows.md) and [numerical-contracts.md](numerical-contracts.md)
open beside this file: the first gives the checklist, the second states the
cross-module numerical rules this module must preserve.

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

Named-form numeric lifting must reject `FormLanguage::Base` arguments before
instantiating `FunctionDerived<RealFunction<T>>`. A trailing deduction-guide
return type needs an immediate-substitution guard around `FunctionDerived`
itself: Clang 19 can instantiate that return type before checking a false
head constraint. Numeric constructors also need the base-type exclusion as
the first constraint, before querying the lifted type. Complex function
coefficients must select the function path without being offered to
`RealFunction`. The named-form complex regression checks function/function,
both mixed function/value directions, and numeric/value construction against
explicitly matched quadrature operators; this is not a continuum-exactness
or refinement test.

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

Matrix-valued P0, P0g, P1 and H1 leaves own their expanded basis values
$\Phi_{a,j}=\phi_a E_j$. Their cached layout depends on both matrix dimensions
and the local basis count, in addition to sampling metadata. Reassigning a
borrowed space can change $2\times3$ into $3\times2$ while preserving its address,
the total component count and the bound sample. These are different matrix
layouts: a formula/sample key or component-count check alone is insufficient.
The owned basis dimensions and size are compared with the current finite
element before a cached value is reused. This is local logical invalidation;
it does not change entity ownership or require communication. Coefficient-vector
resizing after a space definition changes remains a separate caller obligation.

For complex-valued spaces, Trial ∘ Test expressions follow Rodin's
trial-first, test-second sesquilinear convention. See
[numerical-contracts.md](numerical-contracts.md#real-and-complex-form-convention)
for the precise algebra and backend rules.

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

Scalar P0/P0g shape gradients represent the broken derivative
$\nabla_h u_h=0$, not distributional jumps. Their expression contract still
requires a forwarded trial/test leaf and a polynomial order bound; zero is
sufficient for the zero expression. Named diffusion tabulates reference
gradients through `getGradient()` when available, otherwise through the
basis directional derivatives. This accepts constant bases without requiring
a convenience wrapper or changing the existing P1/H1 tabulation path.

Cell binding is an explicit evaluation boundary. Generic local quadrature
rules propagate `setPolytope()` through their expression trees before replacing
owned quadrature. Every shape leaf implements the binding operation explicitly;
the CRTP base forwards directly, including through custom expression wrappers.
Composite nodes forward the operation and coefficient caches
invalidate their point snapshots while retaining owned allocation for reuse. Native P1 Grad/Div/Jacobian reset prepared physical
results on every cell bind, retain allocated buffers, and reuse affine simplex
results across samples only while the cell binding matches. Direct point binds
without a matching cell binding always refresh. Copies start unbound. Reference
basis tables may survive bindings because they depend on the reference formula,
not physical geometry. Custom cell integrators establish the same lifecycle.

GridFunction retains DOF/scratch allocation reuse but recomputes pushed-forward
basis values on each direct integration-point evaluation. Caller-owned
FunctionBase::Cache retains the resulting field value across basis indices for
the current point binding. Shared fields acquire no mutable cell binding.
DeformationMap owns one explicitly prepared moved point: point binding locates
once, parameterless access reuses it, cell binding clears it, and direct point
queries always locate afresh. No transformation token, UUID or revision is used.
Binding setters, invalidation, and direct evaluation mutate this owned working
object and are non-const. Const coefficient evaluation owns a mutable evaluator;
the prepared moved-point accessor remains const.

`FiniteElement.h` defines the element side: `getBasis(local)`,
`getNode(local)` (reference nodes), counts per `Polytope::Type`.

`ReferenceTabulation<Scalar>` owns private dimensions and arrays without
cache state. `resize()` establishes all extents together; writable sample views
fill values and first derivatives without exposing resizable storage. Derivative
orders above one are excluded by the type contract.
`ReferenceTabulationCache<Element, Scalar>` owns the bounded thread-local store, logical formula/geometry keys, eviction and publication.
H1 supplies public `tabulate` overloads for caller-owned values, gradients, or
both, plus a whole-formula overload that sizes and fills a numerical table.
The cache delegates to that public operation without friendship. Elements
without bulk tabulation retain the generic direct-basis fill. Failed fills
leave the selected entry unpublished; neither elements nor consumers own a
second cache.

Scalar H1 uses this operation to share pyramid and tetrahedral modal values
and derivatives within each reference sample. Every nodal contraction keeps
its original mode order and accumulator type. Pyramid values use scalar
accumulators; tetrahedral values and all derivatives use real accumulators
before scalar promotion. Complex derivatives delegate to the real public basis
and must not accumulate in complex arithmetic instead. Scratch is point-local,
and empty sample collections require no Vandermonde initialization. Other geometries
retain direct basis evaluation.
Public pointwise basis calls, mapped quantities, quadrature orders, interpolation
and projection contracts are unchanged. Warm indexed consumers already reuse
the completed table; faster rebuilding is not a universal assembly speedup.

Public scalar H1 `tabulate` shares modal derivatives within a pointwise
basis sweep on tetrahedra and pyramids. Modal basis classes own modal evaluation;
H1 owns the inverse Vandermonde contraction. Other geometries retain public
basis directional derivatives. For immediate field consumption,
`element.tabulate<Geometry>(point)` returns an owned `Evaluation<Geometry>`.
It prepares modal derivatives and obtains the cached inverse once, then
`getGradient(local)` performs one nodal contraction. Scalar field gradients
consume that result immediately, preserving the fused coefficient accumulation
without an intermediate array of all nodal gradients. Vector Jacobians retain
the existing caller-owned dense gradient output, whose measured performance
is preferable for the small tetrahedral cases. Every
reference derivative accumulates in real arithmetic before scalar promotion. Each
inverse weight is loaded once per mode and reused for the three coordinates;
each coordinate retains the original mode summation order.
Evaluation copies own their modal rows; input coordinates and other evaluations
can change without invalidating them. The borrowed inverse has static lifetime.
The caller-owned dense value/gradient and whole-formula overloads remain available.
Their output extents and coordinates are checked before writing. Modal rows retain
the fixed four-entry layout, with an unused initialized value column for gradient
queries, so no dynamic allocation is introduced by a point evaluation.
The field-gradient consumer owns the cell's logical `(dimension, index)` pair
for that invocation and reads current coefficients through the space's public
`getGlobalIndex` operation. It does not replace that operation with a borrowed
`getDOFs` array: distributed spaces can return translated thread-local closure
storage with a different lifetime and construction cost. Associated integration
points retain the existing reference-table path. A pointwise-gradient improvement
must therefore be attributed separately from warm assembly that already uses
reference association. Neither operation introduces communication.

Matrix-valued H1 shape leaves may borrow an already prepared scalar reference
table while filling their owned component basis. Reference association is
required; formula provenance alone does not authorize indexed substitution.
The existing `findTabulation` operation is used without creating a new cache.
Missing or evicted entries retain direct scalar basis evaluation at the actual
point, rather than constructing values and derivatives over the whole formula.
The borrow ends within that binding call, before any subsequent insertion can
evict it. This policy can benefit a reference table already prepared by another
owner; it does not promise a table or a speedup for a standalone pointwise sweep.

## GridFunction

`GridFunction(fes)` = coefficient vector bound to a space.
- Data access: `getData()`/`setData()` (Math::Vector), `getValue(Point)` /
  `getValue(IntegrationPoint)` for FES-agnostic evaluation (this is the
  sanctioned way — works at any order/basis).
- `operator=` from any FunctionBase expression interpolates into the space by
  applying DOF functionals; `projectOnBoundary` variants exist. This is not an
  `L2` projection; when the distinction matters, assemble a mass-matrix
  `Problem` as described in
  [numerical-contracts.md](numerical-contracts.md).

`IntegrationPoint(point, formula, index)` records the originating integration
sample, not proof that the supplied point has that formula's reference
coordinates. Construct from `IntegrationPoint(polytopeQuadrature, index)` to
establish reference association with an owned mapped sample. Reference tables
and coordinate-dependent reuse use `getReferenceQuadratureFormula()`; it
returns null after formula assignment or for an unproven binding. The ordinary
formula accessor remains available for sampling provenance and coefficient
snapshots. Remapping into another cell must not promote origin metadata into
reference association. Grid-function derivatives on facets retain the existing
direct-point contract: evaluate in an adjacent full-dimensional cell, not with
a tangential facet-element table. These operations remain local and
non-collective. Geometry changes still require mapped-point reconstruction;
logical transformation identities invalidate physical basis/located-point
reuse, but do not update old geometric buffers.
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

Unweighted H1 mass and vector-Jacobian binding accept an empty vector component
range. Their local operator then has shape $0\times0$: mapped quadrature retains
its usual lifetime, but scalar reference tables and scalar contractions are not
needed. Equal empty DOF counts do not establish equal scalar basis counts for
mixed polynomial degrees. A triangular shortcut must never infer such equality
before handling the empty range. Matrix-valued H1 constructors separately reject
zero rows or columns; P0/P0g vector constructors require positive component
counts. These are distinct space contracts, not a backend-dependent policy.

For nonempty unweighted H1 Jacobian forms, component separation gives
$A_{bm+c,am+e}=\delta_{ce}S_{ba}$, where $S$ is the scalar trial/test stiffness
block and $m$ is the scalar component count of the value range. Accumulate the
first component block in the original quadrature order, then replicate its
entries into the remaining component-diagonal blocks. Cross-component zeros
remain initialized; no second scalar matrix or cache is needed. Flattened
matrix components obey the same rule. General component-coupling coefficients
do not establish this factorization and retain their separate weighted kernels.

P1 mass and directional-derivative binding likewise separate formula selection,
mapped-sample ownership and local matrix sizing from scalar reference preparation.
An empty component range returns before division, scalar preparation or coefficient
evaluation. A changed geometry/order retires the previous ready flag; the new
scalar collection is published only after filling succeeds. Convection obtains
trial derivatives and test values directly from scalar P1 elements in their
respective scalar fields, rather than decoding component-zero vector bases.
The stored reference data are independent of the current component count and
remain separate from the physical gradients recomputed on each bind.

H1 directional binding contracts each physical trial gradient with the direction
once per sample, before the test-basis loop. The scalar working buffer is
thread-local and grow-only, not a reference-data cache. A coefficient callback
may bind a different integrator of the same type: no buffer-data pointer or
physical-gradient dimension metadata may survive that callback. Afterwards,
every active scalar is replaced before use, physical gradients have invocation-local
ambient dimensions, and scalar table references are acquired anew. The regression
uses exact matrix equality for nested segment/hexahedron bindings in both ambient
dimension directions. This permits different integrator owners; it does not make
recursive mutation of the same owner a supported operation.

Weighted H1 vector-Jacobian binding has the same coefficient-callback boundary.
Trial/test gradient arrays are grow-only thread-local working storage, not a
physical-derivative cache. A nested same-type integrator can change their active
ambient dimensions even when no allocation occurs. The outer binding restores
every active gradient's logical dimension after each coefficient callback and
obtains reference tables afterwards, before filling and contracting gradients.
Equal degrees share one array; mixed degrees retain separate arrays. No pointer
to an array element or table entry is retained across the callback. The regression
compares nested and non-nested operators exactly in both segment/hexahedron
dimension directions, with scalar and nonsymmetric matrix coefficients. Ordinary
quadrature order, component coupling and backend communication are unchanged.

For a scalar weight $\nu$ and $m$ field components, the weighted operator separates
as $A_{bm+e,am+c}=\delta_{ec}S_{ba}$. Its first component block accumulates $S$ in
the original quadrature order. At the final logical sample, completed entries
are copied into the other component-diagonal blocks before ordinary transpose
mirroring. A zero final contribution does not imply a zero accumulated entry:
earlier contributions must still be copied. The existing exact-zero shortcut is
therefore retained only before finalization or at index zero, when no earlier
contributions exist. No norm threshold or coordinate matching is used. Empty
component ranges require neither coefficient evaluation nor scalar tabulation.
Matrix-valued coefficients remain in their fully coupled contraction path.

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

Native single-field and mixed problems own the body and linear system while
borrowing trial/test functions. Copy/move assignment replaces the body, system,
bindings and assembly state but retains the destination problem identity.
Complete-problem assembler assignment similarly retains destination identity;
the base has no execution state, and concrete assemblers transfer their own
policy, including the OpenMP thread setting. This does not enable assignment
of `FormLanguage::Base` or arbitrary expression nodes. Mixed cumulative offsets
are valid from construction, derived from the ordered space sizes, and refreshed
by assembly. No unused per-block form tuple is constructed. These native lifetime
contracts do not establish PETSc handle-copy or MPI collective semantics.

`ProblemBase` and its concrete single/mixed-field problems also support
`+=` and `-=` for terms supported by the existing `ProblemBody` arithmetic:

```cpp
Problem problem(u, v);
problem += Integral(Grad(u), Grad(v));
problem -= Integral(f, v);
problem += DirichletBC(u, g);
```

Each operation owns its expression through clone-on-build and invalidates the
previous assembly through virtual body assignment. `getBody()` is read-only.
Specifically, assignment sets `m_assembled = false`; the next `solve()` rebuilds
the system, whereas `getLinearSystem()` alone does not rebuild it.
Both signs support preassembled `BilinearForm` and `LinearForm` snapshots.
`problem += load` adds the linear form to the residual and therefore subtracts
its vector from the right-hand side; `problem -= load` adds it to the right-hand
side. Changing or destroying the original assembled form does not change its
owned snapshot in the problem.
Referenced coefficient changes still require explicit reassembly. Subtraction
adds a negative form contribution, not removal by UUID. Boundary-condition
subtraction and other operations absent from the body algebra remain unsupported.
Compound operators return a `ProblemBase&`, including through a concrete problem.
Non-variational adapters retain source compatibility and reject body access unless
they provide it.
