# Design philosophy and code patterns

Read this before writing any Rodin code. It is inferred from the code itself;
when in doubt, imitate the nearest existing header, not generic C++ habit.

## The prime directive: the program is the mathematics

User code must read like the weak formulation it implements:

```cpp
H1 vh(std::integral_constant<size_t, 2>{}, mesh);
TrialFunction u(vh);
TestFunction  v(vh);
Problem poisson(u, v);
poisson = Integral(Grad(u), Grad(v))
        - Integral(f, v)
        + DirichletBC(u, g).on(boundaryAttr);
CG(poisson).solve();
```

Every API decision serves that line. Names are the mathematical names
(`Dot`, `Grad`, `Div`, `Jump`, `Average`, `FaceNormal`, `Frobenius`), never
implementation names. If a proposed API cannot be written as the formula it
computes, the design is wrong, not the formula.

## The underlying algebraic concept

The library models **a typed, graded algebra of variational expressions over
a graded incidence complex**. Both halves matter.

### 1. The form language is a typed expression algebra

- Every expression is a node deriving from `FormLanguage::Base`
  (`src/Rodin/FormLanguage/Base.h`), which supplies exactly two services:
  identity (`Identifiable`, UUID) and polymorphic deep copy (`Copyable`,
  `virtual Base* copy() const noexcept`). Nothing else lives in the root.
- **The C++ type of an expression is its parse tree.** Operators are class
  templates parameterized by their operand types:
  `Dot<FunctionBase<L>, FunctionBase<R>>`,
  `Integral<Dot<ShapeFunctionBase<...>, ShapeFunctionBase<...>>>`. There is
  no runtime AST; the compiler holds the tree.
- **Value semantics via clone-on-build.** A node's constructor deep-copies
  its operands (`m_lhs(lhs.copy())`). Expressions therefore own their whole
  subtree; no dangling references to temporaries, no shared mutable state.
  Copy assignment is deleted at the root (`Base`) — nodes are built,
  copy/move-constructed, and cloned; they are not copy-assigned.
- **Two-level polymorphism, deliberately split.** Evaluation (the hot path)
  is CRTP-static: `FunctionBase<Derived>::operator()(p)` casts to
  `Derived::getValue(p)`. Virtual dispatch is reserved for ownership and
  cold paths (`copy()`, `getName()`). Do not introduce virtual calls into
  per-quadrature-point code.

### 2. The algebra is graded, and the grading lives in the type system

An integrand's *grade* is how many free (trial/test) slots it carries:

- grade 0 — `FunctionBase<D>`: known data (coefficients, solutions, `f`).
- grade 1 — `ShapeFunctionBase<D, FES, TestSpace>`: linear-form integrand.
- grade 2 — trial × test pair: bilinear-form integrand.

Products combine grades and the compiler enforces the arithmetic:
Function·Function → Function; Function·ShapeFunction → ShapeFunction of the
same space; Trial·Test → bilinear entry. `Integral(Grad(u), Grad(v))` is a
`BilinearFormIntegrator`; `Integral(f, v)` is a `LinearFormIntegrator`; an
ill-graded form fails to compile. The Trial/Test distinction is a
compile-time enum template parameter (`ShapeFunctionSpaceType`), not a
runtime flag.

Orthogonally, expressions are graded by **range shape** (Boolean / Real /
Complex / Vector / Matrix), propagated through
`FormLanguage::Traits<T>::RangeType` and scalar promotion computed at the
type level (`FormLanguage::Mult<LHSScalar, RHSScalar>::Type` handles e.g.
Real×Complex).

### 3. Traits are the type-level computation layer

`FormLanguage::Traits<T>` (open specialization, one per node family)
publishes associated types: `ScalarType`, `RangeType`, `FESType`,
`ElementType`, `SpaceType`. The primary template is intentionally
incomplete — a missing specialization is a compile error, never a wrong
default. When you add a node, you add its `Traits` specialization in the
same header, above the class.

### 4. Geometry is a graded incidence complex, separate from coordinates

- A mesh is a family of polytopes identified by `(dimension, index)`, with
  **incidence relations d → d′** stored in `Connectivity` and computed
  *explicitly, on demand*: `mesh.getConnectivity().compute(1, 2)`. Nothing
  computes incidence behind your back; code states its topological
  requirements up front (see the `RODIN_GEOMETRY_REQUIRE_*` macros in
  Mesh.h).
- `Polytope::Key` is the vertex multiset with symmetric hash/equality —
  topological identity is independent of vertex ordering, and keys are
  fixed-size stack arrays (no heap).
- **Topology and geometry are separate layers.** Coordinates and curvature
  attach to topology via per-polytope transformations (affine by default,
  `ParametricTransformation<FE>` for curved/isoparametric cells). Higher
  layers repeat the split: the level set owns topology, classifiers own
  discrete attributes, mesh optimizers own geometry — a layer never decides
  what belongs to another.
- Finite element spaces are defined by *distributing DOFs over the incidence
  complex* by dimension (vertices, then edge interiors, then face/cell
  interiors as degree grows). High-order H1 uses Fekete/Dubiner bases —
  never assume "coefficient = nodal value" or a `node*vdim+component` DOF
  layout; get geometry from `Geometry::Point(cell, rc)` and field values
  from the `GridFunction` itself.

## Minimal and sufficient

The codebase is built from small pieces, each doing one thing, composed
rather than extended:

- **One header = one algebraic concept.** `Dot.h`, `Grad.h`, `Sine.h`,
  `Jump.h`, `EQ.h`, ... The `Variational/` directory listing *is* the
  vocabulary of the language. A new operation is a new header with one class
  template, its `Traits`, its specializations, and its deduction guides.
- **Specialize one name; never multiply names.** `Dot` has specializations
  for every operand grading; `CG<Math::LinearSystem<...>>` specializes per
  system type; `H1<K, Range, Mesh>` per range. CTAD deduction guides let
  users write `Dot(a, b)`, `CG(problem)`, `Integral(Grad(u), Grad(v))` with
  no explicit template arguments. Prefer a new specialization of an existing
  concept over a new class name.
- **Vary behavior by template parameter, not by subclass or flag** (scalar
  type, range type, FES, mesh/context — every FES family and solver is
  parameterized this way).
- **Free functions for procedures; classes only where state is owned.**
  Multi-stage algorithms are staged free functions (e.g. `solveWNGIR`),
  each stage independently testable, composed by the caller.
- **New solvers are expressed IN the form language, not beside it.** WNGIR
  is built from `Problem` + `Integral`/`FaceIntegral` + custom integrators
  that slot into the existing assembly — a sentence in the language, not a
  parallel framework. If you need something new, add the minimal node (a
  term, an integrator), not a subsystem.
- **New capability defaults off and preserves behavior.** Optimizer knobs
  ship with neutral defaults (`setSmoothingPasses(0)`); performance work
  must be math-identical against a baseline; there is exactly one canonical
  path per operation — no duplicate "fast paths" (a bolted-on fast path once
  cost an entire PR; see doc/agents/conventions.md).

## No redundancy, no fragmentation

- One vocabulary: house aliases from `Rodin/Types.h` (`Real`, `Integer`,
  `Index`, `IndexVector`, `FlatSet`, `FlatMap`, `UnorderedMap`, `Optional`,
  `StringView`) — never raw `std::`/`boost::` spellings for these. Small
  dimension-generic objects use `Math::SpatialVector`/`SpatialMatrix`
  (fixed stack storage, runtime dims) so 2D/3D share one code path.
- One place per concern: PETSc calls only under `src/Rodin/PETSc/`;
  third-party integrations (MMG, MFEM, Scotch) additive in their own
  directories; the core never depends on them.
- Don't fake genericity: intrinsically-2D constructs stay 2D constructs
  (triangle-specific targets); don't template them into pretend generality.
- Don't split what belongs together: a node's Traits, class, specializations
  and deduction guides live in one header; tests mirror `src/Rodin/` one
  test executable per module component.

## Surface style (imitate exactly)

- Boost Software License header block at the top of every file; include
  guards `RODIN_<NAMESPACE>_<FILE>_H` (no `#pragma once`).
- Members `m_camelCase`; accessors `getX()`/`setX(...)`, setters return a
  reference for chaining (`solver.setTolerance(1e-10).setMaxIterations(...)`).
  Builder-style qualifiers read like math: `.on(attr)`, `.traceOf(attr)`.
- Leaf specializations are `final`. Constructors: default/copy/move spelled
  out (copy ctor deep-copies children via `copy()`); copy assignment
  deleted or omitted.
- Doxygen with real mathematics: `@f$ ... @f$` formulas, a "Mathematical
  Foundation"/usage section, `@defgroup <Name>Specializations` per
  specialization family. Documentation states the formula the code computes.
- Errors: `Rodin::Alert` exceptions (streamable, composable) at API
  boundaries for user mistakes; `assert` for internal invariants (remember
  it vanishes under `-DNDEBUG`); `assert(ierr == PETSC_SUCCESS)` after PETSc
  calls.
- Header-heavy: templates live in headers; `.cpp` files exist only for the
  rare non-template code. Each module has a `ForwardDecls.h`; include it
  instead of heavyweight headers where a declaration suffices.
