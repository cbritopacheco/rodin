# Core vocabulary and infrastructure

The bottom layer: everything else is written in these terms.

## Types.h — the house vocabulary

Always use these aliases, never the raw `std::`/`boost::`/Eigen spellings:

- Scalars: `Real` (= Double), `Float`, `Integer`, `Boolean`, `Complex`,
  `Index` (= size_t).
- Containers: `List`, `Deque`, `Stack`, `FlatSet`, `FlatMap`, `Map`,
  `UnorderedSet`, `UnorderedMap` (boost::container / boost::unordered
  under the hood), `IndexSet` (= FlatSet<Index>), `IndexVector`,
  `IndexMap`, `Optional`, `StringView`.
- Flat containers are the default for small keyed data (cache-friendly,
  ordered); unordered for large hot maps.

## Object identity and lifetime interfaces

Three tiny root interfaces compose into `FormLanguage::Base`:

- `Copyable` — `virtual Copyable* copy() const noexcept = 0`; polymorphic
  deep clone, caller owns the result. This is how expression trees own
  their children.
- `Moveable` — polymorphic move counterpart.
- `Identifiable` — per-instance UUID, used for tracking form-language
  objects (e.g. matching a `TrialFunction` across a `Problem`).

`Cast.h` defines the `Cast` template for explicit conversion operations.

## Alert — errors and warnings

Streamable, composable diagnostics:

```cpp
Alert::MemberFunctionException(*this, __func__)
    << "explanation " << someValue << Alert::Raise;   // throws
```

- `Exception` (throws on `Raise`), `Warning`, `Info`; decorators for
  terminal output (Bold, Color, NewLine, Prefix, Identifier, Notation).
- Specialized bases carry context automatically: `ClassException`,
  `MemberFunctionException`, `NamespacedException`;
  `Variational/Exceptions/` adds domain ones
  (TrialFunctionMismatch, TestFunctionMismatch, UndeterminedTraceDomain).
- Use Alert for *user-facing/API-misuse* errors; use `assert` for internal
  invariants (compiled out under `-DNDEBUG`);
  `assert(ierr == PETSC_SUCCESS)` after PETSc calls (see conventions.md).
- `RODIN_SILENCE_WARNINGS` / `RODIN_SILENCE_EXCEPTIONS` configure macros
  can mute output.

## Threads

`Threads/Mutex.h`, `Threads/Unsafe.h`, `Threads/Mutable.h` — mutex wrapper,
"this is deliberately not thread-safe" marker, and guarded mutable state
(e.g. lazily computed caches inside const objects). `RODIN_THREAD_SAFE` is
a configure option. FormLanguage objects are documented NOT thread-safe;
parallel assembly gives each thread its own iteration state instead of
locking nodes.

## Context — where computation runs

`Context::Local` (single machine) and `Context::MPI` (distributed) are tag
types. Meshes, connectivity, and assembly are class templates over the
context (`Geometry::Mesh<Context::Local>` is the default `Mesh`).
Distributed behavior lives in the mirrored `MPI/` module tree
(MPI/Geometry, MPI/Assembly, MPI/Variational, MPI/IO), not in `#ifdef`s
inside core classes.

## Utility — metaprogramming toolkit

Small single-purpose template utilities used across the library:
`ForConstexpr` (compile-time loops over index sequences), `ParameterPack`,
`IntegerSequence`, `Product`/`Zip`/`Repeat` (type-level combinatorics),
`IsSpecialization`, `IsOneOf`, `IsCompleteType`, `HasTypeMember`/
`HasValueMember` (detection idiom), `Overloaded` (visitor lambdas),
`UnwrapReference`, `OptionalReference`, `DependentValue`, `Make`, `Extract`,
`STLCasts`. Check here before writing a new trait — it probably exists.

## Math

- `Math::Vector<Scalar>` / `Math::Matrix<Scalar>` — dynamic Eigen dense
  types; `Math::SparseMatrix<Scalar>` — Eigen sparse.
- `Math::SpatialVector<Scalar>` / `Math::SpatialMatrix<Scalar>` — fixed
  capacity 3 (stack, no heap) with *runtime* rows/cols. The types for
  dimension-generic (1/2/3D) geometry code. Caveats: no no-arg
  `Identity()`, no `Zero()`, no `M/scalar`, no Eigen comma-init — size
  explicitly (`Identity(2,2)`), write `(1/s)*M`.
- `Math::LinearSystem<Matrix, Vector>` — the Ax=b abstraction the solvers
  and Problem operate on; supports DOF elimination (essential BCs).
- `Rad`/`Deg` — CRTP unit types (`Math::Unit`) with conversions; angles in
  APIs are typed, not raw Real.
- `NewtonRaphson`, `RungeKutta` — small scalar/ODE steppers (used by e.g.
  the Heart 0D model), distinct from `Solver::NewtonSolver`.
- `Math::Traits` — `IsEigenObject` etc., used by the form language to
  absorb Eigen expressions.

## Configuration (Configure.h.in)

Feature macros stamped by CMake: `RODIN_USE_MPI`, `RODIN_USE_OPENMP`,
`RODIN_USE_UMFPACK`, `RODIN_USE_SPQR`, `RODIN_USE_CHOLMOD`,
`RODIN_WITH_PY`, `RODIN_WITH_PLOT`, `RODIN_THREAD_SAFE`,
`RODIN_SILENCE_WARNINGS`, `RODIN_SILENCE_EXCEPTIONS`, plus version macros
and `RODIN_RESOURCES_DIR`/`RODIN_THIRD_PARTY_DIR` paths and
`RODIN_MAXIMUM_POLYTOPE_VERTICES` (= 8, sizes `Polytope::Key`). Gate
optional-dependency code on these, never on ad-hoc macros.
