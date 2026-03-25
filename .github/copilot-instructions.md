# GitHub Copilot Instructions for Rodin

## What Rodin is

Rodin is a modular C++20 finite element framework. It is designed as a set of orthogonal subsystems that must remain composable:

- **Geometry** handles mesh entities, topology, transformations, and evaluation points.
- **QF** handles quadrature formulas.
- **Variational** provides the high-level form language and problem composition.
- **Assembly** handles low-level assembly logic.
- **Solver** provides solver abstractions.
- Optional modules (PETSc, MPI, etc.) extend core behavior without replacing it.

When adding code, preserve this separation. New functionality should fit into the existing architecture instead of bypassing it.

## How to think before writing code

Before implementing anything, ask:

1. Which existing Rodin module should own this responsibility?
2. Is this a local evaluation abstraction, a global assembly abstraction, or a user-facing composition abstraction?
3. Am I extending an existing Rodin pattern, or introducing a foreign mini-framework?
4. Can the user still build problems explicitly using Rodin's existing style?
5. Am I hard-coding a prototype assumption as if it were generic?

Do not optimize for the shortest patch if it conflicts with Rodin's structure.

## Core architectural rules

### 1) Reuse existing abstractions first

If Rodin already has a good abstraction, build on it instead of duplicating it.

Examples:
- Use `Geometry::Point` for geometric quadrature-point context.
- Use existing FE space / basis abstractions instead of hard-coding one element type.
- Use `Variational::Problem` composition style where appropriate.
- Prefer existing field-oriented abstractions over raw vector-only APIs unless unavoidable.

### 2) Keep local and global responsibilities separate

- Geometry / FE / quadrature decide where and how evaluation happens.
- Constitutive laws compute local material response.
- Integrators / assembly accumulate local contributions globally.
- Problem objects / user code compose those pieces.

Do not mix these layers.

### 3) Prefer composable building blocks over monolithic managers

Rodin's style is compositional: trial/test functions, integrals, boundary conditions, explicit `Problem` construction. New modules should preserve this style and avoid introducing "manager owns everything" patterns unless an existing module in the same problem class already does this.

### 4) Do not confuse a working prototype with a generic backend

A prototype that only works for one element, one quadrature rule, one constitutive law, or one geometric assumption must not be represented as a generic module. If the module is intended to be generic, remove prototype assumptions in code, not only in documentation.

## Distinguish core layers explicitly

When reasoning about finite element mechanics and PDE modules, keep these concerns distinct:

- **Geometry** (entities, mappings, points)
- **FE evaluation** (shape functions, gradients, quadrature evaluations)
- **Local constitutive data/laws** (pointwise state and response)
- **Global assembly** (integration and accumulation into global operators)

This separation is critical for maintainability and reusability.

## Preserve the "user builds the problem" philosophy

Prefer APIs that keep user-level problem composition explicit. Rodin users should still be able to assemble and solve by composing trial/test functions, integrals, boundary terms, and solvers directly. Avoid hiding problem definition behind opaque orchestration objects.

The canonical user-facing pattern looks like:

```cpp
Mesh mesh;
mesh = mesh.UniformGrid(Polytope::Type::Triangle, {16, 16});
mesh.getConnectivity().compute(1, 2);

P1 Vh(mesh);
TrialFunction u(Vh);
TestFunction  v(Vh);

Problem problem(u, v);
problem = Integral(Grad(u), Grad(v))
        - Integral(f, v)
        + DirichletBC(u, Zero());

Solver::SparseLU solver;
problem.solve(solver);
```

New modules must integrate naturally into this composition style.

## Module structure overview

```
src/Rodin/
├── Alert/              # Logging and error reporting (compiled)
├── Assembly/           # FE assembly engine (compiled)
├── Context/            # Execution contexts: Local, MPI (compiled)
├── FormLanguage/       # Expression template DSL base (compiled)
├── Geometry/           # Mesh, polytopes, points, connectivity (compiled)
├── IO/                 # I/O: MFEM, MEDIT, XDMF, HDF5 (compiled)
├── Math/               # Linear algebra, vectors, matrices (compiled)
├── QF/                 # Quadrature formulas (compiled)
├── Serialization/      # Serialization support (compiled)
├── Solver/             # Linear solvers interface (header-only)
├── Utility/            # General utilities (compiled)
├── Variational/        # FE spaces and weak formulations (compiled)
│   ├── P0/             # Piecewise constant
│   ├── P0g/            # P0 with gradient enrichment
│   ├── P1/             # Piecewise linear
│   └── H1/             # Arbitrary-order Lagrange
├── Advection/          # Lagrangian advection (header-only)
├── Eikonal/            # Fast Marching Method (header-only)
├── LinearElasticity/   # Linear elasticity integrators (header-only)
├── Solid/              # Hyperelastic solid mechanics (header-only)
├── MMG/                # Mesh adaptation via Mmg (compiled)
├── PETSc/              # PETSc solvers and distributed assembly (compiled)
├── MPI/                # Distributed mesh and parallel assembly (compiled)
└── Scotch/             # Mesh partitioning (conditional)
```

Physics/model modules (`Advection`, `Eikonal`, `LinearElasticity`, `Solid`) are header-only INTERFACE libraries. They follow the pattern: `add_library(RodinX INTERFACE)`, `add_library(Rodin::X ALIAS)`, link `Rodin::Geometry` + `Rodin::Variational`, and install via the `RodinTargets` export set.

## Key abstractions to know

### Geometry::Point

Represents a spatial point on a mesh polytope, bundling reference coordinates, physical coordinates, and Jacobian data. It is the standard quadrature-point context everywhere in Rodin.

Key API:
- `x()`, `y()`, `z()` — Cartesian physical coordinates
- `getPhysicalCoordinates()` / `getReferenceCoordinates()` — full coordinate vectors
- `getPolytope()` — the owning mesh element
- `getJacobian()` / `getJacobianDeterminant()` / `getJacobianInverse()` — geometric transformation

Build by composition over `Geometry::Point` when extending (e.g., `ConstitutivePoint` in the Solid module), not by duplicating coordinate storage.

### FormLanguage expression template system

All form-language objects (functions, operators, integrands) inherit from `FormLanguage::Base`, which provides:
- Polymorphic cloning via a virtual `copy()` method (returns heap-allocated copy)
- Identity via `Identifiable` (UUID-based identity tracking)
- Lifetime management for temporary objects via `object()` helpers

Operators like `Grad`, `Div`, `Dot` are CRTP-based expression templates. They store operands via `std::reference_wrapper` (non-owning) or `std::unique_ptr` (owning), enabling lazy evaluation during assembly.

When adding new operators, follow this pattern:
1. Define a `Traits` specialization in `FormLanguage/Traits.h` for associated types.
2. Create a CRTP base (e.g., `GradBase<Operand, Derived>`) inheriting from an appropriate `FunctionBase`.
3. Store operands with `std::reference_wrapper<const OperandType>`.
4. Add CTAD (Class Template Argument Deduction) guides if needed.

### FormLanguage::Traits

The `Traits<T>` primary template is unspecialized. Each form-language type provides a specialization defining its associated types (`ScalarType`, `RangeType`, `FESType`, etc.). Always specialize `Traits` for new types.

### ForwardDecls.h pattern

Each module has a `ForwardDecls.h` header containing forward declarations, type aliases, and enums. Include only `ForwardDecls.h` (not full headers) when a forward declaration suffices, to minimize compile-time coupling.

### FE space hierarchy

```
FiniteElementSpaceBase (abstract)
├── P0    — piecewise constant, 1 DOF per cell
├── P0g   — P0 + gradient enrichment
├── P1    — piecewise linear, 1 DOF per vertex
└── H1<k> — arbitrary-order Lagrange (k = 0..6+)
```

Each space defines: `getSize()` (total DOFs), `getVectorDimension()`, `getMesh()`, `getGlobalIndex(d, idx, local)`.

### Solver interface

`Solver::LinearSolverBase<LinearSystem>` is the abstract interface. Concrete solvers (CG, GMRES, SparseLU, UMFPACK, CHOLMOD, SPQR) and PETSc extensions (KSP, SNES) override `solve(LinearSystem&)`.

### Alert system (error handling)

Rodin uses a stream-based alert system, not raw exceptions:

```cpp
Alert::MemberFunctionException(*this, __func__)
  << "Incidence " << d << " -> " << dp << " has not been computed."
  << Alert::Raise;
```

Alert classes: `Info` (blue), `Warning` (yellow), `Success` (green), `Exception` (red, throws). Always use `Alert::Exception` (or a derived exception class) with the `<< Alert::Raise` pattern instead of `throw`.

## How to structure new modules

Organize by responsibility, not convenience. Typical layering for physics modules:

- Kinematics/state quantities
- Constitutive input/context data
- Invariant evaluators
- Constitutive laws
- FE integrators
- Derived fields / post-processed fields

Avoid law-specific APIs in generic core components and avoid grouping unrelated logic under one module.

## Expectations for generic code

If something is documented as generic, it should generally be generic across relevant axes:

- FE space / element type
- Quadrature rule
- Geometry / polytope type
- Backend where relevant
- Constitutive law where relevant

Do not hard-code narrow assumptions (single element type, centroid-only quadrature, single-point assumptions, fixed nodal ordering) in APIs intended to be reusable.

## Code style conventions

Rodin follows consistent code style throughout. Match these patterns exactly:

### Formatting

- **Brace style**: Allman (opening brace on its own line).
- **Indentation**: 2 spaces, no tabs.
- **Line length**: soft limit around 100–110 characters.

### Naming

- **Classes and types**: `PascalCase` — `GridFunction`, `BilinearForm`, `RealFunction`.
- **Methods and functions**: `camelCase` with `get`/`set` prefixes — `getSize()`, `getValue()`, `setName()`.
- **Member variables**: `m_` prefix — `m_u`, `m_fes`, `m_value`.
- **Static member variables**: `s_` prefix — `s_id`, `s_nodes`.
- **Template parameters**: `PascalCase` — `Derived`, `ScalarType`, `FESType`.
- **Enums**: `enum class` exclusively, values in `PascalCase` — `Polytope::Type::Triangle`.
- **Namespaces**: `PascalCase` nested — `Rodin::Geometry`, `Rodin::Variational::P1`.

Favor explicit, domain-correct names over shorthand in public APIs. Prefer `DeformationGradient` over `GradU`, `getShearModulus` over `getMu`, `FirstPiolaKirchhoffStress` over `P`. Avoid cryptic abbreviations in new public APIs.

### Type aliases

Use `using` declarations, never `typedef`:

```cpp
using ScalarType = Real;
using Vertices = std::array<Index, RODIN_MAXIMUM_POLYTOPE_VERTICES>;
```

### Headers and includes

- **File extensions**: `.h` for headers, `.hpp` for heavy template implementation files (included at the end of the corresponding `.h`), `.cpp` for non-template source.
- **Header guards**: `#ifndef RODIN_MODULE_FILENAME_H` / `#define` / `#endif` (not `#pragma once`).
- **Include order**: standard library → third-party (Eigen, Boost) → Rodin headers. Use angle brackets for Rodin includes: `#include <Rodin/Geometry.h>`.
- **Forward declarations**: prefer including `ForwardDecls.h` over full headers when a forward declaration suffices.

### Const correctness

Apply `const` rigorously: `const` on methods that do not mutate, `const` references for input parameters, `const_iterator` variants for containers.

### Return types

Use explicit return types. Avoid `auto` for function return types in declarations. No trailing return types (`-> T`). Reserve `auto` for local variable type deduction and lambda storage.

### Ownership and pointers

- **`std::unique_ptr`**: owning pointers (e.g., cloned operands in expression templates).
- **`std::reference_wrapper`**: non-owning references to external objects (e.g., operands, mesh, FE space).
- **`std::shared_ptr`**: rare, used only for internal lifetime management (e.g., `FormLanguage::Base::m_objs`).
- Raw pointers only at API boundaries with C libraries (e.g., MMG, PETSc).

### Thread safety

`FormLanguage::Base` is explicitly not thread-safe. Use `static thread_local` variables for per-thread output caching in evaluation methods:

```cpp
static thread_local RangeType s_out;
```

Do not introduce mutexes or locks in core evaluation paths.

### Operator overloading

Operator overloads return `*this` by reference for method chaining. `operator=` on `GridFunction` delegates to `project()`. Form composition uses `operator+=`, `operator-=`.

### Fluent API style

Boundary conditions and configuration use fluent chaining:

```cpp
DirichletBC(u, g).on(boundaryAttribute);
```

### Virtual methods

Use `override` on all overrides, `= 0` for pure virtuals, `virtual` destructors (`= default`). Use `final` on leaf classes.

### Doxygen documentation

- Use `/** */` for multi-line documentation blocks, `///` or `///<` for inline/trailing comments.
- File-level: `@file`, `@brief`.
- Classes: `@brief`, `@tparam`.
- Methods: `@param[in]`, `@param[out]`, `@returns`.
- Math: `@f$ ... @f$` for inline LaTeX.

### Optional and common types

`Optional<T>` is an alias for `std::optional<T>` (defined in `Types.h`). Other common aliases: `Index` (`size_t`), `Real` (`double`), `Attribute` (`size_t`), `StringView` (`std::string_view`).

## How to extend Rodin correctly

Use this order:

1. Define the local mathematical object/state.
2. Define the local operator/law.
3. Define the FE integrator / global assembly piece.
4. Define user-facing convenience APIs and examples.

Do not jump straight to a high-level wrapper if lower abstractions are not correct.

## Testing strategy in Rodin

Tests should validate both behavior and abstraction boundaries.

- Use `tests/unit/` for local API and class behavior.
- Use `tests/manufactured/` for numerical correctness.

Add targeted tests for genericity claims (e.g., non-centroid quadrature, alternate FE spaces, heterogeneous attributes) when relevant.

Prefer targeted tests during development; broaden scope only after local confidence.

## Build and validation workflow (optimized)

### Fast local workflow

1. Update submodules if needed:
   ```bash
   git submodule update --init --recursive
   ```
2. Configure an out-of-source build:
   ```bash
   cmake -S . -B build \
     -DCMAKE_BUILD_TYPE=Debug \
     -DRODIN_BUILD_SRC=ON \
     -DRODIN_BUILD_UNIT_TESTS=ON \
     -DRODIN_BUILD_MANUFACTURED_TESTS=ON \
     -DRODIN_BUILD_EXAMPLES=OFF \
     -DRODIN_BUILD_DOC=OFF
   ```
3. Build incrementally:
   ```bash
   cmake --build build -j2
   ```
4. Run the smallest relevant test subset first, using the same `unit` / `manufactured` / `slow` CTest labels used in CI:
   ```bash
   ctest --test-dir build/tests -L unit -LE slow --output-on-failure
   ctest --test-dir build/tests -L manufactured -LE slow --output-on-failure
   ```
5. If public headers / exported targets changed, validate install/downstream usage:
   ```bash
   bash tests/installation/test_installation.sh
   ```

### CI-derived hints (actual workflows)

- Typical CI environment variables:
  - `MAKEFLAGS=-j2`
  - `OMP_NUM_THREADS=2`
  - `CTEST_PARALLEL_LEVEL` set to `1` or `2` depending on suite.
- CI test selection relies on labels:
  - `unit`, `manufactured`, and `slow`.
- Build matrix uses:
  - Ubuntu: gcc-12 for Build workflow
  - Ubuntu: gcc-14 for Tests workflow
  - Ubuntu: gcc-10 for Benchmarks workflow
  - macOS: Homebrew toolchain with OpenMP and MPI variants

## Dependencies (from current workflows)

### Ubuntu CI packages

Core packages installed in Build/Tests/Benchmarks/Copilot setup workflows:

- `libboost1.74-all-dev` (or `libboost-all-dev` in Installation workflow)
- `libsuitesparse-dev`
- `libeigen3-dev`
- `libscotch-dev`
- `libmetis-dev`
- `libhdf5-dev`
- `libomp-dev`
- `petsc-dev`
- `mpich`
- `lcov`

The Tests and Coverage workflow jobs additionally install `gcc-14 g++-14`.

### macOS CI packages

Common Homebrew dependencies across Build/Installation workflows:

- `boost` (and sometimes `boost-mpi`)
- `suitesparse`
- `eigen`
- `scotch`
- `metis`
- `hdf5-mpi` (Build workflow) or `hdf5` (Installation workflow)
- `libomp`
- `petsc`
- `open-mpi`
- `lcov`

## CMake options frequently used in CI

- `RODIN_BUILD_SRC=ON`
- `RODIN_BUILD_EXAMPLES=ON/OFF`
- `RODIN_BUILD_UNIT_TESTS=ON/OFF`
- `RODIN_BUILD_MANUFACTURED_TESTS=ON/OFF`
- `RODIN_BUILD_BENCHMARKS=ON/OFF`
- `RODIN_BUILD_DOC=ON/OFF`
- `RODIN_USE_MCSS=ON/OFF`
- `RODIN_MULTITHREADED=ON/OFF`
- `RODIN_USE_MPI=ON/OFF`
- `RODIN_USE_PETSC=ON/OFF`
- `RODIN_USE_ASAN=ON/OFF`
- `RODIN_USE_UBSAN=ON/OFF`
- `RODIN_CODE_COVERAGE=ON/OFF`

## What not to do

- Do not invent a side architecture parallel to Geometry, Variational, Assembly, or Solver.
- Do not add special-case APIs in core modules for one law/application.
- Do not present narrow prototype code as generic infrastructure.
- Do not hide explicit Rodin composition behind opaque manager objects.
- Do not duplicate geometry/FE concepts already present in the repository.
- Do not document capabilities broader than the implementation.
- Do not use `typedef`; use `using` declarations.
- Do not use `#pragma once`; use `#ifndef` guards following `RODIN_MODULE_FILENAME_H`.
- Do not use raw `throw`; use the `Alert::Exception` stream pattern with `<< Alert::Raise`.
- Do not introduce mutex/lock patterns in core evaluation code; use `thread_local` caching instead.

## What good Copilot contributions look like in Rodin

A good change usually has these properties:

- Code is placed in the correct module.
- Public naming is explicit and domain-correct.
- Responsibilities are cleanly separated.
- Existing abstractions are reused.
- User-facing style remains compositional and recognizable as Rodin.
- Tests validate both behavior and intended generality.
- Code follows the style conventions documented above (Allman braces, 2-space indent, `m_` members, `enum class`, explicit return types, `const` correctness).
- New types provide `FormLanguage::Traits` specializations and `ForwardDecls.h` entries.
- Header-only modules use the INTERFACE library CMake pattern.

## Final heuristic

When in doubt: extend Rodin by adding a small number of strong abstractions that compose with the existing framework, instead of adding a large amount of narrowly working code.
