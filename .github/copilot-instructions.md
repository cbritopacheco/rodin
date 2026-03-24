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
4. Can the user still build problems explicitly using Rodin’s existing style?
5. Am I hard-coding a prototype assumption as if it were generic?

Do not optimize for the shortest patch if it conflicts with Rodin’s structure.

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

Rodin’s style is compositional: trial/test functions, integrals, boundary conditions, explicit `Problem` construction. New modules should preserve this style and avoid introducing “manager owns everything” patterns unless an existing module in the same problem class already does this.

### 4) Do not confuse a working prototype with a generic backend

A prototype that only works for one element, one quadrature rule, one constitutive law, or one geometric assumption must not be represented as a generic module. If the module is intended to be generic, remove prototype assumptions in code, not only in documentation.

## Distinguish core layers explicitly

When reasoning about finite element mechanics and PDE modules, keep these concerns distinct:

- **Geometry** (entities, mappings, points)
- **FE evaluation** (shape functions, gradients, quadrature evaluations)
- **Local constitutive data/laws** (pointwise state and response)
- **Global assembly** (integration and accumulation into global operators)

This separation is critical for maintainability and reusability.

## Preserve the “user builds the problem” philosophy

Prefer APIs that keep user-level problem composition explicit. Rodin users should still be able to assemble and solve by composing trial/test functions, integrals, boundary terms, and solvers directly. Avoid hiding problem definition behind opaque orchestration objects.

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

## Naming rules for public APIs

Favor explicit, domain-correct names over shorthand.

Prefer names like:
- `KinematicState`
- `DeformationGradient`
- `RightCauchyGreenTensor`
- `FirstPiolaKirchhoffStress`
- `GreenLagrangeStrain`
- `MaterialTangent`
- `getShearModulus`
- `getLameFirstParameter`

Avoid cryptic public names like `P`, `getMu`, `getLambda`, or ambiguous abbreviations unless Rodin already standardizes them internally.

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
4. Run smallest relevant test subset first:
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

Additional test/coverage tooling in some jobs includes `gcc-14 g++-14`.

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

## What good Copilot contributions look like in Rodin

A good change usually has these properties:

- Code is placed in the correct module.
- Public naming is explicit and domain-correct.
- Responsibilities are cleanly separated.
- Existing abstractions are reused.
- User-facing style remains compositional and recognizable as Rodin.
- Tests validate both behavior and intended generality.

## Final heuristic

When in doubt: extend Rodin by adding a small number of strong abstractions that compose with the existing framework, instead of adding a large amount of narrowly working code.
