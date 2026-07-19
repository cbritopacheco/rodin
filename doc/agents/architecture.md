# Repository architecture

The stable shape of the codebase (~140k lines of C++20, header-heavy).
Details per domain live in the level-3 docs (see README.md); this file is the
map.

## Layout

```
src/Rodin/          The library
tests/unit/         GoogleTest, mirrors src/Rodin module tree
tests/manufactured/ Convergence + targeted regressions (incl. cross-space and PETSc)
tests/benchmarks/   Google Benchmark (target RodinBenchmarks)
examples/           One directory per theme; each example is a CLI executable
third-party/        Git submodules (eigen, googletest, ...)
py/                 pybind11 Python bindings (scikit-build / pyproject.toml)
plugins/            rodin_blender — Blender integration
resources/          Meshes and data used by examples/tests (medit, mfem, gmsh)
doc/agents/         This knowledge base
.claude/skills/     Procedural skills (plain markdown, any agent can read)
```

## The layer stack

Modules stack strictly; lower layers never include higher ones:

```
1. Core vocabulary   Types.h, Copyable/Moveable/Identifiable, Cast, Array,
                     Tuple/Pair/FlatSet, Alert, Threads, Utility, Context
2. Math              Eigen-based dense/sparse types, SpatialVector/Matrix,
                     LinearSystem, angles/units, root-finding, RK schemes
3. Geometry + QF     Mesh = graded incidence complex + attached
                     transformations; quadrature formulas on reference
                     polytopes
4. FormLanguage      Base node class, Traits, polymorphic List
5. Variational       The weak-form DSL: functions, operators, FES,
                     GridFunction, integrals, BCs, Problem
6. Assembly + Solver Backends that evaluate integrators over meshes;
                     linear/nonlinear solver wrappers
7. Domain modules    Solid, Heart, Distance/Eikonal/Advection/Hilbert,
                     Adaptation, Location
8. Integrations      PETSc, MMG, MFEM, Scotch, MPI, IO, Serialization
```

## Module inventory (src/Rodin/)

Sized by header count so you know where the mass is: Variational 127,
Geometry 34, PETSc 32, MPI 24, Solid 23, Utility 22, Solver 22, Alert 20,
Math 17, Heart 14, Adaptation 13, Serialization/IO 11 each, MMG 10, QF 7,
Assembly 7, FormLanguage 6, Distance 6, Test 6.

- **Root headers** — `Types.h` (house aliases), `Copyable.h`/`Moveable.h`
  (polymorphic clone/move interfaces), `Identifiable.h` (UUID identity),
  `Cast.h`, `Array.h`, `Tuple.h`/`Pair.h`/`FlatSet.h`, `Threads.h`
  (Mutex/Unsafe/Mutable wrappers), `Configure.h.in` (feature macros).
  → core.md
- **Alert/** — streamable exception/warning system
  (`Exception() << "..." << Raise;`), colored terminal output. → core.md
- **Context/** — execution contexts: `Context::Local` and `Context::MPI`.
  Meshes, connectivities, and assemblies are templated on the context; the
  MPI stack mirrors the local one. → core.md
- **Utility/** — template metaprogramming toolkit (ForConstexpr,
  ParameterPack, Zip, Product, IsSpecialization, Overloaded, ...). → core.md
- **Math/** — `Vector`/`Matrix` (dynamic, Eigen), `SparseMatrix`,
  `SpatialVector`/`SpatialMatrix` (fixed capacity 3, runtime dims, no heap),
  `LinearSystem` (Ax=b + DOF elimination), `Rad`/`Deg` unit types,
  NewtonRaphson, RungeKutta, Constants, Traits. → core.md
- **Geometry/** — Mesh/SubMesh, Polytope + symmetric-hash `Key`,
  Connectivity (incidence d→d′), transformations (Identity/Parametric),
  `Point` (reference↔physical), PointCloud, iterators, mesh algebra
  (`skin`/`trim`/`keep`/`trace`, `UniformGrid`), Builder, partitioners
  (Greedy/BalancedCompact + Sharder/Shard), `MinSTCut` (s-t cut
  classifier). → geometry.md
- **QF/** — quadrature on reference polytopes: GaussLegendre, GaussLobatto,
  GrundmannMoller (simplex), Centroid; `PolytopeQuadratureFormula`
  dispatcher. → geometry.md
- **FormLanguage/** — `Base` (clone + identity root), `Traits` (open
  type-computation), `List` (polymorphic node container). → variational.md
- **Variational/** — the language: ~60 operator headers (arithmetic,
  calculus, transcendental, comparison, boolean, complex, DG traces),
  function bases (Real/Complex/Vector/Matrix/Boolean), shape functions,
  FES families (P0, P0g global-constant, P1, high-order H1 with
  Fekete/GLL/Dubiner bases), GridFunction, integral families
  (Integral/BoundaryIntegral/FaceIntegral/InterfaceIntegral), DirichletBC
  (value + identification flavours), PeriodicBC, Problem/ProblemBody,
  DenseProblem/SparseProblem, Potential (integral equations), TraceOperator,
  Flow. → variational.md
- **Assembly/** — Sequential and OpenMP mesh-iteration backends,
  `BilinearFormAssemblyInput` (FES + local & global integrator lists),
  `ConstraintMap` (value/identification constraint expansion), `Default.h`
  strategy selection. → solvers-assembly.md
- **Solver/** — ~20 `LinearSolver` wrappers over Eigen (CG, BiCGSTAB, GMRES,
  DGMRES, MINRES, IDRS, LeastSquaresCG, SparseLU, SimplicialLLT/LDLT,
  HouseholderQR, PartialPivLU), SuiteSparse (UMFPack, CHOLMOD, SPQR),
  AppleAccelerate; `NewtonSolver` (constant damping, no line search).
  → solvers-assembly.md
- **Solid/** — finite-strain solid mechanics: Kinematics (KinematicState,
  Invariants), Constitutive laws (Hooke, NeoHookean, MooneyRivlin,
  SaintVenantKirchhoff, HolzapfelOgden, ActiveFiberLaw/ActiveContraction),
  Integrators (InternalVirtualWork residual/tangent façade), Fields (stress/
  strain postprocessing), Local/ConstitutivePoint (tag-based input
  injection), Linear/ (linear elasticity integrators). → physics.md
- **Heart/** — CCMLC2014 0D reduced ventricular model (stepper + reduced
  Holzapfel passive-energy law). → physics.md
- **Distance/, Eikonal/, Advection/, Hilbert/** — the level-set toolkit:
  distance models (Eikonal, Poisson, SignedPoisson, Rvachev,
  SpaldingTucker), FMM fast marching, Lagrangian advection, H¹
  extension-regularization (`H1a`). → physics.md
- **Adaptation/** — moving-interface mesh adaptation: the WNGIR module
  (Welsch natural-gradient interface registration) plus analytic-function
  adapters and a per-cell geometry cache. Related work (TMOP-style target
  matrix optimization, native remeshers) lives on other branches — verify
  presence before referencing it. → physics.md
- **Location/** — AABB bounding-volume-hierarchy point locator.
- **PETSc/** — mirror of Math/Solver/Assembly/Variational/IO over PETSc
  objects. → petsc.md
- **MMG/, Scotch/, MPI/, IO/, Serialization/, MFEM/** — integrations and
  persistence. MFEM/ is currently an empty stub (the MFEM *file format*
  lives in IO/). → integrations.md
- **Test/** — library-side test utilities (Random functions, Utility)
  used by the test suites. → testing.md

## Cross-cutting invariants

- **Topology vs geometry separation** at every level: incidence complex vs
  transformations (Geometry); level-set/classifier vs mesh optimizer
  (Adaptation). A layer never decides what belongs to another.
- **Context mirroring**: distributed functionality lives in parallel module
  trees (`MPI/Geometry`, `MPI/Assembly`, ...) rather than `#ifdef`s in core.
- **Backends are additive**: Eigen is the native path; PETSc/MMG/Scotch
  attach from their own directories; core never includes them.
