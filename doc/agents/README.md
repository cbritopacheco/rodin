# Agent knowledge base — index

Hierarchical: read top-down, stop at the depth your task needs.

## Level 1 — always

- [philosophy.md](philosophy.md) — the design philosophy: the typed graded
  expression algebra the library models, node patterns, minimality rules,
  surface style. **Read before writing any code.**
- [conventions.md](conventions.md) — hard rules: scope discipline, error
  handling, design preferences, testing conventions.

## Level 2 — orientation

- [architecture.md](architecture.md) — the module map: every directory under
  `src/Rodin/`, what it is, how the layers stack.

## Level 3 — per-domain detail (open the one your task touches)

- [core.md](core.md) — root vocabulary and infrastructure: `Types.h`
  aliases, Copyable/Moveable/Identifiable, Alert, Threads, Utility
  metaprogramming, Context, configuration macros.
- [geometry.md](geometry.md) — the mesh model: polytopes, connectivity,
  transformations, Point evaluation, mesh algebra (skin/trim/keep/trace),
  SubMesh, builders, partitioning, classification utilities.
- [variational.md](variational.md) — the form language in depth: operator
  taxonomy, shape-function grading, finite element spaces (P0, P0g, P1,
  high-order H1), GridFunction, integrals, boundary conditions, Problem.
- [solvers-assembly.md](solvers-assembly.md) — assembly backends and inputs,
  constraint expansion, the linear-solver catalogue, NewtonSolver.
- [physics.md](physics.md) — domain modules: Solid mechanics, the Heart 0D
  model, the level-set toolkit (Distance/Eikonal/Advection/Hilbert),
  Adaptation.
- [integrations.md](integrations.md) — MMG, Scotch, MPI, IO formats,
  Serialization, Python bindings, Blender plugin.
- [petsc.md](petsc.md) — PETSc backend: idioms, CI version skew, known traps.
- [testing.md](testing.md) — test-suite layout, test patterns
  (FD-consistency, manufactured solutions, benchmarks), CI workflows.

## Level 4 — foundations (theory ↔ code)

The mathematics the library implements, each anchored to the classes that
embody it. Open when a task requires understanding *why* the code is
shaped the way it is, or when writing new formulations:

- [theory/variational-formulation.md](theory/variational-formulation.md) —
  weak forms, Galerkin projection, conformity, well-posedness, essential
  vs natural BCs, constraints (elimination/identification/multipliers),
  DG facet terms, Newton linearization.
- [theory/finite-elements.md](theory/finite-elements.md) — the Ciarlet
  triple, reference→physical pullback, DOF distribution over the
  incidence complex, nodal vs modal bases (Fekete/GLL/Dubiner),
  interpolation vs projection, quadrature exactness.
- [theory/meshes-and-geometry.md](theory/meshes-and-geometry.md) — the
  incidence-complex mesh model, charts and metric factors, element
  quality, discrete vs continuous mesh improvement, level sets and
  implicit domains, partitioning.
- [theory/shape-optimization.md](theory/shape-optimization.md) — shape
  derivatives and the Hadamard structure, boundary vs distributed
  gradients, the Hilbertian descent framework, level-set evolution and
  its maintenance obligations, body-fitted vs ersatz discretization.
- [theory/continuum-mechanics.md](theory/continuum-mechanics.md) —
  finite-strain kinematics, hyperelasticity and consistent tangents,
  active contraction and internal variables, 0D reduced models, frame
  conventions.

Procedural how-tos (build, run examples) live in `.claude/skills/*/SKILL.md`
— plain markdown, readable by any agent.
