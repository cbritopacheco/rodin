# Domain modules

The physics/application layers built on the form language. Structural
facts only — parameter tuning and experiment history do not belong in the
repo.

For Solid or level-set terms that assemble residuals, tangents, constraints,
or projection-like transfers, keep [workflows.md](workflows.md) and
[numerical-contracts.md](numerical-contracts.md) open with this file.

## Solid (src/Rodin/Solid/) — finite-strain solid mechanics

Layered like the mathematics:

- **Kinematics/** — `KinematicState` (F, C, invariants at a point),
  `Invariants` (isotropic I1..I3 + anisotropic I4f etc. for fiber
  families).
- **Constitutive/** — one law per header, all deriving from
  `HyperElasticLaw`: `Hooke` (linear), `NeoHookean`, `MooneyRivlin`,
  `SaintVenantKirchhoff`, `HolzapfelOgden` (one fiber family),
  `ActiveFiberLaw` (1D Hill–Maxwell active contraction law with internal
  state γ, β), `ActiveContraction` (wraps any passive law + active fiber
  law; static path takes activation via tags, dynamic path runs a local
  update).
- **Local/** — `ConstitutivePoint`: the *tag-based input injection*
  mechanism. Integrators stamp tags (CellIndex, QuadraturePointIndex,
  TimeStep, ElectricalActivation, previous internal state, ...) onto the
  point; a user-supplied input lambda reads tags and feeds the law. This
  keeps laws decoupled from where their inputs come from. `FiberKinematics`
  carries preferred directions.
- **Integrators/** — `InternalVirtualWork` façade =
  `InternalVirtualWorkResidual` + `InternalVirtualWorkTangent` (the
  nonlinear-elasticity pair for NewtonSolver); `Linear/
  LinearElasticityIntegral` for the linear case.
- **Fields/** — postprocessing functions: `FirstPiolaKirchhoffStress`,
  `CauchyStress`, `GreenLagrangeStrain`.

House rule (conventions.md): internal variables (active extension,
γ/β-type state) are **first-class DOFs** on discontinuous spaces solved by
the global Newton — not per-quadrature Schur condensation.

## Heart (src/Rodin/Heart/)

`CCMLC2014` — 0D reduced ventricular model: a stepper over a reduced
Holzapfel-like passive-energy law (`HolzapfelReducedLaw`, `PassiveEnergy`
with reduced invariants and first/second derivatives). Pairs with the
Solid active law for lumped-parameter studies (examples/Heart,
examples/Models).

## Level-set toolkit

- **Distance/** — distance/redistancing models behind a common `Base`:
  `Eikonal` (PDE-based), `Poisson` and `SignedPoisson` (Poisson-equation
  approximations), `Rvachev` and `SpaldingTucker` (normalizations of an
  existing level set).
- **Eikonal/FMM** — fast marching method on simplicial meshes (the
  workhorse redistancer; "fmm" in example flags).
- **Advection/Lagrangian** — Lagrangian variational advection of scalar
  fields (level-set transport).
- **Hilbert/H1a** — the Hilbertian H¹ extension-regularization of shape
  derivatives (velocity extension for shape optimization; the `--ell`
  regularization length in the shape-opt examples).

Advected level sets must be periodically redistanced — a carried φ drifts
regardless of fit quality.

## Adaptation (src/Rodin/Adaptation/)

Moving-interface mesh adaptation under the layering invariant
(philosophy.md): the level set owns topology, classification owns discrete
attributes (`Geometry::MinSTCut` is the classifier primitive), geometry
fitting owns node positions and never decides topology.

On this branch the module is **WNGIR**: Welsch natural-gradient interface
fitting. The unique model is \(M=F+D\) with affine quadratic hinges.
\(F\) is normalized half-squared fitting curvature without the level-set Hessian;
\(D\) penalizes deviations of current deviatoric strain and divergence from
their global current-volume averages. These averages are not elementwise.
Coherent affine motion is free, including translation, infinitesimal rotation,
uniform dilation, shear and anisotropic stretching. Shared finite element DOFs
couple elements. With both distribution coefficients positive, distribution
controls \(H^1\) modulo global affine fields on a fixed connected Lipschitz domain;
it is not a full-space norm or an unconditional uniform refinement bound.
The defaults are `model.fit = 1`, `model.distribution.deviatoric = 1e-4`,
`model.distribution.divergence = 1e-2` and `model.hinge = 10`.
D retains the mesh factor \(h\); there is no shape or constraint Hessian in the metric.
The fitting energy/force remain robust Welsch. Minimizing this integral is not
equivalent to minimizing the maximum geometric residual. Power-loss experiments
are isolated diagnostics, not a production option or schedule.
Surface integration uses at least order 12 and geometric sampling at least 14
(vertices are included separately). These are safeguards for non-polynomial
targets, not exact-integration or Hausdorff guarantees; bulk integration and
physical displacement sampling retain their own order rules.
Inner Newton uses a frozen merit
line search and a stationarity residual relative to the fixed fitting force.
Directional Newton scales the predictor and frozen quadratic form before the
hinge solve. The optional physical-motion bound is disabled by default. Nonpositive
robust directional curvature falls back to positive weighted fitting curvature.
The outer Armijo line search only backtracks, enforcing the actual \(j\) and \(Q\)
budgets. The geometric stopping target is sampled \(D_\infty\), not robust RMS.
There is no mass completion, diagonal shift, logarithmic barrier, nonlinear-hinge
alternative, gauge or separate inertia factorization. Unresolved modes remain
subject to the selected linear backend; no uniqueness correction is imposed.

`WNGIR.h` is the public include; the implementation lives in `Adaptation/WNGIR/`,
without a `Detail` layer. All method-specific types belong to
`Rodin::Adaptation::WNGIR`, without redundant class prefixes.
`Problem.h` retains the outer orchestration and metric
Problem. `HingeProblem.h` assembles the tangent and negative stationarity
residual for Rodin's `Solver::NewtonSolver`; its step policy retains the frozen
inner merit search. `LinearSolver.h` adapts the retained
linear backends to the native solver interface. The local Eigen backend supports CG,
SparseLU, and optional MUMPS solving the same globally centered operator.
`Adapt.h` owns a local vector P1 space and `WNGIR::Problem`, applying accepted
quality-valid displacements to affine simplicial meshes without changing topology.
`WNGIR::Problem(u, v)` follows supplied storage and currently supports only Eigen.
`WNGIR::Adapt(mesh)` defaults to Eigen for local meshes; MPI requires PETSc and
construction is disabled until the distributed specialization is implemented.
The sparse metric is symmetrized before factorization and residual evaluation,
so triangular direct solvers and the true-residual test use the same operator.
Metric and hinge Problems and metric forms are retained, but deformation-dependent forms are
reassembled per outer iteration. MUMPS retains symbolic analysis while the
augmented sparsity pattern is unchanged and numeric factors for identical systems.
`AnalyticFunctionAdapters.h` lifts analytic lambdas into FunctionBase;
`CellGeomCache` caches per-cell geometry. Historical experimental comparisons
remain in experiments/wngir_calibration; they are not current implementations.

Standing principle regardless of branch: interface-fitting constraints
inside variational solves are smooth penalties, never hard projections
(conventions.md). All derivative code is gated by FD-consistency tests
(testing.md).
