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
\(D\) is pulled-back pointwise deviatoric current
strain: \(h\kappa_D\int j\,\operatorname{dev}\epsilon(v):
\operatorname{dev}\epsilon(z)\), where \(\epsilon(v)=\operatorname{sym}(\nabla v F^{-1})\).
It does not penalize local infinitesimal rotations or isotropic strain. There is
no global dilation coupling. Higher-order spaces can contain additional conformal
kernel modes not covered by the similarity gauge; this term alone is not an
\(H^1\) norm. Their independent coefficients
are kappaF and kappaD, both defaulting to one. D retains the mesh
factor \(h\); there is no shape or constraint Hessian in the metric.
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
hinge solve, bounded by physical predictor motion divided by \(h\). Nonpositive
robust directional curvature falls back to positive weighted fitting curvature.
The outer Armijo line search only backtracks, enforcing the actual \(j\) and \(Q\)
budgets. The geometric stopping target is sampled \(D_\infty\), not robust RMS.
There is no mass completion, diagonal shift, logarithmic barrier, nonlinear-hinge
alternative or separate inertia factorization. Each linear solve selects zero
along unresolved similarity modes using a tiny restricted eigensolve. This gauge
is not part of the objective and releases modes when active hinges resolve them;
it is not a global coercivity certificate.

`WNGIR.h` is the public include; parameters/report and form-language
coefficients live beside `WNGIRSolver.h`. The local Eigen backend supports CG,
SparseLU, and optional MUMPS solving the same pointwise deviatoric operator.
The sparse metric is symmetrized before factorization and residual evaluation,
so triangular direct solvers and the true-residual test use the same operator.
One Problem and metric forms are retained, but deformation-dependent forms are
reassembled per outer iteration. MUMPS retains symbolic analysis while the
augmented sparsity pattern is unchanged and numeric factors for identical systems.
`AnalyticFunctionAdapters.h` lifts analytic lambdas into FunctionBase;
`CellGeomCache` caches per-cell geometry. Historical experimental comparisons
remain in experiments/wngir_calibration; they are not current implementations.

Standing principle regardless of branch: interface-fitting constraints
inside variational solves are smooth penalties, never hard projections
(conventions.md). All derivative code is gated by FD-consistency tests
(testing.md).
