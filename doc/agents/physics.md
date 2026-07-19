# Domain modules

The physics/application layers built on the form language. Structural
facts only — parameter tuning and experiment history do not belong in the
repo.

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

On this branch the module is **WNGIR** — Welsch natural-gradient interface
registration: `WNGIR.h` public include; `WNGIRParameters`/`WNGIRReport`,
backend-independent `WNGIRSolver`, metric/force/admissibility integrator
headers (`WNGIRSurfaceObservationMetric`, `WNGIRSurfaceForce`,
`WNGIRAdmissibility*`), and a normal-offset initializer. It is the default
displacement model for fitting a mesh interface to a level set;
dimension-generic, built in the form language; the 2D examples use the
native Eigen path (PETSc-backed 3D solver variants live on the
module/Adaptation lineage). `AnalyticFunctionAdapters.h` lifts analytic lambdas
into `FunctionBase`; `CellGeomCache` caches per-cell geometry.

Standing principle regardless of branch: interface-fitting constraints
inside variational solves are smooth penalties, never hard projections
(conventions.md). All derivative code is gated by FD-consistency tests
(testing.md).
