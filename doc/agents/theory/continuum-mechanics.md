# Continuum mechanics — theory ↔ code

The mathematics behind `src/Rodin/Solid/` and `Heart/`. Companion to
[../physics.md](../physics.md).

## Kinematics

Motion $x = X + u(X)$; deformation gradient $F = I + \nabla_X u$; right
Cauchy–Green tensor $C = F^\top F$; Green–Lagrange strain
$E = \tfrac12(C - I)$. Isotropic invariants
$I_1 = \operatorname{tr} C$, $I_2$, $I_3 = \det C = J^2$; a fiber family
with direction $f_0$ adds anisotropic invariants, chiefly
$I_{4f} = f_0\cdot C f_0$ (squared fiber stretch).

Code: `Solid::KinematicState` carries $(F, C, \dots)$ at a constitutive
point; `Kinematics/Invariants.h` computes the isotropic + anisotropic
sets; `Fields/GreenLagrangeStrain` exposes $E$ as a postprocessing
function. Small-strain theory linearizes to
$\varepsilon = \tfrac12(\nabla u + \nabla u^\top)$ — that regime lives in
`Solid/Linear/` and is a *different bilinear form*, not a special case
switch.

## Hyperelasticity

A strain-energy density $W(F)$ (frame-indifferent ⇒ $W = \hat W(C)$)
generates stress by differentiation:
first Piola–Kirchhoff $P = \partial W/\partial F$; Cauchy stress
$\sigma = J^{-1} P F^\top$ (`Fields/FirstPiolaKirchhoffStress`,
`Fields/CauchyStress` are exactly these formulas). The catalogue in
`Constitutive/` is the standard ladder of $W$:

- `Hooke` — quadratic in $\varepsilon$ (linear elasticity);
- `SaintVenantKirchhoff` — quadratic in $E$ (large rotation, small
  strain; loses ellipticity in strong compression — a known theory
  limitation, not a bug);
- `NeoHookean`, `MooneyRivlin` — compressible rubber-type laws in
  $I_1, I_2, J$;
- `HolzapfelOgden` — anisotropic soft tissue: isotropic matrix +
  exponential fiber term in $I_{4f}$ (one fiber family here);
- `ActiveFiberLaw` / `ActiveContraction` — see below.

All derive from `HyperElasticLaw`: a law's contract is to return, at a
constitutive point, the stress and the consistent material tangent
$\mathbb C = \partial P/\partial F$ (or equivalents).

## The weak form and its linearization

Static equilibrium in the reference configuration, weak form (total
Lagrangian):

$$R(u; v) = \int_{\Omega_0} P(F) : \nabla_X v \; dX - \ell_{ext}(v) = 0.$$

Newton's method needs the directional derivative

$$R'(u)(\delta u, v) = \int_{\Omega_0} \nabla_X\delta u\;
\frac{\partial P}{\partial F}\; : \nabla_X v\; dX,$$

which contains both the *material* stiffness (from
$\partial^2 W/\partial C^2$) and the *geometric* (initial-stress)
stiffness — using only the material part is the classic
loses-quadratic-convergence error. Code:
`Integrators/InternalVirtualWorkResidual` and
`...Tangent` are precisely $R$ and $R'$, packaged by the
`InternalVirtualWork` façade for `Solver::NewtonSolver`; the
FD-consistency tests (testing.md) enforce tangent = derivative of
residual.

The `ConstitutivePoint` tag-injection design (physics.md) is the seam
between the *field-level* integrators (which know meshes and quadrature)
and the *point-level* laws (which know only kinematic state + inputs):
integrators stamp CellIndex/QuadraturePointIndex/time-step/activation
tags; a user lambda maps tags to the law's inputs. Laws stay pure
functions of local state — the property that makes them unit-testable
against FD at a single point.

## Active contraction and internal variables

Muscle-type models (Hill–Maxwell / active-strain-rate family) augment the
passive law with a 1D fiber element carrying *internal state*: active
extension $e_c$ and stiffness/stress-like variables ($\gamma, \beta$ here)
evolving by ODEs driven by electrical activation $|u|$ and mechanical
sliding. Two structural facts with design consequences:

- **Internal variables are first-class DOFs** in this codebase
  (conventions.md): they live as `GridFunction`s on discontinuous spaces
  and enter the *global* coupled Newton residual, rather than being
  eliminated per quadrature point (Schur condensation). The mixed
  "condense locally, correct the tangent" alternative produced
  confirmed-wrong tangents and is rejected.
- **Quasi-statics can latch.** In a quasi-static setting, after
  activation ceases, the contracted configuration is itself an
  equilibrium ($\delta = \Delta e_c = 0$ ⇒ no decay path), so the model
  never relaxes — resolved either by a spontaneous decay rate parameter
  (non-physical regularization) or properly by dynamics: with inertia
  ($\rho\ddot u$ term, Newmark/generalized-α integration) the unbalanced
  elastic force drives the sliding that detaches cross-bridges. When a
  Solid example "won't relax", check the formulation regime before the
  law.

## Reduced (0D) models

`Heart/CCMLC2014` is the dimensional-reduction companion: the ventricle
as a lumped mechanical model — reduced invariants standing in for the 3D
kinematics, a reduced Holzapfel-type passive energy
(`HolzapfelReducedLaw`, with first and second derivatives — the 0D
"stress" and "tangent"), stepped through a cardiac cycle by an ODE
integrator (`Math::RungeKutta` territory, not the FEM Newton). 0D models
are the calibration/validation harness for the 3D laws: parameters should
round-trip between `Solid::ActiveFiberLaw` and the 0D cycle before
anyone debugs a 3D run.

## Frame conventions to keep straight

- Everything Solid assembles is **referential** (total Lagrangian):
  integrals over $\Omega_0$, operators are $\nabla_X$, stress is $P$.
  Cauchy stress is a *postprocess* (`Fields/CauchyStress`), not an
  assembly quantity.
- Fiber directions (`Local/FiberKinematics`) are material vectors
  ($f_0$, in the reference configuration); their pushforward $F f_0$
  happens inside laws/fields, never in user setup.
- Units are the user's contract (the laws are unit-agnostic); mixing a
  Pa-scale modulus with a kPa-scale activation is a silent modeling
  error the code cannot catch.
