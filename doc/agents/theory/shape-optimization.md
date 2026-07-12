# Shape and topology optimization — theory ↔ code

Rodin's raison d'être (the README says so). The mathematics here explains
`Hilbert/`, `Distance/`, `Eikonal/`, `Advection/`, the MMG integration,
`Adaptation/`, and the entire examples/ShapeOptimization,
examples/BoundaryOptimization, examples/DensityOptimization families.

## The problem class

Minimize $J(\Omega)$ over admissible shapes, typically
$J(\Omega) = \int_\Omega j(u_\Omega) + \ell\,|\Omega|$ where $u_\Omega$
solves a PDE on $\Omega$ (state equation) and $\ell$ is a volume
penalization (the `--ell` of the cantilever examples is exactly this
Lagrange multiplier). Compliance minimization —
$J = \int_\Omega Ae(u):e(u)$ with linear elasticity as the state — is the
canonical benchmark (cantilever examples).

## Shape derivatives and the Hadamard structure

Perturb the domain by a vector field: $\Omega_t = (I + t\theta)(\Omega)$.
The shape derivative is
$J'(\Omega)(\theta) = \lim_{t\to0} (J(\Omega_t) - J(\Omega))/t$. The
Hadamard structure theorem: for smooth data the derivative depends only on
the normal trace on the boundary,

$$J'(\Omega)(\theta) = \int_{\partial\Omega} g\; \theta\cdot n\, ds,$$

with a scalar density $g$ (for compliance with traction-free optimized
boundary: $g = Ae(u):e(u) - \ell$). Computing $g$ for PDE-constrained $J$
uses the adjoint method; compliance is self-adjoint, which is why the
examples get away without a separate adjoint solve.

Two discrete realizations of the same derivative — know which one a code
path uses:

- **Boundary (Hadamard) form**: assemble $\int_\Gamma g\,\theta\cdot n$ on
  interface facets. Concentrated on the discrete interface; pointwise
  noisy on curved/graded meshes (the density is quadratic in strain and
  multiplies inverse-metric factors with no $\det J$ cancellation).
- **Distributed (volume/Eshelby) form**: integrate
  $\int_\Omega E : \nabla\theta$ with the Eshelby-type tensor; integration
  by parts (using the state equation) recovers the boundary form
  analytically. Discretely it is a *smoother, more robust* gradient — the
  preferred form here. Correctness arbiter for any change to such a
  density: the **directional FD test** — compare
  $\int E : \nabla\theta$ against a central finite difference of the
  discrete objective under mesh perturbation $\pm t\theta$ for a
  localized $\theta$. An Eshelby density that fails this test leaves a
  bulk residual that does not vanish at optimality (symptom: the descent
  keeps flagging regions the interface already reached).

## From derivative to descent: the Hilbertian framework

$J'(\Omega)$ is a linear functional; to descend you need a *velocity
field*, i.e. a Riesz representative in a chosen Hilbert space:

$$\text{find } \theta \in H:\quad b(\theta, \xi) = -J'(\Omega)(\xi)
\quad \forall\xi \in H,$$

with $b$ an $H^1$-type inner product (mass + $\mathrm{ell}^2$·stiffness).
This is `Hilbert::H1a` — extension-regularization in one step: the
boundary data spreads into the bulk and is smoothed, with the length
scale set by the regularization parameter. Consequences you can predict
from the theory: larger length ⇒ smoother, more global descent
directions, fewer oscillations, slower fine-scale detail; too small ⇒
mesh-frequency noise in the velocity. (This is the shape-gradient `--ell`
distinct from the volume multiplier — two different ells in the
examples' CLI, both theory-mandated.)

Descent then must be *safeguarded*: $J(\Omega_{t})$ is nonconvex, and an
unguarded step can run past the basin (observed failure mode: the shape
floods the whole domain). An Armijo line search on the true objective is
part of the algorithm, not an option.

## Level-set representation of the evolving shape

Represent $\Omega = \{\phi < 0\}$ and move it with the normal velocity
$V = \theta\cdot n$ via the Hamilton–Jacobi equation
$\partial_t\phi + V|\nabla\phi| = 0$ (or transport
$\partial_t\phi + \theta\cdot\nabla\phi = 0$). This handles topology
changes (holes merging/splitting) for free — that is the entire point of
level-set shape optimization vs boundary-tracking. The maintenance
obligations follow from the theory:

- **Redistancing**: transport degrades $|\nabla\phi|$; periodically
  restore the signed-distance property (`Eikonal::FMM`,
  `Distance::*`). Skipping it corrupts every downstream quantity that
  assumes $|\nabla\phi| \approx 1$ (normals, band widths, CFL). A subtle
  discrete trap that has actually bitten: a redistancer that takes the
  *sign* from a stale cell classification can flip regions the interface
  just crossed — creating spurious zero crossings. When fit and
  classification move the interface, sign and distance must come from the
  same, current geometry.
- **CFL-type step control**: the interface should move a bounded number
  of cells per step (`--dt-max` style caps).

The loop of a level-set shape optimization run (all the
LevelSet\*Cantilever examples):
solve state → evaluate objective + shape derivative → Riesz/extend
(`H1a`) → line-search the step → advect φ → redistance → re-discretize
the domain for the next state solve.

## Getting a computable domain from φ

Two standard discretizations of $\{\phi<0\}$, both present:

1. **Body-fitted remeshing** (MMG `LevelSetDiscretizer` / `Adapt`; native
   cutter + optimizers in `Adaptation/`): cut the zero set into the mesh,
   producing interface-conforming cells and an interface attribute.
   Accurate integration on $\Omega$ and $\Gamma$; the cost is mesh
   surgery each step and the resulting element-quality management
   (meshes-and-geometry.md).
2. **Ersatz/density approaches** (examples/DensityOptimization): keep a
   fixed mesh, fill the void phase with weak material (or interpolate
   material by a density) — SIMP-style topology optimization. No
   remeshing, but a modeling error controlled by the ersatz contrast.

The house layering rule when body-fitting a moving interface
(philosophy.md): φ owns topology; a classifier turns φ into discrete cell
attributes; mesh optimization owns geometry/quality and *never* decides
topology. Interface-fitting constraints inside variational mesh solvers
are smooth penalties, never hard projections (conventions.md).

## What to measure

Shape-optimization work is judged by run-level evidence, per the testing
conventions: monotone (line-searched) objective histories, volume
tracking the multiplier's target, interface fit metrics per stage, mesh
validity (0 inverted), and — for any new derivative — the directional FD
test above. "The shape looks right" is not evidence.
