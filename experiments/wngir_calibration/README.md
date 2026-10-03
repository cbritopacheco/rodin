# Canonical WNGIR Campaigns

WNGIR now has one model: $M=F+D$ (Fitting, Distribution), affine quadratic hinges,
directional Newton, frozen inner-merit backtracking, and actual outer j/Q
and fitting-energy Armijo checks. No shape or constraint Hessian is assembled
in the canonical metric.
There is no selectable full-Hessian or logarithmic model, coefficient-space
completion, mass term, or nonlinear-hinge path. PSD does not imply invertibility.

- F is kappa_f times the normalized Hessian of half the squared level-set residual with
  the level-set Hessian omitted. It is not robust-weighted.
- $D$ is $h\kappa_D$ times pointwise deviatoric current-configuration strain.
  Local infinitesimal rotations and isotropic strain have zero cost.
- The fitting energy and force remain robust Welsch.
- Hinge widths are 0.1 times the identity quality margins. The default
  model-decrease-scaled hinge weight is mu_hat=90.
- Defaults are 30 outer / 15 inner corrections and $\kappa_F=\kappa_D=1$.
  There is no shared bulk coefficient. These are choices, not a new calibration result.
  Historical shape-curvature runs are retained separately, not reproduced by
  the canonical metric away from identity.

There is no separate inertia audit or automatic metric repair. Linear residual,
direction, inner merit and actual outer quality/energy checks remain in place.
Fitting must resolve the similarity modes that D leaves free;
for example, a planar interface cannot identify tangential translation.
Both metric terms are PSD on admissible states. Unresolved conformal modes
can still prevent coercivity; testing similarity modes does not exhaust the
kernel for every finite-element space. No general invertibility claim is made.
The distribution integrators tabulate frozen current strain once per basis and
quadrature point, rather than reevaluating the inverse deformation for each
basis pair. MUMPS retains symbolic analysis while the sparsity pattern matches,
and numeric factors while entries match exactly. An entirely inactive hinge
model uses the predictor directly without assembling and solving a zero correction.
Geometry traces record inactive-hinge skips and analysis/factorization counts.

## Remaining Parameters

| Role | Parameters | Defaults |
| --- | --- | --- |
| Metric | kappa_f, kappa_d | 1, 1 |
| Hinge strength | mu_hat; kappa_j, kappa_q | 90; 1, 1 |
| Guard width | quality_guard | 0.1 of identity margins |
| Quality budget | j_safe, q_max; actual line-search Jacobian floor | 0.01, 10; 0.01 |
| Robust fitting | robust_scale | automatic |
| Iteration limits | outer, inner; CG per linear solve | 30, 15; 1000 |
| Inner stopping | stationarity residual relative / absolute | 1e-3 / 1e-12 |
| Geometry target | geometric_sup_tolerance | 0 selects h^(p+1), hence h^2 for P1 |

Legacy active-fit tolerance paths and predictor fallbacks are removed. The primary
calibration controls are the two metric weights and mu_hat; the quality budget
is prescribed rather than inferred from fit.

Termination is uniform:
- Success requires the full-interface sampled D_inf target and the actual sampled j/Q budget.
- Inner convergence requires ||Mv-f+DB(v)||_2 <= atol + rtol*scale, where
  scale=max(initial inner residual norm, ||f||_2). The accepted inner iterate is
  reassembled and checked, including after the last permitted correction.
- Accepted physical displacement is evaluated from the FE field at cell vertices
  and validation quadrature, not from modal coefficients. Small steps or small
  relative energy changes must persist for the same five-iteration window.
- Stagnation, caps, invalid diagnostics and failed solves are not target success.
- Linear solves must meet the requested residual; no hidden acceptance floor remains.

D_inf here is explicitly a sampled normalized level-set residual |phi|/|grad phi|.
Interface vertices supplement quadrature maxima. This is not a certified
Hausdorff distance or proof of geometric order. Invalid samples invalidate the
whole diagnostic irrespective of the target value. The automatic h^(p+1) target
is a screening budget with constant one, not a measured approximation constant.

Historical campaigns retain their frozen binaries and stopping protocols.
They are not pooled with canonical results. New schemas/model identifiers
reject resuming those campaigns.

## Model

For A_k = I + grad u_k, let e_k(v) = sym(grad v A_k^-1), j_k = det A_k,
V_k = integral j_k and g_k = grad phi(x+u_k). With the frozen gradient-scale
normalization N, the metric is

```text
F[v,z] = kappa_f N integral_interface (g_k.v)(g_k.z)
D[v,z] = h kappa_d integral_volume j_k dev(e_k(v)):dev(e_k(z))
```

The unique metric is fitting plus pointwise distribution. Local isotropic
strain and infinitesimal rotations have zero distribution cost. Conforming
degrees of freedom enforce compatibility. No shape/constraint Hessian or
global dilation subtraction remains. Higher-order conformal modes require
separate kernel checks; similarity gauges alone are not a coercivity proof.

The frozen affine slacks are s_J = j_k-j_safe+Dj_k[grad v] and
s_Q = Q_max-Q_k-DQ_k[grad v]. With delta_J=guard*(1-j_safe),
delta_Q=guard*(Q_max-1), the inner objective is

```text
0.5 M[v,v] - f_k[v]
+ 0.5 mu_k integral_volume sum_{a=J,Q} kappa_a (1-s_a(v)/delta_a)_+^2
mu_k = mu_hat * (0.5 max(0,f_k[predictor])) / reference_domain_volume
f_k = -DE_W(u_k)
```

Negative affine slack has finite penalty; no fraction-to-boundary step is
used. Inner merit backtracking handles active-set changes. The nonlinear
outer update still requires actual Jacobian and distortion admissibility
and sufficient decrease of E_W. Metrics do not add a quality force to E_W.

The 2D runner uses MUMPS and independent --kappa-f/--kappa-d grids,
each defaulting to 1e-4,1e-3,1e-2,0.1,1. With mu_hat=0.1,1,10,100,1000,
all five resolutions and eleven lobe counts give 6,875 fitting-distribution cases.
Shape-curvature flags and weights are removed.
The 3D runner uses fixed --kappa-f and --kappa-d/--mu-hat grids.
Executable defaults remain all one. Executables use
--wngir-kappa-f and --wngir-kappa-d; the old metric flags
are rejected. Hinges retain mu_hat, kappa_j, kappa_q and the quality guard.
Use a new output directory: model identifiers and schemas deliberately reject
resuming the historical logarithmic campaigns. Historical calibration data,
frozen binaries, logs, and investigation reports are preserved but are not
evidence of calibration of this canonical model. Retired experimental runners
and experimental CMake targets have been removed.

Both runners default to 30/15. The broad 3D launcher explicitly selects 20/10.
Neither runner is launched automatically by this change.

Each case saves a lossless trace in `<out-dir>/iteration_logs/`, with its case
identity and command in the first two lines. New CSV records identify
kappa_f, kappa_d and mu_hat and cannot resume older coefficient schemas.
Inner rows contain stationarity residual, its tolerance and relative residual, plus
Newton correction and iterate norms as diagnostics,
step factor, linear iteration count/error, and convergence status; failed linear
solves and infeasible corrections also produce rows. Geometry rows contain
complete-interface RMS/maximum distance, normal error, Jacobian and distortion,
cumulative Newton corrections, and elapsed solver wall time including setup.
They cover the initial state, each accepted outer state, and the terminal state.
The `outer` index on geometry rows counts accepted outer steps; inner rows use
zero-based outer attempt indices. Inner work does not itself move the accepted mesh.

Tracing evaluates geometry additionally and therefore incurs measurement cost.
It does not change the fitting energy, damping, or stopping conditions. Timings
include this logging/validation cost and must be compared with like-for-like runs.

To compute first-hit counts after selecting independent per-lobe reference
coefficients, supply a JSON map from lobe index to `C_infty`:

```sh
python3 experiments/wngir_calibration/summarize_iteration_logs.py \
  /path/to/campaign/iteration_logs --coefficients /path/to/reference_coefficients.json \
  --order=2 --output=/path/to/campaign/first_hits.csv
```

The P1 target is `D_infty <= C_infty(lobes) * h^2`, with sampled Jacobian greater
than 0.01 and distortion strictly below 10, matching the runners' admissibility
profile. Coefficients are required, not inferred from candidate terminal errors.
The summary reports inner attempts (including failed linear solves), their mean,
median, and maximum per attempted outer solve, damping/failure counts, and first-hit outer
count, cumulative completed Newton corrections, and time. Missing hits remain
missing rather than being replaced by the iteration cap. A sampled target hit at
one resolution does not establish an error order or certify a Hausdorff bound.

## Historical Cleanup Regression

The 2026-10-02 PSD-only/residual-stop cleanup passed 35 scoped C++ tests,
five manufactured assembly tests and 12 runner/trace tests. The 2D and 3D
reconstruction and sweep targets build. Disabling interface agreement and
automatic-target validity checks made their regressions fail; replacing d/4
by d/2 failed the independent projected-curvature oracle in P1/P2 and 2D/3D.
An exact P2 quadratic interior-maximum test also distinguishes physical-field
measurement from a coefficient maximum.

End-to-end response checks used single-thread MUMPS in
`/tmp/wngir-canonical-smoke.kLPd9f`. At n=10, four lobes, coefficients
(F,S,D,mu_hat)=(1,1e-4,1e-3,0.1), the 2D three-step check reached the sampled
h^2 target (D_inf=0.010402491178703134); the 3D one-step check retained quality
but did not reach it (D_inf=0.039010266465819246). These are output/termination
checks, not a calibration comparison. The attempted 3D n=5 check had no
classified interface (zero inside cells and zero facets), so it exited
`empty-interface` without fitting; this is a classifier-resolution limitation,
not fast optimization convergence. None of these runs is pooled with the live campaign.

Before separating the coefficients, the canonical executables were compared with the frozen full-baseline-affine
MUMPS runs, using the original explicit 20 outer / 10 inner caps and identical
stopping settings (not the production 30/15 defaults).
The table is historical evidence for the old default balance (1,1e-4,1e-4),
not for the new all-one defaults.

| Case | Accepted outer steps | Inner corrections | Final sampled D_inf | Maximum Q |
| --- | ---: | ---: | ---: | ---: |
| 2D n=50, lobes=10 | 5 | 6 | 7.46403537e-4 | 4.69870469 |
| 3D n=12, lobes=4 | 17 | 29 | 5.00368357e-3 | 9.99997640 |
| 3D n=20, lobes=4 | 12 | 30 | 2.29516328e-3 | 9.99988930 |

All accepted-state sampled RMS/sup distances, j/Q, fitting energies and inner
counts agree within 7e-13 absolute error; outer counts and exit reasons match.
The 2D case hits its historical empirical target. The 3D runs have no target:
n=12 stops on energy stagnation and n=20 on line-search failure near the quality
budget. This cleanup does not remove those limitations or establish a 3D error
order. Verification logs are in /tmp/wngir-canonical-verification/.

Verification also passed 40 scoped C++ unit tests, six manufactured assembly
checks and six campaign/trace tests. P1/P2 reconstruction targets build in 2D
and 3D; one-step P2 smoke solves accepted valid updates in both dimensions.
The 2D cantilever builds. The 3D cantilever target is unavailable in this build
because RODIN_USE_PETSC is off, so that integration was not compiled here.
Inner totals in the table include work from a final rejected outer attempt.

The independent F/S/D change passed 21 metric/solver unit tests, six manufactured
assembly checks and seven campaign/trace tests. Both P1 reconstruction targets
build. Mapping the old weights to (1,1e-4,1e-4) reproduces the seven logged states
of the 2D n=50, lobes=10 reference within 3e-14 on fit, geometry, quality and inner
counts. A one-step 3D n=8, lobes=4 smoke solve with the new all-one defaults
accepted a valid update (minimum j=1.196, maximum Q=1.004); it is not a calibration
or convergence result. The removed shared-bulk CLI flag is explicitly rejected.
