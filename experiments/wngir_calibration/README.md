# Canonical WNGIR Campaigns

WNGIR now has one model: M = F + S + D (Fitting, Shape, Distribution), affine quadratic hinges,
directional Newton, frozen inner-merit backtracking, and actual outer j/Q
and fitting-energy checks. There is no selectable logarithmic model, PSD
clipping by default, coefficient-space completion, mass term, or nonlinear-hinge path.
The experimental `--wngir-positive-shape-curvature=1` projects the local
Hessian of (d/4)(Q-1) onto its PSD spectrum before assembly. It adds no force
and does not guarantee invertibility of the total metric.

- F is kappa_f times the normalized Hessian of half the squared level-set residual with
  the level-set Hessian omitted. It is not robust-weighted.
- S is h*kappa_s times D2[(d/4)(Q-1)], without clipping.
- D is h*kappa_d times current-configuration symmetric strain,
  with global uniform dilation projected out. It permits rigid motions.
- The fitting energy and force remain robust Welsch.
- Hinge widths are 0.1 times the identity quality margins. The default
  model-decrease-scaled hinge weight is mu_hat=90.
- Defaults are 30 outer / 15 inner corrections and kappa_f=kappa_s=kappa_d=1.
  There is no shared bulk coefficient. These are choices, not a new calibration result.
  The previous default metric is reproduced by (kappa_f,kappa_s,kappa_d)=(1,1e-4,1e-4).
  More generally, old weights map to (kappa_obs,kappa_bulk*kappa_c,kappa_bulk*kappa_reg).

There is no separate inertia audit or automatic metric repair. Linear residual,
direction, inner merit and actual outer quality/energy checks remain in place.
Fitting must resolve the similarity modes that D leaves free;
for example, a planar interface cannot identify tangential translation.
Full shape curvature is not globally PSD. No general coercivity claim is made.
The distribution integrators tabulate frozen current strain once per basis and
quadrature point, rather than reevaluating the inverse deformation for each
basis pair. MUMPS retains symbolic analysis while the sparsity pattern matches,
and numeric factors while entries match exactly. An entirely inactive hinge
model uses the predictor directly without assembling and solving a zero correction.
Geometry traces record inactive-hinge skips and analysis/factorization counts.

## Remaining Parameters

| Role | Parameters | Defaults |
| --- | --- | --- |
| Metric | kappa_f, kappa_s, kappa_d | 1, 1, 1 |
| Hinge strength | mu_hat; kappa_j, kappa_q | 90; 1, 1 |
| Guard width | quality_guard | 0.1 of identity margins |
| Quality budget | j_safe, q_max; actual line-search Jacobian floor | 0.01, 10; 0.01 |
| Robust fitting | robust_scale | automatic |
| Iteration limits | outer, inner; CG per linear solve | 30, 15; 1000 |
| Inner stopping | relative Newton correction | 1e-3 |
| Geometry target | geometric_sup_tolerance | 0 (explicit target disabled) |

Line-search, stagnation, legacy active-fit tolerances, diagnostics and backend
controls remain unchanged. The primary calibration controls are the three metric
weights and mu_hat; the quality budget is prescribed rather than inferred from fit.

## Model

For A_k = I + grad u_k, let e_k(v) = sym(grad v A_k^-1), j_k = det A_k,
V_k = integral j_k and g_k = grad phi(x+u_k). With the frozen gradient-scale
normalization N, the metric is

```text
F[v,z] = kappa_f N integral_interface (g_k.v)(g_k.z)
S[v,z] = h kappa_s integral_volume (d/4) D2Q(A_k)[grad v,grad z]
D[v,z] = h kappa_d [
  integral_volume j_k e_k(v):e_k(z)
  - (integral_volume j_k tr e_k(v))(integral_volume j_k tr e_k(z))/(d V_k)]
```

Only global dilation is removed, not each element's volumetric strain.
The negative rank-one term is implemented exactly by a sparse augmented
direct solve or the corresponding CG operator. It is part of D, not an
extra completion or penalty.

The frozen affine slacks are s_J = j_k-j_safe+Dj_k[grad v] and
s_Q = Q_max-Q_k-DQ_k[grad v]. With delta_J=guard*(1-j_safe),
delta_Q=guard*(Q_max-1), the inner objective is

```text
0.5 M[v,v] - f_k[v]
+ 0.5 mu_k integral_volume sum_{a=J,Q} kappa_a (1-s_a(v)/delta_a)_+^2
mu_k = mu_hat * (0.5 f_k[predictor]) / reference_domain_volume
f_k = -DE_W(u_k)
```

Negative affine slack has finite penalty; no fraction-to-boundary step is
used. Inner merit backtracking handles active-set changes. The nonlinear
outer update still requires actual Jacobian and distortion admissibility
and sufficient decrease of E_W. Metrics do not add a quality force to E_W.

The 2D runner uses MUMPS and independent --kappa-f/--kappa-s/--kappa-d grids,
each defaulting to 1e-4,1e-3,1e-2,0.1,1. With mu_hat=0.1,1,10,100,1000,
all five resolutions and eleven lobe counts give 34,375 cases per shape variant.
`--shape-curvature=full,psd` compares both (68,750 cases), interleaved for each
coefficient tuple. The schema, resume keys, logs and manifest distinguish the
variants; the stopped fixed-F campaign is preserved separately.
The 3D runner still uses fixed --kappa-f and --kappa-s/--kappa-d grids.
Executable defaults remain all one. Executables use
--wngir-kappa-f, --wngir-kappa-s and --wngir-kappa-d; the old metric flags
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
kappa_f, kappa_s and kappa_d and cannot resume older coefficient schemas.
Inner rows contain Newton correction and iterate norms, relative correction,
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

## Cleanup Regression

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
