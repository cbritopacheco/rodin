# Canonical fitting-distribution mini study

The production metric is now

\[
M_u(v,z)=\frac{\kappa_F}{G_*^2}\int_{\Gamma_0}
(\nabla\phi(T_u)\cdot v)(\nabla\phi(T_u)\cdot z)
+h_0\kappa_D\int_{\Omega_0}j_u
\operatorname{dev}\operatorname{sym}(\nabla v A_u^{-1}):
\operatorname{dev}\operatorname{sym}(\nabla z A_u^{-1}).
\]

No shape or constraint Hessian is assembled. The robust Welsch objective,
pre-hinge directional Newton, affine quadratic hinges, stationarity residual,
similarity gauges and outer Armijo/actual geometry checks are unchanged.
The hinge reaction modifies the increment, not necessarily a separate quality
restoration step. Actual nonlinear quality is enforced by acceptance checks.

## Protocol and preservation

The distribution comparison was stopped at the user's request with 5,738
completed records out of 16,000. Its frozen binaries and records are retained
in `tmp/wngir-distribution-grid-20261003`; `stopped.json` records the interrupted
attempt. Do not append canonical results to this historical dataset.

The new study is in `tmp/wngir-fd-mini-20261003`, with 27 completed attempts,
binary hashes, exact commands and full traces. The nine geometries are 2D
resolution 10 with lobes 0, 2, 4, 6, 8; 3D resolution 10 with lobes 0, 2;
and 2D/3D resolution 20 with lobes 0. Controls are

\[
(\kappa_F,\kappa_D,\widehat\mu)
\in\{(10,0.01,0.1),(10,0.1,0.1),(10,0.01,10)\}.
\]

The exact study executables are preserved in that directory's `bin/`
subdirectory before subsequent formatting-only example rebuilds. Their
hashes match the original manifest; the commands retain their original
build-tree paths.

Each result has an exact geometry/fitting/distribution/hinge match in the
stopped elementwise campaign, with historical shape weight 0.01. Both use
30 outer / 15 inner caps, sampled target equal to background spacing squared,
quality ceiling 10, Jacobian floor 0.01, four MUMPS/OpenMP/OpenBLAS threads,
inner relative tolerance 0.001, and true linear residual tolerance 1e-9.
There are no warmups or repeats. Historical and new CPU endpoint flags are
retained; contaminated timings are not evidence of speedup. Endpoint sampling
cannot exclude intermediate interference. This is a numerical screen, not
a clean performance benchmark or a certified geometric-order study.

## Matched results

Outer means and maxima are among target hits; inner totals include all nine
attempts. Every configuration hits eight targets, with no quality violations.

| Weights | Metric | Hits | Mean outer | Maximum outer | Total inner |
| --- | --- | ---: | ---: | ---: | ---: |
| \((10,0.01,0.1)\) | Historical fitting + shape + distribution | 8/9 | 2.000 | 3 | 2 |
| \((10,0.01,0.1)\) | Canonical fitting + distribution | 8/9 | 2.375 | 5 | 5 |
| \((10,0.1,0.1)\) | Historical fitting + shape + distribution | 8/9 | 3.000 | 7 | 0 |
| \((10,0.1,0.1)\) | Canonical fitting + distribution | 8/9 | 2.625 | 6 | 0 |
| \((10,0.01,10)\) | Historical fitting + shape + distribution | 8/9 | 1.875 | 3 | 3 |
| \((10,0.01,10)\) | Canonical fitting + distribution | 8/9 | 2.250 | 5 | 8 |

The unresolved geometry is the coarse eight-lobe 2D target in every setting.
It stalls with distortion approximately 1.555, well below the budget. Shape
removal therefore does not resolve that obstruction. No approximation-limit
claim follows from this result.

At distribution weight 0.01 and hinge weight 0.1, the 3D resolution-20 sphere
takes five rather than three outer iterations. Final sampled geometric error
is 0.002054, maximum distortion 9.362 and minimum Jacobian 0.0207; this is a
valid target hit, not a quality failure. Increasing distribution weight to
0.1 reaches the same target in two iterations, with distortion 6.412 and
minimum Jacobian 0.0393, but doubles the outer count on the coarse four- and
six-lobe 2D targets from three to six. This is geometry-dependent redistribution,
not evidence for a universally superior coefficient.

## Conclusions

Removing shape preserves observed target coverage and removes its assembly and
local eigensolves. It does not establish a universal reduction in outer counts.
The two remaining metric weights and hinge strength should be calibrated by
target hits, worst and aggregate first-hit counts, then uncontaminated cost.
Quality budget consumption is allowed, not an additional minimization target.
At the time of this mini study, production weights were one and hinge strength was 90; the mini study
does not establish new universal defaults. Pointwise distribution is PSD but
not a continuous H1 norm modulo only similarities, particularly in 2D and
higher-order spaces. The paper's old global-dilation coercivity proof has been
replaced by a fixed finite-dimensional quotient statement.
