# Adaptive affine-hinge calibration

Status: proposed, not launched. Historical datasets and stopped experiments
must remain separate from results of the canonical adaptive implementation.

## Model

Use the robust Welsch energy, quadratic fitting tensor, globally centered
distribution, directionally scaled predictor, adaptive affine quadratic
hinges and outer Armijo acceptance with actual sampled quality checks.
There are no alternative hinge models in the production grid.

The adaptive measures for Jacobian and distortion are independent. Each
preserves the mapped cell mass, mixing equal weights and normalized nonlinear
guard penetration at the current and full-predictor geometries. The mixing
fraction is fixed at one half; risks are capped at 100. Freeze both measures
throughout the inner solve. Do not calibrate these policy constants yet.

## Proposed Grid

| Setting | Proposed values |
| --- | --- |
| Fitting stiffness \(\kappa_F\) | \(1\), fixed reference scale |
| Deviatoric variation \(\kappa_{\rm dev}\) | \(10^{-5},10^{-4},10^{-3},10^{-2}\) |
| Divergence variation \(\kappa_{\rm div}\) | \(0,10^{-4},10^{-3},10^{-2}\) |
| Hinge strength \(\widehat\mu\) | \(0.1,1,10,100,1000\) |
| Spaces and matching geometry degree | \(P_1,P_2\) |
| 2D grid points per coordinate \(n\) | \(8,16,32\) |
| 2D lobes \(\ell\) | \(0,4,6,8\) |
| 3D grid points per coordinate \(n\) | \(4,6,8\) |
| 3D lobes \(\ell\) | \(0,4,6\) |
| Outer / inner correction budgets | \(30/15\) |
| Inner relative / absolute residual | \(10^{-3}/10^{-12}\) |
| True relative linear residual | \(10^{-6}\) |
| Relative Jacobian floor / distortion budget | \(0.01/10\) |
| Guard / Jacobian row weight / distortion row weight | \(0.1/1/1\) |
| Geometric target | \(h_0^{p+1}\), with fixed background \(h_0=1/(n-1)\) |
| Geometry scale / amplitude | \(R_0=0.25,\ a=0.05\) |
| Center / phase | \((0.5,\ldots,0.5)/0\) |
| Predictor motion cap | Unrestricted |
| Surface / volume / quality / geometry orders, affine \(P_1\) | \(8/2/2/32\) |
| Surface / volume / quality / geometry orders, \(P_2\) | \(12/8/16/32\) |
| Solver / execution | MUMPS; one heavy process at a time |
| Warmups / repeats | None |

This is 80 parameter sets. Matched P1/P2 cases give 1,920 attempts in 2D
and 1,440 in 3D, hence 3,360 in total. Fixing fitting stiffness retains the
current preferred normalization and focuses the grid on distribution-to-fit
ratios and hinge recovery, rather than multiplying the screening cost by
another common stiffness scale.

Run 2D first, then 3D. Within each phase increase resolution, then lobes;
test both degrees consecutively for each parameter set and geometry.
Use four OpenMP/OpenBLAS threads only when CPU availability permits; record
the actual thread environment and do not overlap other heavy work.
Reject a case before launch if its memory estimate exceeds the available
budget; record the omission explicitly. No silent substitution of coarser
meshes or reduced quadrature is permitted.

## Responses And Ranking

Log initial and every accepted outer state, inner residuals and merit steps,
the final state and exit reason. Retain exact commands, binary and source
hashes, quadrature settings, runtime, peak resident memory and CPU contention.
Record failed, timed-out and unstarted attempts separately.

Primary responses are sampled maximum distance and its dimensionless score

\[
C_p = D_\infty/h_0^{p+1},
\]

the first target-hit outer iteration and elapsed cost, Welsch energy, accepted
outer count, per-outer inner counts (median, maximum and total), sampled
maximum distortion and minimum Jacobian. Keep best and final geometry
distinct: energy descent need not decrease maximum distance monotonically.

Rank parameter sets first by target coverage within the iteration budgets,
then by worst-case and median final \(C_p\) over all matched cases, then by
iterations and uncontaminated cost. Report a Pareto comparison if these
criteria conflict. Quality-budget usage is allowed; lower distortion breaks
ties but is not the primary objective. Never omit cap/failure cases from
aggregate rankings. Compare P1 and P2 using both absolute distance and
degree-dependent scores; their target hit rates alone are not directly
comparable.

This screen does not establish an asymptotic geometric order. Validate the
selected sets separately on finer meshes and denser quality/geometric
sampling, including topologically mismatched domains, before claiming
transferability. No such follow-up is automatically launched.
