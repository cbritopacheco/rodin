# Adaptive affine-hinge calibration

Implementation update after this relaunch: production now uses a cached
closed-form reference covering with the same subdivision-based point counts.
Future launches through the checkout runner record that policy. The active
campaign below remains frozen to its uniform lattice; do not pool its results
with closed-form covering runs or relabel its witnesses.

Relaunch configuration dated 2026-10-10: canonical uniform witnesses.
Replacement directory: `tmp/swift-barycentric-campaign-relaunch-20261010`.
Stopped replacement directory: `tmp/swift-barycentric-free-campaign-20261010`.
Read the replacement's `progress.json` for execution status and `launch.json` for the runner
PID. Campaign solves are serial, but do not wait for other agents' work.
External CPU activity is recorded for timing interpretation, not used as a gate.
Historical datasets and stopped experiments
must remain separate from results of the canonical adaptive implementation.

The current checkout uses hierarchical parameter flags, for example
`--model-fit`, `--convergence-iterations-outer` and `--sampling-subdivision`.
The active relaunch retains its frozen executable and runner with the previous
spellings. Their numerical settings are unchanged; do not replace either file
mid-campaign. Future launches use the renamed flags.

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
| Deviatoric variation \(\kappa_{\rm dev}\) | \(0,10^{-5},10^{-4},10^{-3},10^{-2},0.1,1\) |
| Divergence variation \(\kappa_{\rm div}\) | \(0,10^{-5},10^{-4},10^{-3},10^{-2},0.1,1\) |
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
| Exterior boundary | Free; no Dirichlet condition |
| Surface / volume integration orders, affine \(P_1\) | \(8/2\) |
| Surface / volume integration orders, \(P_2\) | \(12/8\) |
| Uniform quality-lattice subdivisions, \(P_1/P_2\) | \(2/16\) |
| Triangle witnesses, \(P_1/P_2\) | \(6/153\) |
| Tetrahedron witnesses, \(P_1/P_2\) | \(10/969\) |
| Geometric sampling order | \(32\) plus facet vertices |
| Solver / execution | MUMPS; one heavy process at a time |
| Warmups / repeats | None |

The distribution grids are identical but varied independently, not restricted
to equal coefficient pairs. Zero includes distribution ablations; the
double-zero case can be singular and must be recorded even if its solve fails.

With unit hinge row weights this is 245 parameter sets. Matched P1/P2 cases
give 5,880 attempts in 2D and 4,410 in 3D, hence 10,290 in total.
Fixing fitting stiffness retains the
current preferred normalization and focuses the grid on distribution-to-fit
ratios and hinge recovery, rather than multiplying the screening cost by
another common stiffness scale.

## Proposed Hinge-Weight Extension

The model depends on the products \(\widehat\mu\kappa_j\) and
\(\widehat\mu\kappa_Q\), not three independent hinge strengths: the predictor,
adaptive risks and their normalization do not depend on these row weights.
Multiplying both row weights by the same positive factor and dividing
\(\widehat\mu\) by that factor leaves the frozen inner model unchanged.
Consequently a full independent sweep wastes common-scale comparisons.

The approved extension fixes the distortion row weight
and varies the Jacobian-to-distortion ratio:

| Setting | Proposed values |
| --- | --- |
| \(\kappa_j\) | \(0.1,1,10\) |
| \(\kappa_Q\) | \(1\), normalization |
| \(\widehat\mu\) | \(0.1,1,10,100,1000\), unchanged |

This yields 735 parameter sets and 30,870 attempts over the same geometry
grid: 17,640 in 2D and 13,230 in 3D. An independent three-value sweep for both
row weights would instead require 92,610 attempts, with overlapping effective
strength pairs. Both choices retain the actual sampled quality bounds; they
change penalty balance, not the prescribed quality budget.

## Execution

The fixed-wall campaign in `tmp/swift-adaptive-campaign-20261010` was stopped
at the user's request after 9,395 recorded attempts. Its frozen driver imposed
zero exterior displacement. Those results remain historical fixed-wall data;
they must not be pooled with a free-boundary campaign or resumed to evaluate
the modified driver. A new campaign requires a fresh build and manifest.
The free-boundary barycentric replacement was also stopped at the user's
request after 2,077 recorded attempts, including the interrupted final attempt;
1,622 target hits were recorded. Its frozen sources and binaries remain
unchanged. Production quality sampling now uses only uniform reference
lattices, with adaptive weights unchanged and no exterior Dirichlet condition.
The table above describes the current production policy. The relaunch starts
fresh and does not pool records from either stopped campaign. The completed
64-run Lobatto/covering comparison remains a separate experimental dataset;
its unrestricted covering point set is not part of this calibration.

Run 2D first, then 3D. Within each phase increase resolution, then lobes;
test both degrees consecutively for each parameter set and geometry.
Use four OpenMP/OpenBLAS threads and record the actual thread environment.
Only one campaign solve runs at a time; unrelated work may run concurrently.
The frozen runner uses four threads, a single campaign-wide lock and one
child solver process at a time. Each case has a 30-minute wall-clock timeout
and an 8 GiB live process-tree RSS limit, polled every two seconds. Native
`time -l` also records the completed process's peak RSS; polling can miss a
short-lived memory peak before it exits. Resource-limited cases are retained
as failures, not silently omitted or replaced by smaller cases.
External CPU detection is sampled and recognizes scientific executables and
compilers; it cannot detect every unrelated CPU consumer. Contended timings
are flagged, not grounds for delaying a solve. Controller restarts append only
missing attempts after checking that CSV and JSONL attempt sets agree.
Full traces are compressed per attempt; transient mesh outputs are removed
after parsing to bound disk use. No silent substitution of coarser meshes or
reduced quadrature is permitted.

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
