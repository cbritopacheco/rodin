# PETSc coverage-CI timeout investigation — 2026-10-24

The [measurement data](petsc-ci-timings.json) retain CPU/elapsed times,
PETSc event counts and the temporary test-stage timings used below.

## Scope and comparison

The reference is PR #368 at `3f38d00f560741507d8b92e978e646480bad133a`,
[Tests run 40028067601](https://github.com/cbritopacheco/rodin/actions/runs/40028067601).
Runner CTest elapsed times include process startup and coverage shutdown.
A watchdog termination is a lower bound, not a completed runtime or a numerical
failure. No mesh, degree, quadrature order, solver tolerance, assertion,
configuration, MPI rank or production source is changed by these fixes.

Local comparisons use Ubuntu 24.04/aarch64 in Docker on an Apple M1 Pro,
GCC 14.2, PETSc 3.19.6, `RelWithDebInfo` and atomic GCC coverage, with
`RODIN_MULTITHREADED=OFF`, `OMP_NUM_THREADS=2`, `OPENBLAS_NUM_THREADS=1` and
`OMP_WAIT_POLICY=PASSIVE`. The Linux container has an eight-CPU quota and a
6 GiB memory limit. Runner jobs use Ubuntu 24.04/x86-64 and the corresponding
GCC/PETSc build. Architecture, host load and runner placement remain confounders;
a wall-time ratio for one test is not a universal machine-speed constant.
These are coverage-build timings. No uninstrumented Release comparison was
performed, so they neither quantify coverage overhead nor establish a
production-code performance regression.

The full P test translation unit exceeded the local 6 GiB compilation cap.
Temporary benchmark-source copies therefore instantiate only local nonlinear
contexts and distributed Stokes contexts. Their measured test bodies, geometry
sets and numerical parameters are unchanged; the repository retains all
registrations. This compilation-scope difference is another comparison limit.

Only one local test process/context is executed at a time. Compilation is
completed before timing. Temporary copies of compiler inputs and coverage-note
files use RAM to avoid Docker shared-filesystem overhead; optimization and
coverage instrumentation are unchanged. Temporary Stokes test-helper timers
separate setup, assembly plus solve, residual validation and error norms. The
timers do not change any numerical operation, and are not part of the patch.
PETSc `-log_view` provides nested solver/assembly-event timings.

## Runner observations

The following are complete **local, multithreading-OFF** CTest groups unless
explicitly marked as timed out. P/HP nonlinear groups contain four independent
GoogleTest cases, not one linear solve.

| Study | Tetrahedron | Pyramid | Hexahedron | Wedge |
| --- | ---: | ---: | ---: | ---: |
| P nonlinear Poisson | 766.05 s | >1800.07 s | 121.27 s | 201.61 s |
| HP nonlinear Poisson | 768.56 s | >1800.10 s | 140.35 s | 216.46 s |
| Stokes boundary P2/P1 rates | 188.47 s | 644.99 s | — | — |
| Stokes boundary P3/P2 rates | 641.51 s | >1800.11 s | — | — |

The timed-out P Pyramid group had already completed its degree history in
544.600 s and its tangent checks in 457.008 s. The HP counterparts completed
in 875.967 s and 112.277 s. Those two unthreaded groups were killed during the amplitude-four
incorrect-reaction control; the final independent sensitivity case was never
reached. Raising a timeout does not establish that the remaining cases pass.

The geometry effect is substantial but is **not confined to Pyramid**:

- Threaded Stokes P3/P2 rate histories on Tetrahedron and Hexahedron also exceed
  1800 s. The threaded Tetrahedron P2/P1 history takes 474.67 s, versus 188.47 s
  without threading on a different runner. The P nonlinear Pyramid degree
  history similarly takes 1254.004 s with threading versus 544.600 s without.
  This is an observation, not an
  isolated OpenMP-overhead measurement.
- Curved ReactionDiffusion boundary groups on Tetrahedron and Hexahedron exceed
  1800 s. MPI-2 Tetrahedron completes in 1755.47 s; other contexts hit the limit.
- The H nonlinear Pyramid P1 history alone completes in 455.938 s; the complete
  six-case group is killed by its old 600 s allowance.
- Runner times vary substantially between contexts. Stokes Pyramid P3/P2
  finishes in 1470.34 s with MPI-1, while local and MPI-2/3/4 contexts hit 1800 s.
  These data do not justify multiplying every timeout by one fixed factor.

Representative job logs: [nonlinear P/HP](https://github.com/cbritopacheco/rodin/actions/runs/40028067601/job/120682357443),
[Stokes Pyramid](https://github.com/cbritopacheco/rodin/actions/runs/40028067601/job/120684440311),
[curved ReactionDiffusion](https://github.com/cbritopacheco/rodin/actions/runs/40028067601/job/120687780624).

## Matched local measurements

Elapsed ratios below are runner/local, with the same original test bodies and
parameters. Each row is one serial run; these are workload-specific observations,
not a statistically established hardware multiplier. Local CPU columns sum
user and system time; PETSc profiling and temporary test-stage timers are enabled.

| Study and context | Local elapsed | Local CPU | Runner elapsed | Runner/local |
| --- | ---: | ---: | ---: | ---: |
| P nonlinear Hexahedron, all four cases, local | 41.42 s | 40.02 s | 121.27 s | 2.93 |
| P nonlinear Tetrahedron, all four cases, local | 320.48 s | 319.35 s | 766.05 s | 2.39 |
| P nonlinear Pyramid, degree history only, local (body) | 154.80 s | 154.77 s | 544.60 s | 3.52 |
| Stokes boundary P3/P2 Pyramid, complete history, MPI-1 | 886.40 s | 885.83 s | 1470.34 s | 1.66 |

The Pyramid row compares GoogleTest body times on both machines; its local
process elapsed time is 155.38 s. Its CPU value covers the full process.
The other rows compare complete CTest-group workloads including process
lifecycle. The full four-case Pyramid group has not been completed locally.

For the Hexahedron group PETSc reports 40.50 s total, 20.91 s in `SNESSolve`,
and only 0.017 s in `KSPSolve` across 31 linear solves. Function and Jacobian
evaluation events take 7.71 s and 15.64 s respectively. The latter include
explicit tangent checks outside `SNESSolve`; nested event times must not be
added as independent phases. Independent test timers account separately for
initial assembly, residual validation and field-error integration. The linear
solver is not the dominant cost in this measured workload.

For Tetrahedron the total is 319.8 s, with 184.66 s in `SNESSolve`,
148.70 s in function evaluation and 60.46 s in Jacobian evaluation. KSP
again takes only 0.018 s. The local error-norm checks are comparatively small: a
degree-four order-16 solve spends about 1.64 s on norms, versus 22.54 s in
Newton and 4.00 s in the independent residual assembly. The same sample has
729 global unknowns on both Tetrahedron and Hexahedron; the measured difference
is not explained by global linear-system size alone. These events isolate the
assembly/evaluation path but do not identify its individual C++ hotspots.

The Pyramid history has 154.8 s of PETSc elapsed time: 100.55 s in Newton,
including 71.61 s in function evaluation and 28.93 s in Jacobian evaluation.
KSP takes 0.008 s. The test timers measure approximately 29.0 s of initial
assembly, 17.7 s of independent residual assembly and 7.5 s of error norms;
about 95% of this history is assembly/Newton evaluation rather than norms or
linear solves. All four degrees converge in three Newton iterations. Degree
four alone accounts for about 77% of the history and has 1241 global unknowns,
versus 729 for the Tetrahedron/Hexahedron samples. This is observed high-order
assembly cost, not evidence of nonlinear iteration divergence.

The Stokes MPI-1 history passes all three grids in 886.40 s (885.83 s CPU).
Its finest grid has 3072 cells, 88347 velocity unknowns and 9009 pressure
unknowns, for a mixed dimension of 97356. The finest-grid assembly-plus-solve
takes 613.46 s and its norms take 170.23 s; that grid accounts for about 89% of
the complete history. Peak child RSS is approximately 2.12 GiB, below the
local memory cap.

Across the Stokes history, setup/factorization (`PCSetUp`) takes 40.54 s and
`KSPSolve` 0.58 s; numerical LU (39.75 s) and symbolic LU (0.79 s) are nested
inside PC setup and must not be counted again. The test's assembly-plus-solve
timers sum to 690.18 s. Subtracting PC setup and solve leaves approximately
649.1 s of assembly, or 73% of the history. Norms sum to 195.05 s, about 22%;
linear setup/solve contributes about 5%. This is principally assembly and
independent integration cost, not MUMPS factorization alone.

The measured runner/local factors span 1.66–3.52. Local CPU time stays close
to elapsed time in these runs. The watchdog changes use the actual CI
terminations and completed CI cases rather than a blanket hardware multiplier.
Longer CI allowances are feasible for the identified CPU-intensive groups;
aggregate selections additionally require the execution splits described below.

## Work performed by the tests

The source establishes several reasons for large work counts independently of
machine speed. The nonlinear degree studies assemble with quadrature order 16
and integrate errors at order 18; sensitivity repeats the solve/integration at
orders 18/20. Every sample constructs a fresh problem, and Newton callbacks
reassemble nonlinear residuals and tangents. Incorrect-reaction and tangent
controls add solves/evaluations beyond the convergence histories.

The current quadrature dispatch uses a three-direction Gauss fallback for
Pyramid above its tabulated degree 10, and Tetrahedron above degree 15.
Consequently order 16 uses 729 volume points per cell, order 18 uses 1000, and
order 20 uses 1331. These exactness choices are retained. Hexahedron also uses
729 points at order 16, but its uniform grid has one element per grid cube;
Pyramid and Tetrahedron each have six. This alone makes equal grid dimensions
an unequal element-count comparison against Hexahedron.

For a scalar degree-three element, Pyramid has 30 basis functions and
Tetrahedron 20; at degree four the counts are 55 and 35. Dense local pair work
therefore grows by factors 2.25 and approximately 2.47 before accounting for
basis-evaluation and sparse-solve costs. These are work-count comparisons,
not measured runtime explanations by themselves.

The Stokes P3/P2 test retains the complete three-grid Taylor–Hood history at
3, 5 and 9 grid points per axis. Its finest Pyramid grid has 3072 cells. Each
cell has 90 local velocity components plus 14 pressure functions, compared
with 60 plus 10 on Tetrahedron. The finest order-16 integration visits
2,239,488 volume points per form; several matrix blocks and independent field,
divergence and pressure-mean checks are assembled/integrated. MUMPS is a direct
factorization, so global DOF count and matrix fill also matter. Profiling is
needed to distinguish those costs from Rodin quadrature and assembly.

## Scheduling changes and safeguards

Individual scheduling allowances are increased only in the identified families:

| CTest group | Previous | Updated |
| --- | ---: | ---: |
| P/HP nonlinear Pyramid, complete four cases | 1800 s | 10800 s |
| H nonlinear Pyramid, complete six cases | 600 s | 3600 s |
| Stokes boundary P3/P2 Pyramid | 1800 s | 7200 s |
| Stokes boundary P3/P2 other 3D geometries | 1800 s | 3600 s |
| Curved ReactionDiffusion boundary Tetrahedron/Hexahedron | 1800 s | 7200 s |
| Curved nonlinear Pyramid, each complete rate history | 1800 s | 10800 s |

These values are conservative watchdog policies, not demonstrated upper bounds.
Other process budgets, including the existing four-hour curved-reaction
Pyramid allowance, remain unchanged. The numerical stopping criteria are
independent of these wall-clock budgets.

Job-level timeouts require a separate correction. The real PETSc bulk lane is
killed at two hours after 100 of 145 registrations; the complex curved-boundary
lane is killed at four hours after 188 of 245 registrations, with tests still
passing. Both are partitioned into local and MPI-1/2/3/4 contexts. The curved
nonlinear Pyramid lane is partitioned by context and individual named rate
history, plus the complementary Controls group; each history keeps every mesh
level. The real bulk job allowance is six hours, and curved ReactionDiffusion
lanes receive six hours rather than four. No additional compilation is created.

There are 114 additional execution partitions, for 1354 requests over the same
46 compatible configurations. Registration tests execute the real Bash
selectors against configured test names and require exact, nonoverlapping
coverage of every original test/rank. The configured-manifest replay confirms the same 698 hierarchy build producers
and executable target sets after the split; the expected campaign has 2180
jobs, including 64 combined publishers. The generated workflows retain shared
build dependencies, fresh per-lane counters, mandatory report receipts, combined
Codecov publication, the compilation baseline and the final 95% gate. Complete
coverage remains pending until the replacement CI campaign finishes.

## Associated CI failures and validation

The curved-elasticity selector previously expected the geometry immediately
following the local/rank prefix. `RATE_CASES` instead registers names such as
`MixedTractionP3Rates_Tetrahedron` and `Controls_Tetrahedron`. The selector now
accepts the explicit case component. It keeps local and MPI selections distinct;
configured-name regression checks cover all curved physics families, including
Controls, exactly once per threading configuration.

APT fallback previously replaced direct Ubuntu URLs but left the runner's
`mirror+file:/etc/apt/apt-mirrors.txt` URI unchanged, so all fallback attempts
continued using the failing Azure endpoint. Ubuntu-only mirror-list references
are now resolved to the selected fallback archive. The pinned [runner image source](https://github.com/actions/runner-images/blob/ubuntu24/20261012.334/images/ubuntu/scripts/build/configure-apt-sources.sh)
confirms the three official Ubuntu endpoints and priority annotations handled
by this detection. Source suites, components, signing keys, exact package
versions and third-party definitions are retained.
The new regression reproduces the failure before the change and checks all
fallback attempts afterwards.

The serial Linux infrastructure campaign passes 86 checks, followed by the
additional actual-registration selector check (87 unique checks in total). The
four numerical invocations above also pass, comprising ten GoogleTest cases.
The complete production test matrix remains the responsibility of CI. Workflow regeneration and
`actionlint` pass; the pre-existing `peaceiris/actions-gh-pages@v3` legacy-runtime
diagnostic is excluded from the latter check. The configured-manifest comparison
preserves executable sets and compile flags. Other workers' resource files are
not included in the patch.

For a representative reproduction, configure the same dedicated build recipe
from `ci_builds.definition("PETScDegreeConvergence", {"multithreaded": "OFF",
"suite": "scalar_mixed", "physics": "NonlinearPoisson"})` and build the P/HP
nonlinear targets plus `RodinConvergenceHPETScStokesBoundary`. Set `CC=gcc-14`
and `CXX=g++-14`, use PETSc 3.19.6, and complete compilation before timing.
Execute one CTest registration with `ctest --test-dir BUILD/tests -R '^NAME$'
--verbose --no-tests=error -j 1`. Preserve the original process timeout for a
baseline, or record explicitly that a terminated run supplies only a lower
bound. Run longer diagnostic studies directly with the same GoogleTest filter;
do not change the numerical parameters. PETSc `-log_view :FILE` records nested
events. Python `resource.getrusage(resource.RUSAGE_CHILDREN)` before/after the
single subprocess supplies user/system CPU time separately from elapsed time.
