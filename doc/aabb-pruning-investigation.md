# AABB projection pruning investigation

## Decision

Projection pruning is opt-in: `locator.setProjectionPruning(true)` enables it.
The default avoids additional index construction and storage. This follows the
new-capability policy in `doc/agents/testing.md` and is supported by the measured
tradeoff: large query gains on overlapping cells, but extra storage and setup
with little benefit on regular tensor grids. Repeated point-location workloads
on simplex, curved, or sheared meshes are good candidates for explicit enabling.

This is a synthetic workload study, not a claim about the distribution of
queries or memory budgets in an application. The correctness fixes and Newton
seed retries apply in both modes. Those retries can make unpruned curved queries
substantially slower than the earlier inversion that incorrectly stopped at an
exterior root. Keeping pruning off favors setup and storage; it is not a claim
that unpruned lookup became faster.

## Method

Release Clang 19 / MPICH, PETSc off, on the local macOS arm64 host. All eight
geometries were exercised: Point, Segment, Triangle, Quadrilateral, Tetrahedron,
Hexahedron, Pyramid, and Wedge. The matrix has 190 mesh configurations and
30,836 paired query samples (61,672 diagnostic query executions).

Nonpoint maps use `RealH1Element<K>`, K = 1, 2, 4. The integer-lattice uniform
mesh is normalized to a unit domain, then its nodal coordinates are mapped by
`x_last += a*x_first^K`, with a = 0, 1, 4. In dimension >= 2 this is an injective
triangular shear; the one-dimensional version is monotone. Thus a changes
orientation and curvature without introducing inverted cells. The K=1 cases
use the H1 representation, not the mesh's default `RealP1Element` implementation.
Absolute evaluation costs need not match those of the default element.

Resolution counts vertices per axis: 16/64/256 in 1D, 4/8/16 in 2D, and 3/5/8 in
3D. At most 64 deterministic spatially distributed queries per class are used:

- Interior: a reference point between the centroid and vertex zero, excluding
  the centroid itself so that a solve is required.
- Shared: barycentres of internal reference faces (shared vertices in 1D).
- Miss: points just above the upper domain face, displaced by 1e-4 divided by
  the resolution before applying the shear. Unwarped misses can be rejected at
  the root box; warped misses can enter the tree despite being outside the mesh.

Both modes must return identical membership and cell indices, with reference
coordinates within 1e-8, before timing. The same validation passed in the diagnostic
run. Boundary-cell selection is therefore checked, not merely membership.

Timings use the ordinary executable without counters. Seven interleaved query
groups cover the whole sample set; mode order alternates by group. Cheap groups
repeat for at least 0.01 seconds of elapsed time; expensive groups run at least
once. The reported query cost is the sample-weighted mean CPU time per query,
using `std::clock`, excluding fixture generation and index construction. Build
cost includes locator construction, the first query that builds the index, and
destruction; the process-wide basis cache is warm. Builds were paired and repeated
at least three times, with their median reported. The initial configurations used
seven construction repetitions; the final harness uses three to bound the cost of
large P4 pyramid builds. Existing `FirstP4` benchmarks cover cold conversion caches.

An initial elapsed-time pass was discarded: regular quadrilateral queries showed
large apparent differences despite identical operation counts. A CPU-time repeat
reduced that case to approximately 268 ns off and 257 ns on. Scheduling pauses
can contaminate short elapsed-time comparisons. Small percentage changes in the
remaining results are not treated as evidence of a speedup.

Counters are injected into a generated header for a separate diagnostic executable.
The production header and ordinary timing binary contain no counter increments.
Counts include work until the first accepted cell; candidates never visited after
that acceptance are not counted. Iterations count entries into the Newton loop,
including convergence checks that need no Jacobian. Retries count alternate seeds
after each candidate's centroid seed. Point location uses a distance check, so its
Newton-candidate count is zero.

Memory is retained vector capacity for nodes, entries, boxes, projections and
projection ranges. It excludes mesh storage, the shared basis cache, temporary
construction buffers, allocator metadata and the locator's small fixed fields;
it is not peak RSS.

## Query cost and work

Largest mesh of each type, K=4, a=4. Costs are CPU microseconds per interior
query. Counts are means per query. Differences between geometries also reflect
cell count, element evaluation and mesh decomposition; only on/off within a row
is a controlled comparison.

| Geometry | Cells | CPU off / on (us) | Newton candidates off / on | Index off / on (KiB) |
| --- | ---: | ---: | ---: | ---: |
| Triangle | 450 | 256.0 / 22.3 | 4.12 / 1.31 | 38.9 / 109.9 |
| Quadrilateral | 225 | 18.1 / 2.7 | 1.64 / 1.00 | 19.6 / 55.1 |
| Tetrahedron | 2058 | 9208.8 / 1398.3 | 8.73 / 1.78 | 177.0 / 721.2 |
| Hexahedron | 343 | 1156.8 / 144.7 | 1.92 / 1.08 | 29.6 / 163.0 |
| Pyramid | 2058 | 72607.1 / 4990.3 | 4.73 / 1.31 | 177.0 / 721.2 |
| Wedge | 686 | 2195.2 / 189.9 | 3.11 / 1.16 | 59.1 / 197.9 |

For the largest strongly sheared P4 tetrahedral mesh (2,058 cells), all requested
operation counts are shown below. Transforms and Jacobians count query work only.

| Query | Pruning | Candidates | Transforms | Jacobians | Iterations | Seed retries | CPU (us) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| interior | off | 8.73 | 445.38 | 178.97 | 178.97 | 61.88 | 9208.8 |
| interior | on | 1.78 | 103.41 | 24.94 | 24.94 | 6.25 | 1398.3 |
| shared | off | 9.73 | 746.23 | 243.50 | 243.50 | 69.88 | 12994.5 |
| shared | on | 1.58 | 19.34 | 12.56 | 12.56 | 4.62 | 623.2 |
| miss | off | 13.78 | 594.14 | 319.16 | 319.16 | 110.25 | 16150.0 |
| miss | on | 1.97 | 53.19 | 35.56 | 35.56 | 15.75 | 1744.8 |

Shared faces can be more expensive than interior hits without pruning because
several overlapping boxes are inverted before finding a containing cell. Nearby
misses remain expensive even with pruning: surviving candidates may exhaust all
seeds. Physical distance to the mesh is therefore a poor predictor of cost.

The slow tail matters. In the tetrahedral interior case, the 95th-percentile
Newton-candidate count falls from 22 to 3 (maximum 40 to 4). Transform calls fall
from 1,281 to 737 at the 95th percentile (maximum 4,714 to 1,068). Pruning reduces
the tail but does not eliminate expensive inversions on surviving candidates.

Point and Segment have no projection planes in either mode and their exact work
counts and index storage match. Regular tensor interior queries typically reach
one candidate already; pruning cannot remove that inversion. Regular tetrahedral
interior queries on the largest mesh fall from 3.70 candidates to 1.00.
Three additional 0.03-second repetitions of the regular H1 P1 hexahedral case
show no material query benefit: median interior costs are approximately 596 ns
off and 598 ns on, and shared-face costs are 610 ns off and 615 ns on.

## Construction and memory

On the largest regular H1 P1 tetrahedral fixture, construction rises from
0.713 ms to 2.400 ms while interior queries fall from 3.136 us to 0.926 us.
That is approximately 764 interior queries to amortize construction in this
particular fixture. On the strongly sheared P4 tetrahedral fixture, construction
rises from 200.596 ms to 291.008 ms while interior queries fall from 9.209 ms to
1.398 ms: approximately 12 such queries amortize construction. Neither count is
an application-wide threshold.

The P4 tetrahedral index grows from 181,296 to 738,512 bytes (4.07 times). The P4
hexahedral index grows from 30,344 to 166,904 bytes (5.50 times). These ratios are
for the locator's retained vectors, not the whole mesh or process.

The current exact-nonzero test for projection direction components can retain
nearly axis-aligned normals caused by floating-point noise. On the regular H1 P1
fixtures, quadrilateral index storage grows from 20,024 to 27,720 bytes, and
hexahedral storage grows from 30,344 to 52,216 bytes, despite unchanged inversion
counts. Filtering numerically redundant directions is a concrete follow-up; any
threshold should be a documented policy constant, and skipping such a filter
must preserve membership results. It is not included in this default decision.

## Growth with mesh size and distortion

For K=4, a=4 tetrahedra, mean interior Newton candidates are:

| Cells | Off | On |
| --- | ---: | ---: |
| 48 | 4.75 | 2.10 |
| 384 | 7.14 | 1.67 |
| 2058 | 8.73 | 1.78 |

Pruning keeps this candidate count near one or two despite increasing resolution.
Tree traversal still grows with the mesh. Samples change with resolution, and
nonlinear seed behavior is not monotone; this table is not a universal complexity
law.

At 2,058 cells and K=4, increasing shear gives:

| Shear a | Candidates off / on | Transform calls off / on |
| --- | ---: | ---: |
| 0 | 3.58 / 1.00 | 48.41 / 2.00 |
| 1 | 4.70 / 1.11 | 154.30 / 20.53 |
| 4 | 8.73 / 1.78 | 445.38 / 103.41 |

The primary matrix changes both representation order and the power of the shear.
It cannot alone attribute a time increase to element order. A supplementary
comparison holds the quadratic physical map fixed while representing it at P2
and P4; its results are recorded below.

For the same quadratic shear (a=1) on 2,058 tetrahedra, representation
order alone raises the evaluation cost. With pruning enabled, both interior
cases perform exactly one inversion, three transformations and two Jacobians
per query, with no seed retries:

| Representation | CPU per interior query, off / on (us) | Transformations on | Jacobians on |
| --- | ---: | ---: | ---: |
| P2 | 389.81 / 9.23 | 3 | 2 |
| P4 | 4286.83 / 98.39 | 3 | 2 |

The pruned P4 query is about 10.7 times more expensive than P2 despite identical
inversion work counts. The extra cost is in evaluating the higher-order
transformation/Jacobian representation, not in additional tree candidates or
Newton iterations. This is evidence for investigating geometry evaluation and
reusing evaluations next; it is not a reason to reduce the seed budget, which
protects the exterior-root regression. The unpruned candidate sequences are
slightly different because degree-dependent hull recovery and floating-point
iteration behavior differ even when the physical polynomial is the same.


## Geometry evaluation follow-up

The next comparison uses `4b6925135` as the evaluation baseline. A minimal
change to `ParametricTransformation::jacobian` evaluates each scalar reference
basis derivative once per basis/axis, then multiplies that value by every
physical coordinate coefficient. Previously the derivative call was inside the
physical-coordinate loop. Accumulation order and basis formulas are unchanged;
no cache, reduced seed budget, or changed convergence threshold is introduced.

### Profile and isolated benchmarks

Instruments Time Profiler sampled the original P4 tetrahedron Jacobian benchmark
for 12 seconds. Of 10,950 samples containing the evaluation benchmark, 10,677
(97.5%) contained the scalar basis derivative evaluator. The launched process
was stopped at the recording time limit; this trace identifies the hotspot and
is not used for timing claims. The call-count regression independently observed
three evaluations per basis/axis for a 3D physical map before the change, and
one afterward.

The 96 new `Parametric/Transform`, `Parametric/Jacobian`, and `ParametricNative`
Google Benchmark cases cover all eight geometries at P1, P2, P3 and P4, with both
native physical dimensions and three-dimensional embeddings.
The tables are warmed outside the timed loop. CPU times below are medians of
three alternating baseline/current runs, each with a 0.08-second minimum per
case, using the same final benchmark source and the original evaluation header
for the baseline executable. Each row is a single evaluation, rather than a point-location query.

| Geometry | P4 Jacobian before (us) | After (us) | Speedup |
| --- | ---: | ---: | ---: |
| Segment | 0.117 | 0.049 | 2.39x |
| Triangle | 4.406 | 1.612 | 2.73x |
| Quadrilateral | 1.409 | 0.558 | 2.53x |
| Tetrahedron | 44.366 | 15.516 | 2.86x |
| Hexahedron | 29.134 | 9.492 | 3.07x |
| Pyramid | 542.069 | 188.208 | 2.88x |
| Wedge | 27.662 | 9.907 | 2.79x |

The segment and surface results here include the benefit of embedding in 3D;
ordinary nonembedded cases have fewer physical components to share. P1 and P2
also improved for every nonpoint geometry in the initial three-repetition pass.
The initial chronological pass showed a systematic 1--4% increase even in
unchanged transform evaluations. Alternating executables removed that apparent
regression: P4 transform median ratios after/before ranged from 0.996 to 1.004.
This control demonstrates drift between the earlier measurement blocks.

Point has no reference axes or derivative work. Its isolated P4 zero-column
Jacobian call changed from about 4.32 to 4.55 ns in the alternating runs.
Disassembly shows a different generated prologue before the zero-axis exit,
including an additional move and stack store; the timing change is consistent
with changed compiler bookkeeping, rather than derivative evaluation. AABB point
queries do not call this Jacobian. This small cost is retained explicitly rather
than claiming that every benchmark improved.

The first derivative-hoisting implementation regressed the ordinary P2 segment
AABB hit by about 5% in alternating runs (277 to 291 ns). Disassembly exposed
compiler-generated setup for a general physical-component loop. Restricting the
loop by the existing `SpatialMatrix::MaxSize` makes its actual storage bound
visible, avoiding that unbounded loop machinery without adding another
algorithm. The repeated segment AABB measurement then gave 289 to 281 ns.
Native P4 segment Jacobian evaluation is unchanged (47.75 to 47.39 ns).
Native P1 segment Jacobian evaluation retains a small cost (11.14 to 11.87 ns),
consistent with the changed call bookkeeping seen in the point case; there is
no derivative-call reduction in a one-dimensional physical map. The native
surface P4 Jacobians improve from 2.917 to 1.536 us for triangles and from 0.973
to 0.550 us for quadrilaterals.

### End-to-end controlled workload

The same quadratic shear (`a=1`), P4 representation and 2,058 tetrahedra were
measured with the ordinary counter-free workload, using 0.03-second timing
blocks. Each query class contains 64 samples.

| Query | Off before / after (us) | On before / after (us) |
| --- | ---: | ---: |
| Interior | 4315.55 / 2093.34 | 98.86 / 44.09 |
| Shared boundary | 7726.17 / 3808.22 | 100.50 / 45.29 |
| Nearby miss | 6049.06 / 2533.12 | 955.66 / 404.52 |

The diagnostic executables produced 384 identical per-query records across
both pruning settings: candidates, transforms, Jacobians, Newton loop entries,
seed retries, box candidates, projection rejections and retained index bytes
all matched. The seven transformation tests and 51 AABB tests pass. New exact
scalar-baseline comparisons cover P1/P2/P3/P4 on every geometry, at the centroid
and every reference vertex (including the pyramid apex), with all supported
physical embeddings. The derivative-call regression was observed failing before
and passing after the change.

The existing curved P2 AABB hit benchmarks were also run before and after,
with three repetitions and a 0.06-second minimum (CPU medians):

| Geometry | Default off before / after (us) | Pruned before / after (us) |
| --- | ---: | ---: |
| Triangle | 7.304 / 4.497 | 1.540 / 1.030 |
| Quadrilateral | 2.780 / 2.121 | 0.776 / 0.620 |
| Tetrahedron | 249.841 / 91.864 | 9.714 / 3.647 |
| Hexahedron | 14.542 / 7.212 | 5.191 / 2.551 |
| Pyramid | 705.104 / 271.837 | 53.023 / 20.355 |
| Wedge | 49.198 / 22.511 | 6.384 / 2.987 |

These gains depend on the fixture's share of Jacobian work; the embedded 3D
microbenchmark speedups should not be applied directly to every locator query.
A nonembedded segment has only one physical coordinate and therefore no repeated
physical-component derivative work to remove.

The remaining large cost is repeated modal evaluation across nodal basis
functions and reference derivative axes. For example, each tetrahedron nodal
basis derivative traverses every Dubiner mode and computes all three modal
gradient components, then uses only the selected component. Reusing those modal
values and gradients within an element evaluation is a separate optimization
candidate. Its acceptance must preserve accumulation order, apex behavior,
generic finite-element compatibility, and Newton iteration records; the present
change does not introduce an alternative geometry evaluation path.


## Reproduction

Configure a Release Ninja build with `CMAKE_EXPORT_COMPILE_COMMANDS=ON` and tests
and benchmarks enabled. Paths below use `build` and a scratch output directory.
The optional fourth positional argument selects `degree/shear/resolution`; the
fifth fixes the shear power independently of the representation order.

```sh
cmake --build build --target RodinAABBWorkload RodinBenchmarks
build/tests/benchmarks/RodinAABBWorkload 0.01 > timing.csv
python3 dev/profile_aabb.py --build build --output /tmp/aabb-profile
/tmp/aabb-profile/RodinAABBWorkloadDiagnostic > counts.csv
build/tests/benchmarks/RodinAABBWorkload 0.03 Tetrahedron 2/1/8 2
build/tests/benchmarks/RodinAABBWorkload 0.03 Tetrahedron 4/1/8 2
/tmp/aabb-profile/RodinAABBWorkloadDiagnostic 0.03 Tetrahedron 2/1/8 2
/tmp/aabb-profile/RodinAABBWorkloadDiagnostic 0.03 Tetrahedron 4/1/8 2
build/tests/benchmarks/RodinBenchmarks --benchmark_filter=AABB
```

Google Benchmark has paired regular/curved hit and construction cases; names
containing `Pruned` or `WithProjections` explicitly enable pruning. Other cases
use the default. The standalone diagnostic copy deliberately fails generation
when a source anchor changes, preventing silently incomplete counters.

The existing AABB correctness regressions remain enabled in both modes. The
opt-in regression checks that default construction avoids the Jacobian used for
projection normals, explicit enabling performs it, and disabling restores the
original dependency-call behavior without changing the returned point.


To reproduce the isolated evaluation comparison, save the baseline benchmark
executable before changing production evaluation, then alternate it with the
current executable:

```sh
build/tests/benchmarks/RodinBenchmarks \
  --benchmark_filter='^Parametric' \
  --benchmark_min_time=0.10s --benchmark_repetitions=3
build/tests/benchmarks/RodinAABBWorkload 0.03 Tetrahedron 4/1/8 2
```


### P3 coverage extension

The AABB curved-enclosure regression already exercises P3 on every nonpoint
geometry, with both pruning settings. The later geometry-evaluation follow-up
extends exact scalar-baseline comparisons and derivative-call regressions to
every degree from P1 through P4. Isolated evaluation benchmarks now include P3
in native and embedded dimensions, and the paired workload accepts P3 fixtures.
The historical scaling matrix and before/after tables above retain their
original P1/P2/P4 sampling; they do not claim a measured P3 speedup. The current
full workload enumerates 253 mesh configurations instead of the original 190.

The expanded derivative-call regression was observed failing against the
original evaluator and passing with the fix, including P3. All seven
transformation tests and 51 AABB tests pass. All 24 P3 evaluation benchmark
cases executed with three repetitions, and the cubic tetrahedron paired
workload completed with matching pruning-on/off membership and references.
