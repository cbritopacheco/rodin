# Large interpolation workloads and cache-key audit

## Workloads and correctness

`GridFunctionLarge.cpp` adds 64 cases to `RodinGridFunctionBenchmarks`: scalar
and vector P1/P2, triangles/tetrahedra, four workloads, and cached quadrature
versus direct pointwise expansion. Connectivity and geometry are immutable
throughout each benchmark. The cache remains the existing single-entry,
per-thread cache; this experiment does not introduce a mesh-sized basis cache.

| Mesh | Grid vertices per direction | Cells | P1 scalar DOFs | P2 scalar DOFs | P1 samples/sweep | P2 samples/sweep |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Triangles | 129 | 32768 | 16641 | 66049 | 98304 | 196608 |
| Tetrahedra | 17 | 24576 | 4913 | 35937 | 98304 | 270336 |

Vectors have two components in 2D and three in 3D; their DOF counts scale
accordingly. Meshes, fields, spaces, formulas and mapped points persist across
repetitions. Point coordinates are warmed outside timing. All fields interpolate
the analytic linear function `1 + x + 2*y (+ 3*z)`; vector component `d` is
`(d+1)` times that function. Before timing, every sample/component is checked
on a quadrature miss, repeated hit and direct pointwise evaluation (absolute
error tolerance `1e-10`). Each timed sweep checks its accumulated checksum
(relative tolerance `1e-9`). Output and checksum enter `DoNotOptimize`.

| Workload | Evaluation order | Nominal basis-hit fraction |
| --- | --- | ---: |
| SlowSweep | Every quadrature sample of every cell once | 0% |
| Mixed75 | Each sample four consecutive times, then advance | 75% |
| FastBlocks | Each sample 32 consecutive times, then advance | 96.875% |
| PureHit | One prewarmed sample repeatedly, for the same evaluation count as FastBlocks | 100% |

The first warmed evaluation is excluded from these nominal fractions. SlowSweep
rebuilds the basis entry for each sample and the DOF entry when advancing cells.
Mixed75/FastBlocks traverse the complete large field. PureHit isolates the
hot path of that field; it does not simulate a traversal with all entries
already cached. Direct evaluation still uses the DOF cache, but calls the
space's expansion evaluator without the quadrature basis cache. All times
include output contraction and checksum accumulation. Neither initial
projection nor field re-interpolation is timed here; those operations remain
covered by [the interpolation/re-interpolation suite](GridFunctionInterpolation.md).

## Protocol and alternatives

Host: Apple M1 Pro, 10 cores; Clang 19.1.7; C++20, `-O3 -DNDEBUG`;
Google Benchmark `v1.9.4-34-g4b7b129e`, Release library. Single-threaded,
`RODIN_MULTITHREADED=OFF`, `RODIN_USE_PETSC=OFF`. The host was not dedicated;
main-run one-minute load averages were approximately 3.7–7.9. No compilation
ran alongside measurements. These are observations, not CI thresholds.

Four binaries use identical sources, dependencies and compiler options, changing
only `GridFunction.h`. All include the move-assignment correction and quadrature
identity base class, including the coordinate variant:

* **Coordinates:** header from `34adfb3db`, plus the move-assignment correction.
* **Identity:** baseline identity key from `3371d4e20`, plus the correction.
* **Reduced:** remove owner/formula address fields and repeated basis-context
  checks; retrieve the element from the DOF cache on a miss.
* **Reduced + lookup:** same reduced key, but retrieve the element through the
  space on a basis miss, as production does.

Runs 1–6 reverse variant order: reduced, identity, coordinates, coordinates,
identity, reduced. Each covers all 64 cases with three randomized repetitions
and minimum time 50 ms. Runs 7 and 9 repeat the full matrix for reduced + lookup;
run 8 confirms the identity P2Vector 3D SlowSweep/Mixed75 cases with five
repetitions. The table uses six observations per case and variant, excludes
run 8, and reports median CPU **ns per evaluation**. Direct is the pointwise
median from the identity binary, so its cache/direct comparison uses the same
executable. Pointwise timings in every binary are retained as controls in the
raw data; differing pointwise timings cannot be attributed to the quadrature
comparison. Whole-implementation comparisons include changed code/data layout,
not just the cost of individual comparisons.

| Field | Dim | Workload | Coordinates | Identity | Reduced | Reduced + lookup | Direct |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| P1Scalar | 2 | SlowSweep | 48.3 | 46.6 | 40.3 | 40.1 | 42.8 |
| P1Scalar | 2 | Mixed75 | 20.6 | 18.8 | 14.5 | 15.7 | 25.9 |
| P1Scalar | 2 | FastBlocks | 13.6 | 12.1 | 7.0 | 9.2 | 22.2 |
| P1Scalar | 2 | PureHit | 12.4 | 10.9 | 6.3 | 7.9 | 21.3 |
| P1Scalar | 3 | SlowSweep | 55.9 | 50.2 | 44.5 | 46.1 | 44.1 |
| P1Scalar | 3 | Mixed75 | 24.1 | 20.9 | 16.5 | 18.1 | 33.0 |
| P1Scalar | 3 | FastBlocks | 15.1 | 12.8 | 8.2 | 9.8 | 29.2 |
| P1Scalar | 3 | PureHit | 13.6 | 11.3 | 6.8 | 8.5 | 28.1 |
| P1Vector | 2 | SlowSweep | 68.3 | 65.3 | 60.3 | 58.6 | 42.0 |
| P1Vector | 2 | Mixed75 | 35.0 | 33.0 | 29.5 | 29.3 | 29.2 |
| P1Vector | 2 | FastBlocks | 23.7 | 21.9 | 18.5 | 18.6 | 26.1 |
| P1Vector | 2 | PureHit | 22.1 | 20.2 | 17.0 | 17.0 | 25.5 |
| P1Vector | 3 | SlowSweep | 121.7 | 104.6 | 98.8 | 98.0 | 49.3 |
| P1Vector | 3 | Mixed75 | 55.7 | 49.7 | 45.6 | 45.2 | 40.8 |
| P1Vector | 3 | FastBlocks | 33.3 | 31.4 | 28.0 | 27.8 | 38.2 |
| P1Vector | 3 | PureHit | 30.3 | 28.7 | 25.4 | 25.3 | 37.4 |
| P2Scalar | 2 | SlowSweep | 85.5 | 82.0 | 75.9 | 77.8 | 76.4 |
| P2Scalar | 2 | Mixed75 | 33.8 | 31.1 | 24.9 | 26.1 | 69.3 |
| P2Scalar | 2 | FastBlocks | 17.9 | 16.2 | 10.3 | 10.4 | 65.7 |
| P2Scalar | 2 | PureHit | 15.4 | 13.5 | 8.1 | 8.0 | 64.8 |
| P2Scalar | 3 | SlowSweep | 276.2 | 296.7 | 267.0 | 266.8 | 276.1 |
| P2Scalar | 3 | Mixed75 | 82.4 | 86.8 | 74.3 | 74.4 | 267.4 |
| P2Scalar | 3 | FastBlocks | 25.2 | 24.5 | 17.6 | 17.8 | 263.4 |
| P2Scalar | 3 | PureHit | 17.2 | 15.6 | 9.7 | 10.0 | 255.5 |
| P2Vector | 2 | SlowSweep | 179.3 | 178.3 | 170.6 | 172.4 | 80.5 |
| P2Vector | 2 | Mixed75 | 67.6 | 65.3 | 62.0 | 63.3 | 74.5 |
| P2Vector | 2 | FastBlocks | 34.1 | 32.2 | 28.7 | 29.1 | 70.8 |
| P2Vector | 2 | PureHit | 29.1 | 27.8 | 24.0 | 24.3 | 70.2 |
| P2Vector | 3 | SlowSweep | 904.8 | 840.6 | 949.6 | 964.3 | 280.3 |
| P2Vector | 3 | Mixed75 | 263.6 | 250.7 | 275.4 | 280.1 | 272.1 |
| P2Vector | 3 | FastBlocks | 80.7 | 76.4 | 76.0 | 77.8 | 269.5 |
| P2Vector | 3 | PureHit | 55.6 | 52.4 | 48.4 | 49.1 | 265.2 |

For P2Vector 3D, baseline identity caching costs about **3.0 times** direct expansion
on an all-miss traversal (841 versus 280 ns), gives a modest gain at 75% hits
(251 versus 272 ns), and improves about **3.5 times** at 96.875% hits
(76 versus 269 ns), or **5.1 times** for pure hits (52 versus 265 ns).
P1Vector 3D remains slower than direct at 75% hits. A cache is therefore useful
when its entries are reused; it is not a universal interpolation speedup.
The coordinate variant is usually slower on hits, but is faster in the
P2Scalar 3D all-miss measurements. No blanket coordinate-cost attribution is
made for that case.

The vector slow-path disadvantage has a concrete implementation explanation:
[`H1Element::evaluate`](../../src/Rodin/Variational/H1/H1Element.h) evaluates
one scalar basis value per node and reuses it across components. The quadrature
cache refill evaluates every vector basis separately and stores full vector
values. It repeats scalar-basis work per component, followed by coefficient
contraction. Reducing vector-refill duplication is a possible future optimization;
it is not implemented or claimed by these measurements.

### Investigation of the reduced-key miss regression

Reduced improves many hits but makes P2Vector 3D SlowSweep about 13% slower
than the baseline identity variant. Fresh element lookup does not resolve it. Runs 10–15 isolate
only this case, with five repetitions of at least 100 ms, reversed order,
and no random interleaving:

| Run | Variant | Median ns/evaluation |
| --- | --- | ---: |
| 10 | full-identity | 861.1 |
| 11 | reduced-identity | 937.0 |
| 12 | reduced-fresh | 929.0 |
| 13 | reduced-fresh | 922.3 |
| 14 | reduced-identity | 937.9 |
| 15 | full-identity | 848.6 |

The slowdown persists without the full matrix's allocation/workload context.
Clang's basis-cache routine shrinks from 195 instructions in the baseline identity variant to
167 (reduced) or 176 (reduced + lookup), while the scalar-basis evaluation
kernel has the same 130-instruction sequence in all three binaries. Fewer
comparisons/instructions alone do not establish lower runtime. Code/data
placement and surrounding generated code differ; the precise hardware cause
was not isolated. The existing identity checks are retained to avoid introducing this
measured miss regression, and the missing mesh-address guard is added for
correctness. Experimental patches are checked in for reproduction.

## Mesh-guard confirmation (before the coordinate correction)

Runs 16–17 measure the implementation with the mesh-address guard, before
the later exact-coordinate correction,
using the same complete 64-case matrix, three randomized repetitions per run,
and minimum time 50 ms. Every case passes its numerical oracle. These results
are kept in the separate `final-matrix` protocol. Representative P2Vector 3D
medians (six observations, ns/evaluation) are:

| Workload | Corrected identity | Direct expansion |
| --- | ---: | ---: |
| SlowSweep | 843.8 | 281.3 |
| Mixed75 | 246.5 | 275.5 |
| FastBlocks | 77.0 | 267.6 |
| PureHit | 51.9 | 264.1 |

The miss-heavy disadvantage and reuse benefit persist in the corrected code.
The original coordinate and reduced-key comparisons above refer to their
explicitly identified baseline binaries; they do not isolate the cost of
adding the mesh guard. The new guard fixes incorrect DOF reuse when space
assignment changes meshes, while this benchmark keeps space/mesh bindings fixed.

## Minimality under immutable connectivity

This audit assumes that the formula/index uniquely binds the point coordinates.
Mapped face points violate that assumption in the current API; the later
[coordinate correctness study](GridFunctionCoordinateKey.md) adds coordinates
to the basis entry key. The table below describes the earlier restricted audit.

For the current valid point/formula binding and built-in P1/H1 elements, the
following validation retains the necessary contexts in the current cache
representation (with fixed object lifetimes and valid space definitions):

| Entry | Context validated |
| --- | --- |
| DOFs | field lifetime identity, space address, mesh address, element address, polytope dimension/index |
| Basis, after DOF validation | basis-valid flag, formula lifetime/version identity, sample index |

The unique field identity makes its address redundant. The unique formula
identity (refreshed on assignment) makes its address redundant. DOF validation
already invalidates the basis entry when context changes, so repeating those
context checks in basis lookup is redundant. Formula identity/index identifies
the reference coordinates under the IntegrationPoint contract. Coefficient
values do not belong in this key: each evaluation contracts current coefficients.

Immutable **connectivity alone** does not make the space/element checks redundant.
Move assignment can rebind a field while preserving its field identity, even
when both meshes are immutable and their finite-element pointer is shared.
Space assignment can change vector dimension and element on the same immutable
mesh. It can also rebind the same space object to another immutable mesh without
changing its element pointer. The existing production key missed this case;
`SpaceMeshRebindingRefreshesCache` fails before the mesh-address guard is added
(observed error 1.56) and passes with it. The mesh-guard implementation includes
the mesh address in DOF validation, which invalidates basis values as well.
`MoveAssignmentRefreshesSpaceCache` and
`ElementDefinitionChangeRefreshesCaches` cover these independently. Suppressing
both guards in the reduced variant makes both regressions fail; restored guards
pass. This is an audit of the existing representation, not a proof of an
absolute minimum for arbitrary mutable custom elements/pushforwards. A stronger
immutable-space contract or an explicit space revision/lifetime token would
allow a different key. Mesh geometry changes are outside these benchmarks.

This key audit does not cover destroying/reconstructing a mesh behind a live
space, in-place DOF reordering, or mutable custom pushforwards. Supporting those
operations needs an explicit revision contract, rather than more address checks.

The move regression exposed a separate overload issue: forwarding a derived
GridFunction to `Parent::operator=` selected interpolation assignment. Casting
to the base rvalue selects its move assignment and correctly rebinds the space
before moving coefficients. This fix is included in production.

## Reproduce and inspect

```sh
cmake --build build --target RodinGridFunctionBenchmarks -j 4
build/tests/benchmarks/RodinGridFunctionBenchmarks \
  --benchmark_filter=^LargeInterpolation/ \
  --benchmark_min_time=0.05s --benchmark_repetitions=3 \
  --benchmark_enable_random_interleaving=true \
  --benchmark_out=large.json --benchmark_out_format=json
```

Use the dependency/configuration options described in the small-suite report.
Save each built executable before changing the header, and run them with no
concurrent builds. To test the smaller key, apply
[reduced-identity](variants/GridFunctionLarge-reduced-identity.patch) or
[reduced-fresh](variants/GridFunctionLarge-reduced-fresh.patch) with
`git apply`, build/save the executable, and reverse the patch before trying
the other variant. These patches reproduce the measured experimental headers,
including omission of the mesh-address and exact-coordinate guards: use them only for the fixed
space/mesh bindings in these benchmarks. For coordinates, extract the header
at `34adfb3db` and
apply the same base-rvalue move fix; retain the current quadrature base class.

[Raw per-repetition measurements](results/GridFunctionLarge-2026-09-30.csv)
include protocol, variant, iteration count, sweep CPU/wall times, evaluations,
normalized CPU time, error and mesh sizes. There are 1960 timed observations;
all passed numerical checks. Maximum checksum relative error was `2.09e-10`,
from sequential summation in multi-million-evaluation pure-hit sweeps. Do not
combine matrix, final-matrix, targeted-confirmation and isolated protocols into
one median.

Final regression verification: 127 tests passed across GridFunction, NamedForm,
NamedFormConsistency, GaussLegendre, QuadratureExactness and QF umbrella suites;
one existing geometry-specific elasticity test was skipped. The 48-case small
suite also passes after the production corrections.

The subsequent correctness-first coordinate comparison study is recorded in
[Exact coordinate cache keys](GridFunctionCoordinateKey.md). Its results refer
to a new paired protocol and should not be pooled with the earlier runs.
