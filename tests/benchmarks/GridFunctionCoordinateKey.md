# Exact coordinate cache keys: correctness and cost

## Decision and scope

GridFunction basis caching now validates exact reference-coordinate representations
in addition to the existing identity/context key. `Math::SpatialPoint` stores
the snapshot inline, with a maximum of three components. No coordinate-key
allocation or hashing occurs. This restores correctness for mapped quadrature
points without introducing a separate mapped-point evaluation API or bypassing
basis reuse whenever a point is mapped.

The comparison is of **two or three coordinates at the evaluated point**, not
all vertices or all mesh coordinates. The accepted small cost preserves basis
reuse. Formula lifetime/version identity, mesh/space/element guards and current
coefficient contraction remain in place.

This study assumes fixed mesh geometry and valid space definitions, like the
previous large-mesh measurements. Arbitrary mesh mutation and custom mutable
pushforwards need their own revision contract; coordinates do not certify those
unrelated dependencies.

## Why exact coordinates are needed and accurate

`IntegrationPoint` can retain a face quadrature formula/index after the point
is mapped into a cell, as in `Average` and `Jump`. Formula/index is provenance,
not a unique identifier for the cell reference coordinates. Those operators'
particular traversal order can force misses already, but does not establish a
general uniqueness guarantee for the public IntegrationPoint API.

A diagnostic used P1 on a two-triangle mesh, `u = 1 + x + 10*y`, and two mapped
edge midpoints on the same cell with the same segment formula/index. The
identity-only implementation returned **1.5 instead of 6** at the second point;
the earlier coordinate-key implementation returned 6. This demonstrates an
actual coordinate-key alias, not a tolerance-level interpolation discrepancy.

The production correction compares each Real's bytes with `std::memcmp`, along
with coordinate dimension. It compares individual scalars, not object padding.
For the current binary64 Real, up to 24 coordinate bytes participate. A hit
requires identical coordinate representations; there is no epsilon, rounding,
quantization or collision-prone hash. Even adjacent representable values and
positive/negative zero are distinct. Slightly different inverse-map results
can produce extra misses; they cannot merge different coordinate representations
into one cache hit. The check introduces no approximation of basis values.
Nonfinite coordinates are not certified as meaningful finite-element inputs.

The copied snapshot is independent of the Point object's lifetime. Misses
recompute the existing basis/pushforward code at the actual Point; hits contract
current coefficients. Formula assignment/reconstruction and field/space/mesh
rebinding regressions remain covered.

## Correctness gates

* `MappedQuadratureSamplesRefreshBasisCache` alternates actual quadrature samples
  embedded on distinct reference-cell faces, while holding formula/index and
  cell fixed. Scalar/vector P1/P2 are checked in 2D/3D against direct expansion,
  including repeated hits and returns to the first face.
* `CoordinateKeyUsesExactRepresentations` uses the existing instrumented range
  pushforward to verify that a one-ULP coordinate change and a signed-zero
  change each rebuild the basis, while an identical repeat reuses it. This
  distinguishes invalidation even when the final rounded field value agrees.

Both tests pass with the correction. Removing the coordinate predicate makes
both tests fail, so passing output checks alone is not the only evidence.
The final six suites pass **129 tests**, with one existing geometry-specific
elasticity skip. The comparison is implemented in the cache owner, with no
changes to Average/Jump traversal or finite-element arithmetic.

## Paired benchmark protocol

Host: Apple M1 Pro, 10 cores; Clang 19.1.7; C++20, `-O3 -DNDEBUG`;
Google Benchmark `v1.9.4-34-g4b7b129e`, Release library, single-threaded.
The machine is not dedicated to benchmarking. No compilation ran concurrently
with measurement. Results describe this build/host, not universal thresholds.

The identity executable uses `7ef37bd41`, including the mesh guard and move
assignment fix. The exact-coordinate executable differs only by the coordinate
snapshot and predicate in GridFunction.h. Thus the earlier coordinate-vector
implementation is not substituted for the current code in this experiment.

Runs 1–4 use identity/exact/exact/identity order. Each runs all **112 cases**
(64 large, 48 small), with three randomized repetitions and minimum time 30 ms.
Large fixtures, workload hit ratios and analytic oracles are described in
[GridFunctionLarge.md](GridFunctionLarge.md); re-interpolation operations in
[GridFunctionInterpolation.md](GridFunctionInterpolation.md). Both binaries
pass all ordinary quadrature benchmark oracles. The identity-only binary fails
the mapped-point correctness gate and is therefore a performance reference,
not an acceptable general solution.

The following tables pool six CPU observations per case/variant. Large times
are CPU time per sweep divided by evaluations per sweep; small times are per
operation. Do not combine this protocol with older 64-case-only runs: fixture
initialization and measurement context differ.

### Hit-path overhead (ns/evaluation)

| Field | Dim | Identity | Exact coordinates | Added ns | Added % |
| --- | --- | ---: | ---: | ---: | ---: |
| P1Scalar | 2 | 12.91 | 14.15 | 1.24 | 9.6% |
| P1Scalar | 3 | 13.22 | 15.06 | 1.84 | 13.9% |
| P1Vector | 2 | 20.67 | 21.66 | 0.99 | 4.8% |
| P1Vector | 3 | 28.89 | 30.72 | 1.83 | 6.3% |
| P2Scalar | 2 | 14.19 | 15.81 | 1.62 | 11.4% |
| P2Scalar | 3 | 16.05 | 17.24 | 1.19 | 7.4% |
| P2Vector | 2 | 28.19 | 29.69 | 1.51 | 5.3% |
| P2Vector | 3 | 53.02 | 54.83 | 1.80 | 3.4% |

Most hit paths add approximately **1–2 ns**, with relative overhead about
**3–14%** in this very small operation. This is the whole implementation
change (predicate, snapshot/layout and resulting code), not a claim to time
bare floating-point comparisons. The absolute overhead and relative overhead
should both be retained when judging performance.

### Miss, mixed and fast paths (ns/evaluation)

| Field | Dim | Workload | Identity | Exact coordinates | Direct (exact binary) |
| --- | --- | --- | ---: | ---: | ---: |
| P1Scalar | 2 | SlowSweep | 52.5 | 50.1 | 43.0 |
| P1Scalar | 2 | Mixed75 | 21.6 | 22.8 | 27.5 |
| P1Scalar | 2 | FastBlocks | 14.8 | 16.1 | 23.6 |
| P1Scalar | 2 | PureHit | 12.9 | 14.2 | 23.4 |
| P1Scalar | 3 | SlowSweep | 53.7 | 55.6 | 44.1 |
| P1Scalar | 3 | Mixed75 | 23.2 | 25.6 | 34.1 |
| P1Scalar | 3 | FastBlocks | 15.3 | 16.7 | 31.3 |
| P1Scalar | 3 | PureHit | 13.2 | 15.1 | 30.3 |
| P1Vector | 2 | SlowSweep | 66.5 | 64.0 | 38.9 |
| P1Vector | 2 | Mixed75 | 33.8 | 34.7 | 29.0 |
| P1Vector | 2 | FastBlocks | 22.4 | 23.3 | 26.1 |
| P1Vector | 2 | PureHit | 20.7 | 21.7 | 26.5 |
| P1Vector | 3 | SlowSweep | 108.3 | 114.8 | 46.9 |
| P1Vector | 3 | Mixed75 | 50.6 | 53.9 | 41.1 |
| P1Vector | 3 | FastBlocks | 31.6 | 34.2 | 38.0 |
| P1Vector | 3 | PureHit | 28.9 | 30.7 | 38.9 |
| P2Scalar | 2 | SlowSweep | 83.3 | 84.2 | 76.9 |
| P2Scalar | 2 | Mixed75 | 32.1 | 33.9 | 69.4 |
| P2Scalar | 2 | FastBlocks | 16.6 | 18.2 | 65.6 |
| P2Scalar | 2 | PureHit | 14.2 | 15.8 | 65.2 |
| P2Scalar | 3 | SlowSweep | 302.8 | 278.9 | 287.7 |
| P2Scalar | 3 | Mixed75 | 87.8 | 83.6 | 280.7 |
| P2Scalar | 3 | FastBlocks | 25.1 | 25.6 | 274.9 |
| P2Scalar | 3 | PureHit | 16.0 | 17.2 | 268.4 |
| P2Vector | 2 | SlowSweep | 175.9 | 177.0 | 80.6 |
| P2Vector | 2 | Mixed75 | 66.5 | 66.7 | 73.7 |
| P2Vector | 2 | FastBlocks | 33.1 | 34.2 | 70.8 |
| P2Vector | 2 | PureHit | 28.2 | 29.7 | 70.5 |
| P2Vector | 3 | SlowSweep | 946.0 | 886.2 | 286.7 |
| P2Vector | 3 | Mixed75 | 282.3 | 264.6 | 283.6 |
| P2Vector | 3 | FastBlocks | 81.1 | 80.6 | 272.3 |
| P2Vector | 3 | PureHit | 53.0 | 54.8 | 269.6 |

Exact-coordinate caching remains substantially faster than direct expansion
when basis values are reused. Its vector all-miss path still costs much more
than direct expansion; the earlier explanation about repeated scalar-basis
work in vector refill still applies. Coordinates do not fix that separate cost.
No miss-path speedup is attributed to the coordinate comparison.

### Whole-field interpolation and re-interpolation (ns/operation)

I/E means identity/exact-coordinate. Pointwise controls and all 48 small
cases are retained in the raw data.

| Field | Dim | Reproject I/E | Cross-space assignment I/E | First interpolation I/E |
| --- | --- | ---: | ---: | ---: |
| P1Scalar | 2 | 3053.2 / 3051.7 | 3232.6 / 3208.0 | 3116.0 / 3183.3 |
| P1Scalar | 3 | 41926.7 / 42035.0 | 40892.1 / 41224.3 | 42041.4 / 41600.4 |
| P1Vector | 2 | 5830.4 / 5744.5 | 4974.5 / 4879.2 | 5917.6 / 5876.0 |
| P1Vector | 3 | 129424.1 / 127215.8 | 99672.0 / 100034.9 | 129057.9 / 127504.3 |
| P2Scalar | 2 | 4849.1 / 4805.4 | 9360.8 / 9244.3 | 4915.0 / 4953.9 |
| P2Scalar | 3 | 83468.7 / 84003.0 | 450367.3 / 446412.6 | 84044.5 / 84559.9 |
| P2Vector | 2 | 10243.1 / 9863.3 | 18902.8 / 18866.3 | 10304.9 / 10140.2 |
| P2Vector | 3 | 313581.6 / 288882.8 | 1422429.2 / 1367670.3 | 307421.7 / 291326.5 |

These whole-field differences vary with case and context. No universal
whole-field speedup or negligible-cost claim is made from these short runs.

## Isolated comparison and snapshot controls

The full matrix showed unexpectedly lower miss times with the exact key in
3D P2. To investigate, a third binary preserves the new coordinate snapshot
and miss-path writes, but removes its hit predicate. The compiler removes the
now-unused comparison. This **snapshot-only** variant is a diagnostic ablation;
it has the same mapped-point correctness defect as identity-only and is not
eligible for production.

Runs 5–10 execute only P2Vector 3D SlowSweep (quadrature and pointwise) and
PureHit (quadrature), with five repetitions of at least 100 ms, fixed case
order, and reversed variant order identity/snapshot/exact/exact/snapshot/identity.
Median CPU ns/evaluation per run:

| Run | Variant | All misses | Pure hits | Pointwise miss control |
| --- | --- | ---: | ---: | ---: |
| 5 | identity | 870.3 | 52.9 | 280.6 |
| 6 | snapshot-only | 910.4 | 54.0 | 291.3 |
| 7 | exact | 907.9 | 55.8 | 296.4 |
| 8 | exact | 917.5 | 55.5 | 295.9 |
| 9 | snapshot-only | 902.5 | 54.0 | 288.5 |
| 10 | identity | 910.0 | 53.6 | 285.3 |

The full-matrix miss improvement does **not** reproduce in isolation. The
snapshot-only control is close to exact on misses. Pointwise control timings
also differ despite bypassing the coordinate predicate. Therefore matrix
miss differences cannot be interpreted as the predicate's intrinsic cost or
as a universal improvement. Context and generated code/data layout contribute;
the exact microarchitectural mechanism is not established.

On pure hits, isolated pooled medians are approximately **53.3 ns identity,
54.0 ns snapshot-only, and 55.7 ns exact**. This puts the complete correction's
observed cost at about **2.4 ns**, and the predicate's ablation difference at
about **1.6 ns**. Removing the predicate also changes generated code, so the
ablation is an estimate rather than a standalone instruction latency claim.

Bounded disassembly of the H1 P2 vector basis-cache routine shows 196 instructions
for identity, 222 for snapshot-only and 243 for exact in this build. The
component comparisons are inlined; there is no external memcmp/bcmp call in
that routine. This verifies that the source-level library call does not imply
an out-of-line function call per coordinate. Instruction count alone is not
a runtime model.

## Reproduce

Build `RodinGridFunctionBenchmarks` in Release using the dependency configuration
from the interpolation report. Save the identity executable built at
`7ef37bd41`, then build/save the current coordinate-key executable. Run each
without concurrent builds:

```sh
build/tests/benchmarks/RodinGridFunctionBenchmarks \
  --benchmark_min_time=0.03s --benchmark_repetitions=3 \
  --benchmark_enable_random_interleaving=true \
  --benchmark_out=coordinate-key.json --benchmark_out_format=json
```

For the isolated protocol, filter
`LargeInterpolation/P2Vector/(SlowSweep/(Quadrature|Pointwise)|PureHit/Quadrature)/dimension:3`,
use minimum time `0.1s`, five repetitions, and omit random interleaving.
The snapshot-only ablation removes `|| !sameReferenceCoordinates` from the
basis-miss condition while retaining coordinate storage/writes; it must never
be used as a correctness-preserving implementation.

All **1434** benchmark observations passed their applicable numerical oracle. Maximum large
checksum relative error was `2.09e-10`; maximum small absolute error `1.79e-14`.
Mapped correctness is verified separately by the regression gates above, and
must not be inferred from native-quadrature performance cases alone.
