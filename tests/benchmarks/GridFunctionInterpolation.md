# Grid-function interpolation benchmarks

`RodinGridFunctionBenchmarks` measures scalar and vector P1 and H1 P2 fields
on triangular (2D) and tetrahedral (3D) meshes. Each mesh has four vertices
per coordinate direction. The benchmark names identify the field, operation,
and dimension. Meshes, spaces, source fields and quadrature points are set up
outside timing. Every case checks its final result against an analytic linear
field and reports `absolute_error` and the number of DOFs. A wrong numerical
result fails the case through `SkipWithError`; timings alone are not acceptance.

| Operation | Timed work |
| --- | --- |
| QuadratureHit | Repeated evaluation of the same field, cell and sample |
| QuadratureMiss | Alternate two sample indices in the same cell, refreshing basis values |
| Pointwise | Direct geometric-point evaluation, bypassing the quadrature basis cache |
| Reproject | Interpolate changing analytic data into an existing field, then evaluate |
| CrossSpaceAssignment | Assign alternating source fields from a distinct space of the same order on the same mesh, then evaluate |
| FirstInterpolation | Construct a field, interpolate analytic data, evaluate and destroy the field |

`Reproject` and `CrossSpaceAssignment` alternate two different fields, so the
compiler cannot collapse repeated identical inputs. P1/P2 here describes the
approximation order; cross-space assignment uses distinct spaces of equal
order. The unit regression suite additionally exercises repeated P1-to-P2
assignment. Construction of meshes, spaces, and quadrature formulas is excluded
from the timings, including the quadrature identity's construction-time atomic
increment. FirstInterpolation does include grid-function allocation and its
existing identity increment. These measurements do not estimate the quadrature
identity's construction cost or the speed of an entire PDE solve.

## Run

Enable `RODIN_BUILD_BENCHMARKS=ON` in a Release build and build only the dedicated
target:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRODIN_BUILD_BENCHMARKS=ON
cmake --build build --target RodinGridFunctionBenchmarks -j 4
build/tests/benchmarks/RodinGridFunctionBenchmarks \
  --benchmark_filter=^Interpolation/ \
  --benchmark_min_time=0.05s --benchmark_repetitions=5 \
  --benchmark_enable_random_interleaving=true \
  --benchmark_out=interpolation.json --benchmark_out_format=json
```

The normal project dependency options still apply. The recorded local build
used `RODIN_MULTITHREADED=OFF`, `RODIN_USE_PETSC=OFF`, and
`RODIN_INSTALL_RESOURCES=OFF`, with MacPorts' MPICH/Clang wrapper to satisfy the
installed HDF5's MPI header requirement. No compilation ran during measurement.

## Coordinate comparison experiment — 2026-09-30

Host: Apple M1 Pro, 10 cores; Clang 19.1.7; C++20, `-O3 -DNDEBUG`;
Google Benchmark `v1.9.4-34-g4b7b129e`, Release library. Measurements were
single-threaded. The machine was not dedicated to benchmarking: one-minute
load averages during the main runs were 4.5–5.3. Results are observations on
this host, not portable thresholds or CI performance assertions.

Two executables were built from identical benchmark sources, dependency code,
and compiler options. Only `GridFunction.h` changed: the coordinate variant
used the version at `34adfb3db`, while the identity variant used `b7b11cdd8`.
In particular, both executables included the new quadrature identity base
class, so construction changes cannot confound the cache-key comparison.
The experiment measures the complete cache implementations, including the
coordinate vector and its updates; it is not a timing of bare comparisons.

Run order was identity, coordinates, coordinates, identity. Each run covered
all 48 cases, with five randomized repetitions and at least 50 ms per timed
sample. The table reports medians across ten observations per variant.
“Coordinate overhead” is `(coordinates / identity - 1) * 100`, rather than a
claim that the full application improves by that percentage.

| Field | Dimension | Identity (ns) | Coordinates (ns) | Difference (ns) | Coordinate overhead |
| --- | --- | --- | --- | --- | --- |
| P1Scalar | 2 | 13.00 | 14.59 | 1.59 | 12.3% |
| P1Scalar | 3 | 12.93 | 15.89 | 2.97 | 23.0% |
| P1Vector | 2 | 20.06 | 21.24 | 1.17 | 5.8% |
| P1Vector | 3 | 28.02 | 30.14 | 2.12 | 7.6% |
| P2Scalar | 2 | 15.57 | 19.54 | 3.96 | 25.5% |
| P2Scalar | 3 | 16.97 | 22.72 | 5.75 | 33.9% |
| P2Vector | 2 | 27.16 | 28.68 | 1.52 | 5.6% |
| P2Vector | 3 | 51.34 | 55.33 | 3.98 | 7.8% |

The coordinate cache adds about 1.2–5.8 ns in these repeated-hit cases. Its
relative cost is visible because the whole hit path is only about 13–55 ns.
This supports removing the coordinate comparison from that hot path; it does
not establish that coordinate comparisons dominate interpolation generally.
Basis-cache misses show smaller or more variable relative differences as
basis evaluation begins to dominate.

For whole-field re-interpolation and first interpolation, the differences
were small, varied by field and dimension, and sometimes favored coordinates.
No causal whole-field speedup is claimed. A longer follow-up on the scalar P1
2D cases used ten repetitions of at least 200 ms: Reproject measured 3163.85 ns
with identities and 3165.22 ns with coordinates, a difference of 0.04%.
The repeated-hit case still measured 13.60 versus 14.80 ns (8.9% coordinate
overhead).

### Control investigation

Pointwise evaluation never executes the changed quadrature cache check.
In the main runs, scalar P1 2D pointwise medians were 36.32 ns (identity) and
36.09 ns (coordinates). One mixed-workload follow-up measured 27.14 ns for
the identical identity executable, with individual samples spanning 25.9–35.7
ns, versus 35.58 ns for coordinates. That result could not be attributed to
coordinate comparison, because the operation bypasses it.

Four isolated pointwise runs, in identity/coordinates/identity/coordinates
order, returned medians 36.60, 36.97, 36.61, and 37.01 ns. The discrepancy did
not reproduce in isolation. This identifies the earlier speedup attribution
as dependent on measurement context rather than evidence about the changed
check; the specific hardware/scheduling mechanism was not isolated. The
pointwise case is retained as a control, and no pointwise speedup is claimed.
Repeated randomized runs, reversed variant order, and controls are necessary
before interpreting nanosecond differences.

All measured cases passed their numerical oracle. The largest absolute error
across the main runs was 8.15e-15. Raw per-repetition CPU and wall times, errors,
iteration counts, protocol and variant are recorded in
[the measurement CSV](results/GridFunctionInterpolation-2026-09-30.csv).
Runs 1–4 are the full matrix, 5–6 the longer four-case confirmation, and 7–10
the isolated pointwise controls. Do not combine protocols into one median.

## Correctness gates

`RodinVariationalGridFunctionTest` covers reconstructed grid-function storage,
reconstructed and assigned quadrature formulas, changing coefficients under a
warmed cache, repeated P1-to-P2 assignment, scalar/vector re-interpolation, and
alternating quadrature and generic points. Removing the quadrature identity
comparison makes both formula storage/assignment regressions fail. Multiplying
cached basis values by two makes all four new interpolation/evaluation tests
fail. Restoring the implementation makes them pass. The benchmarks retain
analytic-value checks independently of these unit tests.

For large immutable-mesh traversals, controlled miss/hit mixtures, direct
expansion comparisons, and a cache-key audit, see
[the large-workload report](GridFunctionLarge.md).

The subsequent correctness-first coordinate comparison study is recorded in
[Exact coordinate cache keys](GridFunctionCoordinateKey.md). Its results refer
to a new paired protocol and should not be pooled with the earlier runs.
