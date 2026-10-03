# Spatial value benchmarks

The spatial arithmetic benchmarks are part of the existing `RodinBenchmarks`
executable and benchmark publication workflow. Build an optimized configuration
with `RODIN_BUILD_BENCHMARKS=ON`, then run:

```sh
build/tests/benchmarks/RodinBenchmarks \
  --benchmark_filter='^Spatial(Vector|Matrix|Tensor[34])/' \
  --benchmark_min_time=0.1s --benchmark_repetitions=3 \
  --benchmark_out=spatial-algebra.json --benchmark_out_format=json
```

Each result measures one construction, copy, addition, scalar multiplication,
inner product, or contraction. Operand initialization, counters and independent
reference checks run outside the timed loop. Compiler barriers expose input
objects and materialize outputs on each iteration. Timings include these
barriers, so very small operations should be compared under the same compiler,
optimization flags and hardware. `entries` reports the active input entries;
`object_bytes` reports the complete object size, including unused capacity and
metadata. Throughput counts active entries (multiply-accumulate terms for
`MatMat`), not floating-point operations.

Both real and complex values are registered. Vectors cover sizes zero through
three. Matrices cover the empty shape and all nine nonempty shapes, with all
27 compatible matrix-product dimension combinations. Rank-three and rank-four
tensors cover the empty shape, uniform axis extents one through three, and two
rectangular shapes each. These tensor cases are representative, not an exhaustive
axis-extent sweep or coverage of every possible rank.

`Contract` measures matrix-vector, rank-three tensor-vector, or rank-four
tensor-matrix multiplication. Products use ordinary linear contractions, checked
against independent sums before timing. `Dot` measures each type's existing
method: complex vector `dot` is bilinear, while complex matrix and tensor `dot`
conjugate the second operand. Their complex timings therefore measure different
arithmetic contracts.

`SpatialTensor` owns a fixed-capacity `std::array` of entries. Its capacity is
three to the power of its compile-time rank under the current configuration;
changing runtime extents does not allocate storage. A local tensor stores that
array on the stack; a tensor embedded in another object shares that object's
storage duration. Rank-three and rank-four real tensors reserve 27 and 81 scalar
entries, respectively. The object size counter also includes their extents and
active-entry count.

The separate `RodinSpatialMatrixBenchmarks` executable measures finite-element
assembly against nine scalar assemblies. Those assembly costs should not be
inferred from these value-operation timings.
