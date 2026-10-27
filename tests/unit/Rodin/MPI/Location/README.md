# MPI AABB contract and validation

Include `Rodin/MPI/Location.h` (or `Rodin/MPI.h`) and construct
`Location::AABB locator(mesh)` for an MPI mesh or MPI submesh.

Construction and queries are noncollective. The wrapper selects owned entities
at each requested dimension and composes `AABB<LocalMesh>` over those indices.
Ghost/shared support geometry remains available, but ghost/shared entities are
never search candidates, including in exhaustive fallback. Results refer to the
MPI mesh and retain shard-local indices; use `getGlobalIndex(d, i)` explicitly
when a distributed identifier is needed. A miss is local. Two ranks may each
find an owned incident cell at a shared spatial boundary.

Tolerance uses the full shard's vertex bounding-box diagonal, and can vary with
partitioning/overlap. Reconstruct after geometry, topology, ownership or
reconciliation changes. The mesh must outlive the locator and returned points.
Built-in MPI shard transport registers real P1 and real H1 geometry orders 1–4;
applications must register other polymorphic transformation types with Boost.
Point geometry has no nontrivial polynomial order.

The tests cover:

- Sparse owned indices and overlapping ghost/shared/owned geometric images.
- Owned cells using nonowned vertices, and dimension-specific ownership.
- Shared partition boundaries and ghost-only interior misses with fallback.
- Construction and queries on one rank only, and unequal query counts.
- Empty parent shards and MPI submeshes with collective logical dimensions.
- All eight geometries, native/embedded configurations and H1 orders 1–4.
- Preservation of coordinates/Jacobians through parent sharding and MPI transfer,
  including cached default P1 maps as well as attached H1 maps.
- Lifted MPI mesh identity, local/distributed index maps and P0 field evaluation.

Run `RodinMPILocationAABBTest_np*` with CTest on 1–4 ranks. The tests have a
bounded timeout so accidentally introduced collectives fail instead of hanging.

## Benchmarks

`RodinMPIAABBBenchmarks` compares the full shard, its owned candidate subset and
the MPI wrapper. It includes build, interior, boundary and nearby exterior miss
cases across all geometries, orders 1–4, mesh resolutions, curvature and overlap.
Mesh setup is excluded; build includes first lookup. Synthetic cyclic partitions
and full-overlap cases stress ownership filtering rather than MPI communication.
Benchmark-only barriers and reductions report the slowest rank; locator calls
remain noncollective. Timings use adaptive repetition for cheap queries/builds.

For example:

```sh
mpiexec -n 2 ./RodinMPIAABBBenchmarks \
  --benchmark_filter='MPIAABB/Triangle/P3/.*/Interior/5/1/2' \
  --benchmark_out=timings.json
```

Other ranks write suffixed output files when `--benchmark_out=...` is supplied.
For operation counts, use `dev/profile_aabb.py --target RodinMPIAABBBenchmarks`
with `--build` and `--output`. This generates an instrumented executable outside
the source tree. Its counters report Newton candidates, transformations,
Jacobians, iterations, retries and active tree storage. Candidate identifier
storage is reported separately. Diagnostic timings are not performance evidence.
