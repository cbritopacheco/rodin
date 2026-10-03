# Spatial matrix assembly: CPU investigation

Measurements on 3 October 2026 reproduce an assembly overhead for native matrix
fields on the small meshes below. The additional cost is concentrated in local
entry traversal and triplet generation. Scalar tabulation is already shared,
and sparse conversion has nearly equal CPU cost for the two representations.
No production assembly optimization is included in this investigation.

## Measurement contract

The `RodinSpatialMatrixBenchmarks` target uses `MeasureProcessCPUTime()`.
The `CPU` column and JSON `cpu_time` count CPU consumed by the benchmark process,
including worker threads. They exclude scheduling delays and CPU consumed by
other processes. Shared caches, memory bandwidth and CPU frequency remain possible
sources of variability on a busy machine; seven randomized, interleaved repetitions
are summarized by their median.

Configuration: Clang 19, `RelWithDebInfo` (`-O2 -g -DNDEBUG`), local Eigen assembly
on a ten-core MacBook Pro. OpenMP is compiled in, but sequential cases explicitly
use `Assembly::Sequential`. `OMP_NUM_THREADS=1` and `OMP_WAIT_POLICY=PASSIVE`
are set for the sequential baseline. Triangle meshes have 98 cells; tetrahedral
meshes have 48 cells. Orders P1, P2 and P3 are tested.

One native assembly of a real three-by-three matrix field is compared with nine
repeated assemblies of the same scalar form. This keeps the nine components'
mathematics and total DOFs comparable while retaining warmed caches. This is an
uncoupled mass/stiffness comparison, not Oscar's full viscoelastic application
or a flattened vector implementation.

Before timing, the native sparse operator is compared with the scalar block
operator, using H1's documented component ordering:

$$
  A_{\mathrm{matrix}}=A_{\mathrm{scalar}}\otimes I_9.
$$

The maximum observed relative Frobenius discrepancy was approximately
$1.35\times10^{-16}$. Triplet conversion and the timed complete assembly are
also checked against a sequential reference. An exact sparse nonzero-count check
was rejected: cancellation in a P2 triangle stiffness case can retain different
near-zero entries despite equivalent numerical operators. Numerical comparisons
therefore use a relative tolerance of $10^{-12}$, rather than assuming identical
storage patterns.

## Complete sequential assembly

Every row is the CPU time for one nine-component problem, in milliseconds.

| Order | Dimension | Form | Native matrix CPU ms | Nine scalar CPU ms | Matrix/scalar |
| --- | --- | --- | ---: | ---: | ---: |
| P1 | 2 | Mass | 0.169 | 0.126 | 1.33 |
| P1 | 2 | Stiffness | 0.156 | 0.134 | 1.16 |
| P1 | 3 | Mass | 0.151 | 0.099 | 1.52 |
| P1 | 3 | Stiffness | 0.144 | 0.100 | 1.44 |
| P2 | 2 | Mass | 0.683 | 0.467 | 1.46 |
| P2 | 2 | Stiffness | 0.667 | 0.459 | 1.45 |
| P2 | 3 | Mass | 0.966 | 0.664 | 1.45 |
| P2 | 3 | Stiffness | 0.896 | 0.592 | 1.51 |
| P3 | 2 | Mass | 2.045 | 1.365 | 1.50 |
| P3 | 2 | Stiffness | 1.959 | 1.322 | 1.48 |
| P3 | 3 | Mass | 4.546 | 3.439 | 1.32 |
| P3 | 3 | Stiffness | 4.257 | 3.484 | 1.22 |

## Stage attribution

The diagnostic stages distinguish element binding/local matrix construction,
scanning a cached local matrix, triplet generation and sparse conversion. Each
uses the same forms and quadrature orders as complete assembly. The cached scan
visits one cell; other stages traverse the complete mesh. Isolated triplet
generation retains vector capacity, whereas complete assembly creates a fresh
triplet vector. Stage timings should therefore not be added as an exact
reconstruction of complete assembly.

For 3D P3 stiffness:

| Stage | Native CPU ms | Nine scalar CPU ms | Matrix/scalar |
| --- | ---: | ---: | ---: |
| Full | 4.2575 | 3.4838 | 1.22 |
| Triplets | 3.0879 | 2.3517 | 1.31 |
| Bind | 1.3194 | 1.9927 | 0.66 |
| Scan | 0.1169 | 0.0130 | 9.02 |
| Sparse | 1.0585 | 1.0385 | 1.02 |

Across the sequential sweep, cached scans cost about nine times as much for the
matrix representation. Element binding/local construction is faster than the
nine scalar bindings; sparse conversion is close to parity. Triplet generation
is slower for the matrix representation.

For a scalar element with $n$ basis functions, the nine-component matrix element
has $9n$ local DOFs. The optimized plain mass and gradient integrators allocate
and zero a dense $(9n)\times(9n)$ local matrix, fill equal-component entries,
and mirror its triangular part where applicable. The generic assembler then
visits all $81n^2$ entries and rejects numerical zeros. The nine scalar
assemblies visit $9n^2$ entries in total. Thus $8/9$ of the matrix component
pairs are structurally zero in these uncoupled forms, but still incur local
storage and traversal work.

The source paths are the dense local matrices and component replication in
[`H1/QuadratureRule.h`](../src/Rodin/Variational/H1/QuadratureRule.h), and the
nested local-entry loops preceding the zero check in
[`Assembly/Sequential.h`](../src/Rodin/Assembly/Sequential.h).

A twelve-second Instruments Time Profiler recording of the 3D P3 native stiffness
case independently supports the split. Waiting threads are excluded. After
excluding benchmark initialization and retaining stacks inside the sparse
assembly call, 10.957 seconds of weighted CPU samples remain. About 74.45% include
triplet generation and 25.55% include Eigen's sparse conversion. Element binding
accounts for about 30.30%, cached integrator entry access for 12.82%, and inlined
DOF-array size/access locations for 14.87%. These are inclusive stack attributions
and overlap; they are not independent percentages to sum.

A temporary prototype cached local DOF loop bounds once per cell. Full native
assembly changes ranged from approximately a 2% regression to a 3% improvement
relative to the first sweep, with no consistent gain. The prototype still visits
all zero component pairs and materializes the same dense local matrix. It was
removed; the production assembler is unchanged. The profile's inlined metadata
attribution alone is insufficient evidence that loop-bound caching is useful.

## OpenMP CPU costs

The same operator checks pass for one and four workers. Total process CPU time
includes all of them; it is a work/CPU-efficiency measure rather than parallel
elapsed-time speedup. P3 examples:

| Workers | Dimension | Form | Native CPU ms | Nine scalar CPU ms | Matrix/scalar |
| --- | --- | --- | ---: | ---: | ---: |
| 1 | 2 | Mass | 2.050 | 1.415 | 1.45 |
| 1 | 2 | Stiffness | 1.980 | 1.351 | 1.46 |
| 1 | 3 | Mass | 4.622 | 3.446 | 1.34 |
| 1 | 3 | Stiffness | 4.178 | 3.529 | 1.18 |
| 4 | 2 | Mass | 2.301 | 2.266 | 1.02 |
| 4 | 2 | Stiffness | 2.219 | 2.302 | 0.96 |
| 4 | 3 | Mass | 4.914 | 4.629 | 1.06 |
| 4 | 3 | Stiffness | 4.549 | 4.626 | 0.98 |

## Next optimization and remaining scope

The measured target is the common local-entry iteration. Structural component
coupling information from the integrator could avoid visiting zero component
pairs, and compact local storage could avoid materializing those entries.
Any such change must preserve full coupling for rank-four coefficients and
other forms whose component blocks are dense. Inferring structural zeros from a
floating-point threshold would change the numerical model and is inappropriate.

The current cases use real coefficients, affine triangles/tetrahedra and
uncoupled forms. These results identify a Rodin overhead but do not establish
Oscar's dominant application cost. His exact constitutive, transport and
coupling forms, quadrature orders, mesh and backend still need an application
profile. Tensor-value arithmetic benchmarks alone do not reproduce these costs.

## Reproduction

Build `RodinSpatialMatrixBenchmarks` in an optimized configuration with
`RODIN_BUILD_BENCHMARKS=ON`. For the sequential CPU sweep:

```sh
OMP_NUM_THREADS=1 OMP_WAIT_POLICY=PASSIVE \
  build/tests/benchmarks/RodinSpatialMatrixBenchmarks \
  --benchmark_filter='/Sequential/' \
  --benchmark_min_time=0.05s --benchmark_repetitions=7 \
  --benchmark_enable_random_interleaving=true \
  --benchmark_out=spatial-assembly-cpu.json --benchmark_out_format=json
```

For OpenMP, select `/OpenMP/` and set the worker count to one or four. Read the
JSON `cpu_time` field and median aggregates. The reporter also emits elapsed time;
that column was not used for this comparison.

Local validation: all 120 sequential stage cases, and all 24 OpenMP cases at each
worker count, completed seven repetitions without a numerical-reference error.
Benchmark initialization, operator checks and work counters are outside timing.
