# Named forms and P2 evaluation: review and measurements

## Review outcome

PR #330 adds specialized assembled operators, not string metadata: mass,
diffusion, Helmholtz and linear elasticity. No new numerical defect was found
in this review of the local sequential paths. The exact coordinate guard is
retained; its separate paired ablation is recorded in
[GridFunctionCoordinateKey.md](GridFunctionCoordinateKey.md).

The added consistency coverage checks P2 vector mass, P2 elasticity over all
six applicable cell geometries, copied/cloned/moved forms evaluated after the
original is destroyed, temporary function coefficients, repeated updates to a
P2 GridFunction coefficient, and P2 composition with a load and Dirichlet
conditions (matrix and RHS). Existing tests cover seven geometries, mixed
orders, restrictions, scatter invalidation and backend interfaces.

Ownership is documented on the four classes. Trial/test functions, their
spaces and mesh are borrowed. Function expressions are cloned; external
captures remain the caller's responsibility. GridFunction coefficients retain
a live field reference, including in copied forms. Constructors assemble
immediately; changing coefficient data requires reassembly to update the matrix.

All 102 tests in the three affected executables pass (35 consistency, 32
named-form and 35 GridFunction tests); one existing Segment elasticity case
is skipped. Doubling the named mass kernel's factor made all three new tests
fail and the benchmark reject the matrix before timing. Restoring the factor
restored passing tests. No assembler or interpolation arithmetic was changed
by this review. The three new tests also pass with Apple Clang 21
AddressSanitizer, including stack-use-after-return detection. Only the test
translation unit (and its header kernels) was instrumented; linked library
objects were the existing Release build. The MacPorts Clang 19 runtime stalled
in sanitizer initialization before main; sampling identified recursive allocator
initialization. Recompiling and linking with the Apple runtime removed that
startup problem.

Two existing assembled-form algebra limits were found by compiling candidate
compositions: `A - B` and `(A + B) + C` for named forms fail to compile.
The searches included assembled-form overloads in `Sum.h`, `Minus.h`,
`UnaryMinus.h`, `ProblemBody.h` and the whole Variational tree. These limits
also affect general BilinearFormBase expressions and predate this PR. Signed
coefficients and supported two-form composition are tested. Extending this
algebra is separate work; unrestricted integral-style composition is not claimed.

## Protocol and numerical gates

Apple M1 Pro, Clang 19.1.7, C++20, Release `-O3 -DNDEBUG`, local Eigen assembly,
multithreading disabled and PETSc disabled. This does not certify distributed,
PETSc or genuinely concurrent assembly. Two runs, three randomized repetitions
per case, 20 ms minimum for assembly and 30 ms for field operations. No
compilation or other benchmark ran concurrently. The host is not dedicated.

`RodinNamedFormBenchmarks` contains 160 cases: P1/P2, 2D/3D, five operators,
constant/field coefficients, named/integral paths, warm/first assembly. The
2D mesh has 512 triangles; the 3D mesh has 384 tetrahedra. Coefficient fields
use the same order as the solution and interpolate `1 + sum((d+1)*x_d^2)`;
P2 reproduces this polynomial. Constants have value 2. Each complete sparse
matrix is compared to Integral before timing and to a frozen Integral matrix
after timing, with relative tolerance 1e-11. Across 960 timed observations,
the maximum matrix-relative difference was 3.57e-16.

Warm means reassembly on an existing form; named local kernels are freshly
created by each assembly, so their reference tables are prepared again.
First includes construction, automatic named assembly (explicit integral
assembly), and destruction. Mesh quadrature is already warm from the oracle;
these are not fully cold mesh runs. First retains the last matrix by swap for
the final numerical check; the same retention cost applies to both paths.
Setup, oracle checks, mesh/space construction and field interpolation are
outside the timed loop. Assembly times include kernel preparation, local
matrices and scatter. They are not bare coefficient evaluation timings.

## P2 complete assembly

Median CPU milliseconds, pooling six observations. Ratio is Integral / Named;
values below one mean the named path is slower. Mesh sizes differ by dimension,
so times should only be compared within one dimension.

| Operator | Coefficient | Dim | Named warm ms | Integral warm ms | Ratio | Named first ms | Integral first ms |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Mass | Constant | 2 | 0.096 | 0.261 | 2.72 | 0.295 | 0.262 |
| Mass | Constant | 3 | 0.164 | 0.633 | 3.86 | 0.756 | 0.625 |
| Mass | FieldCoefficient | 2 | 0.541 | 0.782 | 1.45 | 0.745 | 0.767 |
| Mass | FieldCoefficient | 3 | 2.929 | 3.607 | 1.23 | 3.538 | 3.595 |
| Diffusion | Constant | 2 | 0.181 | 0.276 | 1.52 | 0.390 | 0.277 |
| Diffusion | Constant | 3 | 0.454 | 0.599 | 1.32 | 1.055 | 0.611 |
| Diffusion | FieldCoefficient | 2 | 0.580 | 0.578 | 1.00 | 0.784 | 0.577 |
| Diffusion | FieldCoefficient | 3 | 2.209 | 2.188 | 0.99 | 2.860 | 2.218 |
| Helmholtz | Constant | 2 | 0.268 | 0.528 | 1.97 | 0.471 | 0.532 |
| Helmholtz | Constant | 3 | 0.589 | 1.233 | 2.09 | 1.210 | 1.226 |
| Helmholtz | FieldCoefficient | 2 | 1.112 | 1.349 | 1.21 | 1.305 | 1.341 |
| Helmholtz | FieldCoefficient | 3 | 4.962 | 5.869 | 1.18 | 5.480 | 5.830 |
| VectorMass | Constant | 2 | 0.183 | 0.490 | 2.68 | 1.061 | 0.492 |
| VectorMass | Constant | 3 | 1.367 | 1.856 | 1.36 | 7.006 | 1.834 |
| VectorMass | FieldCoefficient | 2 | 0.731 | 1.041 | 1.42 | 1.641 | 1.051 |
| VectorMass | FieldCoefficient | 3 | 5.230 | 5.064 | 0.97 | 11.034 | 4.970 |
| Elasticity | Constant | 2 | 2.686 | 18.890 | 7.03 | 3.591 | 18.846 |
| Elasticity | Constant | 3 | 37.391 | 314.948 | 8.42 | 42.902 | 314.546 |
| Elasticity | FieldCoefficient | 2 | 5.291 | 28.411 | 5.37 | 6.307 | 28.657 |
| Elasticity | FieldCoefficient | 3 | 79.575 | 636.275 | 8.00 | 86.152 | 633.385 |

P2 named elasticity is about 7–8.4x faster for constant coefficients and
5.4–8x faster for field coefficients. Named mass and Helmholtz benefit less
with field coefficients. Field-weighted diffusion is effectively tied;
field-weighted 3D vector mass is about 3% slower in these runs. The benchmark
therefore does not establish a universal named-form speedup.

First assembly has an explicit setup tradeoff: constant-coefficient 3D P2
vector mass takes 7.006 ms for named construction/assembly versus 1.834 ms
for Integral, despite faster named reassembly. The ScatterMap first builds
triplets and the full topological pattern, then binary-searches that pattern
for every local row/column pair to construct its reusable index array.
In this fixture that is 345600 scatter slots. A follow-up counter run confirms
153225 stored named entries versus 51075 Integral entries: the named pattern
retains cross-component zeros, whereas Integral filters zero entries. The
extra initial index construction and larger pattern are concrete setup work;
this cost must be included when deciding whether repeated assembly amortizes
it. The two-case storage check also reproduced the roughly 3.8x first-assembly
ratio. Its raw records are in
[the storage CSV](results/named-storage-2026-10-01.csv).

A P2 scalar field increases the requested mass quadrature order from 4 to 6
(6 to 11 triangle points; 11 to 23 tetrahedron points). Diffusion increases
from order 2 to 4 (3 to 6 triangle points; 4 to 11 tetrahedron points).
Elasticity follows its existing integrator's order convention: order 4 for
constants, 6 for P2 fields. Helmholtz uses each term's own quadrature order.
Both named and Integral paths evaluate a scalar coefficient once per quadrature
sample. More samples, field evaluation and kernel work all contribute to the
constant/field difference; the difference is not a pure interpolation ablation.

The vector mass comparison has a concrete structural difference: Integral's
H1 weighted mass specialization tabulates scalar bases and accumulates only
matching component blocks. Named mass builds dense reference matrices and
adds all entries, including cross-component zeros. A P2 tetrahedron has 30
vector basis entries, so its reference matrices have 900 entries rather than
300 nonzero block entries. Increasing to 23 quadrature samples raises this
dense table to 165600 bytes and its accumulation to 20700 entries per cell.
This explains additional work in that path; attributing the exact 3% timing
difference would require a separate kernel experiment. No such optimization
was bundled into this review.

## P2 field evaluation

The existing interpolation benchmark was rerun on the retained exact-coordinate
implementation: 56 P2 cases, six observations each, no numerical failures.
Large immutable meshes have 32768 triangles or 24576 tetrahedra. Sample-wise
analytic checks precede timing; timed checksums provide a second numerical
gate. See [GridFunctionLarge.md](GridFunctionLarge.md) for workload definitions.
The largest relative checksum error was 2.09e-10;
the largest absolute error in the smaller interpolation cases was 9.27e-15.

Median CPU ns/evaluation, pooling six observations:

| Field | Dim | Workload | Quadrature cache | Direct pointwise |
| --- | ---: | --- | ---: | ---: |
| P2Scalar | 2 | SlowSweep | 82.1 | 74.3 |
| P2Scalar | 2 | Mixed75 | 33.5 | 68.3 |
| P2Scalar | 2 | FastBlocks | 18.2 | 64.9 |
| P2Scalar | 2 | PureHit | 15.6 | 64.1 |
| P2Scalar | 3 | SlowSweep | 283.6 | 274.9 |
| P2Scalar | 3 | Mixed75 | 84.4 | 279.1 |
| P2Scalar | 3 | FastBlocks | 25.3 | 271.2 |
| P2Scalar | 3 | PureHit | 17.3 | 261.5 |
| P2Vector | 2 | SlowSweep | 174.4 | 80.2 |
| P2Vector | 2 | Mixed75 | 66.1 | 74.1 |
| P2Vector | 2 | FastBlocks | 34.6 | 70.3 |
| P2Vector | 2 | PureHit | 28.9 | 70.5 |
| P2Vector | 3 | SlowSweep | 894.4 | 288.7 |
| P2Vector | 3 | Mixed75 | 257.7 | 277.6 |
| P2Vector | 3 | FastBlocks | 79.1 | 270.2 |
| P2Vector | 3 | PureHit | 53.7 | 265.2 |

The scalar 3D all-miss cost is about 284 ns/evaluation. The vector 3D all-miss
cost is about 894 ns versus 289 ns for direct expansion; with 75% basis hits
it drops to 258 ns and with pure hits to 54 ns. Cache refill evaluates each
vector basis separately, repeating scalar basis evaluation per component.
`H1Element::evaluate` instead evaluates each scalar basis once and reuses it
across components. This is the measured P2 vector miss-path concern; retaining
exact coordinate validation is compatible with addressing it separately.

Small-fixture operation costs (microseconds per complete operation):

| Field | Dim | Re-interpolate | Cross-space assignment | First interpolation |
| --- | ---: | ---: | ---: | ---: |
| P2Scalar | 2 | 4.72 | 9.18 | 4.80 |
| P2Scalar | 3 | 82.55 | 438.81 | 82.36 |
| P2Vector | 2 | 9.68 | 18.38 | 9.85 |
| P2Vector | 3 | 283.33 | 1331.30 | 286.74 |

These are nodal interpolation operations, not L2 projection. Fixture sizes and
cross-space behavior are specified in
[GridFunctionInterpolation.md](GridFunctionInterpolation.md); these numbers
cannot be extrapolated to a large full-field interpolation by treating them as
single-point costs.

## Reproduction

```sh
cmake --build BUILD --target RodinNamedFormBenchmarks RodinGridFunctionBenchmarks
BUILD/tests/benchmarks/RodinNamedFormBenchmarks --benchmark_min_time=0.02s --benchmark_repetitions=3 --benchmark_enable_random_interleaving=true
BUILD/tests/benchmarks/RodinGridFunctionBenchmarks --benchmark_filter='.*P2.*' --benchmark_min_time=0.03s --benchmark_repetitions=3 --benchmark_enable_random_interleaving=true
```

Run each command twice, sequentially. Raw iteration records are in
[the assembly CSV](results/named-evaluation-2026-10-01.csv) and
[the P2 evaluation CSV](results/p2-evaluation-2026-10-01.csv). P1 controls,
first assembly costs and every repeat are retained there. The older
`NamedForms.cpp` benchmarks remain available, but their nonzero-count checks
alone are not numerical equivalence oracles.
