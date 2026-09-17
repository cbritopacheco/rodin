# Integrator benchmarks

The `RodinBenchmarks` target contains two complementary integrator benchmark
levels under the `Integrator/` prefix.

- `Integrator/Kernel/` measures element binding, quadrature, basis evaluation,
  coefficient evaluation, and local matrix or vector entry generation. It does
  not include sparse assembly.
- `Integrator/Assembly/` measures the same variational expressions through the
  standard Rodin form and assembly interfaces. It therefore includes mesh
  traversal, local-to-global scatter, and sparse-operator construction.
  `Cold` cases construct a fresh form and operator on every iteration, whereas
  `Reassembly` and boundary cases perform one unmeasured warm-up assembly and
  measure reuse of the existing form and assembly storage.

The distinction is required when assessing a `QuadratureRule` optimization. A
local-kernel improvement may be hidden by sparse assembly, while an assembly
improvement may be unrelated to the integrator itself.

## Coverage

The benchmark matrix includes the following variables that materially affect
integrator cost:

- finite-element family (`P1` and `H1`), polynomial degree, and mixed trial and
  test degrees;
- scalar and vector ranges;
- two- and three-dimensional meshes;
- affine and degree-two parametric transformations;
- segment, triangle, quadrilateral, tetrahedron, hexahedron, wedge, and pyramid
  geometries;
- constant, analytic coordinate-dependent, matrix-valued, and finite-element
  coefficients;
- volume, boundary, local, and global integrators;
- local basis count, quadrature point count, and resulting local entry count;
- complete assembly, whose threading is controlled by the normal Rodin and
  `OMP_NUM_THREADS` configuration.

Every local specialization in `Variational/P1/QuadratureRule.h` and
`Variational/H1/QuadratureRule.h` is represented by one or more of the
following expressions:

- basis and source loads;
- mass, weighted mass in both supported expression orderings, anisotropic
  vector mass, and finite-element-coefficient mass;
- gradient-gradient and weighted gradient-gradient forms;
- Jacobian-Jacobian and weighted Jacobian-Jacobian forms;
- divergence-pressure forms in both trial-test orientations;
- Jacobian-advection forms;
- the P1 global potential form.

Generic controls include divergence-divergence and symmetric-gradient
elasticity, flux-gradient loads, and P0 mass and source forms. Boundary mass
and source assembly provide non-volume controls. Reaction-diffusion provides a
common multi-integrator assembly case, while analytic and finite-element
transport fields distinguish expression structure from coefficient-evaluation
cost.

## Running

Build and list the repertoire with:

```sh
cmake --build build -j --target RodinBenchmarks
build/tests/benchmarks/RodinBenchmarks \
  --benchmark_list_tests | grep '^Integrator/'
```

Run only local kernels or complete assembly with:

```sh
build/tests/benchmarks/RodinBenchmarks \
  --benchmark_filter='Integrator/Kernel/.*'

OMP_NUM_THREADS=8 build/tests/benchmarks/RodinBenchmarks \
  --benchmark_filter='Integrator/Assembly/.*'
```

For a reproducible comparison, use the same build type, compiler, thread
count, CPU frequency policy, and benchmark filter on both revisions, and emit
Google Benchmark JSON:

```sh
build/tests/benchmarks/RodinBenchmarks \
  --benchmark_filter='Integrator/.*' \
  --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true \
  --benchmark_out=integrators.json \
  --benchmark_out_format=json
```

Compare two saved runs with the Google Benchmark comparison tool vendored in
the repository:

```sh
python3 third-party/benchmark/tools/compare.py benchmarks \
  integrators-before.json integrators-after.json
```

Kernel cases report processed local entries per second. Assembly cases report
processed cells per second and the assembled nonzero count. Comparisons across
different polynomial degrees should use both elapsed time and throughput:
elapsed time describes application cost, whereas entry throughput describes
the efficiency of the integrator as its mathematically required local work
increases.
