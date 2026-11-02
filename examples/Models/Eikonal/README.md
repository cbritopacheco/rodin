# Arrival times with fast marching

`FMM.cpp` solves a point-source problem on the unit square with constant
positive speeds 1 and 2. It uses 16, 32 and 64 intervals per axis, placing an
exact mesh vertex at the center. The equation and analytic reference are

\[
F|\nabla T|=1,\qquad T(x)=\frac{|x-(0.5,0.5)|}{F}.
\]

Build `Eikonal_FMM` and run it in a scratch directory. No input data are
required. `FMMErrors.csv` records maximum and mean nodal errors for every
case. `FMMArrival.xdmf` and its HDF5 companions contain the finest unit-speed
computed and analytic nodal fields. All cases run sequentially.

The executable checks finite nonnegative values, a unique source, decreasing
maximum error with refinement and the loose regression envelope described in
the source. The observed unit-speed maximum errors are 0.0474507, 0.0299216
and 0.0182986. Speed 2 halves these errors. These measurements apply to the
generated triangular grids and central point source.

Render the actual field and contours with ParaView.

```sh
pvpython examples/Models/Eikonal/visualize_paraview.py \
  /path/to/results doc/Images/FMMArrival.png
```

Titles use the bundled Manrope font, mathematical labels use TeX, and the
background is transparent. Contours clarify the front without dense mesh edges.

## Native-cell fast iteration

`Rodin/Eikonal/FIM.h` supplies a local real P1 arrival-time solver on native
first-order straight cells, with optional OpenMP execution:

```cpp
#include <Rodin/Eikonal/FIM.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// mesh and seeds are an existing local mesh and source-vertex container.
const size_t dimension = mesh.getDimension();
mesh.getConnectivity().compute(dimension, 0);
mesh.getConnectivity().compute(0, dimension);
P1 space(mesh);
GridFunction arrival(space);

const auto speed = [](const Point&) { return Real(1); };
Eikonal::FIM solver(arrival, speed);
solver.seed(seeds).solve();
```

Scalar coefficients specify positive speeds; matrix coefficients supply ambient
symmetric positive definite squared-velocity tensors. For example:

```cpp
Math::SpatialMatrix<Real> tensor(2, 2);
tensor(0, 0) = 4;
tensor(0, 1) = tensor(1, 0) = 1;
tensor(1, 1) = 1;
MatrixFunction coefficient(tensor);
Eikonal::FIM solver(arrival, coefficient);
solver.seed(seeds).solve();
```

The equations are

\[
c|\nabla T|=1,
\qquad
\sqrt{\nabla T^T M\nabla T}=1.
\]

Every source receives time zero; seeding resets other times to positive
infinity. Invalid or empty source lists raise; duplicate indices are harmless.
Disconnected components without sources remain unreachable. The field and mesh
must outlive the solver. Coefficient expressions and callables are owned by
value; borrowed data captured by a callable must remain valid.

### Cellwise coefficients and straight travel costs

Before a pending solve begins iteration, the coefficient is evaluated once at
each cell's reference centroid, using a cell-associated `Geometry::Point`.
Scalar speeds are validated and tensor metrics are factored once. These fixed
cellwise inputs are reused by every local search, active-list revisit, and
OpenMP worker. No coefficient evaluation occurs inside the face minimization.

For a cell \(K\) with native map \(X_K\) and reference centroid \(\widehat x_K\),
the prescribed input discretization is

\[
c_K=c(X_K(\widehat x_K)),
\qquad M_K=M(X_K(\widehat x_K)).
\]

This is an explicit piecewise-constant approximation for spatially varying
coefficients. Discontinuous coefficients are evaluated on the owning cell
trace. It does not integrate within-cell coefficient variation along a path.
Reseed after changing geometry or borrowed coefficient data; the next pending
solve rebuilds all cell samples. To replace the owned coefficient expression,
construct a new solver. Calling `solve()` again with no pending active vertices
leaves the converged field unchanged.

A native boundary entity \(F\), not containing target vertex \(x_i\), supplies
an interpolated entry time and physical entry point. The local candidate is

\[
\min_{q\in\widehat F}
\left[
  \sum_j N_j(q)T_j+
  \sqrt{d(q)^T M_K^{-1}d(q)}
\right],
\qquad d(q)=\sum_jN_j(q)x_j-x_i.
\]

The scalar travel cost is \(|d|/c_K\). Physical vertices are evaluated through
`Geometry::Point`, retaining the supplied first-order chart. Native point,
segment, triangle, and quadrilateral traces are used directly; no cell is
subdivided into simplices. Point, edge, and triangular trace minima are computed
analytically, including constrained boundary minima. Candidates may precede
some supporting times; FIM imposes no marching-order cutoff. Quadrilateral
traces retain multistart projected search of their native bilinear interpolant.
There is no path quadrature or integration-order setting.

Point, segment, triangle, quadrilateral, tetrahedron, hexahedron, wedge, and
pyramid reference types are covered, including mixed meshes. The input geometry
must be native first order, straight, convex, and have planar faces; embedded
surface cells must be planar. Higher-order geometric charts are rejected,
including higher-order representations of straight cells. First-order
multilinear geometry alone does not certify validity, convexity, or planar
faces: these are caller preconditions. Curved cells and high-order arrival-time
spaces are unsupported.

On embedded planar cells, tensor dimensions match the ambient embedding. The
restriction of the ambient inverse-velocity metric is used on tangential
displacements. A tangential velocity tensor may be completed with a decoupled
positive normal block.

### Iteration and OpenMP

`setTolerance` controls absolute accepted time decreases. `setLocalTolerance`
controls the reference-coordinate resolution of the multistart projected
quadrilateral search; it has no effect on analytic point, edge, or triangular
updates. `setMaxIterations` limits active-list rounds; exhaustion raises and
retains the last iterate. `getIterations` reports the last solve's round count.
Face search remains approximate, without a global-minimum certificate for
arbitrary native quadrilateral traces. Search resolution, spatial refinement,
and coefficient discretization must be assessed for the intended problem.

Configure Rodin with `RODIN_USE_OPENMP=ON`, then opt in:

```cpp
Eikonal::FIM solver(arrival, coefficient);
solver.setThreadCount(4).seed(seeds).solve();
```

The default worker count is one. Builds without OpenMP reject requests for more
workers. Dependency levels preserve the serial Gauss-Seidel order between
vertices sharing a cell, and next-round vertices are collected in the original
order. Serial and parallel executions therefore retain identical fields and
iteration counts for the same sampled inputs.

Coefficient callbacks and tensor validation run serially before iteration.
Workers only read the prepared cell coefficients and factors. No worker-local
coefficient copies are needed. Worker exceptions from local updates are rethrown
after joining. Level construction and barriers have a cost, so additional
workers need not improve small problems. GPU, MPI, and PETSc implementations
are not supplied.

Build and run `RodinFIMTest` and `RodinManufacturedEikonalFIMTest`. Tests cover
every reference type, native face-interior minima, scalar/tensor reduction,
embedded planar cells, mixed material traces, exactly one coefficient sample
per cell, reseeding after coefficient changes, rejection of higher-order
geometry, and exact serial/parallel agreement. Manufactured checks cover plane
waves, constant-tensor point-source refinement, and refinement of a cellwise
variable-speed approximation.

The point-source error is integrated over a fixed physical domain. Nodal
averages alone can increase on very coarse meshes because the changing sample
set overweights the exact direct-source stencil.

### Tensor FMM and matched comparisons

FMM accepts the same scalar-speed or tensor coefficient interface:

```cpp
Eikonal::FMM solver(arrival, coefficient);
solver.seed(seeds).solve();
```

The tensor overload freezes the coefficient at each cell centroid. It requires
native first-order affine simplices nonobtuse in the inverse-velocity metric;
unsupported cells raise. Shared `Eikonal::Metric` validation checks matrix
dimensions, symmetry, and positive definiteness. FMM uses metric Gram entries
in analytic simplex updates. The original scalar-speed path is unchanged.

For an upwind simplex update with Gram matrix \(G\) and supporting times \(u\),

\[
t\mathbf 1-u=G\mu,\qquad\mu\ge0.
\]

Nonnegative Gram entries ensure that supporting times precede the candidate.
FIM revisits values and does not impose this FMM metric-angle condition.
Anisotropy, native geometry, and parallel execution are separate reasons to
assess the methods for a particular problem.

Enable `RODIN_BUILD_BENCHMARKS` and `RODIN_USE_OPENMP`, then run the registered
matched comparison:

```sh
cmake --build build --target RodinEikonalPropagationBenchmarks -j
./build/tests/benchmarks/RodinEikonalPropagationBenchmarks \
  --benchmark_filter='Eikonal/(FMM|FIM)/2/32/2/' \
  --benchmark_min_time=0.03s --benchmark_repetitions=5 \
  --benchmark_out=propagation.json --benchmark_out_format=json
```

Case arguments are dimension, reference-grid intervals, velocity ratio, and
requested workers. A rotation and stretch make the affine simplex mesh
nonobtuse in the constant tensor metric. Different velocity ratios change the
physical stretch while fixing the reference travel-time domain. Within each
case the physical mesh, coefficients, and sources are identical for both
methods. Timing includes coefficient preparation, seeding, and iteration;
setup, reference solves, and error checks are excluded.

Counters report cells, vertices, coefficient
evaluations per solve (counted in untimed reference solves), analytic maximum
error, difference from FMM, and iteration count. The benchmark checks exactly
one FIM coefficient evaluation per cell and exact worker-count invariance.
Both methods use analytic local updates on these simplex meshes. FMM evaluates
and factors tensors during local updates, whereas FIM caches one factor per
cell. Runtime differences therefore include coefficient preparation, local
update implementation, and queue versus active-list scheduling. These cases do
not measure the quadrilateral search used on other native geometries.
