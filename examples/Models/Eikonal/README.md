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

`Rodin/Eikonal/FIM.h` supplies a local fast iterative solver with optional OpenMP and
an output `GridFunction<P1<Real, LocalMesh>>`. Its calling pattern matches
fast marching:

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

A scalar coefficient is a positive propagation speed. A matrix-valued Rodin
function or callable supplies the symmetric positive definite squared-velocity
tensor. For a constant tensor in two dimensions, replace the solver construction
with:

```cpp
Math::SpatialMatrix<Real> tensor(2, 2);
tensor(0, 0) = 4;
tensor(0, 1) = tensor(1, 0) = 1;
tensor(1, 1) = 1;
MatrixFunction coefficient(tensor);
Eikonal::FIM solver(arrival, coefficient);
solver.seed(seeds).solve();
```

The scalar and tensor equations are respectively

\[
c|\nabla T|=1,
\qquad
\sqrt{\nabla T^T M\nabla T}=1.
\]

Every source receives time zero; seeding resets all other times to positive
infinity. Duplicate source indices are harmless, invalid or empty source lists
raise, and disconnected components without sources remain unreachable. Reseed
after changing geometry, coefficient data, or numerical controls. The output
field and its mesh must outlive the solver. Function expressions and callables
are owned by value; borrowed data captured by a callable must remain valid.

### Native geometry and curved paths

The implementation covers point, segment, triangle, quadrilateral, tetrahedron,
hexahedron, wedge, and pyramid reference types, including mixed meshes. It uses
the existing native P1 traces on boundary points, edges, triangles, and
quadrilaterals. No cell is subdivided into simplices. The arrival-time space
remains first order even when the geometric transformation has higher order.

For a native cell map \(X_K\), a target reference vertex \(r_i\), and a boundary
entry \(r_F\), the admissible path is

\[
x(s)=X_K\bigl(r_i+s(r_F-r_i)\bigr),\qquad 0\le s\le1.
\]

The local candidate is the native face interpolation of the entry time plus

\[
\int_0^1\frac{|\dot x(s)|}{c(x(s))}\,ds
\quad\text{or}\quad
\int_0^1\sqrt{\dot x(s)^T M(x(s))^{-1}\dot x(s)}\,ds.
\]

Geometry and coefficients are evaluated on the owning cell trace through
`Geometry::Point`; this retains discontinuous material traces and the cell's
actual Jacobian. Reference cells are convex, so the mapped reference segments
stay within their cells. On a curved embedded mesh these are mapped paths,
rather than physical endpoint chords. Valid regular maps and conforming
interfaces are preconditions.

A tensor on an embedded mesh has the ambient embedding size and defines an
ambient path metric whose restriction is used on tangents. A tangential velocity
tensor can be completed by a decoupled positive normal block. Coupled tangent
and normal entries instead define a different induced metric and must be chosen
intentionally.

### Approximation controls and validation

`setTolerance` controls absolute accepted time decreases. `setLocalTolerance`
controls reference-coordinate search resolution; `setIntegrationOrder` controls
positive segment quadrature along each path. `setMaxIterations` limits
active-list rounds; exhaustion raises while retaining the last iterate.
`getIterations` reports the last solve's round count.

Native face minimization uses a deterministic multistart projected pattern
search. This is an approximate semi-Lagrangian discretization, not a certified
global optimizer on arbitrary curved maps or heterogeneous coefficients.
Reducing the search resolution does not certify discovery of every nonconvex
minimum. Neither quadrature order nor active-list convergence alone establishes
arrival-time accuracy. Check path quadrature, search resolution, and spatial
refinement for the intended geometry and anisotropy. There is no GPU,
MPI, PETSc, or high-order arrival-space specialization in this implementation.

Build and run `RodinFIMTest` and `RodinManufacturedEikonalFIMTest`. The tests
exercise every reference type, native face-interior minima, scalar/tensor
reduction, quadratic cell maps, embedded curved paths, mixed cells with
material jumps, invalid inputs, reseeding, and convergence of constant-tensor
point-source fields. The point-source check integrates field error over the
fixed physical domain. A nodal average can increase on the coarsest meshes
because the changing sample set overweights the source-cell stencil, where
direct source-to-vertex paths are exact.

### OpenMP execution

Configure Rodin with `RODIN_USE_OPENMP=ON`, then opt in before solving:

```cpp
Eikonal::FIM solver(arrival, coefficient);
solver.setThreadCount(4).seed(seeds).solve();
```

The default worker count is one. A build without OpenMP rejects a request for
more than one worker. Dependency levels preserve the serial active-list order
between all vertices sharing a cell. Updates within a level have no co-cell
read/write conflicts and may run concurrently. Next-round vertices are collected
in the original order after workers join. Deterministic coefficients therefore
produce the same values and iteration counts as the serial algorithm; this is
checked by exact field comparisons, including curved tensor cases.

Each worker owns a copy of the coefficient expression. External data captured
by a callback must allow concurrent reads and must remain fixed during a solve.
Worker exceptions are rethrown after joining. Level construction, barriers, and
coefficient copying have a cost, so additional workers need not improve small
problems. This implementation preserves Gauss-Seidel dependencies rather than
changing the iteration to simultaneous Jacobi updates.

### Tensor FMM and matched comparisons

FMM also accepts the squared-velocity tensor through the same coefficient API:

```cpp
Eikonal::FMM solver(arrival, coefficient);
solver.seed(seeds).solve();
```

The tensor overload freezes the coefficient at each cell centroid. It requires
native first-order affine simplices, nonobtuse in the inverse-velocity metric;
unsupported cells raise. Matrix dimensions, symmetry, and positive definiteness
are checked by the shared `Eikonal::Metric` implementation. FMM uses metric Gram
entries in its analytic simplex updates and metric lengths in edge updates.
The scalar-speed path retains its original arithmetic and behavior.

The metric angle condition is sufficient for causal acceptance. For an upwind
simplex update, let the Gram matrix be \(G\) and the supporting-time vector be
\(u\). An admissible update satisfies

\[
t\mathbf 1-u=G\mu,\qquad \mu\ge0.
\]

Nonnegative Gram entries ensure that every supporting time precedes the
candidate. FIM can revisit values and does not impose this FMM mesh condition.
Its approximate local face search still needs its own accuracy assessment.
Anisotropy alone is therefore not a reason to choose FIM; native/curved geometry,
mesh restrictions, and parallel execution are separate considerations.

Build `RodinEikonalPropagationBenchmarks` with `RODIN_BUILD_BENCHMARKS=ON` to
compare seeded propagation on matched constant-tensor triangle and tetrahedron
meshes. For example:

```sh
cmake --build build --target RodinEikonalPropagationBenchmarks -j
./build/tests/benchmarks/RodinEikonalPropagationBenchmarks \
  --benchmark_filter='Eikonal/(FMM|FIM)/2/32/2/' \
  --benchmark_min_time=0.01s --benchmark_repetitions=5 \
  --benchmark_out=propagation.json --benchmark_out_format=json
```

Case arguments are dimension, reference-grid intervals, velocity ratio, and
requested workers. A rotation and stretch make the mesh nonobtuse in the tensor
metric. Different velocity ratios change the physical domain stretch while
holding the reference travel-time domain fixed. Within each case, FMM and FIM
use identical physical meshes, coefficients, and sources. Timing includes
seeding and solving; setup, reference solves, and error checks are excluded.
FIM quadrature order and search resolution, cells, vertices, analytic maximum
errors, differences from FMM, and iteration counts are reported as counters.

This compares the complete current implementations. FMM uses analytic simplex
updates, whereas native FIM performs numerical face minimization and path
quadrature. Runtime differences include that local-update work as well as queue
versus active-list scheduling. FMM can be substantially faster on these simple
causal meshes, even when OpenMP accelerates FIM. The comparison is not evidence
that FIM is intrinsically slower or that either method is more accurate on
arbitrary anisotropic or curved domains.
