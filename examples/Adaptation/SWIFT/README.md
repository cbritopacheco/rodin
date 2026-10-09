# Reconstruction

These self-contained examples demonstrate the Strain-distributed Welsch
Implicit-interface Fitting Technique (SWIFT). A classified mesh interface is
fitted to a circle in two dimensions or a sphere in three dimensions. The
mesh topology is retained; the outer box boundary is held fixed.

## Examples

| Source | Executable | Usage |
|--------|------------|-------|
| `ReconstructionP1.cpp` | `SWIFT_ReconstructionP1` | `SWIFT::Adapt` owns the solve and applies valid vertex displacement. |
| `ReconstructionP2.cpp` | `SWIFT_ReconstructionP2` | `SWIFT::Problem` writes displacement; the application updates quadratic transformations. |
| `ReconstructionP3.cpp` | `SWIFT_ReconstructionP3` | `SWIFT::Problem` writes displacement; the application updates cubic transformations. |

All three examples support both spatial dimensions. The polynomial degree
controls the displacement space and resulting geometric transformations,
not the target geometry.

## Build

Run the following commands from the repository root. Initialize the submodules
for a fresh checkout and configure a local build with examples enabled:

```sh
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DRODIN_BUILD_EXAMPLES=ON -DRODIN_USE_PETSC=OFF \
  -DRODIN_INSTALL_RESOURCES=OFF
cmake --build build --target SWIFT_ReconstructionP1 SWIFT_ReconstructionP2 SWIFT_ReconstructionP3 -j1
```

The standard Rodin dependencies, including Boost, Eigen and HDF5, must be
available to CMake. These examples do not require downloaded mesh resources
or PETSc: they construct local meshes and use Eigen-backed displacement
storage. See the repository's [build workflow](../../../doc/agents/workflows.md)
for general build procedures.

For optional MUMPS support, add `-DRODIN_USE_MUMPS=ON` when configuring, with
MUMPS and its dependencies installed. SWIFT selects MUMPS when compiled in;
otherwise its default direct solver is SparseLU. Add `-DRODIN_USE_OPENMP=ON`
for an OpenMP-enabled build. Thread environment variables do not enable
features absent from the configured build.

An existing configured build can be reused by replacing `build` in the
commands. Only the three named targets need to be built.

## Run

The executable syntax is:

```text
SWIFT_ReconstructionP1 [grid-points-per-axis] [dimension]
SWIFT_ReconstructionP2 [grid-points-per-axis] [dimension]
SWIFT_ReconstructionP3 [grid-points-per-axis] [dimension]
```

| Argument | Meaning | Default | Accepted values |
|----------|---------|---------|-----------------|
| First | Number of grid points per axis, not the number of cells | `16` | Integer at least `4` |
| Second | Spatial dimension | `2` | `2` or `3` |

For example, run from a scratch directory so generated files remain outside
the source tree. The following commands assume the shell initially starts at
the repository root:

```sh
build_dir="$(pwd)/build"
mkdir -p /tmp/rodin-swift
cd /tmp/rodin-swift

# Reconstruct a circle with linear, quadratic and cubic displacement.
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP1" 16 2
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP2" 16 2
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP3" 16 2

# Reconstruct a sphere on a smaller three-dimensional grid.
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP1" 8 3
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP2" 8 3
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP3" 8 3
```

Three-dimensional meshes and higher-order spaces require more memory and
factorization work. Start with a small grid. Runs of the same degree in the
same working directory overwrite that degree's output; use separate working
directories to retain several resolutions or dimensions.

Only the two positional arguments are parsed. Model and solver controls are
changed through the C++ parameter object, not command-line flags.

## Geometry and Workflow

The background domain is the unit square or cube, triangulated or
tetrahedralized from a uniform grid. Its fixed reference spacing is

\[
h=\frac{1}{n-1}.
\]

The target is the zero level set of the signed distance

\[
\phi(x)=\lVert x-c\rVert-R,
\qquad c=(0.5,\ldots,0.5),\qquad R=0.25.
\]

Cells are classified by the sign of the level set at their centroids. Facets
separating the two cell attributes form the interface to fit. This classified
interface is only an approximation to the target; no intersection cutting or
remeshing is performed.

| Attribute | Value | Role |
|-----------|-------|------|
| `Inside` | `1` | Cells with negative centroid level-set value |
| `Outside` | `2` | Remaining cells |
| `Interface` | `10` | Facets separating inside and outside cells |
| `Boundary` | `20` | Exterior box facets, with zero displacement increments |

Each example prepares mesh connectivity, constructs the level set and its
analytic gradient, marks the facets, solves the fit, checks the quality
report, and writes the resulting geometry.

## API Usage and Ownership

The linear example uses `SWIFT::Adapt`, which owns the displacement space and
problem and applies a quality-valid displacement to the supplied mesh:

```cpp
// mesh, phi, gradient and facet attributes are prepared beforehand.
Adaptation::SWIFT::Adapt fitting(mesh);
Adaptation::SWIFT::Parameters parameters;
parameters.model.h = h; // Keep the background reference scale fixed.
fitting.setParameters(parameters).setInterfaceAttribute(Interface);

// Keep the vector alive while the boundary expression is used.
const Math::Vector<Real> zero = Math::Vector<Real>::Zero(dimension);
auto& u = fitting.getTrialFunction();
fitting.getProblem() += DirichletBC(u, VectorFunction(zero)).on(Boundary);
const auto report = fitting.execute(phi, gradient); // Updates mesh when valid.
```

The quadratic and cubic examples use `SWIFT::Problem` directly. For example:

```cpp
// Choose degree 2 for quadratic displacement; use 3 for cubic displacement.
H1 space(std::integral_constant<size_t, 2>{}, mesh, dimension);
TrialFunction u(space);
TestFunction v(space);
Adaptation::SWIFT::Problem fitting(u, v);
Adaptation::SWIFT::Parameters parameters;
parameters.model.h = h;
fitting.setParameters(parameters).setInterfaceAttribute(Interface);
const Math::Vector<Real> zero = Math::Vector<Real>::Zero(dimension);
fitting += DirichletBC(u, VectorFunction(zero)).on(Boundary);
const auto report = fitting.solve(phi, gradient);
const auto& displacement = u.getSolution(); // mesh has not been modified.
```

The mesh, space, trial and test functions must outlive the problem. The
supplied gradient must differentiate the level set; its Hessian is not
required. Homogeneous Dirichlet conditions constrain each increment.

After a valid solve, the higher-order examples create a separate mesh and
evaluate the displacement at the original mesh's geometry nodes. They update
both vertices and all positive-dimensional entity transformations with
`ParametricTransformation`. Moving corner vertices alone would discard the
quadratic or cubic geometry. The original mesh remains the reference mesh for
displacement evaluation throughout this step.

## Model and Parameters

The canonical displacement metric is

\[
M=F+D.
\]

Fitting observes motion along the target gradient. Distribution penalizes
variation of current strain about its **global domain average**, with separate
deviatoric and divergence weights. It is not centered independently in each
element. Compatible global affine motions are free in the distribution term;
shared finite-element degrees of freedom couple neighboring elements.

The force and outer merit use robust Welsch fitting energy. Directional Newton
scales the quadratic model before affine quadratic quality hinges recover the
direction. Armijo backtracking accepts a trial only when energy decrease and
the actual sampled Jacobian/distortion checks both pass. Quality is not a
separate shape term in the metric.

The examples change only `model.h`; all other values below are production
defaults. Edit the parameter object before `setParameters()` to change them.

| C++ parameter | What it controls | Default |
|---------------|------------------|---------|
| `model.h` | Fixed background reference spacing | Required; set to \(1/(n-1)\) |
| `model.fit` | Target-normal fitting stiffness, not a force multiplier | `1` |
| `model.distribution.deviatoric` | Globally centered deviatoric-strain stiffness | `1e-4` |
| `model.distribution.divergence` | Globally centered divergence stiffness, normalized by dimension | `1e-2` |
| `model.hinge` | Hinge strength relative to predicted fitting improvement | `10` |
| `model.distortion` | Spendable relative-distortion budget | `10` |
| `model.jacobian` | Relative Jacobian floor | `1e-2` |
| `model.qualityGuard` | Fraction of identity margins used for hinge activation | `0.1` |
| `model.jacobianWeight`, `model.distortionWeight` | Relative weights of the two hinge rows | `1` each |
| `model.robustScale` | Fixed Welsch scale when positive | `0`: automatic |
| `globalization.directionalNewton` | Predictor/model directional scaling | `true` |
| `globalization.maxStepOverH` | Optional predictor-motion cap divided by reference spacing | `0`: unrestricted |
| `globalization.armijo` | Sufficient-decrease coefficient | `1e-4` |
| `linear.solver` | Local linear backend | MUMPS if compiled in; otherwise SparseLU |
| `linear.threads` | Backend thread policy | `0`: unchanged |
| `trace` | Iteration diagnostics | `false` |
| `traceQualityWitness` | Limiting-cell predicted/actual quality diagnostics | `false` |

For example, enable tracing or select an available backend through the API:

```cpp
parameters.trace = true;
parameters.linear.solver = Adaptation::SWIFT::Parameters::LinearSolver::SparseLU;
```

### Convergence

| C++ parameter under `convergence` | Meaning | Default |
|----------------------------------|---------|---------|
| `tolerance.geometric` | Sampled maximum-error success target | `0`: automatic \(h^{p+1}\) for displacement degree \(p\) |
| `tolerance.innerRelative` | Inner stationarity residual relative to fitting force | `1e-3` |
| `tolerance.innerAbsolute` | Absolute inner stationarity allowance | `1e-12` |
| `tolerance.linearRelative` | Linear residual tolerance for all backends | `1e-6` |
| `tolerance.energy` | Relative energy-change stagnation threshold | `1e-8` |
| `tolerance.step` | Absolute accepted-motion stagnation threshold | `0` |
| `tolerance.stepOverH` | Accepted-motion/reference-spacing stagnation threshold | `5e-4` |
| `iterations.outer` | Outer fitting cap | `30` |
| `iterations.inner` | Newton-correction cap per outer iteration | `15` |
| `iterations.linear` | Iteration cap per CG solve; does not cap direct solves | `1000` |
| `iterations.backtracks` | Outer trial-halving cap | `32` |
| `iterations.stagnation` | Consecutive small steps or energy changes before a best-effort exit | `5` |

A positive geometric tolerance overrides the automatic target. Zero does not
disable geometric stopping. Stagnation and iteration-limit exits are distinct
from geometric success.

### Quadrature and Validation

| Displacement/geometry | Surface integration | Volume integration | Quality sampling | Geometric sampling |
|-----------------------|---------------------|--------------------|------------------|--------------------|
| Affine simplicial \(P_1\) | `8` | `2` | `2` plus vertices | `32` plus facet vertices |
| These \(P_2/P_3\) examples | `12` | `8` | `16` plus vertices | `32` plus facet vertices |

`quadrature.order` overrides common integration order; `quadrature.surface`
and `quadrature.volume` override their respective integrations independently.
`quadrature.quality` and `quadrature.validation` control independent quality
and geometric sampling. All default to `0`, selecting the automatic policy.
The higher-order policies are provisional, not exactness guarantees for
non-polynomial level sets or nonlinear quality expressions.

## Diagnostics and Output

Each run prints one summary with the following fields:

| Field | Interpretation |
|-------|----------------|
| `exit` | Solver exit reason; distinguishes geometric success from best effort |
| `energy` | Final robust fitting energy |
| `D_inf` | Maximum sampled gradient-normalized level-set discrepancy |
| `C` | Sampled geometric discrepancy divided by \(h^{p+1}\) |
| `outer` | Reported outer fitting iteration count |
| `inner` | Accumulated inner Newton corrections across the solve |
| `min_j` | Minimum sampled relative Jacobian |
| `max_Q` | Maximum sampled relative distortion |

For this signed-distance target, the normalized level-set discrepancy measures
distance to the circle or sphere at sampled points. It is not a certified
continuous supremum or two-sided Hausdorff distance. A single target hit does
not establish a refinement order.

The program returns `1` for invalid argument ranges or an unsatisfied quality
budget. A quality-valid best-effort result is written and returns `0`, even if
the geometric target was missed; inspect `exit` and `D_inf` to distinguish it
from a target hit. Other invalid inputs or runtime errors may raise exceptions.

| Degree | Mesh description | Mesh data |
|--------|------------------|-----------|
| \(P_1\) | `swift/ReconstructionP1.xdmf` | `swift/ReconstructionP1.mesh.h5` |
| \(P_2\) | `swift/ReconstructionP2.xdmf` | `swift/ReconstructionP2.mesh.h5` |
| \(P_3\) | `swift/ReconstructionP3.xdmf` | `swift/ReconstructionP3.mesh.h5` |

Paths are relative to the working directory. Open the XDMF file in an
XDMF-compatible viewer such as ParaView and keep its HDF5 companion alongside
it. The examples write mesh geometry, not a campaign CSV or trajectory archive.

## Further Documentation

- [SWIFT guide source](../../../doc/Guides/SWIFT.dox): mathematical model,
  solver stages, supported contexts, and extension points.
- [Parameter definitions](../../../src/Rodin/Adaptation/SWIFT/Parameters.h):
  complete parameter API and automatic quadrature policies.
- [Problem implementation](../../../src/Rodin/Adaptation/SWIFT/Problem.h):
  direct variational-problem interface.
- [Adapt workflow](../../../src/Rodin/Adaptation/SWIFT/Adapt.h): owned linear
  displacement and in-place mesh adaptation.

Campaign parsing, response extraction and trajectory drivers remain in
[`experiments/swift_calibration`](../../../experiments/swift_calibration),
outside these examples and excluded from the default build.
