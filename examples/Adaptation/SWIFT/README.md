# Reconstruction

These self-contained examples demonstrate the Strain-distributed Welsch
Implicit-interface Fitting Technique (SWIFT). A classified mesh interface is
fitted to a circle or lobed circle in two dimensions, or a sphere or lobed
sphere in three dimensions. The
mesh topology is retained; the outer box boundary is free to move.

## Examples

| Source | Executable | Displacement and geometry degree |
|--------|------------|----------------------------------|
| `Reconstruction.cpp` | `SWIFT_ReconstructionP1` | \(P_1\) |
| `Reconstruction.cpp` | `SWIFT_ReconstructionP2` | \(P_2\) |
| `Reconstruction.cpp` | `SWIFT_ReconstructionP3` | \(P_3\) |

All three examples support both spatial dimensions. The polynomial degree
controls the displacement space and resulting geometric transformations,
not the target geometry. CMake compiles the same source three times with a
different `RODIN_SWIFT_RECONSTRUCTION_DEGREE`. Parsing, classification,
`SWIFT::Problem`, geometry application and diagnostics are shared; no separate
degree-specific driver is maintained.

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
SWIFT_ReconstructionP1 [--name=value | --name value]...
SWIFT_ReconstructionP2 [--name=value | --name value]...
SWIFT_ReconstructionP3 [--name=value | --name value]...
```

| Flag | Meaning | Default |
|------|---------|---------|
| `--n` | Grid points per axis, not cells; integer at least `4` | `16` |
| `--dimension` | Spatial dimension: `2` or `3` | `2` |
| `--lobes` | Angular frequency; nonnegative and integer in 2D | `0`: circle/sphere |
| `--amp` | Radial perturbation amplitude, ignored when lobes is zero | `0.05` |
| `--R0` | Base radius; must exceed amplitude | `0.25` |
| `--phase` | Angular phase in radians; a rotation about the third axis in 3D | `0` |
| `--cx`, `--cy`, `--cz` | Target-center coordinates; third coordinate used only in 3D | `0.5` each |
| `--output` | Output stem without extension | `swift/ReconstructionP<degree>` |
| `--help` | Print every flag, its meaning and current value, then exit | Off |

Both `--name=value` and `--name value` are accepted. Boolean flags accept `0`
or `1`, or enable their option when supplied without a value. Unknown flags,
missing values, nonfinite/malformed numbers and unsupported solvers are
rejected. The previous positional syntax is replaced by named flags.

For example, run from a scratch directory so generated files remain outside
the source tree. The following commands assume the shell initially starts at
the repository root:

```sh
build_dir="$(pwd)/build"
mkdir -p /tmp/rodin-swift
cd /tmp/rodin-swift

# Reconstruct a circle with linear, quadratic and cubic displacement.
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP1" --n=16 --dimension=2
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP2" --n=16 --dimension=2
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP3" --n=16 --dimension=2

# Reconstruct a sphere on a smaller three-dimensional grid.
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP1" --n=8 --dimension=3
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP2" --n=8 --dimension=3
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP3" --n=8 --dimension=3

# Fit the same four-lobe target with explicit model weights and work budgets.
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP2" \
  --n=16 --dimension=2 --lobes=4 --amp=0.05 --R0=0.25 \
  --swift-fit=1 --swift-distribution-deviatoric=1e-4 \
  --swift-distribution-divergence=1e-2 --swift-hinge=10 \
  --swift-outer-iterations=30 --swift-inner-iterations=15 \
  --swift-linear-solver=sparse-lu --swift-trace \
  --output=results/lobed-p2
```

Three-dimensional meshes and higher-order spaces require more memory and
factorization work. Start with a small grid. Runs of the same degree in the
same working directory overwrite that degree's output; use separate working
directories to retain several resolutions or dimensions.

Replace only the executable name to repeat a configuration with another
degree. Model defaults are identical, while automatic quadrature and geometric
target policies still depend on degree. These illustrative examples retain
centroid classification and impose no exterior Dirichlet condition; they do
not reproduce the historical calibration drivers' classifier or boundary setup.

## Geometry and Workflow

The background domain is the unit square or cube, triangulated or
tetrahedralized from a uniform grid. Its fixed reference spacing is

\[
h=\frac{1}{n-1}.
\]

With `--lobes=0`, the target is the zero level set of the signed distance

\[
\phi(x)=\lVert x-c\rVert-R,
\qquad c=(0.5,\ldots,0.5),\qquad R=0.25.
\]

For a positive lobe frequency in 2D, polar coordinates about the specified
center define

\[
\phi(x)=r-R_0-a\cos\bigl(\ell(\theta-\theta_0)\bigr).
\]

In 3D, the axis-balanced target uses a unit radial direction rotated about the
third axis by the negative phase, denoted by \(\widehat n\), and

\[
\phi(x)=r-R_0-\frac{a}{3}\sum_{i=1}^{3}\cos(\ell\widehat n_i).
\]

Both supply analytic gradients. In 3D, `--lobes` controls directional
frequency, not an exact count of visible lobes. Lobed targets are implicit
radial functions rather than exact signed distances. Keep the target inside
the background box for a closed reconstruction.

Cells are classified by the sign of the level set at their centroids. Facets
separating the two cell attributes form the interface to fit. This classified
interface is only an approximation to the target; no intersection cutting or
remeshing is performed.

| Attribute | Value | Role |
|-----------|-------|------|
| `Inside` | `1` | Cells with negative centroid level-set value |
| `Outside` | `2` | Remaining cells |
| `Interface` | `10` | Facets separating inside and outside cells |
| `Boundary` | `20` | Exterior box facets; label only, with no prescribed displacement |

Each example prepares mesh connectivity, constructs the level set and its
analytic gradient, marks the facets, solves the fit, checks the quality
report, and writes the resulting geometry.

## API Usage and Ownership

Every executable uses `SWIFT::Problem` directly. The shared source selects
the displacement degree at compilation:

```cpp
// Degree is 1, 2 or 3, selected by the executable's CMake target.
H1 space(std::integral_constant<size_t, Degree>{}, mesh, dimension);
TrialFunction u(space);
TestFunction v(space);
Adaptation::SWIFT::Problem fitting(u, v);
// Named flags populate the same hierarchical parameter object for every degree.
fitting.setParameters(options.parameters).setInterfaceAttribute(Interface);
const auto report = fitting.solve(phi, gradient);
const auto& displacement = u.getSolution(); // mesh has not been modified.
```

The mesh, space, trial and test functions must outlive the problem. The
supplied gradient must differentiate the level set; its Hessian is not
required. No Dirichlet condition is imposed by this example: exterior mesh
facets move with the computed displacement. Applications may explicitly add
boundary conditions to the problem when their geometry requires them.

After a valid solve, all degrees create a separate mesh and
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

The reference scale `model.h` is derived from `--n` and held fixed. All other
values below start from production defaults and can be overridden by the
named flags. The complete model/control mapping is:

| Flag | C++ parameter |
|------|---------------|
| `--swift-fit` | `model.fit` |
| `--swift-distribution-deviatoric` | `model.distribution.deviatoric` |
| `--swift-distribution-divergence` | `model.distribution.divergence` |
| `--swift-hinge` | `model.hinge` |
| `--swift-jacobian` | `model.jacobian` |
| `--swift-distortion` | `model.distortion` |
| `--swift-quality-guard` | `model.qualityGuard` |
| `--swift-jacobian-weight` | `model.jacobianWeight` |
| `--swift-distortion-weight` | `model.distortionWeight` |
| `--swift-robust-scale` | `model.robustScale` |
| `--swift-directional-newton` | `globalization.directionalNewton` |
| `--swift-max-step-over-h` | `globalization.maxStepOverH` |
| `--swift-armijo` | `globalization.armijo` |
| `--swift-linear-solver` | `linear.solver`: `cg`, `sparse-lu`, or compiled-in `mumps` |
| `--swift-linear-threads` | `linear.threads` |
| `--swift-trace`, `--trace` | `trace` |
| `--swift-quality-witness` | `traceQualityWitness` |

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

| Flag | C++ parameter under `convergence` |
|------|----------------------------------|
| `--swift-geometric-tolerance` | `tolerance.geometric` |
| `--swift-inner-relative-tolerance` | `tolerance.innerRelative` |
| `--swift-inner-absolute-tolerance` | `tolerance.innerAbsolute` |
| `--swift-linear-relative-tolerance` | `tolerance.linearRelative` |
| `--swift-energy-tolerance` | `tolerance.energy` |
| `--swift-step-tolerance` | `tolerance.step` |
| `--swift-step-over-h-tolerance` | `tolerance.stepOverH` |
| `--swift-outer-iterations` | `iterations.outer` |
| `--swift-inner-iterations` | `iterations.inner` |
| `--swift-linear-iterations` | `iterations.linear` |
| `--swift-backtracks` | `iterations.backtracks` |
| `--swift-stagnation-iterations` | `iterations.stagnation` |

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

| Flag | C++ parameter |
|------|---------------|
| `--quad-order` | `quadrature.order` |
| `--surface-quadrature-order` | `quadrature.surface` |
| `--volume-quadrature-order` | `quadrature.volume` |
| `--quality-validation-order` | `quadrature.quality` |
| `--geometric-validation-order` | `quadrature.validation` |

| Displacement/geometry | Surface integration | Volume integration | Quality sampling | Geometric sampling |
|-----------------------|---------------------|--------------------|------------------|--------------------|
| Affine simplicial \(P_1\) | `8` | `2` | `2` lattice subdivisions | `32` plus facet vertices |
| These \(P_2/P_3\) examples | `12` | `8` | `16` lattice subdivisions | `32` plus facet vertices |

`quadrature.order` overrides common integration order; `quadrature.surface`
and `quadrature.volume` override their respective integrations independently.
`quadrature.quality` selects the shared inner-hinge and actual-quality witnesses:
uniform barycentric lattices on simplices, including vertices. Tensor cells use
Cartesian grids, wedges use triangle-times-segment grids, and pyramids use
shrinking square layers. No supplemental points are added. The setting counts
subdivisions per reference edge, not polynomial degree. Equal positive reference
weights sum to reference volume; mapped weights determine each cell's discrete
mass. Separate adaptive Jacobian and distortion measures mix equal mass with
normalized nonlinear guard penetration
at the current and full-predictor geometries. These positive measures preserve
cell mass and remain frozen throughout the inner solve, with the same weights
in its energy, residual and tangent. The equal fraction is one half, risks are
capped at 100, and an inverted predictor receives maximal distortion risk.
Actual quality checks take unweighted extrema over all witnesses. These sets
are not degree-exact integration rules or continuous quality certificates.
The default simplex counts are 6/10 witnesses on triangles/tetrahedra for
affine \(P_1\), and 153/969 for \(P_2/P_3\).
`quadrature.validation` controls independent geometric sampling.
All default to `0`, selecting the automatic policy.
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
| `target`, `target_hit` | Actual geometric threshold and whether it was reached |
| `quality_ok` | Whether the sampled quality budget is satisfied |
| `outer` | Reported outer fitting iteration count |
| `inner` | Accumulated inner Newton corrections across the solve |
| `inner_max`, `inner_last` | Maximum and final corrections per outer step |
| `inner_residual` | Final inner stationarity residual |
| `min_j` | Minimum sampled relative Jacobian |
| `max_Q` | Maximum sampled relative distortion |

For this signed-distance target, the normalized level-set discrepancy measures
distance to the circle or sphere at sampled points. It is not a certified
continuous supremum or two-sided Hausdorff distance. A single target hit does
not establish a refinement order.

The program returns `1` for invalid arguments, caught solver/runtime errors or an unsatisfied quality
budget. A quality-valid best-effort result is written and returns `0`, even if
the geometric target was missed; inspect `exit` and `D_inf` to distinguish it
from a target hit.

| Degree | Mesh description | Mesh data |
|--------|------------------|-----------|
| \(P_1\) | `swift/ReconstructionP1.xdmf` | `swift/ReconstructionP1.{background,moved}.mesh.h5` |
| \(P_2\) | `swift/ReconstructionP2.xdmf` | `swift/ReconstructionP2.{background,moved}.mesh.h5` |
| \(P_3\) | `swift/ReconstructionP3.xdmf` | `swift/ReconstructionP3.{background,moved}.mesh.h5` |

Default paths are relative to the working directory; `--output` replaces the
stem. Open the XDMF file in an XDMF-compatible viewer such as ParaView and keep
all its HDF5 companions alongside it, including the field files. It contains
two named blocks: `background` is the unchanged reference mesh, and `moved` is
the fitted isoparametric mesh. Use **Extract Block** to select either one,
then **Threshold** on `cell_label` (or `Attribute`) with value `1` for the
inside domain or `2` for the outside domain. Classification is retained from
the background; it is not recomputed after fitting.

| Field on both blocks | Location | Interpretation |
|----------------------|----------|----------------|
| `Attribute`, `cell_label` | Cell | Original inside/outside classification |
| `displacement` | Node | Accumulated displacement from the background configuration |
| `phi` | Node | Target level set interpolated on the respective configuration |
| `j` | Cell | Relative deformation Jacobian at the cell centroid |
| `q_rel` | Cell | Relative deformation distortion at the cell centroid |

The displacement on `moved` is the transported background field, not an
additional deformation to apply. Both blocks carry the final deformation
diagnostics. The cell-centered quality fields are centroid samples, not
cellwise extrema or the solver's validation certificate. The examples do not
write a campaign CSV or trajectory archive.

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
