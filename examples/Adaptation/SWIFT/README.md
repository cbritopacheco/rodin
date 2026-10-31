# Interface reconstruction

These self-contained examples demonstrate the Strain-distributed Welsch
Implicit-interface Fitting Technique (SWIFT). A classified mesh interface is
fitted to a circle or lobed circle in two dimensions, or a sphere or lobed
sphere in three dimensions. The
mesh topology is retained; the outer box boundary is free to move.

`InterfaceFitting.cpp` is the introductory P1 tutorial. It creates a unit-square
mesh, defines a circle and its gradient, classifies an internal interface,
solves for a vector P1 displacement and writes reference/fitted mesh blocks.
The executable requires both a geometric target hit and a valid sampled
quality budget. The exterior boundary is free to move.

```sh
cmake --build build --target SWIFT_InterfaceFitting
mkdir -p results/interface-fitting
cd results/interface-fitting
../../build/examples/Adaptation/SWIFT/SWIFT_InterfaceFitting
```

The generated **Fitting a mesh interface with SWIFT** guide explains each code
stage and shows a short P2 extension. Render its computed P1 panels with the following command.

```sh
pvpython examples/Adaptation/SWIFT/visualize_background.py \
  results/interface-fitting doc/Images/SWIFTP1Tutorial.png \
  --reconstruction InterfaceFitting --degree 1
```

`BackgroundFitting.cpp` adds a step-by-step MMG–MinSTCut–SWIFT example with
the outer square held fixed. MMG prepares a graded background before
classification; SWIFT then retains its topology while fitting the circle.
The generated documentation's **Adapting and fitting an implicit interface**
tutorial explains the formulas, ownership and measured results.

```sh
cmake --build build --target SWIFT_BackgroundFitting
mkdir -p results/background-fitting
cd results/background-fitting
../../build/examples/Adaptation/SWIFT/SWIFT_BackgroundFitting
```

It writes three named mesh blocks to `BackgroundFitting.xdmf`, HDF5 companion
files and `BackgroundFitting.csv`. The program requires a geometric target hit,
a valid sampled quality budget, positive signed triangle Jacobians, unit area,
fixed exterior vertices and unchanged counts during fitting. An independent
segment-to-circle calculation checks the entire polygonal interface against
the geometric tolerance. Outputs belong in a scratch directory.

Render the actual mesh comparison from the repository root.

```sh
pvpython examples/Adaptation/SWIFT/visualize_background.py \
  results/background-fitting doc/Images/SWIFTBackgroundFitting.png
```

The figure uses Manrope, the Rodin palette and a transparent background. The
navy circle is an analytic reference over the computed meshes.

## Examples

| Source | Executable | Displacement and geometry degree |
|--------|------------|----------------------------------|
| `Reconstruction.cpp` | `SWIFT_ReconstructionP1` | \(P_1\) |
| `Reconstruction.cpp` | `SWIFT_ReconstructionP2` | \(P_2\) |
| `Reconstruction.cpp` | `SWIFT_ReconstructionP3` | \(P_3\) |

All three examples support both spatial dimensions. The polynomial degree
selects the displacement space and resulting geometric transformations. CMake compiles the same source three times with a
different `RODIN_SWIFT_RECONSTRUCTION_DEGREE`. Parsing, classification, `SWIFT::Problem`, geometry application and
diagnostics share one implementation.

## Build

Run the following commands from the repository root. Initialize the submodules
for a fresh checkout and configure a local build with examples enabled.

```sh
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DRODIN_BUILD_EXAMPLES=ON -DRODIN_USE_PETSC=OFF \
  -DRODIN_INSTALL_RESOURCES=OFF
cmake --build build --target SWIFT_ReconstructionP1 SWIFT_ReconstructionP2 SWIFT_ReconstructionP3 -j1
```

The standard Rodin dependencies, including Boost, Eigen and HDF5, must be
available to CMake. These examples do not require downloaded mesh resources
or PETSc because they construct local meshes and use Eigen-backed displacement
storage. See the repository's [build workflow](../../../doc/agents/workflows.md)
for general build procedures.

For optional MUMPS support, add `-DRODIN_USE_MUMPS=ON` when configuring, with
MUMPS and its dependencies installed. SWIFT selects MUMPS when compiled in;
otherwise its default direct solver is SparseLU. Add `-DRODIN_USE_OPENMP=ON`
for an OpenMP-enabled build. Thread environment variables control the features enabled at configuration.

An existing configured build can be reused by replacing `build` in the
commands. Only the three named targets need to be built.

## Run

The executable accepts the following syntax.

```text
SWIFT_ReconstructionP1 [--name=value | --name value]...
SWIFT_ReconstructionP2 [--name=value | --name value]...
SWIFT_ReconstructionP3 [--name=value | --name value]...
```

| Flag | Meaning | Default |
|------|---------|---------|
| `--n` | Grid points per axis; integer at least `4` | `16` |
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
the repository root.

```sh
build_dir="$(pwd)/build"
mkdir -p /tmp/rodin-swift
cd /tmp/rodin-swift

# Reconstruct a circle with linear, quadratic and cubic displacement.
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP1" --n=16 --dimension=2 --classification-smoothing=0.12
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP2" --n=16 --dimension=2 --classification-smoothing=0.12
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP3" --n=16 --dimension=2 --classification-smoothing=0.12

# Reconstruct a sphere on a smaller three-dimensional grid.
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP1" --n=8 --dimension=3 --classification-smoothing=0.028
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP2" --n=8 --dimension=3 --classification-smoothing=0.028
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP3" --n=8 --dimension=3 --classification-smoothing=0.028

# Fit the same four-lobe target with explicit model weights and work budgets.
"$build_dir/examples/Adaptation/SWIFT/SWIFT_ReconstructionP2" \
  --n=16 --dimension=2 --lobes=4 --amp=0.05 --R0=0.25 \
  --classification-smoothing=0.12 \
  --model-fit=1 --model-distribution-deviatoric=1e-4 \
  --model-distribution-divergence=1e-4 --model-hinge-jacobian=10 --model-hinge-distortion=100 \
  --convergence-iterations-outer=30 --convergence-iterations-inner=15 \
  --linear-solver=sparse-lu --trace \
  --output=results/lobed-p2
```

Three-dimensional meshes and higher-order spaces require more memory and
factorization work. Start with a small grid. Runs of the same degree in the
same working directory overwrite that degree's output; use separate working
directories to retain several resolutions or dimensions.

Replace only the executable name to repeat a configuration with another
degree. Model defaults are identical, while automatic quadrature and geometric
target policies still depend on degree. All reconstruction examples use
`Geometry::MinSTCut` and impose no exterior Dirichlet condition. Historical
campaign data retain their original classification and boundary setup.

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

Both supply analytic gradients. In 3D, `--lobes` controls the directional frequency of the implicit radial target. Keep the target inside
the background box for a closed reconstruction.

`Geometry::MinSTCut` is the canonical classifier. Its data are cell averages
of a smooth signed phase indicator.

\[
m_K=\frac{1}{|K|}\int_K\tanh\!\left(\frac{\phi(x)}{\varepsilon}\right)dx.
\]

The binary Potts objective balances agreement with these averages against
the measure of separating facets. The returned `inside` and `outside` cell
indices receive attributes; `cut` identifies the interface to fit. There is
no centroid-sign alternative or sign pinning. This classified interface is
only an approximation to the target; no intersection cutting or remeshing is
performed. Classification precedes SWIFT and is independent of displacement
degree.

| Classifier Setting | Flag | Default |
|--------------------|------|---------|
| Fidelity | `--classification-fidelity` | \(1\) |
| Dimensionless facet smoothing | `--classification-smoothing` | \(0.01\) |
| Phase transition scale | Fixed policy | \(\varepsilon=1.25h\) |
| Cell-average quadrature order | Fixed policy | \(2\) |

The reconstruction driver converts the dimensionless smoothing coefficient
\(\widehat s\) into the physical facet weight

\[
s=h\widehat s,\qquad h=\frac{1}{n-1}.
\]

MinSTCut itself performs no mesh scaling and uses the supplied facet function
directly. Increasing \(\widehat s\) relative to fidelity regularizes the
classified interface more strongly. This scaling balances cell and facet costs for mesh-scale label changes. The
  smoothing weight must be chosen to retain the features of interest; larger
  values can remove small regions.

The flag previously supplied the physical weight \(s\). To reproduce an old
run at the same resolution, now pass \(\widehat s=s/h\). For example, the
previous 2D weight \(0.008\) at \(n=16\) becomes `--classification-smoothing=0.12`;
the successful weight \(0.001\) at that resolution becomes
`--classification-smoothing=0.015`. Fidelity defaults to \(1\); dimensionless
smoothing defaults to \(0.01\). MinSTCut's own default facet function remains one.

Quadrature approximates the phase average of the nonpolynomial integrand. Smoothing is a facet
function in the MinSTCut API; these examples use a constant function.

| Attribute | Value | Role |
|-----------|-------|------|
| `Inside` | `1` | Cells in the minimum-cut inside set |
| `Outside` | `2` | Cells in the minimum-cut outside set |
| `Interface` | `10` | Facets separating inside and outside cells |
| `Boundary` | `20` | Exterior box facets; label only, with no prescribed displacement |

Each example prepares mesh connectivity, constructs the level set and its
analytic gradient, marks the facets, solves the fit, checks the quality
report, and writes the resulting geometry.

## API Usage and Ownership

Every executable uses `SWIFT::Problem` directly. The shared source selects
the displacement degree at compilation.

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
solver uses the supplied level set and its exact gradient. No Dirichlet condition is imposed by this example, so exterior mesh
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
deviatoric and divergence weights. The centering uses the full current domain. Compatible global affine motions are free in the distribution term;
shared finite-element degrees of freedom couple neighboring elements.

The force and outer merit use robust Welsch fitting energy. Directional Newton
scales the quadratic model before affine quadratic quality hinges recover the
direction. Armijo backtracking accepts a trial only when energy decrease and
the actual sampled Jacobian/distortion checks both pass. The sampled quality budget is enforced during trial acceptance.

The reference scale `model.h` is derived from `--n` and held fixed. All other
values below start from production defaults and can be overridden by the
named flags. The complete model/control mapping is:

| Flag | C++ parameter |
|------|---------------|
| `--model-fit` | `model.fit` |
| `--model-distribution-deviatoric` | `model.distribution.deviatoric` |
| `--model-distribution-divergence` | `model.distribution.divergence` |
| `--model-jacobian` | `model.jacobian` |
| `--model-distortion` | `model.distortion` |
| `--model-quality-guard` | `model.qualityGuard` |
| `--model-hinge-jacobian` | `model.hinge.jacobian` |
| `--model-hinge-distortion` | `model.hinge.distortion` |
| `--model-robust-scale` | `model.robustScale` |
| `--globalization-directional-newton` | `globalization.directionalNewton` |
| `--globalization-max-step-over-h` | `globalization.maxStepOverH` |
| `--globalization-armijo` | `globalization.armijo` |
| `--linear-solver` | `linear.solver`: `cg`, `sparse-lu`, or compiled-in `mumps` |
| `--linear-threads` | `linear.threads` |
| `--trace` | `trace` |
| `--trace-quality-witness` | `traceQualityWitness` |

| C++ parameter | What it controls | Default |
|---------------|------------------|---------|
| `model.h` | Fixed background reference spacing | Required; set to \(1/(n-1)\) |
| `model.fit` | Target-normal fitting stiffness | `1` |
| `model.distribution.deviatoric` | Globally centered deviatoric-strain stiffness | `1e-4` |
| `model.distribution.divergence` | Globally centered divergence stiffness, normalized by dimension | `1e-4` |
| `model.distortion` | Spendable relative-distortion budget | `10` |
| `model.jacobian` | Relative Jacobian floor | `1e-2` |
| `model.qualityGuard` | Fraction of identity margins used for hinge activation | `0.1` |
| `model.hinge.jacobian` | Jacobian recovery strength | `10` |
| `model.hinge.distortion` | Distortion recovery strength | `100` |
| `model.robustScale` | Fixed Welsch scale when positive | `0`: automatic |
| `globalization.directionalNewton` | Predictor/model directional scaling | `true` |
| `globalization.maxStepOverH` | Optional predictor-motion cap divided by reference spacing | `0`: unrestricted |
| `globalization.armijo` | Sufficient-decrease coefficient | `1e-4` |
| `linear.solver` | Local linear backend | MUMPS if compiled in; otherwise SparseLU |
| `linear.threads` | Backend thread policy | `0`: unchanged |
| `trace` | Iteration diagnostics | `false` |
| `traceQualityWitness` | Limiting-cell predicted/actual quality diagnostics | `false` |

These model and classifier defaults follow the current shared leader of the
ongoing 2D calibration. They remain provisional until the campaign, including
its 3D stage, finishes. Work budgets, tolerances and sampling policies are unchanged.

For example, enable tracing or select an available backend through the API.

```cpp
parameters.trace = true;
parameters.linear.solver = Adaptation::SWIFT::Parameters::LinearSolver::SparseLU;
```

### Convergence

| Flag | C++ parameter under `convergence` |
|------|----------------------------------|
| `--convergence-tolerance-geometric` | `convergence.tolerance.geometric` |
| `--convergence-tolerance-inner-relative` | `convergence.tolerance.innerRelative` |
| `--convergence-tolerance-inner-absolute` | `convergence.tolerance.innerAbsolute` |
| `--convergence-tolerance-linear-relative` | `convergence.tolerance.linearRelative` |
| `--convergence-tolerance-energy` | `convergence.tolerance.energy` |
| `--convergence-tolerance-step` | `convergence.tolerance.step` |
| `--convergence-tolerance-step-over-h` | `convergence.tolerance.stepOverH` |
| `--convergence-iterations-outer` | `convergence.iterations.outer` |
| `--convergence-iterations-inner` | `convergence.iterations.inner` |
| `--convergence-iterations-linear` | `convergence.iterations.linear` |
| `--convergence-iterations-backtracks` | `convergence.iterations.backtracks` |
| `--convergence-iterations-stagnation` | `convergence.iterations.stagnation` |

| C++ parameter under `convergence` | Meaning | Default |
|----------------------------------|---------|---------|
| `convergence.tolerance.geometric` | Sampled maximum-error success target | `0`: automatic \(h^{p+1}\) for displacement degree \(p\) |
| `convergence.tolerance.innerRelative` | Inner stationarity residual relative to fitting force | `1e-3` |
| `convergence.tolerance.innerAbsolute` | Absolute inner stationarity allowance | `1e-12` |
| `convergence.tolerance.linearRelative` | Linear residual tolerance for all backends | `1e-6` |
| `convergence.tolerance.energy` | Relative energy-change stagnation threshold | `1e-8` |
| `convergence.tolerance.step` | Absolute accepted-motion stagnation threshold | `0` |
| `convergence.tolerance.stepOverH` | Accepted-motion/reference-spacing stagnation threshold | `5e-4` |
| `convergence.iterations.outer` | Outer fitting cap | `30` |
| `convergence.iterations.inner` | Newton-correction cap per outer iteration | `15` |
| `convergence.iterations.linear` | Iteration cap per CG solve | `1000` |
| `convergence.iterations.backtracks` | Outer trial-halving cap | `32` |
| `convergence.iterations.stagnation` | Consecutive small steps or energy changes before a best-effort exit | `5` |

A positive geometric tolerance overrides the automatic target. Zero selects the automatic geometric target. Stagnation and iteration-limit exits are distinct
from geometric success.

### Quadrature and Sampling

Parameter flags follow the C++ hierarchy: dots become hyphens and camel-case
words are separated by hyphens. For example, `parameters.sampling.subdivision`
maps to `--sampling-subdivision`, and `parameters.convergence.iterations.outer`
maps to `--convergence-iterations-outer`. Previous flat flag names are rejected.

| Flag | C++ parameter |
|------|---------------|
| `--quadrature-order` | `quadrature.order` |
| `--quadrature-surface` | `quadrature.surface` |
| `--quadrature-volume` | `quadrature.volume` |
| `--sampling-subdivision` | `sampling.subdivision` |
| `--quadrature-validation` | `quadrature.validation` |

| Displacement/geometry | Surface integration | Volume integration | Quality sampling | Geometric sampling |
|-----------------------|---------------------|--------------------|------------------|--------------------|
| Affine simplicial \(P_1\) | `8` | `2` | `2` covering subdivisions | Facet covering, `2` subdivisions |
| These \(P_2/P_3\) examples | `12` | `8` | `16` covering subdivisions | Facet covering, `16` subdivisions |

`quadrature.order` overrides common integration order; `quadrature.surface`
and `quadrature.volume` override their respective integrations independently.
`sampling.subdivision` selects the shared inner-hinge and actual-quality witnesses
from closed-form shifted/contracted reference lattices, without supplemental points.
The setting specifies subdivisions of the generating lattice. Shifting and
  contracting the lattice preserves its point count and places its witnesses in
  the interior.
Except on tetrahedra, sites are \((I+\tfrac12\boldsymbol1)/(m+1)\).
Tetrahedra use \((I+a\boldsymbol1)/(m+3a)\), with \(a=1/3\) for \(m=1\)
and \(a=1/(2\sqrt2)\) otherwise. These structured sets have known reference covering radii. Continuous quality
  between witnesses requires further analysis.
Equal positive reference
weights sum to reference volume; mapped weights determine each cell's discrete
mass. Separate adaptive Jacobian and distortion measures mix equal mass with
normalized nonlinear guard penetration
at the current and full-predictor geometries. These positive measures preserve
cell mass and remain frozen throughout the inner solve, with the same weights
in its energy, residual and tangent. The equal fraction is one half, risks are
capped at 100, and an inverted predictor receives maximal distortion risk.
Actual quality checks take unweighted extrema over all witnesses. These sets define the sample points for nonlinear quality checks.
The default simplex counts are 6/10 witnesses on triangles/tetrahedra for
affine \(P_1\), and 153/969 for \(P_2/P_3\).
The interface maximum uses the same covering construction on facets, with
resolution selected by `sampling.subdivision`. It does not add vertices or
quadrature points. `quadrature.validation` controls the independent
RMS and normal-error integrals and physical-motion diagnostics (automatic
order `32` for interface integrals).
All default to `0`, selecting the automatic policy.
The higher-order policies are provisional integration and sampling choices. Their
  accuracy for nonpolynomial level sets and nonlinear quality expressions is
  assessed by refinement.

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
distance to the circle or sphere at sampled points. A continuous supremum or two-sided Hausdorff bound requires additional geometric
  analysis. Refinement order is measured over a sequence of meshes.

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
two named blocks. `background` is the unchanged reference mesh, and `moved` is
the fitted isoparametric mesh. Use **Extract Block** to select either one,
then **Threshold** on `cell_label` (or `Attribute`) with value `1` for the
inside domain or `2` for the outside domain. The moved mesh retains the background classification.

| Field on both blocks | Location | Interpretation |
|----------------------|----------|----------------|
| `Attribute`, `cell_label` | Cell | Original inside/outside classification |
| `displacement` | Node | Accumulated displacement from the background configuration |
| `phi` | Node | Target level set interpolated on the respective configuration |
| `j` | Cell | Relative deformation Jacobian at the cell centroid |
| `q_rel` | Cell | Relative deformation distortion at the cell centroid |

The displacement on `moved` is the background field transported to the fitted
  geometry. Both blocks carry the final deformation
diagnostics. Cell-centered quality fields record centroid samples; the solver report records
  extrema over its quality witnesses.

## Further Documentation

- [Reference witness generator](Witness/README.md): offline numerical covering
  sets for every Rodin reference geometry, with configurable point counts and
  SVG comparisons. The production solver retains its configured witness sets.

- [SWIFT guide source](../../../doc/Guides/SWIFT.dox): step-by-step P1
  interface-fitting tutorial and a P2 result.
- [Parameter definitions](../../../src/Rodin/Adaptation/SWIFT/Parameters.h):
  complete parameter API and automatic quadrature policies.
- [Problem implementation](../../../src/Rodin/Adaptation/SWIFT/Problem.h):
  direct variational-problem interface.
- [Adapt workflow](../../../src/Rodin/Adaptation/SWIFT/Adapt.h): owned linear
  displacement and in-place mesh adaptation.

Parameter studies use these same reconstruction executables and named controls,
without a separate reconstruction or trajectory implementation.
