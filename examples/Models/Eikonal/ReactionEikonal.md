# Anisotropic reaction–eikonal propagation

`Eikonal_ReactionEikonal` implements the RE+ structure in equation (12) of
Nagel et al., *Comparison of Propagation Models and Forward Calculation Methods
on Cellular, Tissue and Organ Scale Atrial Electrophysiology*, IEEE TBME 70(2),
2023, [doi:10.1109/TBME.2022.3196144](https://doi.org/10.1109/TBME.2022.3196144).
An anisotropic FIM solve schedules a foot current in a Courtemanche
reaction–diffusion simulation on an axis-aligned uniform grid. Monodomain
strand simulations calibrate the squared-velocity tensor.

The example uses a homogeneous two-dimensional sheet, baseline ionic kinetics
and explicitly chosen conductivities. The paper's anatomical mesh, regional
conductivity tuning, ionic remodeling, torso and ECG calculation are outside
this example. Its regional parameters are in separate supplementary tables;
they are not inferred from the main article.

## Equations and units

The eikonal solve computes foot-current onset times satisfying

\[
\nabla T^T M\nabla T=1,\qquad T|_{\Gamma}=0.
\]

The voltage and ionic variables then evolve according to

\[
\beta C_m\partial_t V
=\nabla\cdot(\sigma_m\nabla V)+I_{\mathrm{foot}}-\beta I_{\mathrm{ion}}(V,\eta),
\qquad \partial_t\eta=f(V,\eta),
\qquad n\cdot\sigma_m\nabla V=0.
\]

Lengths are in millimetres, times in milliseconds, voltages in millivolts and
concentrations in millimolar. The normalized diffusion tensor is

\[
D\;[\mathrm{mm^2/ms}]
=1000\frac{\sigma_m\;[\mathrm{S/m}]}
{\beta\;[\mathrm{m^{-1}}]C_m\;[\mathrm{F/m^2}]}.
\]

The chosen parameters are \(\beta=140000\,\mathrm{m^{-1}}\),
\(C_m=0.01\,\mathrm{F/m^2}\), and transverse conductivity
\(\sigma_t=0.1\,\mathrm{S/m}\). For requested velocity ratio \(r\), the
longitudinal conductivity is \(\sigma_l=r^2\sigma_t\). These are illustrative
inputs, not the paper's regional fitted values.

The physical domain is the fixed square \(\Omega=[0,10]^2\,\mathrm{mm^2}\).
`LocalMesh::UniformGrid(Triangle, {n+1,n+1})` supplies regularly spaced,
axis-aligned vertices and two native right triangles per Cartesian grid cell.
The mesh is scaled by \(h=10/n\) mm, with no rotation, shear or stretch.
The isotropic and anisotropic runs therefore have identical physical meshes
and source nodes.

Fiber orientation is prescribed independently as \(f=(0.8,0.6)^T\), with
transverse direction \(s=(-0.6,0.8)^T\) and \(R=[f,s]\). This rotation acts
only on the conductivity and squared-velocity tensors:

\[
\sigma_m=\sigma_t R\operatorname{diag}(r^2,1)R^T,
\qquad M=R\operatorname{diag}(c_l^2,c_t^2)R^T.
\]

The example uses FIM for the activation schedule. At this off-axis fiber
orientation, the standard right triangles are metric-obtuse for the
anisotropic tensor, so they violate the current tensor FMM causal-stencil
contract. FMM is not run on this RE+ setup. The separate
`RodinEikonalPropagationBenchmarks` comparison retains its supported causal
simplex meshes.

## Ionic kinetics and foot current

All 21 baseline Courtemanche–Ramirez–Nattel state equations are included in
[Courtemanche1998.h](Courtemanche1998.h). The cell is first paced at 1 Hz for
the requested number of conditioning beats. The default is ten; this is a
finite conditioning procedure, not a convergence test for a limit cycle.
Setting the beat count to zero uses the published resting state directly.

The normalized positive depolarizing foot current is

\[
j_{\mathrm{foot}}(x,t)
=\frac{A_F}{\tau_F}
 \exp\!\left(\frac{t-T(x)}{\tau_F}\right),
\qquad T(x)\leq t<t_{\mathrm{th}}(x),
\]

and is zero otherwise. Here \(j_{\mathrm{foot}}=I_{\mathrm{foot}}/(\beta C_m)\)
has units mV/ms. The current is permanently disabled for that node after its
voltage first reaches \(V_{\mathrm{th}}=-30\,\mathrm{mV}\), including an
earlier crossing caused by diffusion. Its sign is positive in the voltage
equation; openCARP's outward-current convention assigns the opposite sign.
The exponential is averaged over the time step, including a fractional onset
step. Threshold cutoff is detected at step boundaries.

This follows the exponential-foot formulation in the
[openCARP manual, section 3.2.5](https://opencarp.org/manual/opencarp-manual-latest.pdf)
and the RE construction of
[Neic et al. (2017)](https://doi.org/10.1016/j.jcp.2017.06.020).
The defaults \(A_F=2.464218456421841\,\mathrm{mV}\) and
\(\tau_F=0.43917740479584627\,\mathrm{ms}\) are the published baseline
Courtemanche fit in the
[openCARP foot-current example](https://opencarp.org/documentation/examples/03_eikonal/04_fit_ifoot).
They are demonstration parameters, not a fit to the conditioning protocol or
regional remodeling of this example. Such changes require refitting.

## Calibration and discrete scheme

The squared-velocity tensor is not set equal to conductivity. Two monodomain
strand simulations measure longitudinal and transverse velocities between
probes near physical coordinates 3 and 7 mm:

\[
c_a=\frac{x_{a,2}-x_{a,1}}{t_{a,2}-t_{a,1}},
\qquad M=R\operatorname{diag}(c_l^2,c_t^2)R^T.
\]

Both strands have physical length 10 mm and spacing \(h=10/n\) mm. Each
receives a 20 pA/pF depolarizing pulse for 2 ms on its first 1 mm. Their
principal diffusion coefficients differ according to the requested conductivity
ratio. The sheet source set contains vertices within physical radius 0.5 mm of the
bottom-right corner \((10,0)\) mm. FIM onset times drive the foot current
at every sheet node; the RE+ sheet does not receive the rectangular monodomain stimulus.

Rodin assembles the mass and anisotropic stiffness forms as

```cpp
mass = Integral(u, v);
stiffness = Integral(MatrixFunction(diffusion) * Grad(u), Grad(v));
```

The mass matrix is row-sum lumped. The ionic substep uses frozen-coefficient
Rush–Larsen updates for gates and forward Euler for voltage and concentrations.
Backward Euler then advances diffusion:

\[
(M_L+\Delta t K_D)V^{n+1}=M_LV^*,
\qquad (K_D)_{ij}=\int_\Omega\nabla\phi_i^TD\nabla\phi_j\,dx.
\]

This is a first-order operator split. All ionic variables are explicit nodal
grid functions. The fixed diffusion matrix is factored once and its
factorization reused. A smaller time step and finer spatial mesh are needed
for a numerical-convergence assessment; gating stability does not establish
voltage or concentration accuracy.

A local activation time (LAT) is the midpoint of the step with greatest
positive voltage slope. A node must also attain voltage above -40 mV to be
classified as activated. FIM's \(T\) is a foot-current onset, whereas LAT marks
the subsequent upstroke. Their difference is reported as onset-to-upstroke
lag, not as numerical error. Diffusion can trigger a node before its scheduled
foot current. The default 80 ms run covers activation; use a longer duration
for repolarization studies. Optional animation postprocessing estimates APD90
from exported snapshots as described below. ECG metrics are not calculated.

## Build and run

Build the example with `RODIN_BUILD_EXAMPLES=ON`:

```sh
cmake --build build --target Eikonal_ReactionEikonal -j
mkdir -p /tmp/rodin-replus
cd /tmp/rodin-replus
/path/to/build/examples/Models/Eikonal/Eikonal_ReactionEikonal
```

Optional positional arguments are
`intervals dt_ms end_ms pace_beats velocity_ratio workers`.
Defaults are `40 0.02 80 10 2 1`. The worker count applies to FIM; the ionic
and diffusion time integrator is serial. More than one FIM worker requires
`RODIN_USE_OPENMP=ON`.

Run comparison cases from separate output directories:

```sh
Eikonal_ReactionEikonal 40 0.02 80 10 1 1  # isotropic control
Eikonal_ReactionEikonal 40 0.01 80 10 2 1  # smaller time step
Eikonal_ReactionEikonal 80 0.02 80 10 2 1  # finer mesh, recalibrated velocities
Eikonal_ReactionEikonal 40 0.02 500 10 2 1 # include repolarization
```

`ReactionEikonalVoltage.xdmf` contains the RE+ voltage/calcium time series,
sampled approximately every millisecond.
`ReactionEikonalActivation.xdmf` contains the RE+ LATs,
the FIM foot-onset field and onset-to-upstroke lag. The associated CSV gives
these quantities at every vertex. `ReactionEikonalCell.csv` records the final
conditioning beat. The console reports calibrated velocities, FIM rounds,
and lag statistics outside the stimulated region.
Plot the exported native triangle fields with NumPy, Matplotlib and h5py:

```sh
python3 /path/to/rodin/examples/Models/Eikonal/plot_reaction_eikonal.py . --time 20
```

The script produces `ReactionEikonal.png` and `ReactionEikonal.svg`.

Invalid ionic states or unactivated tissue cause failure. Mesh and time-step
refinement remain necessary for numerical-convergence assessment.

## Observed numerical checks

The current uniform-grid setup is checked with isotropic and anisotropic
500 ms runs, default settings and four FIM workers. Physical grid coordinates,
source sets, finite ionic states and complete activation are checked. The
nodal fields and repolarization estimates for the bottom-right source are:

| Quantity | Isotropic | Anisotropic |
| --- | ---: | ---: |
| Latest RE+ LAT (ms) | 36.49 | 29.91 |
| Activated nodes at 20 ms | 44.1% | 70.8% |
| Latest estimated 90% repolarization (ms) | 340.7 | 329.1 |
| Estimated mean nodal APD90 (ms) | 314.4 | 307.4 |
| FIM active-list rounds | 39 | 41 |

Each run used 1,681 vertices and 3,200 triangles. All six source vertices lie
within 0.5 mm of \((10,0)\) mm. Exact Cartesian coordinates, identical meshes
and source sets, complete activation and complete 90% repolarization were
checked. Default serial and four-worker FIM runs reproduced the CSV fields
and exported voltage snapshots exactly. The larger round counts reflect
propagation from the opposite end of the serial vertex order.


A separate Python evaluation of the original CellML right-hand sides was
integrated with adaptive BDF at relative tolerance \(10^{-9}\) and absolute
tolerance \(10^{-11}\). For one 2 ms, 20 pA/pF stimulus from the published
resting state, its peak voltage was approximately 24.492 mV. The example's
frozen-gate step produced peaks 26.292, 25.382 and 24.935 mV at steps 0.02,
0.01 and 0.005 ms. Voltage at 300 ms approached the reference -68.8379 mV
through -68.9665, -68.9130 and -68.8759 mV. This checks approach to the
source ODE trajectory; it is not validation against experimental data.
The observed changes with mesh and time step prevent treating the default
run as a numerically converged reproduction of the paper.

## Ionic model provenance

The baseline equations and initial values were adapted from the
[Physiome Courtemanche–Ramirez–Nattel 1998 CellML model](https://models.cellml.org/workspace/courtemanche_ramirez_nattel_1998/file/78335420620450b2eab817577c0faf4450bbf8cb/courtemanche_ramirez_nattel_1998.cellml),
revision `78335420620450b2eab817577c0faf4450bbf8cb`, curated by Penny Noble,
Oxford University, and exported by Alan Garny's COR software. The downloaded
CellML file has SHA-256
`10c8a66fa537b17f94c93354a793e90fd493b21356a448134172b66fa7a212fc`.
The model material is available under
[CC BY 3.0](https://creativecommons.org/licenses/by/3.0/).

Changes are translation of the equations to C++, replacement of the embedded
stimulus scheduler by an external normalized depolarizing current, and the
explicit/Rush–Larsen stepping described above. The whole-cell capacitance of
100 pF in the source model is retained for ionic-current normalization; it is
distinct from the membrane capacitance per unit area in the tissue equation.

The original model reference is Courtemanche, Ramirez and Nattel,
*Ionic mechanisms underlying human atrial action potential properties:
insights from a mathematical model*, AJP Heart 275, H301–H321 (1998),
[doi:10.1152/ajpheart.1998.275.1.H301](https://doi.org/10.1152/ajpheart.1998.275.1.H301).
The repository reference is Lloyd et al., *The CellML Model Repository*,
Bioinformatics 24, 2122–2123 (2008),
[doi:10.1093/bioinformatics/btn390](https://doi.org/10.1093/bioinformatics/btn390).

## Animation on identical uniform grids

Run the two cases from separate output directories:

```sh
mkdir -p isotropic anisotropic
(cd isotropic && /path/to/Eikonal_ReactionEikonal 40 0.02 80 10 1 1 > run.log)
(cd anisotropic && /path/to/Eikonal_ReactionEikonal 40 0.02 80 10 2 1 > run.log)
python3 /path/to/rodin/examples/Models/Eikonal/animate_reaction_eikonal.py \
  isotropic anisotropic animation --end-ms 50 --duration 10
```

The renderer requires NumPy, Matplotlib, h5py, and FFmpeg with `libvpx-vp9` and
`libx264`. It verifies the regular axis-aligned square lattice, identical native
coordinates, connectivity, source nodes, output times, time steps and
conditioning counts. Run inputs, fiber directions and calibrated velocities
are read from `ReactionEikonalParameters.csv`.

`REPlusComparison.webm` retains a transparent background; the MP4 version uses
a white background for player compatibility. Both are encoded at 30 fps.
The first 50 ms of simulated time are shown over approximately ten seconds,
with physical time labelled explicitly. Voltage is interpolated linearly
between exported snapshots; both panels use the same -85 to 35 mV color scale,
camera bounds and physical coordinate scale. The arrow marks the anisotropic
fiber direction. The activated-node percentage uses the measured RE+ LATs.
A poster PNG and JSON record of inputs, latest activation and playback are
also written.

Preview the fiber axes on the actual uniform mesh with:

```sh
python3 /path/to/rodin/examples/Models/Eikonal/plot_reaction_eikonal.py anisotropic --fibers
```

This writes `ReactionEikonalFibers.png` and `.svg`. Uniform fibers remain at
\(36.87^\circ\) while the grid edges remain aligned with the coordinate axes.

## Repolarization comparison

For the complete action potential, rerun both uniform-grid cases to 500 ms
and enable the repolarization display:

```sh
(cd isotropic && /path/to/Eikonal_ReactionEikonal 40 0.02 500 10 1 1 > run.log)
(cd anisotropic && /path/to/Eikonal_ReactionEikonal 40 0.02 500 10 2 1 > run.log)
python3 /path/to/rodin/examples/Models/Eikonal/animate_reaction_eikonal.py \
  isotropic anisotropic repolarization --end-ms 500 --duration 20 --repolarization
```

The 30 fps animation shows 500 ms over approximately twenty seconds. The
additional indicator counts nodes that have crossed their 90% repolarization
level after their action-potential peak. With initial voltage \(V_i^0\) and
maximum exported voltage \(V_i^{\mathrm{peak}}\), this level is

\[
V_{90,i}=V_i^0+0.1(V_i^{\mathrm{peak}}-V_i^0).
\]

The first downward post-peak crossing defines \(t_{90,i}\). Crossing times
are interpolated linearly between exported snapshots. The estimated duration
is \(\mathrm{APD}_{90,i}=t_{90,i}-t_{\mathrm{LAT},i}\), using the simulation's
steepest-upstroke LAT. Unobserved crossings retain infinity internally and are
excluded from reported APD statistics. The nodal mean is unweighted.
Both the peak and crossing estimates use the approximately 1 ms output cadence;
they do not have the full time integrator's 0.02 ms resolution.

`REPlusRepolarizationTraces.png` and `.svg` compare voltage at the same central
mesh vertex in absolute time and after alignment by local LAT. Markers indicate
the 90% crossings. The aligned traces make the AP-duration difference visible
separately from the earlier activation. The JSON record also contains nodal
APD90 ranges, latest crossing times and the probe coordinates and duration.
Repolarization is generated by the ionic kinetics and retained diffusion; FIM
supplies the activation schedule only.
