# Anisotropic reaction–eikonal propagation

`Eikonal_ReactionEikonal` implements the RE+ structure in equation (12) of
Nagel et al., *Comparison of Propagation Models and Forward Calculation Methods
on Cellular, Tissue and Organ Scale Atrial Electrophysiology*, IEEE TBME 70(2),
2023, [doi:10.1109/TBME.2022.3196144](https://doi.org/10.1109/TBME.2022.3196144).
An anisotropic FIM solve schedules a foot current in a Courtemanche
reaction–diffusion simulation. A separate tensor FMM solve checks the eikonal
field. Monodomain strand simulations calibrate the squared-velocity tensor.

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

The fiber directions are \(f=(0.8,0.6)^T\) and \(s=(-0.6,0.8)^T\), with
\(R=[f,s]\). An equilateral triangular lattice spans the reference parallelogram with
vertices \(0\), \((10,0)\), \((5,5\sqrt{3})\), and
\((15,5\sqrt{3})\) mm. It is mapped by
\(x=R\operatorname{diag}(r,1)\xi\). Both the geometry and diffusion tensor
are transformed consistently. This also provides a useful numerical control:
changing \(r\) changes the physical propagation distances, while preserving
the normalized reference problem. The triangles are strictly acute in the
calibrated inverse-velocity metric, as required by tensor FMM.

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
probes near reference coordinates 3 and 7 mm:

\[
c_a=\frac{x_{a,2}-x_{a,1}}{t_{a,2}-t_{a,1}},
\qquad M=R\operatorname{diag}(c_l^2,c_t^2)R^T.
\]

The strand spacing matches the nominal reference edge length, stretched
by the corresponding principal factor. Each strand
receives a 20 pA/pF depolarizing pulse for 2 ms on its first reference 0.5 mm.
The sheet source set contains vertices within reference radius 0.5 mm of the
corner. FIM onset times drive the foot current at every sheet node; the RE+
sheet does not receive the rectangular monodomain stimulus.

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
for repolarization studies. APD and ECG metrics are not calculated.

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
FIM/FMM foot-onset fields and onset-to-upstroke lag. The associated CSV gives
these quantities at every vertex. `ReactionEikonalCell.csv` records the final
conditioning beat. The console reports calibrated velocities, FIM rounds,
FIM/FMM differences and lag statistics outside the stimulated region.
Plot the exported native triangle fields with NumPy, Matplotlib and h5py:

```sh
python3 /path/to/rodin/examples/Models/Eikonal/plot_reaction_eikonal.py . --time 20
```

The script produces `ReactionEikonal.png` and `ReactionEikonal.svg`.

Invalid ionic states, unactivated tissue or disagreement between FIM and FMM
on this matched stencil cause failure. Agreement with FMM checks the eikonal
calculation; it does not validate the ionic kinetics or physiology.

## Observed numerical checks

The following runs used ten conditioning beats and an 80 ms activation window.
All nodes activated, all ionic states remained finite, and concentrations stayed
positive. Both principal velocities were recalibrated for every run.

| Intervals | Time step (ms) | Velocity ratio | Transverse velocity (mm/ms) | Latest RE+ LAT (ms) | Maximum FIM/FMM difference (ms) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 40 | 0.02 | 1 | 0.361664 | 45.110 | \(1.55\times10^{-12}\) |
| 40 | 0.02 | 2 | 0.361664 | 45.110 | \(1.31\times10^{-12}\) |
| 40 | 0.01 | 2 | 0.364964 | 44.655 | \(3.49\times10^{-12}\) |
| 80 | 0.02 | 2 | 0.400000 | 41.290 | \(4.52\times10^{-12}\) |

Four FIM workers reproduced the default CSV exactly. Comparing corresponding
vertices in the isotropic and stretched cases gave identical RE+ LATs and FIM
onset differences below \(3\times10^{-14}\) ms. The physical coordinates
differ according to the prescribed stretch. In the default case,
\(t_{\mathrm{LAT}}-T\) ranged from approximately -1.808 to 1.610 ms;
negative values are consistent with activation by diffusion before the local
foot current begins.

A separate Python evaluation of the original CellML right-hand sides was
integrated with adaptive BDF at relative tolerance \(10^{-9}\) and absolute
tolerance \(10^{-11}\). For one 2 ms, 20 pA/pF stimulus from the published
resting state, its peak voltage was approximately 24.492 mV. The example's
frozen-gate step produced peaks 26.292, 25.382 and 24.935 mV at steps 0.02,
0.01 and 0.005 ms. Voltage at 300 ms approached the reference -68.8379 mV
through -68.9665, -68.9130 and -68.8759 mV. This checks approach to the
source ODE trajectory; it is not validation against experimental data.
The observed changes with mesh and time step prevent treating the default
run as a numerically converged reproduction of the paper. A separate 500 ms
run with 20 intervals completed through repolarization; final voltages ranged
from approximately -78.98 to -78.55 mV.

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
