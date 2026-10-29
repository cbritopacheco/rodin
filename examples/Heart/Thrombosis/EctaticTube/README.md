# Ectatic tube flow, coagulation, and washout

A tube with an enlarged segment supplies a fixed geometry in which residence
and coagulation can be examined together. The fluid is solved at every time
step and transports thrombin $T$, fibrinogen $G$, fibrin $F$, and an independent
washout tracer $c$. Transport is one-way: concentrations do not change the flow.

## Fluid and species models

Continuous P1 velocity and pressure solve a transient, shear-thinning Oseen
problem with equal-order pressure stabilization, orthogonal-subscale VMS, and
grad-div terms. Carreau–Yasuda stress is
$\sigma=-pI+2\mu(\dot\gamma)D(u)$. The wall is fixed and no-slip. A prescribed
pulsatile inlet pressure has impedance and tangential damping; the outlet is
traction-free with backflow stabilization. There is no ventricular or lumped
outlet circulation model.

With reaction constant $k$ and species diffusivities $D_i$, the transport balances are

$$
\partial_tT+u\cdot\nabla T-D_T\Delta T=0,\qquad
\partial_tG+u\cdot\nabla G-D_G\Delta G=-kGT,\qquad
\partial_tF+u\cdot\nabla F-D_F\Delta F=kGT.
$$

Activated endothelium supplies thrombin through the boundary flux
$D_T\nabla T\cdot n=j_T A$; $A$ depends on completed-cycle wall-shear indices.
Fibrinogen inflow supplies the feed concentration; thrombin and fibrin inflow
values are zero. [`ThrombosisParameters`](../../CoronaryArtery/ThrombosisModel.h)
defines concentrations, reaction rates, diffusion, and activation parameters.
The model does not include platelets, lysis, clot permeability, or mechanical
occlusion.

The washout equation and observable are

$$
\partial_tc+u\cdot\nabla c=0,\qquad c(X,0)=1,\qquad
c|_{\Gamma_{\mathrm{in}}}=0,\qquad
R(t)=\frac{1}{|\Omega|}\int_\Omega c(X,t)\,dX.
$$

Thus $R$ measures the remaining fraction of the original bolus, independently
of thrombin production or fibrinogen conversion.

## Time integration, stabilization, and output

Backward Euler advances the P1 transport fields. SUPG terms test their transient,
advective, diffusive, and reaction residuals; residual-based crosswind diffusion
controls sharp concentration fronts. Numerical stabilization does not by itself
prove positivity or a concentration maximum principle.

Flow warm-up precedes species activation. Cycle-integrated wall shear provides
TAWSS and OSI, used by the wall-source indicator. The geometry is a conduit
with an ectatic segment, not a patient-specific coronary reconstruction.
The generated mesh uses metre coordinates; velocity, pressure, time, and
concentration use SI units, including $\mathrm{mol/m^3}$ for species.
XDMF stores flow, species, and washout fields. Scalar histories report fluxes,
wall-shear indices, concentration content, and the volume-weighted washout $R$.

Requires `RODIN_USE_PETSC=ON`; MPI is not required.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target EctaticTubeThrombosis -j4
cd build
./examples/Heart/Thrombosis/EctaticTube/EctaticTubeThrombosis
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/EctaticTubeThrombosis.mesh`. Generate it from the repository root with:

```sh
python3 examples/Heart/Thrombosis/EctaticTube/make_test_mesh.py resources/examples/Heart/EctaticTubeThrombosis.mesh
```

The mesh can be selected with `-pt_mesh`. The `-pt_flow_cycles` option controls
flow warm-up, `-pt_species_cycles` controls the subsequent transport interval,
and `-pt_dt` selects the time step in seconds.
