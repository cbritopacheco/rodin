# Coronary rigid-wall CFD

Coronary flow is driven by a lumped left ventricle and loaded by intramyocardial
outlet beds. A fixed lumen separates perfusion and rheology from wall deformation.
`CoupledLV0DCoronary3D` owns the ventricular state, outlet compartments, retained
flow/projection problems, and wall-shear observables.

## Flow and circulation model

In the fixed lumen $\Omega_f$, velocity $u$ and pressure $p$ satisfy

$$
\rho_f(\partial_tu+u\cdot\nabla u)-\nabla\cdot\sigma_f=0,\qquad
\nabla\cdot u=0,\qquad
\sigma_f=-pI+2\mu(\dot\gamma)D(u).
$$

The wall has zero velocity. The inlet pressure is supplied by the 0D arterial
state, with normal impedance and tangential damping. Backflow terms stabilize
the pressure boundaries. Carreau–Yasuda viscosity depends on the regularized
strain-rate magnitude; its definition is given in the
[coronary fluid model](../FSI/README.md#component-models).

For each outlet, let $Q$ be outward 3D flux, $p_c$ compartment pressure,
$p_{\mathrm{im}}$ intramyocardial pressure, and $p_{\mathrm{tm}}=p_c-p_{\mathrm{im}}$.
The compartment and proximal relation are

$$
C\dot p_{\mathrm{tm}}=Q-q_v,\qquad
p_{\mathrm{im}}=\chi p_v,\qquad
p_{\mathrm{out}}=p_c+R_a\Phi_a(Q)Q.
$$

A regularized Starling venous law defines $q_v$ from compartment/drain pressures
and lumen state. The viscosity multipliers $\Phi_a,\Phi_v$ come from apparent
viscosities normalized by the calibration viscosity. The scalar compartment
update uses backward Euler. Proximal resistance is represented by a normal
boundary impedance using the opening area; its local velocity form is a
reduced outlet-profile approximation.

## Numerical procedure and output

Continuous quadratic velocity and linear pressure form a Taylor–Hood pair.
Backward Euler, Temam convection, orthogonal-subscale VMS and grad-div terms
are assembled. `FlowMode::Oseen` freezes advection and viscosity; `FlowMode::Newton`
uses a nonlinear residual and tangent, including the viscosity derivative.
A small pressure penalty regularizes the pressure block. Step reduction restores
the flow state after an unsuccessful 3D solve.

XDMF contains velocity, pressure, and wall-shear fields. CSV records
ventricular/outlet states, fluxes, and solver observables. Source mesh coordinates use the configured conversion to metres
(default factor $10^{-2}$); solved quantities use SI units. No wall displacement
or ALE problem is present.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/CFD/CoronaryArtery
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/CoronaryArtery/Output3D.mesh`. See `CoronaryArtery.cpp` for the PETSc options, including `-coronary_dt`, `-coronary_nsteps`, and `-coronary_output_prefix`.

## Outlet closures and calibration

The bed resistances of every outlet are calibrated once, at rest, from the
Murray split of `lcaTargetFlow` and the pressure budget. The viscosity at which
they are defined is set with `-coronary_calibration`:

- `newtonian` (default) defines them at `-coronary_calibration_viscosity`
  ($\mu_N$, 3.5 mPa s).
- `reference` defines them for the reference blood at its resting operating
  point. The reference blood is the healthy Carreau–Yasuda fit (override with
  `-coronary_cal_cy_*`), and the normalisation is
  $\mu_{\mathrm{cal},k}=\mu_{\mathrm{ap,ref}}(4v_k/r_k)$. The reference blood
  then reproduces the prescribed resting flow ($\Phi=1$), and any other blood
  perfuses the same bed.

Each limb may be a series of generations, set with
`-coronary_arteriolar_tree` or `-coronary_venular_tree` as `r:v:f,...` (m,
m/s, share of the limb's pressure drop). Its factor is

$$
\Phi(Q)=\sum_k f_k\,\mu_{\mathrm{ap}}(g_k\lvert Q\rvert/q_0)/\mu_{\mathrm{cal},k},
\qquad g_k = 4v_k/r_k .
$$

This is a resistance-weighted spectrum of resting shear rates. A single
generation recovers the single parallel-tube compartment.

`-coronary_outlet_closure` selects the reduced closure:

- `rmu` is the shear-dependent closure $\Phi(Q)$.
- `rest` freezes $\Phi$ at the resting operating point of the running blood, so
  it coincides with `rmu` at rest.
- `plateau` uses $\mu_\infty$ for every generation.

Both constant closures are multiplied by `-coronary_outlet_constant_scale`.

The calibration can be written with `-coronary_calibration_out` and read back
with `-coronary_calibration_in`. Outlet areas are always measured on the mesh,
so a mesh study can hold the terminal parameters fixed.

Every step also appends the total dissipation, the wall-mean $\lvert\tau_w\rvert$
and the mean inlet and outlet pressures to `CoronaryArtery_fields.csv`.
`run_coronary_closures.sh` runs the closure comparison for the three blood
fits, with optional time-step, mesh, drainage-pressure and multi-generation
studies, and with a secant-matched constant resistance. Its last step,
`compare_coronary_closures.py`, reduces the last cycle of each run.

