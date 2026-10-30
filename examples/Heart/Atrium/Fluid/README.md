# Atrium fluid

The fluid uses continuous linear velocity and pressure on affine tetrahedra.
For lagged velocity $w=u^n$ and viscosity $\mu_n$, the Oseen momentum residual is

$$
R_m=\rho\frac{u^{n+1}-u^n}{\Delta t}
 +\rho(w\cdot\nabla)u^{n+1}
 +\frac{\rho}{2}(\nabla\cdot w)u^{n+1}
 -\nabla\cdot(2\mu_n D(u^{n+1}))+\nabla p^{n+1},
\qquad D(u)=\frac{\nabla u+\nabla u^T}{2}.
$$

The Temam term accounts for a lagged velocity that is not exactly divergence
free. Carreau–Yasuda viscosity is evaluated at the lagged strain rate and is
constant within each affine linear tetrahedron. The elementwise strong viscous
divergence therefore vanishes in the stabilization residual.

The equal-order method adds SUPG and PSPG terms

$$
\sum_K\int_K R_m\cdot
 \left[\tau_s\rho(w\cdot\nabla)v+\tau_p\nabla q\right]dX,
$$

and a grad-div term. The same time, advection, Temam, and pressure contributions
enter both stabilization residuals. With element size $h_K$, the parameters are

$$
\tau_1=\left(\frac{4\mu_n}{\rho h_K^2}
                    +\frac{2|w|}{h_K}\right)^{-1},\qquad
\tau_s=\frac{s_s}{\rho/\Delta t+\rho/\tau_1},\qquad
\tau_p=\frac{s_p}{\rho/\Delta t+\rho/\tau_1}.
$$

Backward Euler advances time. Boundary labels 1, 2, 3, and 5 are pressure-driven
inlets, label 6 is the outlet, and the extracted material-4 interface is the
fixed no-slip wall. Pressure penalties, boundary impedance, and backflow terms
are retained in the flow problem. The default uses a prescribed outlet pressure.
Ventricular coupling, species transport, and orthogonal-subscale corrections are
optional parameters, disabled by default.

All reported opening fluxes use the outward fluid normal: inflow is negative
and outflow is positive. For pressure-penalty coefficient $\epsilon_p$, the
constant pressure test gives the global balance

$$
Q_{\mathrm{in}}+Q_{\mathrm{out}}
=-\epsilon_p\int_{\Omega_f}p\,dX.
$$

This global balance does not enforce pointwise zero divergence. The reported
maximum velocity is $\max_X\|u_h(X)\|_2$, evaluated from nodal speeds; for
an affine P1 tetrahedral field the maximum occurs at a vertex.


## Architecture

`Fluid` owns the fluid submesh, P1 velocity and pressure spaces, the retained
Oseen `Problem`, KSP, and wall-shear and projection Problems. `Output` is nested
as `Fluid::Output` and exports observed fields. The common coordinator supplies
one `RunParameters` configuration and advances the component through
`initialize()`, `step()`, and `finish()`. Physical and stabilization coefficients
belong to `Fluid::Parameters`; common time and output settings are defined by
[the run configuration](../README.md#configuration).

`Velocity`, `Pressure`, `WallShearStress`, `TimeAveragedWallShearStress`, and
`OscillatoryShearIndex` are exported in `Fluid.xdmf`. The latter two are updated
at completed statistics cycles. Optional coagulation transport exports
`Thrombin`, `Fibrinogen`, `Fibrin`, and `ThrombosisActivation`; it is disabled by
default and is separate from solid activation. The 0D ventricular driver and
orthogonal subscales are also disabled by default and available only in fluid
mode.

## Run

```sh
cmake --build build --target Atrium -j4
mpirun -n 4 build/examples/Heart/Atrium/Atrium \
  -atrium_mode fluid \
  -atrium_mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  -atrium_output output/atrium-fluid \
  -atrium_dt 0.001 -atrium_steps 850 -atrium_output_every 1 \
  -atrium_fluid_lv_model false -atrium_fluid_kinetics false
```

The run advances one default pressure period of $0.85\,\mathrm{s}$. The initial
rest state and every advanced state are saved. A single cycle from rest does
not imply a periodic solution. The shared `Diagnostics.csv` records pressure,
flux, speed, cavity volume, solver iterations, and preconditioner statistics at
every step. [Volume rendering](../Visualization/README.md) uses `Fluid.xdmf`.

`-atrium_fluid_max_velocity` sets the numerical divergence guard in metres per
second for fluid and FSI trials. Its default is $1000\,\mathrm{m/s}$, three
orders above a $1\,\mathrm{m/s}$ reference speed. This heuristic guard does not
define physiological admissibility. Nonfinite peak speeds and speeds above the
guard reject the trial; its measured peak and the configured guard are logged.

## Preconditioner refresh

`-atrium_fluid_lag_preconditioner N` refreshes numeric factors on the first
step and every $N$ physical timesteps. Its default is one; larger values require
an iterative KSP. Within FSI, trials retain factors on that schedule; one
refreshes every trial. A failed retained solve retries the same system once
with fresh factors. Symbolic analysis can be retained while sparsity remains
unchanged. This policy takes precedence over `-ksp_reuse_preconditioner` for the
flow KSP and does not reuse numeric projection operators.
