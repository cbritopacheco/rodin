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
fixed no-slip wall, or the Robin interface in [FSI](../FSI/README.md). Pressure
penalties, boundary impedance, and backflow terms are retained in the flow
problem. The default uses a prescribed outlet pressure. Ventricular coupling,
species transport, and orthogonal subscales are optional parameters, disabled
by default and available in fluid and FSI modes.

## Split orthogonal subscales

`-atrium_fluid_vms` replaces SUPG/PSPG by term-by-term orthogonal subscales.
On the moving domain, with $a=u^n-w$ and mesh velocity $w$, the stabilization is

$$
\sum_K\int_K\tau_s\rho^2\left[(\nabla u^{n+1})a-\xi^n\right]\cdot(\nabla v)a
+\tau_c\,(\nabla\cdot u^{n+1})(\nabla\cdot v)
-\sqrt{\tau_c}\,\pi^n\,\nabla\cdot v
+\tau_p\nabla p^{n+1}\cdot\nabla q\,dX,
$$

$$
\xi^n=\Pi_h\!\left[(\nabla u^n)(u^n-w^n)\right]+\frac{1}{\Delta t}\tilde u^n,
\qquad
\pi^n=\Pi_h\!\left[\sqrt{\tau_c}\,\nabla\cdot u^n\right],
\qquad
\tau_c=\frac{s_c\rho h_K^2}{4\tau_1},\qquad \tau_p=\frac{s_p\tau_1}{\rho}.
$$

$\Pi_h$ is the $L^2$ projection onto the P1 space and $\tilde u^n$ the dynamic
convective subscale; `-atrium_fluid_quasi_static_subscales` omits its history
and uses $\tau_s=s_s\tau_1/\rho$. The projections are computed once per
physical step on the previous configuration and are then nodal fields carried
by the mesh. Convection, divergence, and pressure gradient are stabilized
separately; the pressure-gradient term has no projection. With the Brinkman
drag, $\tau_1^{-1}$ includes $\sigma_B/\rho$.

## Lumped ventricle

`-atrium_fluid_lv_model` advances Rodin's lumped left-ventricle model once per
step and applies its pressure at the outlet through a diode resistance: open
while the mean outlet pressure of the previous step exceeds the ventricular
pressure, closed otherwise. The ventricle is driven by its own prescribed
atrial pressure waveform; the three-dimensional flow does not feed back to it.
`-atrium_fluid_pressure_ramp` scales the venous and ventricular pressures
linearly from zero over the given time.

## Coagulation species

`-atrium_fluid_kinetics` transports prothrombin $P$, thrombin $T$, antithrombin
$A$, and fibrinogen $G$ after `-atrium_fluid_flow_cycles` completed cycles
(default three). With $\mathcal L c=\partial_tc|_\chi+(u-w)\cdot\nabla c
-\nabla\cdot(D\nabla c)$,

$$
\mathcal LP=-R,\qquad \mathcal LT=R-R_I,\qquad
\mathcal LA=-R_I,\qquad \mathcal LG=-R_F,
$$

$$
R=k_a\frac{T^2}{K_A+T}\frac{P}{P_\mathrm{ref}},\qquad
R_I=k_\mathrm{on}AT,\qquad R_F=k_F\frac{TG}{K_G+G}.
$$

Nodal concentrations move with the mesh, so the backward-Euler difference is
the time derivative at fixed mesh coordinate and the convecting velocity is
$u-w$. Transport uses SUPG and Codina crosswind diffusion with one quadrature
rule for all terms; it preserves uniform states and therefore the invariant
$P+T-A$. A nodal production–destruction step follows, positive and
budget-exact for any time step. Initiation is the volumetric channel
$(j_0/\delta)\tilde EP/P_\mathrm{ref}$ in an endocardial layer of thickness
$\delta$, where $\tilde E$ extends the wall activation built from TAWSS and OSI;
$\delta=0$ selects the wall flux. The fibrin surrogate $F=G_0-G$ adds the drag
$\sigma_B(F)\,u$ to the momentum equation, implicit in $u$ and lagged in $F$.
`Coagulation.h` lists the constants.

Open pressure boundaries permit flow reversal. With outward normal $n$ and
relative advecting velocity $b=u^n-a$, both inlet and outlet residuals include

$$
\int_{\Gamma_o}\frac{\rho}{2}[-b\cdot n]_+\,u^{n+1}\cdot v\,dS.
$$

The correction acts on the full velocity, including normal and tangential
components. For fixed openings, its sum with the convective boundary-energy
term is $\frac{\rho}{2}\int_{\Gamma_o}[b\cdot n]_+|u^{n+1}|^2\,dS\geq0$
when the backflow coefficient is one. It supplies no prescribed flow and does
not introduce a valve. Opening impedance and tangential damping, when enabled,
remain separate terms.

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
at completed statistics cycles and vanish away from the wall. Optional
coagulation transport exports `Prothrombin`, `Thrombin`, `Antithrombin`,
`Fibrinogen`, `Fibrin`, `ThrombinAntithrombin`, `BrinkmanDrag`,
`InvariantDrift`, `ThrombosisActivation`, and `ActivationLayer`; it is separate
from solid activation. Species solvers use the `-atrium_ptag_` option prefix
and default to GMRES with block Jacobi.

The lagged viscosity and the cell size are constant on an affine P1
tetrahedron. They are evaluated once per cell for each lagged state and
geometry and shared by all stabilization coefficients and by the viscous term.

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
