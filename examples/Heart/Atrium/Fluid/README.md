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

The equal-order pair is stabilized by term-by-term orthogonal subscales, the
default, or by residual-based SUPG and PSPG terms with
`-atrium_fluid_vms false`. The residual-based method adds

$$
\sum_K\int_K R_m\cdot
 \left[\tau_s\rho(w\cdot\nabla)v+\tau_p\nabla q\right]dX,
$$

and a grad-div term. The same time, advection, Temam, drag, and pressure
contributions enter both stabilization residuals. With element size $h_K$,
the parameters are

$$
\tau_1=\left(\frac{4\mu_n}{\rho h_K^2}
                    +\frac{2|w|_K}{h_K}+\frac{\sigma_B}{\rho}\right)^{-1},\qquad
\tau_s=\frac{s_s}{\rho/\Delta t+\rho/\tau_1},\qquad
\tau_p=\frac{s_p}{\rho/\Delta t+\rho/\tau_1}.
$$

The advecting speed $|w|_K$ and the Brinkman coefficient $\sigma_B$ are
sampled at the cell centroid, so $\tau_1$ and every parameter derived from it
are constant on a cell.

Backward Euler advances time. Boundary labels 1, 2, 3, and 5 are pressure-driven
inlets, label 6 is the outlet, and the extracted material-4 interface is the
fixed no-slip wall, or the Robin interface in [FSI](../FSI/README.md). Pressure
penalties, boundary impedance, and backflow terms are retained in the flow
problem. The default uses a prescribed outlet pressure. Ventricular coupling
and species transport are optional parameters, disabled by default and
available in fluid and FSI modes.

## Split orthogonal subscales

Term-by-term orthogonal subscales stabilize convection, divergence, and
pressure gradient separately. On the moving domain, with $a=u^n-w$ and mesh
velocity $w$, the stabilization is

$$
\sum_K\int_K\tau_s\rho^2\left[(\nabla u^{n+1})a-\xi^n
-\frac{1}{\Delta t}\tilde u^n_K\right]\cdot(\nabla v)a
+\tau_c\,(\nabla\cdot u^{n+1})(\nabla\cdot v)
-\sqrt{\tau_c}\,\pi^n\,\nabla\cdot v
+\tau_p\left(\nabla p^{n+1}-\zeta^n\right)\cdot\nabla q\,dX,
$$

$$
\xi^n=\Pi_h\!\left[(\nabla u^n)(u^n-w^n)\right],
\qquad
\pi^n=\Pi_h\!\left[\sqrt{\tau_c}\,\nabla\cdot u^n\right],
\qquad
\zeta^n=\Pi_h\!\left[\nabla p^n\right],
\qquad
\tau_c=\frac{s_c\rho h_K^2}{4\tau_1},\qquad \tau_p=\frac{s_p\tau_1}{\rho}.
$$

$\Pi_h$ is the $L^2$ projection onto the P1 space. The projections are
computed once per physical step on the previous configuration and are then
nodal fields carried by the mesh. `-atrium_fluid_pressure_projection false`
omits $\zeta^n$ and penalizes the whole pressure gradient.

The dynamic convective subscale obeys
$\rho\,\partial_t\tilde u+(\rho/\tau_1)\tilde u
=-\rho\,[(\nabla u)a-\Pi_h((\nabla u)a)]$ and is advanced by backward Euler,

$$
\tilde u^n_K=\tau_s\rho\left[\frac{\tilde u^{n-1}_K}{\Delta t}
-\left((\nabla u^n)a^n-\xi^n\right)(x_K)\right],\qquad
\tau_s=\frac{s_s}{\rho/\Delta t+\rho/\tau_1}.
$$

A subscale is orthogonal to the finite-element space, so its $L^2$ projection
onto that space would discard it. It is stored as one value per cell, the
cell mean of the expression above; the convective residual is affine on a
cell and its mean is the value at the centroid $x_K$.
`-atrium_fluid_quasi_static_subscales` omits the history and uses
$\tau_s=s_s\tau_1/\rho$.

## Lumped ventricle and mitral valve

`-atrium_fluid_lv_model` advances Rodin's lumped left-ventricle model once per
step, before the flow solve, and couples it to the mitral plane
$\Gamma_m$. The ventricular model contains the valve law

$$
Q=K\,(p_\mathrm{at}-p_v),\qquad
K=\begin{cases}K_\mathrm{at},&p_\mathrm{at}>p_v,\\K_p,&\text{otherwise,}\end{cases}
$$

with the open-valve conductance $K_\mathrm{at}$ and the leakage $K_p$. The
ventricle reads the mean pressure of the mitral plane at the last accepted
step as its atrial pressure,

$$
p_\mathrm{at}^{n}=\frac{1}{|\Gamma_m|}\int_{\Gamma_m}p^n\,dS,
$$

and the flow applies the same law as the resistive condition

$$
\sigma n=-\left(p_v^{n+1}+\frac{|\Gamma_m|}{K}\,u^{n+1}\cdot n\right)n
\quad\text{on }\Gamma_m,
$$

with $K$ selected by the test $p_\mathrm{at}^n>p_v^{n+1}$, which is the test
of the ventricular model. For a uniform normal velocity this is
$\bar p-p_v=Q/K$, so the volume leaving the atrium is the volume entering the
ventricle in both valve states. The closed valve is the pointwise penalty
form of the zero-flux constraint: the constant $|\Gamma_m|/K_p$ bounds the
normal velocity, and the mean pressure of the plane plays the role of the
multiplier that reopens the valve. The exchange is explicit in the
ventricular pressure and implicit in the resistance.

`-atrium_fluid_lv_feedback false` drives the ventricle with its prescribed
atrial waveform instead; the flow then does not feed back to it, and the
resistances are `mitralResistanceOpen` and `mitralResistanceClosed`.
`-atrium_fluid_pressure_ramp` scales the venous pressure linearly from zero
over the given time. A coupled ventricle keeps its own pressure during the
ramp, since its valve remains closed while the atrium is below it; the
one-way ventricular pressure is scaled as well.

## Coagulation species

`-atrium_fluid_kinetics` transports prothrombin $P$, thrombin $T$, antithrombin
$A$, fibrinogen $G$, and factor Xa $X$ after `-atrium_fluid_flow_cycles`
completed cycles (default three). With $\mathcal L c=\partial_tc|_\chi
+(u-w)\cdot\nabla c-\nabla\cdot(D\nabla c)$,

$$
\mathcal LP=-R,\qquad \mathcal LT=R-R_I,\qquad
\mathcal LA=-R_I,\qquad \mathcal LG=-R_F,\qquad
\mathcal LX=-k_{dX}X,\qquad \partial_tF|_\chi=R_F,
$$

$$
R=\left(k_a\frac{T^2}{K_A+T}+\kappa X\right)\frac{P}{P_\mathrm{ref}},\qquad
R_I=k_\mathrm{on}AT,\qquad R_F=k_F\frac{TG}{K_G+G}.
$$

Nodal concentrations move with the mesh, so the backward-Euler difference is
the time derivative at fixed mesh coordinate and the convecting velocity is
$u-w$. Transport uses SUPG and Codina crosswind diffusion with one quadrature
rule for all terms; it preserves uniform states and therefore the invariant
$P+T-A$. A nodal production–destruction step follows, positive and
budget-exact for any time step; it conserves $P+T-A$ and $G+F$.

The equal-order velocity is not pointwise divergence free. On the moving
domain the balance of $\int c^2$ contains $\nabla\cdot(u-w)$ from the advection
and $\nabla\cdot w$ from the domain, so a deviation from the inflow state
grows at the rate $|\nabla\cdot u_h|/2$. The transport therefore adds

$$
\frac12\int_{\Omega}(\nabla\cdot u_h)\,(c-c_\mathrm{in})\,v\,dX,
$$

the skew-symmetric form about the inflow concentration $c_\mathrm{in}$. The
inflow state remains an exact solution, and the operator acting on the
deviation is energy neutral. `-atrium_fluid_ptag_skew false` removes the term.

Initiation passes through factor Xa. The endocardium emits $j_XE$, where $E$
is the wall activation built from TAWSS and OSI; in an endocardial layer of
thickness $\delta$ the emission is the volumetric source
$(j_X/\delta)\chi\tilde E$, where $\tilde E$ extends $E$ into the lumen and
$\chi$ normalizes the dose, and $\delta=0$ selects the wall flux
$-D\,\partial_nX=j_XE$. Factor Xa is a catalyst, so $P+T-A$ is unaffected.
Without flow $X$ tends to $j_X\chi\tilde E/(\delta k_{dX})$; the default
$j_X=j_0k_{dX}/\kappa$ then reproduces the direct channel
$(j_0/\delta)\chi\tilde EP/P_\mathrm{ref}$, and the two initiations differ
only where the flow washes $X$ out. `-atrium_fluid_ptag_xa false` selects the
direct channel $P\to T$. `-atrium_fluid_ptag_e_osi false` uses
$E=s(\mathrm{TAWSS})$ without the factor $\min(2\,\mathrm{OSI},1)$.

Fibrin is immobile: $F$ is a nodal field advanced by the reaction step and
carried by the mesh nodes. `-atrium_fluid_ptag_immobile_fibrin false` uses
the advected surrogate $F=G_0-G$. The drag $\sigma_B(F)\,u$ enters the
momentum equation, implicit in $u$ and lagged in $F$. `Coagulation.h` lists
the constants.

| Option | Symbol | Default |
|---|---|---|
| `-atrium_fluid_ptag_j_x` | $j_X$ ($\mathrm{mol\,m^{-2}\,s^{-1}}$) | $j_0k_{dX}/\kappa$ |
| `-atrium_fluid_ptag_kappa_x` | $\kappa$ ($\mathrm{s^{-1}}$) | $3$ |
| `-atrium_fluid_ptag_k_dx` | $k_{dX}$ ($\mathrm{s^{-1}}$) | $10^{-2}$ |
| `-atrium_fluid_ptag_k_on` | $k_\mathrm{on}$ ($\mathrm{m^3\,mol^{-1}\,s^{-1}}$) | $7.1$ |
| `-atrium_fluid_ptag_P0`, `-atrium_fluid_ptag_A0` | $P_0$, $A_0$ ($\mathrm{mol\,m^{-3}}$) | $1.4\times10^{-3}$, $2.8\times10^{-3}$ |
| `-atrium_fluid_ptag_tau_a` | activation shear stress (Pa) | $0.4$ |

The decay time of factor Xa is $1/k_{dX}=100\,\mathrm{s}$, so $X$ approaches
its stagnant level only after a comparable simulated time.

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
`Fibrinogen`, `FactorXa`, `Fibrin`, `ThrombinAntithrombin`, `BrinkmanDrag`,
`InvariantDrift`, `ThrombosisActivation`, and `ActivationLayer`; it is separate
from solid activation. Species solvers use the `-atrium_ptag_` option prefix
and default to GMRES with block Jacobi.

## Assembly

The lagged viscosity, the cell size, and $\tau_1$ are constant on a cell.
They are evaluated once per cell for each lagged state and geometry and
shared by the viscous term, every stabilization coefficient, and the
interface tangent. With cellwise coefficients every polynomial integrand of
the affine-P1 Oseen form has degree at most two, and the four-point rule
integrates it exactly; the fibrin drag is the only non-polynomial coefficient.

The mass operators of the $L^2$ projections and of the wall-shear recovery
depend on the geometry only. They are assembled once per geometry, and later
projections on that geometry assemble the right-hand side alone. On a fixed
domain they are assembled once per run.

The scheme is linear in each step. The advecting velocity, the viscosity,
the stabilization parameters, the three projections, the velocity subscale,
the fibrin of the drag, the valve state, and the ventricular pressure are
taken from the previous step.

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
