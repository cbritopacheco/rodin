# Atrium FSI

`FSI` couples the Demiray–Holzapfel–Gasser–Ogden solid with Kelvin–Voigt
dissipation to the ALE Oseen fluid. By default, active wall contraction drives
an initially resting fluid through openings at zero pressure. With
`-atrium_fluid_lv_model`, the pulmonary venous pressure and the lumped left
ventricle behind the mitral diode load the same coupled problem. The fluid uses
P1/P1 velocity and pressure with either residual-based SUPG/PSPG or
term-by-term orthogonal subscales, and optionally transports the P–T–A–G
coagulation species on the moving mesh.

## Interface model

The solid and fluid submeshes retain their common parent vertex identities.
Their interface is conforming, and retained MPI scatter plans transfer its
fields. Let $n_f$ point out of the current fluid domain, $N_s$ out of the
reference solid, and $J_\Gamma=dS/dS_0$. The target conditions are

$$
d_f=d_s,\qquad u_f=v_\Gamma,\qquad
P_sN_s=J_\Gamma(-\sigma_fn_f),\qquad
\sigma_f=-pI+2\mu D(u_f).
$$

The interface rate uses BDF1,

$$
v_\Gamma=\frac{d_s^{n+1}-d_s^n}{\Delta t}.
$$

Solid bulk dynamics retain Newmark integration. Its endpoint velocity differs
from the BDF1 interface rate at finite timestep. A harmonic extension moves the
fluid interior; opening caps and solid opening cut faces remain fixed. No
additional interface ring is clamped. Without the ventricular model, all fluid
openings have zero imposed pressure and no prescribed flux. With it, the
openings retain the venous pressure, impedances, and mitral diode of the
[fluid model](../Fluid/README.md); the lumped ventricle advances once per
physical step, before the coupling iterations, and its pressure and valve
state remain fixed during them.

For reference-area impedance $\alpha_0>0$, a solid/ALE/fluid evaluation solves

$$
P_sN_s+\alpha_0v_\Gamma=g.
$$

The fluid weak residual contains the Robin and interface convection terms

$$
\int_{\Gamma(t)}\frac{\alpha_0(u_f-v_\Gamma)+T}{J_\Gamma}\cdot w\,dS
-\frac{\rho_f}{2}\int_{\Gamma(t)}[(u_f^n-a)\cdot n_f](u_f\cdot w)\,dS,
$$

where $a$ is the ALE mesh velocity. Division by $J_\Gamma$ converts reference
Robin data to the current fluid surface. The interface velocity is not imposed
as a strong fluid Dirichlet condition.

For analytic Newton, the variational Robin reaction defines the exchange map

$$
H_T(g,T)=\Pi_\Gamma T+\alpha_0(u_f-v_\Gamma),\qquad
H_g(g,T)=H_T(g,T)+\alpha_0u_f.
$$

Here $\Pi_\Gamma$ restricts the transferred traction to the common interface.
Reference traction is exchanged in kPa and solid velocities in mm/s; fluid
assembly converts to Pa and m/s. In the Newton path, recovered velocity gradients serve wall-shear
diagnostics; they do not define the exchanged traction.

The automatic impedance is

$$
\alpha_0=\sqrt{\rho_{s,0}E_m},\qquad
E_m=\frac{9K_m\mu_m}{3K_m+\mu_m},\qquad
\mu_m=a_m,\quad K_m=2\mu_v.
$$

Material coefficients are converted to Pa. The incompressible mixed wall uses
$E_m=3\mu_m$. This matrix-based scale omits
fiber reinforcement and is an algorithmic choice. `-atrium_fsi_robin` supplies
an explicit value in $\mathrm{Pa\,s/m}$.

## Analytic interface Newton

With `-atrium_fsi_newton`, the interface unknown $x=(g,T)$ solves

$$
R(x)=\begin{pmatrix}g-H_g(g,T)\\T-H_T(g,T)\end{pmatrix}=0,\qquad
R'(x_k)\delta x=-R(x_k).
$$

`FSI::solveNewton()` uses SNES with line search and a matrix-free GMRES
operator. Each analytic action solves the solid tangent problem, harmonically
extends the displacement increment, and solves the differentiated fluid
problem. The fluid derivative includes geometry and surface measure, ALE mass
ratios, relative convection, Temam terms, Carreau–Yasuda viscosity and
SUPG/PSPG/grad-div coefficients. Residual and derivative volume forms share
order-four quadrature. No finite differences construct the solver operator.

With split orthogonal subscales, the differentiated residual contains the
streamline term on convection minus its projection, the grad-div term minus
the projected divergence, and the pressure-gradient term. The projections
belong to the previous time level and are nodal fields carried by the mesh;
they are evaluated once per physical step on the previous configuration, so
the tangent contains no derivative of a projection. The Brinkman coefficient
of the coagulation model enters in the same way. The lumped ventricle and the
opening terms act on fixed caps and contribute no geometric derivative.

An evaluation retains its input, residual and solved bulk state. Exact equality
with the most recently evaluated input allows reuse during initial residual
requests, Jacobian setup and final acceptance. A different trial point requires
a fresh evaluation. This cache lives for one physical timestep. Tangent
matrices are assembled at the current Newton point; zero-direction diagnostic
solves are omitted.

`-atrium_fsi_atol` defaults to $10^{-6}$ and applies to the global Euclidean norm
of the two nodal Robin residual blocks, in kPa. This norm is mesh dependent;
it is not a normalized interface error. `-atrium_fsi_solid_atol` defaults to
$10^{-10}$ for the inner solid solve. Native interface SNES/KSP options use
`-atrium_fsi_interface_`, for example
`-atrium_fsi_interface_ksp_monitor`. A trial exceeding the fluid velocity guard is reported to SNES as a function-domain
error, allowing line-search backtracking. Such trials do not publish a valid
evaluated-state cache or commit physical histories. An invalid initial state,
unrecoverable solver error, or Newton failure rejects the step.
`-atrium_fsi_check_tangent` enables residual Taylor checks for verification.

## Budgeted Robin sweeps

Without `-atrium_fsi_newton`, each sweep solves the solid with Robin data
$g=T+\alpha_0u_f$, moves the mesh, solves the fluid, and updates the traction
with the same variational reaction as Newton,

$$
T\leftarrow T+\alpha_0(u_f-v_\Gamma).
$$

Its fixed point satisfies $u_f=v_\Gamma$ and coincides with the Newton solution.
One sweep per step (`-atrium_fsi_iterations 1`) is the loosely coupled
Robin–Robin scheme of Burman, Durst, Fernández, Guzmán and Ruz; further sweeps
are its correction iterations. `-atrium_fsi_recovered_traction` instead
recovers the fluid velocity gradient and projects physical stress to the
reference interface with lumped trace mass; that exchange leaves an interface
slip at its fixed point. Starting at the second completed sweep, the relative
solid displacement-change norm determines early termination:

$$
r_d=\frac{\|d_s^{k+1}-d_s^k\|_2}{\|d_s^{k+1}\|_2}.
$$

At zero denominator, the absolute increment norm in metres is used. Defaults
are `-atrium_fsi_rtol 1e-6` and `-atrium_fsi_iterations 5`; zero tolerance uses
the full budget. Budget exhaustion accepts the last complete sweep and marks
it unconverged. A small budget does not establish interface continuity.

## Ownership and output

`FSI` retains trace scatters, interface data and the harmonic mesh-motion
Problem. `Solid` and `Fluid` own their assembly Problems and histories. Previous
physical states remain fixed through all coupling trials. Mesh-coordinate
updates invalidate geometric caches. Accepted steps commit both components
and write one common `Diagnostics.csv`, including cavity volumes, fluxes,
interface mismatch, coupling status and mesh Jacobian diagnostics.

`Solid.xdmf` uses reference coordinates and `Displacement`; apply **Warp By
Vector** for deformation. `Fluid.xdmf` already contains moving coordinates.

```sh
cmake --build build --target Atrium -j4
mpirun -n 6 build/examples/Heart/Atrium/Atrium \
  -atrium_mode fsi -atrium_fsi_newton \
  -atrium_mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  -atrium_output output/atrium-fsi-newton \
  -atrium_dt 0.001 -atrium_steps 850 -atrium_output_every 1 \
  -atrium_solid_contractility 2000 -atrium_solid_cache_constitutive \
  -atrium_solid_apex 46.39638367,187.44307667,-205.03317767 \
  -atrium_fsi_atol 1e-6 -atrium_fsi_solid_atol 1e-10 \
  -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps \
  -mat_mumps_icntl_20 0 -ksp_diagonal_scale false -ksp_diagonal_scale_fix false
```

This advances one $0.85\,\mathrm{s}$ cycle. The contraction strength is an
explicit example parameter; it does not constitute a physiological calibration.

## Pressure-driven coupling, ventricle, and species

```sh
mpirun -n 6 build/examples/Heart/Atrium/Atrium \
  -atrium_mode fsi -atrium_fsi_newton \
  -atrium_mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  -atrium_output output/atrium-fsi-ptag \
  -atrium_dt 0.001 -atrium_steps 11050 -atrium_output_every 10 \
  -atrium_solid_contractility 2000 -atrium_solid_cache_constitutive \
  -atrium_solid_apex 46.39638367,187.44307667,-205.03317767 \
  -atrium_fluid_vms -atrium_fluid_lv_model -atrium_fluid_pressure_ramp 0.1 \
  -atrium_fluid_kinetics -atrium_fluid_flow_cycles 3 \
  -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps \
  -mat_mumps_icntl_20 0 -ksp_diagonal_scale false -ksp_diagonal_scale_fix false
```

Three flow cycles precede ten species cycles. The species are advanced after
each accepted coupled step with the accepted velocity relative to the mesh,
$u-w$, on the accepted geometry; the fibrin drag acts on the flow of the next
step. The solid reference configuration is unloaded: the opening pressures are
ramped linearly over `-atrium_fluid_pressure_ramp` seconds, and no prestress is
computed.
