# Atrium FSI

Active wall contraction drives an initially resting fluid through open boundaries.
`FSI` couples the Demiray–Holzapfel–Gasser–Ogden solid with Kelvin–Voigt
dissipation to the ALE Oseen fluid. The fluid uses P1/P1 velocity and pressure
with SUPG, PSPG and grad-div stabilization. No ventricular driver, imposed flow,
coagulation transport or orthogonal subscales enter this coupled model.

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
additional interface ring is clamped. All fluid openings have zero imposed
pressure and no prescribed flux.

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

Material coefficients are converted to Pa. This matrix-based scale omits
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
`-atrium_fsi_interface_ksp_monitor`. Newton failure rejects the step.
`-atrium_fsi_check_tangent` enables residual Taylor checks for verification.

## Budgeted Robin sweeps

Without `-atrium_fsi_newton`, the legacy Robin sweeps recover the fluid velocity
gradient and project physical stress to the reference interface with lumped
trace mass. This exchange differs from Newton's variational reaction map.
Starting at the second completed sweep, the relative
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
