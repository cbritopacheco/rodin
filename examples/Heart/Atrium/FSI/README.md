# Atrium FSI

Fluid inertia loads the contracting wall, so independently advancing the two
components does not enforce their interface conditions. `FSI` uses a partitioned
Robin–Robin construction with the fluid impedance and interface convective
correction of the [coronary P1/P1 example](../../CoronaryArtery/FSI/Explicit_P1P1_PETSc_Seq/README.md).
Both subproblems receive Robin data on material label 4.

Let $n_f$ point out of the current fluid domain, $N_s$ out of the reference solid,
and $J_\Gamma=dS/dS_0$. The continuous target is displacement compatibility,
velocity continuity, and traction equilibrium:

$$
d_f=d_s,\qquad u_f=\dot d_s,\qquad
P_sN_s=J_\Gamma(-\sigma_fn_f),\qquad
\sigma_f=-pI+2\mu D(u_f).
$$

The time-discrete interface rate is

$$
v_\Gamma^{n+1}=\frac{d_s^{n+1}-d_s^n}{\Delta t}.
$$

It matches the BDF1 ALE boundary velocity. For the selected average-acceleration
Newmark parameters, it also equals $(v_s^{n+1}+v_s^n)/2$. The solid retains its
Newmark endpoint velocity for bulk inertia and tissue viscosity; that endpoint
velocity is distinct from the interface rate at finite time step.

### Robin data and reference-area traction

Let $T^k$ denote the fluid reaction on the solid per reference area and let
$\alpha_0>0$ be a reference-area Robin impedance. A correction sweep solves

$$
P_s^{k+1}N_s+\alpha_0v_\Gamma^{k+1}=g^k,\qquad
g^k=T^k+\alpha_0u_f^k.
$$

The fluid receives the following boundary contribution to its weak residual:

$$
\int_{\Gamma(t)}\frac{\alpha_0(u_f^{k+1}-v_\Gamma^{k+1})+T^k}{J_\Gamma}
\cdot v\,dS
-\frac{\rho_f}{2}\int_{\Gamma(t)}[(u_f^n-a)\cdot n_f]
(u_f^{k+1}\cdot v)\,dS.
$$

Here $a$ is the ALE mesh velocity. Division by $J_\Gamma$ converts reference-area
Robin data to the moving fluid surface. The second integral is the weakly
consistent interface convective correction. Fluid velocity is not constrained
by an interface Dirichlet condition.

After each fluid solve, the transmission update is

$$
T^{k+1}=T^k+\alpha_0(u_f^{k+1}-v_\Gamma^{k+1}),\qquad
g^{k+1}=T^{k+1}+\alpha_0u_f^{k+1}.
$$

$T^{k+1}$ is the reference fluid variational reaction, including the
interface convective correction; matching P1 traces permit the update in
their nodal space. It follows the discrete Robin stress update of
[Burman et al. (2025), equations (13), (16)–(19)](https://discovery.ucl.ac.uk/id/eprint/10210489/).
There is no relaxation. The sweep structure, stopping test and budget follow
the coronary example, but its recovered-gradient traction is not used: with
that traction the fixed point satisfies equilibrium with the variational
reaction $\lambda_h$ while leaving the slip
$u_f-v_\Gamma=(\lambda_h-T_{\mathrm{rec}})/\alpha_0$, which does not vanish
as the sweeps converge. With the algebraic update the fixed point has
$u_f=v_\Gamma$. Both transmission velocities retain BDF1; the solid Robin
tangent uses $\alpha_0/\Delta t$.

### Displacement stopping test and sweep budget

After a complete solid/ALE/fluid sweep, define

$$
r_d=\frac{\|d_s^{k+1}-d_s^k\|_2}{\|d_s^{k+1}\|_2}.
$$

The norm is the global Euclidean norm of the complete solid displacement
coefficient vector, not a surface-weighted interface norm. If the denominator
is zero, the absolute increment norm in metres is used instead. Starting with
the second sweep, $r_d<\varepsilon_d$ terminates the correction loop. Defaults
are `-atrium_fsi_rtol 1e-6` and `-atrium_fsi_iterations 5`. A zero
tolerance disables early termination and uses the complete sweep budget.

Budget exhaustion accepts the last completed sweep. `fsi_converged` and
`fsi_budget_exhausted` in the common CSV distinguish criterion satisfaction
from acceptance at the budget. Inner solver failures, invalid constitutive
states, and mesh inversion still reject the step. A small budget does not
establish interface continuity; velocity slip, traction-data change, and
geometry mismatch remain diagnostics. Unlike the reference's early exit
before its fluid solve, each accepted Atrium sweep includes ALE and fluid
updates, so both exported geometries correspond to the accepted solid state.

No wall-face ring band is introduced. Existing solid opening cut faces remain
fixed, ALE caps remain fixed, and fluid opening interiors remain open. The
interface is not given an additional strong velocity constraint.

The automatic impedance uses the coronary scaling with a matrix-based isotropic
small-strain modulus:

$$
\alpha_0=\sqrt{\rho_{s,0}E_m},\qquad
E_m=\frac{9K_m\mu_m}{3K_m+\mu_m},\qquad
\mu_m=a,\quad K_m=2\mu_v.
$$

Material coefficients are converted to Pa for this calculation. This scale
omits directional fiber reinforcement and is an algorithmic choice, not a full
anisotropic Young modulus. `-atrium_fsi_robin` supplies an explicit positive
value in $\mathrm{Pa\,s/m}$.

### ALE evolution and problem ownership

A harmonic extension of solid displacement moves the fluid interior; opening
caps remain fixed. With $a=(d_f^{n+1}-d_f^n)/\Delta t$ and cell-volume ratio
$r_K=|K^n|/|K^{n+1}|$, the conservative first-order ALE residual is

$$
R_m^{\mathrm{ALE}}=
\rho_f\frac{u_f^{n+1}-r_Ku_f^n}{\Delta t}
+\rho_f[(u_f^n-a)\cdot\nabla]u_f^{n+1}
+\rho_f\left(\frac12\nabla\cdot u_f^n-\nabla\cdot a\right)u_f^{n+1}
-\nabla\cdot(2\mu_nD(u_f^{n+1}))+\nabla p^{n+1}.
$$

The default stabilization is the orthogonal-subscale VMS of the coronary
example, evaluated with the relative advecting velocity $u_f^n-a$: the
bilinear term $\int\tau_K\rho^2[(\nabla u)(u^n-a)]\cdot[(\nabla v)(u^n-a)]$, the
subtraction of $\Pi_h[(\nabla(u^n-a))(u^n-a)]$ and the dynamic subscale,
grad-div with $\tau_C=\rho h_K^2/(4\tau_1)$ minus $\Pi_h[\sqrt{\tau_C}\nabla\cdot u^n]$,
and PSPG $\int(\tau_1/\rho)\nabla p\cdot\nabla q$; $h_K$ is the longest edge of
the current cell. The projections are recomputed after every mesh motion.
`-atrium_fluid_vms false` selects SUPG and PSPG on the residual above. Previous
fluid velocity and solid internal states remain fixed through coupling trials.
Mesh motion invalidates Rodin's geometric caches through its coordinate API.
The same spaces, Problems, solvers, and MPI trace scatter plan are retained.
Recovered gradients and wall-shear projections are computed once per accepted
physical step and remain output/statistics operations.

The common CSV records the Robin coefficient, traction-data update,
displacement-change metric and tolerance, criterion and budget status,
fluid/BDF1 interface velocity mismatch, fluid/Newmark endpoint mismatch, and
shared-geometry mismatch separately. Fluid `MeshVelocity` gives the BDF1 target
at the wall; solid `Velocity` gives the Newmark endpoint velocity. The two fields
should not be identified at finite time step.

All openings have zero imposed pressure and no prescribed flux. Fluid starts
at rest and active wall contraction drives the flow. Backflow stabilization
remains active. FSI excludes the 0D ventricular driver and species transport. The outer solid surface is unloaded and opening rims are
clamped. No valve or circulatory model is represented.


## Architecture and run

`FSI` binds the existing `Solid` and `Fluid` objects and retains its trace
scatter plans, interface data, harmonic mesh-motion Problem, and solver. It
owns the coupling iteration; each component owns its assembly and physical
state. The shared `Parameters::timeStep` is used for both components. Previous
states remain frozen through coupling trials, and output and time histories
are committed after the displacement criterion is satisfied or the sweep
budget is exhausted.

```sh
cmake --build build --target Atrium -j4
mpirun -n 4 build/examples/Heart/Atrium/Atrium \
  -atrium_mode fsi \
  -atrium_mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  -atrium_output output/atrium-fsi \
  -atrium_dt 0.001 -atrium_steps 24 -atrium_output_every 1 \
  -atrium_solid_cache_constitutive \
  -atrium_fluid_lag_preconditioner 5 \
  -ksp_type gmres -ksp_rtol 1e-12 -ksp_max_it 500 \
  -ksp_norm_type unpreconditioned -pc_type lu \
  -pc_factor_mat_solver_type mumps -mat_mumps_icntl_20 0 \
  -ksp_diagonal_scale false -ksp_diagonal_scale_fix false
```

`Solid.xdmf` stores reference coordinates and `Displacement`; apply **Warp By
Vector** to visualize deformation. `Fluid.xdmf` stores moving coordinates and
must be viewed directly. `Diagnostics.csv` combines both component diagnostics
and coupling status. The run saves 25 states including the initial one.
The [shared configuration](../README.md#configuration) lists the component
option prefixes; `-help` includes displacement-change tolerance and sweep budget.
