# Coronary fluid–structure interaction

A deforming arterial wall changes the fluid domain, while fluid stress loads the
wall. Let $\Omega_{s,0}$ be the reference solid, $\Omega_f(t)$ the current lumen,
and $\Gamma(t)$ their interface. Solid displacement $d_s$ defines
$F_s=I+\nabla_Xd_s$ and $J_s=\det F_s$. An ALE displacement $d_f$ extends the wall
motion into the fluid, with mesh velocity $a$.

## Component models

The fluid stress and ALE balance are

$$
\sigma_f=-pI+2\mu(\dot\gamma)D(u_f),\qquad
\rho_f\left(\left.\partial_tu_f\right|_{\mathrm{ALE}}
+[(u_f-a)\cdot\nabla]u_f\right)-\nabla\cdot\sigma_f=0,\qquad
\nabla\cdot u_f=0.
$$

Here $D(u)=(\nabla u+\nabla u^T)/2$ and
$\dot\gamma=(\dot\gamma_{\mathrm{reg}}^2+2D:D)^{1/2}$. The Carreau–Yasuda law is

$$
\mu(\dot\gamma)=\mu_\infty+(\mu_0-\mu_\infty)
[1+(\lambda_{\mathrm{CY}}\dot\gamma)^m]^{(n-1)/m}.
$$

The reference solid satisfies $\rho_{s,0}\ddot d_s-\mathrm{Div}_XP_s=0$.
The displacement-only Neo-Hookean variants use the stored energy and elastic stress

$$
W=\frac{\mu_s}{2}(\mathrm{tr}(F_s^TF_s)-3)-\mu_s\log J_s
 +\frac{\lambda_s}{2}(\log J_s)^2,\qquad
P_s^{\mathrm{el}}=\mu_sF_s+(\lambda_s\log J_s-\mu_s)F_s^{-T}.
$$

The sequential staggered variants instead use compressible Yeoh elasticity,
with $\bar I_1=J_s^{-2/3}\mathrm{tr}(F_s^TF_s)$:

$$
W_Y=c_1(\bar I_1-3)+c_2(\bar I_1-3)^2+c_3(\bar I_1-3)^3
+\frac{\kappa}{2}(J_s-1)^2,\qquad P_s^{\mathrm{el}}=\frac{\partial W_Y}{\partial F_s}.
$$

The P1/P1 driver multiplies this energy by a smooth transmural factor $m(\xi)$;
the sequential P2/P1 driver uses a homogeneous Yeoh law. The corresponding
stress and tangent use the same multiplier, held fixed in a material-point
linearization.

Some variants add a reference-gradient damping form
$\eta_s\int_{\Omega_{s,0}}\nabla_Xv_s:\nabla_Xw\,dX$.
This specific term differs from the Green–Lagrange strain-rate viscosity in
the [atrium](../../Atrium/Solid/README.md). Pressure prestress, ring clamps,
outlet supports, and heart-motion data are specified per variant.

## Transmission conditions

With $n_f$ outward from the fluid and $N_s$ outward from the reference solid,
the physical interface conditions are

$$
d_f=d_s,\qquad u_f=v_s,\qquad
P_sN_s=J_\Gamma t_f,\qquad t_f=-\sigma_fn_f,
$$

where $J_\Gamma=dS/dS_0$ converts a spatial traction into a reference-area load.
Opposite fluid and solid normals determine the minus sign in $t_f$.

The staggered examples impose Robin data on **both** component problems. In a
correction sweep, the solid and fluid conditions have the forms

$$
P_s^{k+1}N_s+\alpha_s v_s^{k+1}
 =J_\Gamma t_f^k+\alpha_su_f^k,\qquad
\sigma_f^{k+1}n_f+\alpha_fu_f^{k+1}
 =-t_f^k+\alpha_f\widetilde v_s^{k+1}.
$$

The solid Robin coefficient $\alpha_s$ is expressed per reference area; the
fluid coefficient $\alpha_f$ is expressed per current area. The P2/P1 staggered
forms use $\alpha_s=J_\Gamma\alpha_f$, while the P1/P1 form uses its Robin
coefficient directly on the reference solid surface. These are different
normalizations at finite deformation. The automatic fluid scale is
$\alpha_f=\gamma_R\sqrt{\rho_sE_{\mathrm{eq}}}$; an explicit coefficient can be
selected instead. Here $E_{\mathrm{eq}}$ is a representative small-strain Young
modulus. For Yeoh, $\mu_{\mathrm{eq}}=2c_1$ and
$E_{\mathrm{eq}}=9\kappa\mu_{\mathrm{eq}}/(3\kappa+\mu_{\mathrm{eq}})$;
Neo-Hookean variants use their configured Young modulus. The Robin coefficient
is an algorithmic impedance, not an extra elastic law.

These expressions specify the traction/velocity part of transmission. The fluid
weak form also contains the interface convective correction

$$
-\frac{\rho_f}{2}\int_{\Gamma(t)}
[(u_f^n-a)\cdot n_f]\,u_f^{n+1}\cdot v\,dS.
$$

A single sweep is a loose partitioned step and generally leaves interface slip
and a traction imbalance. Further correction sweeps reduce the coupling error;
a configured iteration budget alone does not establish convergence. The ideal
conditions are recovered only with compatible interface kinematics and converged
transmission data. In particular, the sequential drivers can target the BDF1
wall rate $\widetilde v_s=(d_s^{n+1}-d_s^n)/\Delta t$ instead of the Newmark
endpoint velocity. These two velocities differ at finite time step.

The [atrium FSI](../../Atrium/FSI/README.md) also uses
Robin–Robin transmission and the BDF1 interface rate. Its two Robin forms use
one reference-area impedance, and its traction is updated from the discrete
Robin equation rather than reconstructed from output stress gradients. Its
bulk solid material and time integration are specified by the atrium model.

## Problem ownership

Each staggered driver has a `Config` containing mesh, material, time, and
circulation data. Its `InterfaceMap` connects conforming fluid and solid
boundary entities; fields and `Problem` objects represent the solid solve,
fluid solve, ALE extension, and projections. The monolithic driver instead
assembles corrections on one common mesh. These example-local constructions
are distinct from the atrium's `Atrium`, `Solid`, `Fluid`, and `FSI` components.

## Discretization and variants

The staggered fluid problems use conservative two-configuration BDF1 ALE
integration with lagged Oseen transport, a Temam correction, and mesh-divergence
terms. The solid uses Newmark integration. The sequential drivers solve a weighted
harmonic ALE problem,

$$
\int_{\Omega_{f,0}}\omega_K\nabla_Xd_f:\nabla_X\varphi\,dX=0,\qquad
\omega_K=(h_{\mathrm{ref}}/h_K)^p,
$$

with default powers $p=0.75$ for P1/P1 and $p=0.5$ for P2/P1. The MPI staggered
driver uses unit weight. The P2/P1 drivers fix all opening caps; the P1/P1 driver
fixes the inlet and ring band and leaves the ALE outlet caps with their natural
condition. The shared wall follows transferred solid displacement. Mesh
Jacobians must remain positive.

| Example | Fluid / solid spaces | Wall and coupling model |
|---|---|---|
| [Explicit_P1P1_PETSc_Seq](Explicit_P1P1_PETSc_Seq/README.md) | P1/P1 velocity–pressure; P1 wall. | Graded Yeoh; Robin–Robin, equal-order pressure stabilization, intramyocardial R–$\mu$ beds. |
| [Explicit_PETSc_Seq](Explicit_PETSc_Seq/README.md) | P2/P1 velocity–pressure; P1 wall. | Homogeneous Yeoh; Robin–Robin and non-Newtonian intramyocardial Windkessels. |
| [Explicit_PETSc_MPI](Explicit_PETSc_MPI/README.md) | P2/P1 velocity–pressure; P1 wall. | Neo-Hookean; Robin–Robin with Aitken corrections; interface mapping requires one rank. |
| [PETSc_MPI](PETSc_MPI/README.md) | P2/P1 velocity–pressure; P1 displacement on a common mesh. | Neo-Hookean; monolithic ALE/SNES correction system and lumped outlet states. |

Conforming interface geometry is required. Mesh transfer, projection, stabilization,
and a finite correction count are numerical approximations; the continuous
interface equations alone do not establish discrete conservation or stability.
