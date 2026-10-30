# Lumped ventricular mechanics and circulation

`CCMLC2014T` represents the reduced ventricular mechanics of Caruel, Chabiniok,
Moireau, Lecarpentier, and Chapelle. A spherical chamber couples radial inertia,
passive/viscous stress, active contractile variables, valves, and proximal/distal
arterial compliances. No spatial mesh or finite-element space is involved.

## Chamber mechanics and active state

Let $R_0,d_0$ be reference cavity radius and wall thickness, $y$ radial displacement,
$v=\dot y$, $\ell=1+y/R_0$, and $S$ the reduced total stress. The driver reports

$$
V=\frac{4\pi}{3}(R_0+y)^3,\qquad e=\frac{\ell^2-1}{2},\qquad
\dot y=v,\qquad
\rho_0d_0\dot v+\frac{d_0}{R_0}\ell S-p_v\ell^2=0.
$$

Passive, viscous, and active stresses contribute to $S$. Contractile strain
$e_c$, stiffness $k_c$, stress $\tau_c$, and a load-dependent relaxation variable
complete the material state. In particular, their rates include

$$
\dot k_c+r k_c=\zeta k_0u_+,\qquad
\dot\tau_c=k_c\dot e_c+\zeta\sigma_0u_+-r\tau_c,
$$

where $u_+$ is positive electrical drive and $r,\zeta$ are the model's destruction
rate and recruitment factor. The independent strain equation couples these
variables to chamber strain and contractile viscosity; they are not prescribed
active pressure amplitudes.

## Cavity and circulation balance

With positive $Q_{\mathrm{cav}}$ denoting net flux out of the cavity and cavity
capacity $C_v$, the pressure row is

$$
C_v\dot p_v+Q_{\mathrm{cav}}+4\pi(R_0+y)^2v=0.
$$

Mitral and aortic pressure differences select filling/ejection fluxes. Let
$p_{ar},p_d$ be proximal/distal arterial pressures, $q_p$ proximal branch flow,
$q_d$ distal return flow, and $Q_{ao}$ aortic outflow. The two compliance rows are

$$
C_p\dot p_{ar}=Q_{ao}-q_p,\qquad C_d\dot p_d=q_p+q_d,
\qquad Q_{ao}=\max(K_{ar}(p_v-p_{ar}),0).
$$

The model evaluates $q_d$ in the direction from systemic venous pressure toward
$p_d$, so it is normally negative during drainage. Branch-flow laws and valve
conductances are specified by `Model::Input`; these states form the reduced
four-element circulation rather than a resolved vascular network.

## Numerical procedure and output

The dynamic state uses the configured implicit scheme (BDF2 by default with
backward-Euler startup) and an exact dense residual Jacobian. The driver uses
$\Delta t=10^{-3}\,\mathrm{s}$ over three cycles of period
$0.85\,\mathrm{s}$. Electrical drive and atrial pressure are prescribed periodic
functions in the driver. `initialize()` sets the mechanical/active state and
`step()` solves the next nonlinear system.

`ccmlc2014_0d_cycle.csv` records time, radial motion, chamber/arterial pressures,
active variables, volume, and flow. Length, stress, pressure, and time are in SI
units. The spherical reduction omits spatial fiber architecture, regional
wall deformation, valve geometry, and 3D hemodynamics.

No PETSc or MPI dependency is required.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CCMLC2014_0D -j4
cd build
./examples/Heart/Ventricle/CCMLC2014_0D/CCMLC2014_0D
```

Output files are written relative to the working directory.

## Input and configuration

No input mesh is required.
