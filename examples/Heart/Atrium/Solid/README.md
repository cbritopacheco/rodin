# Atrial solid electromechanics

This example couples finite-strain atrial wall motion to independent active
fiber variables. An exponential Demiray matrix with single-family Holzapfel–Gasser–Ogden-type reinforcement and Green–Lagrange strain-rate
viscosity form the passive tissue response. Rodin's `ActiveFiberLaw` supplies
the Hill–Maxwell equations. Displacement uses continuous linear wedge shape
functions; active extension, stiffness, and stress variables use discontinuous
cellwise constants. The coupled Newton problem retains its spaces, forms,
fields, and solvers across time steps.

## Passive material and dissipation

For displacement $u$, reference fiber $f_0$, and reference velocity $v$,

$$
F=I+\nabla_Xu,\quad C=F^TF,\quad J=\det F,\quad
E=\frac{C-I}{2},\quad
\bar I_1=J^{-2/3}\mathrm{tr}C,\quad
\bar I_4=J^{-2/3}f_0\cdot Cf_0.
$$

`DemirayHolzapfelGasserOgden` combines a Demiray matrix and a single-family
HGO-type reinforcement. Its three-dimensional stored energy is

$$
W_p=\frac{a}{2b}\left[e^{b(\bar I_1-3)}-1\right]
 +\frac{a_f}{2b_f}\left[e^{b_f(\bar I_4-1)^2}-1\right]
 +\mu_v(\log J)^2.
$$

The exponential terms use their continuous limits when an exponent is zero.
The volumetric coefficient gives initial bulk modulus $2\mu_v$. The fiber
term acts in extension and compression. Reference energy and stress vanish
at $F=I$.

The viscous dissipation potential and its associated stresses are

$$
\Phi=\frac{\nu}{2}\dot E:\dot E,\qquad
\dot E=\frac{F^T\nabla_Xv+(\nabla_Xv)^TF}{2},\qquad
\Sigma_v=\nu\dot E,\qquad P_v=\nu F\dot E.
$$

Thus $P_v:\nabla_Xv=\nu\dot E:\dot E\geq0$ for $\nu\geq0$.
This is a finite-strain constitutive viscosity. Stored energy excludes the
dissipation potential. The tangent includes changes of both $F$ and the
velocity gradient through the time integrator.

## Active fiber model

The fiber Green–Lagrange strain is $e_f=f_0\cdot Ef_0$. The active second
Piola stress and serial constitutive relation are

$$
\Sigma_a=\sigma_{1D}f_0\otimes f_0,\qquad
\sigma_{1D}=E_s\frac{e_f-e_c}{(1+2e_c)^2},\qquad
\tau_c+\mu_c\dot e_c=
E_s\frac{(e_f-e_c)(1+2e_f)}{(1+2e_c)^3}.
$$

The internal contractile stiffness and stress obey

$$
\dot k_c=n_0k_0u_e^+-(|u_e|+\alpha|\dot e_c|)k_c,
$$

$$
\dot\tau_c=n_0\sigma_0u_e^+-(|u_e|+\alpha|\dot e_c|)\tau_c
 +\dot e_c k_c.
$$

Here $u_e^+=\max(u_e,0)$, $n_0$ is Rodin's length-dependent Starling factor,
$E_s$ is serial stiffness, $\mu_c$ is the contractile dashpot coefficient,
and $\sigma_0$ controls the active stress source. The numerical unknowns are

$$
e_c,\qquad \gamma=\sqrt{k_c},\qquad \beta=\frac{\tau_c}{\gamma}.
$$

All three are independent discontinuous finite-element fields participating
in the global Newton solve with displacement. Their residuals and all
cross derivatives come from `ActiveFiberLaw`; no active field is eliminated
at a quadrature point. The library's small positive stiffness floor regularizes
the stress transformation at zero contractile stiffness. Admissible iterates
have $\gamma>0$, $1+2e_c>0$, and positive deformed wedge Jacobians.

## Mechanical problem and initialization

The first Piola stress and momentum residual are

$$
P=\frac{\partial W_p}{\partial F}+\nu F\dot E+F\Sigma_a,
\qquad
R(u;w)=\int_{\Omega_0}\rho_0\ddot u\cdot w+P:\nabla_Xw\,dX.
$$

Rim faces with labels 1, 2, 3, 5, and 6 have fixed displacement. Inner and
outer wall faces, labeled 4 and 7, are traction-free in the solid-only problem.
Chamber pressure and fluid interface forces belong to the fluid/coupling model.

Initialization solves the passive static problem with no active stress and
zero rate. With these unloaded wall boundaries the reference configuration is
an equilibrium. The active extension is initialized to the passive fiber strain,
$e_c=e_f$, and the contractile stress is zero. Initial velocity and acceleration
are zero. Electrical drive subsequently evolves the internal variables; it
does not impose an instantaneous target contraction at the initial state.

## Time integration and nonlinear solve

Newmark average acceleration uses $\beta_N=1/4$ and $\gamma_N=1/2$:

$$
u^*=u_n+\Delta t\,v_n+\Delta t^2(1/2-\beta_N)\ddot u_n,
\qquad v^*=v_n+\Delta t(1-\gamma_N)\ddot u_n,
$$

$$
\ddot u_{n+1}=\frac{u_{n+1}-u^*}{\beta_N\Delta t^2},\qquad
v_{n+1}=v^*+\gamma_N\Delta t\,\ddot u_{n+1}.
$$

The mechanical residual uses the current deformation and Newmark velocity.
The active law uses its compatible history-dependent update: backward Euler
for stiffness and transformed stress, midpoint active extension and the mean
of current/previous fiber strains in the serial relation. The first-order
active updates limit the temporal order of the coupled method.

Rodin's `InternalVirtualWork`, `Integral`, PETSc `KSP`, and `SNES`
assemble and solve the four-field system. Every nonlinear and coupling term
uses the same order-four quadrature policy. A collective state-domain check
allows SNES line search to reject invalid active stretches and nonpositive wedges
before constitutive assembly. An electrical
upper bound over each wedge supplies only an initial Newton predictor; the
active equations always evaluate their source at quadrature points.

The pulse traversal time is $T=(\max_X s(X)+3w)/v$. By default, `--steps N`
sets $\Delta t=T/N$ and writes $N+1$ states including initialization.
`--dt` instead sets the physical time increment in seconds. `--stop-after`
limits completed intervals without changing the chosen time step. Temporal,
spatial, and quadrature refinement are needed for quantitative dynamics.

## Parameters and units

Coordinates are in mm, time in seconds, and stresses in kPa. `--density`
accepts reference density in $\mathrm{kg/m^3}$ and converts it using

$$
\rho_{\mathrm{internal}}=10^{-9}\rho_{\mathrm{SI}},\qquad
1060\,\mathrm{kg/m^3}=1.06\times10^{-6}\,\mathrm{kPa\,s^2/mm^2}.
$$

The defaults below are demonstration coefficients, not a calibrated atrial dataset.

| Parameter | Default | Option |
|---|---:|---|
| Matrix scale $a$ | $20\,\mathrm{kPa}$ | `--matrix-scale` |
| Matrix exponent $b$ | $5$ | `--matrix-exponent` |
| Fiber scale $a_f$ | $20\,\mathrm{kPa}$ | `--fiber-scale` |
| Fiber exponent $b_f$ | $10$ | `--fiber-exponent` |
| Volume penalty $\mu_v$ | $500\,\mathrm{kPa}$ | `--volume-penalty` |
| Tissue viscosity $\nu$ | $0.2\,\mathrm{kPa\,s}$ | `--viscosity` |
| Serial stiffness $E_s$ | $20\,\mathrm{kPa}$ | `--active-stiffness` |
| Contractile viscosity $\mu_c$ | $3\,\mathrm{kPa\,s}$ | `--active-damping` |
| Cross-bridge stiffness source $k_0$ | $100\,\mathrm{kPa}$ | `--cross-bridge-stiffness` |
| Contractility source $\sigma_0$ | $20\,\mathrm{kPa}$ | `--contractility` |
| Destruction factor $\alpha$ | $0.4$ | `--destruction-rate` |
| Positive electrical peak $A_+$ | $20\,\mathrm{s^{-1}}$ | `--activation` |
| Resting drive magnitude $A_-$ | $10\,\mathrm{s^{-1}}$ | `--relaxation` |
| Reference density | $1060\,\mathrm{kg/m^3}$ | `--density` |

`--contractility` controls stress generation. Active extension is a solved
field and has no prescribed target-strain amplitude option. Setting
`--activation 0` supplies no positive electrical drive for a resting-state check.

## Architecture and geometry

- [`Solid.h`](Solid.h) and [`Solid.cpp`](Solid.cpp): directly owned mechanical
  and active fields, constitutive inputs, retained coupled problem, nonlinear
  solve, passive initialization, and Newmark updates.
- [`AtriumContraction.cpp`](AtriumContraction.cpp): command-line entry point.
- [`ProblemData.h`](ProblemData.h) and [`ProblemData.cpp`](ProblemData.cpp):
  mesh checks, reference fibers, basal plane, and electrical source.
- [`Output.h`](Output.h) and [`Output.cpp`](Output.cpp): observation of solved
  fields, XDMF export, metadata, and mechanical diagnostics.

The reusable `Atrium::Solid` component exposes `initialize()`, `step()`, and
`finish()`. The coordinator partitions the combined fluid–solid parent once with
Scotch, then the component extracts its material-4 MPI submesh. Its mesh,
problems, fields, and solvers are retained across steps. The mesh resources are in `resources/examples/Heart/Atrium/`:

- `AtriumCoarse_FluidSolid_1mm.mesh`: conforming fluid tetrahedra and solid wedges.
- `Coarse.mesh`: source atrial volume and labeled boundary.
- `AtriumWall_Repaired_Surface.mesh`: prepared atrial wall surface.
- `AtriumWall_Repaired_Wedges_1mm.mesh`: one outward layer of six-node wedges,
  with a nodal extrusion distance of $1\,\mathrm{mm}$.

The wall construction uses the material-4 surface. Inlet and outlet regions
are open, with their labels transferred to rim quadrilaterals. The surface
preparation includes a local handle-to-disk repair and MMGS adaptation, so its
topology is a modeling assumption requiring anatomical assessment.

## Plane electrical wave from the anatomical base

Label 6 is the basal cap in the source volume. In the wedge wall, its label
is preserved on the adjacent rim. Inner rim vertices define a least-squares
plane through their mean $b$. Its normal $d$ points toward the wall centroid.

$$
s(X)=d\cdot(X-b),\qquad
 a(X,t)=\exp\left[-\left(\frac{s(X)-vt}{w}\right)^2\right],\qquad
 u_e(X,t)=(A_++A_-)a(X,t)-A_-.
$$

The crest starts at the basal plane and travels along positive $d$. The signed
electrical input permits cross-bridge formation near the wave crest and
relaxation behind it. Points with $s(X)<0$ receive only the initial tail of this
one-direction pulse. A wall-distance wave would be a different prescribed input.
The fitted plane and coordinate range are recorded in geometry metadata;
`WaveCoordinate` is a static nodal inspection field.

## Simple landmark-based fibers

The default left-appendage apex landmark is

$$
X_{\mathrm{LAA}}=(46.39638367,187.44307667,-205.03317767)\ \mathrm{mm}.
$$

Coordinates are used instead of vertex or triangle numbers because MMG
renumbers entities. An apex patch is selected using inner-triangle centroids
within 3 mm of this point, together with the corresponding outer triangles.
A finite patch avoids imposing a point-only Dirichlet condition in the wall
volume. Its radius and center can be changed from the command line.

A single harmonic coordinate is solved with Rodin:

$$
-\Delta\psi=0\quad\text{in }\Omega_{\mathrm{wall}},\qquad
\psi=0\quad\text{on rim 6},\qquad
\psi=1\quad\text{on the appendage patch},
$$

with zero normal flux on the remaining boundary. A temporary boundary label
8 is used for this solve and then restored, preserving all original wall and
rim labels.

For each wedge, the known extrusion vectors supply a unit transmural
normal $m$. The harmonic gradient is projected into the wall tangent
plane and the fiber is constructed as

$$
\ell=\frac{(I-mm^T)\nabla\psi}{\|(I-mm^T)\nabla\psi\|},
\qquad f=\ell\times m.
$$

Thus the fibers are tangent to the wall and follow contours of the
base-to-appendage coordinate. They are stored as one constant unit direction
per wedge. There is no transmural rotation in this one-layer model.

This is a simplified use of the Laplace-gradient frame in Piersanti et al.,
*Modeling cardiac muscle fibers in ventricular and atrial electrophysiology
simulations*, [arXiv:2101.10960v1](https://arxiv.org/abs/2101.10960v1), equations
(3) and (5). It does not implement their regional atrial bundle selection,
pulmonary-vein rules, or a validated anatomical fiber reconstruction.

Critical points are handled explicitly. If the harmonic tangent vanishes,
the projected direction toward the apex is used; at a remaining pole, a
Cartesian tangent is chosen. `FiberFallback` flags these cells: 0 means the
harmonic construction, 1 the geometric-apex fallback, and 2 the pole fallback.
`FiberCoordinate`, `LongitudinalDirection`, `WallNormal`, and
`FiberDirection` are exported for inspection.

## Build and run

From the repository root with a configured build tree:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DRODIN_USE_MPI=ON -DRODIN_USE_PETSC=ON -DRODIN_USE_SCOTCH=ON \
  -DRODIN_CODE_COVERAGE=OFF
cmake --build build --target AtriumContraction -j4
mpirun -n 4 build/examples/Heart/Atrium/Solid/AtriumContraction \
  --mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  --output output/atrium-active/Contraction \
  --base-label 6 --apex 46.39638367 187.44307667 -205.03317767 \
  --apex-radius 3 --steps 24 --activation 20 --relaxation 10 \
  --contractility 20 --viscosity 0.2 --active-damping 3
```

`--fibers-only` exports reference geometry and fiber diagnostics without advancing
mechanics. `--help` lists all model and source parameters. Output geometry checks
use the exact minimum of each six-node wedge's Jacobian determinant over its
reference cell, including possible thickness-interior minima.

## ParaView

Open `<stem>.xdmf` with its companion HDF5 files in the same directory. Apply
**Warp By Vector**, select `Displacement`, and use scale factor 1 for physical
motion in mm. `Velocity` and `Acceleration` use mm/s and mm/s². `Activation`
is the dimensionless Gaussian shape; `ElectricalActivation` is the signed
source in inverse seconds.

`ActiveExtension`, `ActiveGamma`, and `ActiveBeta` are solved cell fields.
`ContractileStiffness` stores $k_c=\gamma^2$ and `ContractileStress` stores
$\tau_c=\gamma\beta$; the latter is an internal contractile stress, distinct
from the serial fiber stress $\sigma_{1D}$. For fibers, apply **Cell Centers**
and line glyphs oriented by `FiberDirection`. `FiberCoordinate`, `WallNormal`,
`LongitudinalDirection`, and `FiberFallback` describe the reference frame.

The CSV records time, displacement/velocity/acceleration maxima, activation,
active-extension extrema, contractile stress, Newton residuals, and minimum
scaled Jacobian. Metadata records geometry, material coefficients, electrical
source parameters, density, and time-integration coefficients.
