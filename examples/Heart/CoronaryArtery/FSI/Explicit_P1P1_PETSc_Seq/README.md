# Coronary FSI with stabilized linear velocity and pressure

This sequential example couples coronary flow to a deformable arterial wall
using the [Robin–Robin splitting](../README.md#transmission-conditions).
The fluid uses continuous P1 velocity and P1 pressure; wall and ALE displacement
also use P1 fields. Separate conforming fluid and solid meshes are related by
an interface map.

## Model and transmission

The [component balances](../README.md#component-models) combine Carreau–Yasuda
flow, compressible graded Yeoh wall elasticity, and inertia.
A harmonic coordinate $\xi$ is zero on the lumen and one on the outer wall.
`GradedYeoh` defines the stored energy by

$$
W(X,F)=m(\xi(X))W_Y(F).
$$

The factors for intima, media, and adventitia are $0.5$, $1$, and $1.5$,
blended smoothly near $\xi=1/3$ and $2/3$. Base constants are
$c_1=8\times10^4$, $c_2=4\times10^5$, $c_3=5\times10^6$, and
$\kappa=1.2\times10^7$, all in Pa. The scale is sampled at constitutive
quadrature points. A single wedge layer spans all three bands, so this grading
does not by itself provide through-thickness spatial resolution. The wall adds the
reference-gradient viscous form

$$
\eta_s\int_{\Omega_{s,0}}\nabla_Xv_s:\nabla_Xw\,dX.
$$

Its default coefficient is $\eta_s=4\times10^3\,\mathrm{Pa\,s}$.
Static follower-pressure increments initialize the wall at a configured fraction
of the initial arterial pressure (default $0.975$); the remaining pressure is
ramped over ten dynamic steps by default. The dynamic load contains
that pressure level and a transferred traction remainder. The inlet and ring
band are clamped; outlet springs/dashpots and a normal 0D heart-motion penalty provide
additional support data. These are boundary loads, not active coronary fibers.

Both the fluid and solid receive Robin conditions. Spatial traction is converted
by $J_\Gamma$ before entering the reference solid problem; this variant applies
the Robin velocity terms directly per reference area. The fluid target velocity
is the BDF1 mesh-boundary rate by default, selected by
`-coronary_mesh_consistent_interface_vel`. The Newmark wall endpoint velocity
is a distinct finite-time-step quantity. The Robin fluid operator includes the
interface convective correction described in the parent model.

## Stabilization and time integration

Equal-order velocity and pressure require pressure stabilization. The weak
problem adds pointwise pressure-gradient stabilization together with convective
orthogonal-subscale VMS and grad-div terms. The subscale projections are $L^2$
projections; local stabilization parameters depend on lagged viscosity and speed.
This pressure-gradient closure is distinct from assuming an inf-sup-stable pair.

The fluid uses conservative BDF1 ALE/Oseen integration. Weighted harmonic
mesh extension fixes the inlet and ring band; the outlet mesh caps have a
natural condition. The wall uses Newmark
with $\beta=1/4$ and $\gamma=1/2$. The default is one coupling sweep; increasing
`-coronary_coupling_iterations` permits corrections but does not by itself certify
the interface tolerance.

The 0D heart supplies inlet pressure and heart-motion data. Outlet R–$\mu$ beds
advance transmural compartment pressure and use apparent-viscosity factors, as
in the [rigid-wall circulation model](../../CFD/README.md#flow-and-circulation-model).

## Geometry, units, and output

Fluid wall label 1 is paired with the solid interface. Fluid inlet 27 and outlets
24, 25, 26, 28, 30, 31 correspond to solid inlet 153 and outlets
150, 151, 152, 154, 155, 156. Label 99 identifies the clamped interface band;
110–120 identify heart-contact patches; the assembled motion penalty acts on
label 110. Source coordinates are converted to metres
with the default factor $10^{-3}$. Stress, time, and velocity use SI units.

XDMF records fluid and solid fields; CSV records pressures, flows, wall motion,
interface slip/power, and coupling observables. The resolved fields comprise
flow, wall motion, harmonic coordinates, and lumped circulation states.

Requires `RODIN_USE_PETSC=ON`; MPI is not required.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq -j4
cd build
./examples/Heart/CoronaryArtery/FSI/Explicit_P1P1_PETSc_Seq/CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default meshes are `resources/examples/Heart/CoronaryArtery/Stenosis80.mesh` and `Prismatic80.mesh` in the same directory.
