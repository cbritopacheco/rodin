# Coronary staggered FSI with PETSc/MPI

This driver constructs PETSc/MPI fluid and solid problems, but its cross-mesh
interface map currently requires a **single MPI rank**. It is not a distributed
interface-transfer implementation. The fluid uses continuous P2 velocity and P1
pressure and P1 wall displacement. The harmonic ALE field uses the P2 velocity
space; mesh coordinates are updated at vertices.

## Model and transmission

The [shared coronary equations](../README.md#component-models) combine
Carreau–Yasuda blood and a compressible Neo-Hookean wall. Both subproblems use
[Robin–Robin data](../README.md#transmission-conditions),
including the fluid interface convective correction. The solid Robin terms and
fluid traction load carry the current-to-reference area ratio $J_\Gamma$.

Static follower-pressure prestress can initialize the wall when
`-coronary_prestress_steps` is positive; it is disabled by default. The dynamic load keeps
this follower-pressure level and adds the projected traction remainder. Solid
opening caps and interface ring band are clamped. Plain harmonic extension
moves the fluid mesh while fixing its opening caps. The fluid Robin target uses
the transferred Newmark solid velocity; its BDF1 mesh velocity is computed
separately from displacement differences.

## Time integration and circulation

The fluid uses conservative BDF1 ALE/Oseen integration, lagged viscosity,
Temam/geometric correction, a pressure penalty, and boundary backflow terms. The wall uses
Newmark with $\beta=0.3025$ and $\gamma=0.6$. The default permits three interface
corrections with bounded Aitken relaxation of the displacement iterate.
`-coronary_coupling_iterations` and `-coronary_coupling_tolerance` configure the
correction loop. A fixed correction budget is not a guarantee of convergence.

The 0D ventricular model supplies inlet pressure. Non-Newtonian proximal and
distal branch laws close intramyocardial outlet compartments, with

$$
C(\dot p_c-\dot p_{\mathrm{im}})=Q-q_d,\qquad
p_{\mathrm{out}}=p_c+\Delta p_p(Q).
$$

The [sequential Taylor–Hood model](../Explicit_PETSc_Seq/README.md) describes
this circulation closure. The driver does not include an active coronary-wall
fiber model or a resolved ventricle.

## Inputs and output

Mesh coordinates are converted to metres by the default factor $10^{-3}$.
Pressure, stress, velocity, and time use SI units. XDMF stores fluid and wall
fields; CSV records ventricular/outlet states, fluxes, wall motion, and interface
observables. The MPI label on the driver describes its algebra/mesh types and
does not remove the single-rank interface-map restriction.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery_FSI_Explicit_PETSc_MPI -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/FSI/Explicit_PETSc_MPI/CoronaryArtery_FSI_Explicit_PETSc_MPI
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default meshes are `resources/examples/Heart/CoronaryArtery/FSIFluid.mesh` and `FSISolid.mesh` in the same directory.
