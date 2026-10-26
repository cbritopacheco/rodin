# Atrium models

The atrium examples contain a distributed fluid model and a separate
[active solid model](Solid/README.md). The fluid model uses a fixed tetrahedral
mesh. The solid model uses a wedge wall and prescribed fibers and activation.
Each can be run independently; fluid–structure coupling is not implemented.

## Fluid model

The fluid problem solves incompressible, non-Newtonian flow with a
Carreau–Yasuda viscosity. For velocity $u$, pressure $p$, density $\rho$, and
symmetric velocity gradient $D(u)$, the interior equations are

$$
\rho\left(\frac{\partial u}{\partial t}+(u\cdot\nabla)u\right)
-\nabla\cdot\left(2\mu(\dot\gamma)D(u)\right)+\nabla p=0,
\qquad \nabla\cdot u\approx 0.
$$

The equal-order velocity–pressure discretization uses PSPG and optional VMS
stabilization. Time integration uses backward Euler with lagged advection and
viscosity. Boundary labels 1, 2, 3, and 5 are pressure-driven inlets; label 6
is the outlet and label 4 is the no-slip wall.

The default outlet pressure is supplied by a lumped left-ventricular model.
That model receives a prescribed atrial-pressure waveform; the computed 3D
atrial flow is not fed back into it. The mitral resistance switches using the
previous outlet pressure. An alternative prescribed-pressure outlet disables
the ventricular time step. Thrombin, fibrinogen, and fibrin transport is
optional; wall-shear-derived quantities are retained for its boundary model.

Coordinates in the source mesh are scaled by `meshScale` before the flow
solve. The default scale converts millimetres to metres. Fluid velocity,
pressure, and density use SI units. The solid example documents its own units.

## Architecture

`Atrium` coordinates the components selected by `Atrium::Parameters`:
`SolidOnly`, `FluidOnly`, or `Independent`. The independent mode advances
both components sequentially with the same time step; it does not exchange
interface motion or traction. The local solid currently requires a single MPI
rank. The fluid-only mode can use the distributed MPI mesh.

`Atrium::Solid` owns the wedge mesh, fibers, material, mechanics problem,
solver, and output. `Atrium::Fluid` owns the distributed volume mesh, fields,
flow and transport problems, solvers, and fluid output. Each component retains
its objects between steps. Fluid forms are constructed during initialization
and reassembled as their coefficients change, using Rodin's mesh-owned
transformation and quadrature caches. Each output class writes its component's
fields; there is no shared FSI output yet.

`Atrium::step()` initializes the selected components on its first call and
advances them once. `getTime()`, `getStepCount()`, `getSolid()`, and
`getFluid()` expose their state; a disabled component is null. `finish()`
finalizes output. Calls involving the fluid are collective on its communicator.
A failed coordinated step cannot be retried, since a component may already
have advanced. Independent mode requires a positive `solid.timeStep` equal to
`fluid.dt`, both in seconds. It does not support the solid fiber-inspection
mode.

The fluid component can also be used directly through `AtriumFluid`. Set
`useVentricularModel` and `solveKinetics` to `false` for a prescribed-outlet
flow solve, and set `prescribedOutletPressure` in pascals. A prescribed outlet
resistance in Pa s/m³ can be supplied with `prescribedOutletResistance`.
When the ventricular model is disabled, the CSV `pv` column records the
prescribed outlet pressure; the unavailable arterial and distal pressures
are written as `nan`.

## Build and run

With PETSc and MPI enabled, build from the repository root:

```sh
cmake --build build --target Atrium -j4
cd build
mpirun -n 1 ./examples/Heart/Atrium/Atrium
```

The default volume is `resources/examples/Heart/Atrium/Coarse.mesh`; hydrate
mesh resources with `git lfs pull`. Select a mesh with `-atrium_mesh`, the
time step with `-atrium_dt`, and the cycle counts with `-atrium_flow_cycles`
and `-atrium_species_cycles`. The prescribed-outlet mode is selected with
`-atrium_lv_model false`, `-atrium_outlet_pressure`, and optionally
`-atrium_outlet_resistance`. Disable transport with `-atrium_kinetics false`.
The executable resolves its default mesh path from the build root.

A program in the same CMake build can link `AtriumModels` and advance the
selected model directly:

```cpp
using Atrium = Rodin::Examples::Heart::Atrium;
Atrium::Parameters parameters;
parameters.mode = Atrium::Mode::FluidOnly;
parameters.fluid.useVentricularModel = false;
parameters.fluid.solveKinetics = false;
parameters.fluid.prescribedOutletPressure = 500.0; // Pa

Atrium atrium(context, parameters);
atrium.step();
// Read atrium.getFluid()->getVelocity() and getPressure().
atrium.finish();
```

The `AtriumModels` target exposes fluid modes when both PETSc and MPI are
enabled. The `AtriumSolid` and `AtriumFluid` targets can also be linked
independently. The solid-only coordinator can be constructed without an MPI
context; see the [solid example](Solid/README.md) for its parameters.
