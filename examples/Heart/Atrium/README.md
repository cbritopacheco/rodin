# Atrium

`Atrium` advances an active solid wall, fixed-domain incompressible flow, or
contraction-driven fluid–structure interaction on a common partitioned mesh.
All example types belong to `Rodin::Examples::Heart::Atrium`.

## Architecture

One executable, `Atrium`, constructs the coordinator and selects its components
with `-atrium_mode`:

| Mode | Components and evolution |
|---|---|
| `solid` | Active solid mechanics with passive initialization and Newmark integration |
| `fluid` | Fixed-domain stabilized P1–P1 Oseen flow; the default mode |
| `fsi` | Solid and fluid with ALE motion and converged Robin–Robin transmission |

`Atrium` owns the common MPI mesh, selected `Solid` and `Fluid` components,
the optional `FSI` coupling object, cavity-volume measurements, and one CSV
writer. `initialize()` prepares the components, `step()` advances one accepted
physical timestep, and `finish()` flushes final states and closes outputs.
These operations are collective. FSI freezes both previous states during
coupling iterations and commits only after interface convergence.

| Location | Responsibility |
|---|---|
| `main.cpp` | Executable entry point and run loop |
| `Parameters.h/.cpp` | Common configuration, grouped component parameters, and PETSc option parsing |
| `RunParameters.h` | Time and output settings supplied to both components |
| `Atrium.h/.cpp` | Mesh ownership, component lifetime, stepping, and shared diagnostics |
| [Solid](Solid/README.md) | Passive, viscous, and active material laws; mechanics Problems and output |
| [Fluid](Fluid/README.md) | Oseen flow, stabilization, retained Problems, solver policy, and output |
| [FSI](FSI/README.md) | Robin–Robin transmission and harmonic ALE motion |
| `CavityVolume.h/.cpp` | Cavity volume from the reference mesh and displaced interface |
| [Visualization](Visualization/README.md) | ParaView volume rendering and scene generation |
| `Tests/` | Fluid, coupled, and configuration checks |

Active contraction is part of the solid constitutive model, rather than a
separate executable or mode. `Solid::Output` and `Fluid::Output` observe their
own fields; the coordinator combines their scalar snapshots into one CSV.

## Mesh ownership

`Parameters::meshPath` selects a combined MEDIT volume mesh in mm, with fluid
material 1 and solid material 4. The coordinator loads it on rank zero,
Scotch-partitions it once, and distributes reconciled MPI shards. Solid and
fluid components extract submeshes from each rank's parent shard. Parent
vertex identities define the shared interface; components do not load or
repartition the source. The parent outlives both components.

Each component retains its finite-element spaces, fields, Problems, linear
systems, solvers, and output. The solid retains passive-equilibrium and mixed
transient mechanics Problems. The fluid retains flow, projection, and
wall-shear Problems. FSI retains interface scatter plans and mesh-motion
assembly. `getMesh()`, `getSolid()`, and `getFluid()` expose the selected objects;
disabled components return null.

## Configuration

`Parameters` contains common run settings and the `solid`, `fluid`, and `fsi`
parameter groups. Both components receive the same time and output settings;
there are no separate solid/fluid timestep or output-path options.

| Common option | Default | Meaning |
|---|---|---|
| `-atrium_mode` | `fluid` | `solid`, `fluid`, or `fsi` |
| `-atrium_mesh` | `resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh` | Combined mesh |
| `-atrium_dt` | `0.001` | Common timestep in seconds |
| `-atrium_steps` | `24` | Physical advances, excluding the initial state |
| `-atrium_output` | `output/atrium` | Run output directory |
| `-atrium_output_every` | `1` | Export every this many advances; zero keeps endpoints and fluid cycle boundaries |

Component options use `-atrium_solid_*`, `-atrium_fluid_*`, and `-atrium_fsi_*`.
For example, `-atrium_solid_contractility 20` sets the active stress source in
kPa, `-atrium_fluid_pin 1200` sets mean inlet pressure in Pa, and
`-atrium_fsi_atol 1e-4` sets the interface residual tolerance in Pa. The apex is
specified as three comma-separated coordinates in mm:
`-atrium_solid_apex 46.39638367,187.44307667,-205.03317767`.
`-help` describes the full option set without loading the mesh.

PETSc solver options retain their standard spelling, including `-ksp_type`,
`-pc_type`, `-snes_monitor`, and `-snes_converged_reason`. The default linear
solver is GMRES with overlapping additive Schwarz, local ILU, and diagonal
scaling with matrix restoration. Component projection solvers use their own
PETSc option prefixes.

`-atrium_steps 24 -atrium_dt 0.001` advances 24 intervals to
$0.024\,\mathrm{s}$ and saves 25 states when the export interval is one.
A full default fluid pressure period uses 850 intervals. Steps control run
length in every mode; optional species warm-up cycles do not override it.

The output directory contains `Solid.xdmf` and/or `Fluid.xdmf`, companion
per-rank HDF5 files, and `Diagnostics.csv`. Solid geometry and model metadata
are written to `Solid.geometry.json`. Diagnostics are flushed at every accepted
step regardless of field export cadence. Solid coordinates are reference
coordinates with displacement in mm; fluid coordinates are physical coordinates
in metres and already include ALE motion in FSI.

## Units and mesh resource

The parent mesh uses mm. The solid retains mm, kPa, and seconds; the fluid
submesh converts coordinates to metres with `fluid.meshScale` and uses SI units.
The conversion changes neither parent ancestry nor interface vertex identities.
`FSI` converts displacement, velocity, and traction explicitly during exchange.

`resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh` contains a
coarse fluid volume sharing the repaired wall's inner surface and a single
outward wedge layer with nodal extrusion distance $1\,\mathrm{mm}$. The wall's
surface preparation includes a local handle-to-disk repair and MMGS adaptation.
The fluid is tetrahedralized against that surface and adapted with MMG3D.
Shared wall triangles and vertices and opening-rim edges are marked as required;
cap interiors and fluid tetrahedra remain free to adapt. The solid coordinates,
wedge connectivity, and shared interface triangulation are preserved. The geometry is a demonstration
model requiring anatomical assessment.

`resources/examples/Heart/Atrium/AtriumFine_FluidSolid_1mm.mesh` provides finer
fluid tetrahedra and opening caps. The shared wall triangulation and opening
rims are prescribed during MMG3D adaptation; the outward solid layer retains
the same interface and nodal thickness. The coarse mesh remains the default.

## Enclosed cavity volume

`CavityVolume` measures the filled reference fluid region and the volume swept
by the moving inner solid wall. The common mesh supplies that fluid region even
in solid-only mode. With cavity-outward normal $N$, the divergence
theorem gives

$$
V_s(t)=V_0+\frac13\int_{\Gamma_{fs,0}}
\left[(X+u)\cdot\mathrm{cof}(I+\nabla_Xu)N-X\cdot N\right]\,dS_0.
$$

Inlet and outlet caps remain fixed, consistent with the prescribed zero rim
displacement. The surface integral uses the displaced face tangents; it does
not extend displacement into the fluid interior. MPI ownership excludes overlap
replicas from both cell and interface integrals.

The coordinator writes **one CSV per run** at initialization and after every
completed step. It combines displacement, activation, active-stress extrema,
Newton convergence, cavity volumes, and the enabled fluid pressure, flux,
wall-shear, and transport diagnostics. Volumes are in mL. All modes write
`<output-directory>/Diagnostics.csv`. Component
writers supply cached diagnostic snapshots to this common writer. Disabled
quantities have empty entries; fluid diagnostics are available after the first
completed flow step. All modes use the same column names.

`getSolidCavityVolume()` and `getFluidDomainVolume()` expose the latest volume
values without MPI collectives. Each component's `getDiagnostics()` similarly
exposes the most recent scalar snapshot.

Fluid mode has a fixed domain, so its volume remains constant. In FSI mode,
the solid-enclosed cavity and moving fluid domain describe the same geometry.
The CSV includes their volume difference, opening and wall fluxes, interface
velocity and geometry mismatches, coupling residual, and minimum ALE Jacobian.

Field states are written and flushed at the common export interval, including
the initial and final states. Fluid cycle boundaries are also exported. Open the XDMF series
in ParaView, apply **Warp By Vector** using **Displacement**, and use
**File > Reload Files > Reload existing file(s)** to discover completed states.
The CSV diagnostics are flushed at each step as well.


## Build and run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DRODIN_USE_MPI=ON -DRODIN_USE_PETSC=ON -DRODIN_USE_SCOTCH=ON \
  -DRODIN_CODE_COVERAGE=OFF
cmake --build build --target Atrium -j4
mpirun -n 4 build/examples/Heart/Atrium/Atrium \
  -atrium_mode solid \
  -atrium_mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  -atrium_output output/atrium-solid \
  -atrium_dt 0.001 -atrium_steps 24 -atrium_output_every 1 \
  -atrium_solid_contractility 20 -atrium_solid_cache_constitutive \
  -snes_monitor -snes_converged_reason
```

Select `fluid` or `fsi` with the same common options. The component READMEs
describe their models and provide pressure-driven fluid and contraction-driven
FSI commands.

Programmatic configuration uses the same parameter groups:

```cpp
using Atrium = Rodin::Examples::Heart::Atrium::Atrium;
Atrium::Parameters parameters;
parameters.mode = Atrium::Mode::FSI;
parameters.timeStep = 0.001;
parameters.steps = 24;
parameters.outputDirectory = "output/atrium-fsi";
parameters.solid.contractility = 20;
parameters.fluid.preconditionerLag = 5;

Atrium simulation(context, parameters);
simulation.initialize();
for (size_t step = 0; step < parameters.steps; ++step)
  simulation.step();
simulation.finish();
```

The numerical checks require MPI and MUMPS:

```sh
python3 examples/Heart/Atrium/Tests/Fluid.py --mpiexec mpiexec
python3 examples/Heart/Atrium/Tests/FSI.py --mpiexec mpiexec
python3 examples/Heart/Atrium/Tests/Parameters.py --mpiexec mpiexec
```

The checks cover resting states, pressure reversal, pressure-reference shifts,
MPI consistency, preconditioner refresh, time refinement, and interface
convergence on a conforming tube. They do not establish spatial convergence on
the anatomical mesh or periodic physiological flow.

## Console diagnostics

A physical timestep starts with its advance index, target simulation time,
and common timestep size in seconds. Timed stages use `|---------` separator
lines; indentation shows the nesting of FSI coupling iterations, solid and fluid
solves, ALE motion, and diagnostic/export operations. Each stage prints `begin`
and a matching `end wall_time=... s (root)` line. A stage interrupted by an
exception prints `failed` instead of `end`.

Solid mechanics separates the initial assembly from Newton, whose elapsed time
includes all residual and tangent reassembly inside SNES. Fluid stages separate
assembly, field splits, KSP/preconditioner setup, retry when needed, wall-shear
recovery, fluxes, species transport when enabled, and field export. Coupled
iterations also time ALE motion, traction exchange, and interface convergence
and relaxation. The outer timestep includes the shared cavity-volume and CSV
work. These are inclusive wall times on rank zero: parents contain child work
and communication waits. They are neither CPU times nor maximum-over-rank
measurements, and nested timings should not be added together. Instrumentation
adds no MPI synchronization.

Indented component messages retain `[Solid/Newton]`, `[Fluid/KSP]`, and
`[FSI/Coupling]` labels for solver iterations, residuals, convergence reasons,
preconditioner status, and interface tolerances. Native PETSc monitors enabled
with `-snes_monitor` or `-ksp_monitor` appear between the corresponding timed
solve separators and retain PETSc's own formatting. Residual norms retain the
solver-dependent scaling. State quantities carry physical units, real values
use six decimal places in scientific notation, and the shared CSV retains full
numeric precision. Fluid phase totals accumulate across coupling trials;
its `elapsed_step` spans the staged step, including the coupled solid work.

Stage markers use the terminal palette: magenta for solid, blue for fluid,
cyan for coupling, and yellow for ALE or retries. Success and failure markers
use green and red. Labels, status words, times, and values always use the
terminal's default foreground; no information depends on a marker's color or
contrast. Timestep and diagnostic headings use the default foreground, with no
fixed white, gray, background, or dim styling. Labels are bold when styling is
available.

Redirected stage logs are plain text. `NO_COLOR=1`, `CLICOLOR=0`, or `TERM=dumb`
disable styling. MPI launchers can forward output through pipes; set
`CLICOLOR_FORCE=1` before `mpiexec` to retain colored markers in that case.
`NO_COLOR=1` takes precedence over forced color. Native PETSc monitors retain
their own formatting.
