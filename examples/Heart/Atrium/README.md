# Atrium

`Atrium` advances a solid wall, an incompressible fluid, or both components
independently. A combined volume mesh contains fluid tetrahedra with material
label 1 and solid wedges with material label 4. Their interface triangles share
vertices in the parent mesh.

## Mesh ownership and architecture

`Atrium::Parameters::meshPath` selects the common mesh. `Atrium` loads it on
rank zero, partitions it once with Scotch, distributes its shards, and retains
the reconciled MPI mesh. Every rank holds its owned cells and the framework's
existing overlap for shared entities.

`Atrium::Solid` and `Atrium::Fluid` extract distributed `SubMesh` objects from
that parent using their configured material labels. The extraction preserves
parent-to-submesh index maps and reconciles ownership inside each selected
region. Components do not load or repartition meshes. The parent must outlive
both components; `Atrium` owns the parent before its optional component members.

Each component retains its finite-element spaces, fields, `Problem` objects,
linear systems, solvers, and output. The solid retains both its mixed transient
problem and passive-equilibrium problem. The fluid retains its Oseen problem,
wall-shear problem, and projection problems. Assembly reuses the fixed mesh and
Rodin's transformation and quadrature machinery across time steps.

`Atrium::Parameters::mode` selects `Solid`, `Fluid`, or `FSI`.
`initialize()` prepares the selected components, `step()` advances each once,
and `finish()` closes their outputs. These operations are collective on the
parent communicator. FSI mode requires equal positive solid and fluid
time increments. It currently advances the two equations independently without
exchanging interface motion or traction. A failed coordinated step cannot be
retried because one component may already have advanced.

`getMesh()` exposes the combined parent. `getSolid()` and `getFluid()` return
the selected components, or null when disabled. Each component's `getMesh()`
exposes its submesh and ancestry maps. The `AtriumModels`, `AtriumSolid`, and
`AtriumFluid` CMake targets expose the coordinator and components respectively.

## Solid model

The [solid component](Solid/README.md) couples displacement to independent
contractile extension, square-root stiffness, and transformed-stress fields.
`DemirayHolzapfelGasserOgden` supplies the passive exponential energy;
`KelvinVoigt` supplies Green–Lagrange strain-rate dissipation. Passive equilibrium
initializes the wall before active dynamics. Newmark time integration and Rodin's
active-fiber update form a mixed PETSc/SNES problem.

## Fluid model

The fluid uses continuous linear velocity and pressure on affine tetrahedra.
For lagged velocity $w=u^n$ and viscosity $\mu_n$, the Oseen momentum residual is

$$
R_m=\rho\frac{u^{n+1}-u^n}{\Delta t}
 +\rho(w\cdot\nabla)u^{n+1}
 +\frac{\rho}{2}(\nabla\cdot w)u^{n+1}
 -\nabla\cdot(2\mu_n D(u^{n+1}))+\nabla p^{n+1},
\qquad D(u)=\frac{\nabla u+\nabla u^T}{2}.
$$

The Temam term accounts for a lagged velocity that is not exactly divergence
free. Carreau–Yasuda viscosity is evaluated at the lagged strain rate and is
constant within each affine linear tetrahedron. The elementwise strong viscous
divergence therefore vanishes in the stabilization residual.

The equal-order method adds SUPG and PSPG terms

$$
\sum_K\int_K R_m\cdot
 \left[\tau_s\rho(w\cdot\nabla)v+\tau_p\nabla q\right]dX,
$$

and a grad-div term. The same time, advection, Temam, and pressure contributions
enter both stabilization residuals. With element size $h_K$, the parameters are

$$
\tau_1=\left(\frac{4\mu_n}{\rho h_K^2}
                    +\frac{2|w|}{h_K}\right)^{-1},\qquad
\tau_s=\frac{s_s}{\rho/\Delta t+\rho/\tau_1},\qquad
\tau_p=\frac{s_p}{\rho/\Delta t+\rho/\tau_1}.
$$

Backward Euler advances time. Boundary labels 1, 2, 3, and 5 are pressure-driven
inlets, label 6 is the outlet, and the extracted material-4 interface is the
fixed no-slip wall. Pressure penalties, boundary impedance, and backflow terms
are retained in the flow problem. The default uses a prescribed outlet pressure.
Ventricular coupling, species transport, and orthogonal-subscale corrections are
optional parameters, disabled by default.

## Units and mesh resource

The parent mesh uses mm. The solid retains mm, kPa, and seconds; the fluid
submesh converts coordinates to metres with `fluid.meshScale` and uses SI units.
The conversion changes neither parent ancestry nor interface vertex identities.
An interface exchange must explicitly convert displacement and traction units.

`resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh` contains a
coarse fluid volume sharing the repaired wall's inner surface and a single
outward wedge layer with nodal extrusion distance $1\,\mathrm{mm}$. The wall's
surface preparation includes a local handle-to-disk repair and MMGS adaptation.
The fluid is tetrahedralized against that surface and coarsened with MMG3D
while preserving its boundary triangulation. The geometry is a demonstration
model requiring anatomical assessment.

## Build and run

From the repository root, configure MPI, PETSc, and Scotch and build the examples:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DRODIN_USE_MPI=ON -DRODIN_USE_PETSC=ON -DRODIN_USE_SCOTCH=ON \
  -DRODIN_CODE_COVERAGE=OFF
cmake --build build --target Atrium AtriumContraction -j4
```

Run the solid for 24 time intervals plus its initial state:

```sh
mpirun -n 4 build/examples/Heart/Atrium/Solid/AtriumContraction \
  --mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  --output output/atrium-solid/Contraction --steps 24 \
  --activation 20 --relaxation 10 --contractility 20 \
  --viscosity 0.2 --active-damping 3
```

Run the fluid independently for 24 time intervals:

```sh
mpirun -n 4 build/examples/Heart/Atrium/Atrium \
  -atrium_mode fluid \
  -atrium_mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  -atrium_output output/atrium-fluid/Flow \
  -atrium_steps 24 -atrium_output_every 1 -atrium_dt 0.001
```

`-atrium_mode fsi` selects both components on the same parent partition.
Use `--dt` for the solid and `-atrium_dt` for the fluid, with equal values in
seconds. Solid options use double dashes; fluid and PETSc solver options use
single dashes. The executables default to GMRES with overlapping additive
Schwarz, local ILU, and diagonal scaling with matrix restoration, configurable
through PETSc options. Output consists of a master XDMF file and per-rank HDF5
pieces.

A program linked to `AtriumModels` can configure the common source and advance
the selected components directly:

```cpp
using Atrium = Rodin::Examples::Heart::Atrium;
Atrium::Parameters parameters;
parameters.meshPath = "resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh";
parameters.mode = Atrium::Mode::FSI;
parameters.solid.timeStep = 0.001;
parameters.fluid.dt = 0.001;

Atrium atrium(context, parameters);
atrium.initialize();
atrium.step();
atrium.finish();
```
