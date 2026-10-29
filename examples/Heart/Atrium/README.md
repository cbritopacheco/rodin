# Atrium

`Atrium` advances a solid wall, an incompressible fluid, or their coupled
dynamics. A combined volume mesh contains fluid tetrahedra with material
label 1 and solid wedges with material label 4. Their interface triangles share
vertices in the parent mesh.

All example types belong to `Rodin::Examples::Heart::Atrium`. The `Atrium`
coordinator owns `Solid` and `Fluid` components. The component aliases on the
coordinator remain available as `Atrium::Solid` and `Atrium::Fluid`.

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
time increments. `FSI` freezes both previous physical states while iterating
interface motion and traction. Both components commit their states only after
interface convergence; a failed step terminates the run.

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

All reported opening fluxes use the outward fluid normal: inflow is negative
and outflow is positive. For pressure-penalty coefficient $\epsilon_p$, the
constant pressure test gives the global balance

$$
Q_{\mathrm{in}}+Q_{\mathrm{out}}
=-\epsilon_p\int_{\Omega_f}p\,dX.
$$

This global balance does not enforce pointwise zero divergence. The reported
maximum velocity is $\max_X\|u_h(X)\|_2$, evaluated from nodal speeds; for
an affine P1 tetrahedral field the maximum occurs at a vertex.

## Contraction-driven FSI

Fluid inertia feeds back on the thin atrial wall, so a single motion/traction
exchange need not produce a consistent timestep. `FSI` uses a partitioned
Robin–Dirichlet iteration. With the fluid-outward normal $n_f$ and fluid stress
$\sigma_f=-pI+2\mu D(u_f)$, the interface conditions are

$$
u_f=v_s,\qquad d_f=d_s,\qquad t_s=-\sigma_f n_f.
$$

Here $d_f$ is the fluid mesh displacement and $t_s$ is the spatial load on the
solid. Its reference-area value is $T_s=J_\Gamma t_s$, where $J_\Gamma$ is the
current-to-reference interface area ratio. The solid trial solves

$$
P_sN_s+\alpha v_s=g^k,\qquad
\widehat g^{k+1}=J_\Gamma(-\sigma_f n_f)+\alpha u_f.
$$

This is Robin–Dirichlet transmission: the fluid velocity is imposed strongly.
It is distinct from the [coronary Robin–Robin operator](../CoronaryArtery/FSI/README.md#transmission-conditions),
which also places a Robin condition on the fluid.

Both Robin terms use reference area. At convergence, velocity continuity
cancels them and recovers traction equilibrium. The coefficient $\alpha$ has
units $\mathrm{Pa\,s/m}$; its automatic value is
$\rho_f V_f^{1/3}/\Delta t$. This impedance scale is a numerical choice,
not a tissue parameter. Aitken relaxation updates
$g^{k+1}=g^k+\omega_k(\widehat g^{k+1}-g^k)$, with bounded factors. Convergence
uses the unrelaxed nodal Robin residual, measured in Pa. Fluid traction is
projected onto the reference interface from pressure and the recovered
velocity gradient; this is a stress-based discrete transfer.

The fluid receives the solid velocity as an essential condition on label 4.
A harmonic extension of solid displacement moves the fluid interior; opening
caps remain fixed. With mesh velocity $a=(d_f^{n+1}-d_f^n)/\Delta t$ and cell
volume ratio $r_K=|K^n|/|K^{n+1}|$, the conservative first-order ALE residual is

$$
R_m^{\mathrm{ALE}}=
\rho\frac{u_f^{n+1}-r_Ku_f^n}{\Delta t}
+\rho[(u_f^n-a)\cdot\nabla]u_f^{n+1}
+\rho\left(\frac12\nabla\cdot u_f^n-\nabla\cdot a\right)u_f^{n+1}
-\nabla\cdot(2\mu_nD(u_f^{n+1}))+\nabla p^{n+1}.
$$

SUPG and PSPG use this residual and the relative advecting velocity. Previous
velocity and solid internal states remain fixed through the coupling trials.
The fluid boundary receives the Newmark solid velocity, whereas the ALE mesh
velocity uses a backward displacement difference. These rates differ at finite
time step; exact time-discrete agreement is not imposed.

Mesh motion invalidates Rodin's geometric caches through the mesh coordinate
API. The coupled driver retains the same spaces, Problems, solvers, and MPI
interface scatter plans. Gradient and traction projections are reassembled.

All openings have zero imposed pressure and no prescribed flux. Fluid starts
at rest and active wall contraction drives the flow. Backflow stabilization
remains active. FSI currently excludes the 0D ventricular driver, species
transport, and orthogonal subscales. The outer solid surface is unloaded and
opening rims are clamped. No valve or circulatory model is represented.

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
wall-shear, and transport diagnostics. Volumes are in mL. Solid and `FSI` modes
use `<solid-output>.csv`; fluid mode uses its configured `csvPath`. Component
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

Solid states are written and flushed at each time step. Open the XDMF series
in ParaView, apply **Warp By Vector** using **Displacement**, and use
**File > Reload Files > Reload existing file(s)** to discover completed states.
The CSV diagnostics are flushed at each step as well.

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
  -atrium_steps 24 -atrium_output_every 1 -atrium_dt 0.001 \
  -atrium_lv_model false -atrium_kinetics false
```

Run contraction-driven FSI on the common coarse mesh:

```sh
mpirun -n 4 build/examples/Heart/Atrium/Atrium \
  -atrium_mode fsi \
  -atrium_mesh resources/examples/Heart/Atrium/AtriumCoarse_FluidSolid_1mm.mesh \
  --output output/atrium-fsi/Solid -atrium_output output/atrium-fsi/Flow \
  --dt 0.001 -atrium_dt 0.001 -atrium_steps 24 \
  --cache-constitutive -atrium_output_every 1 \
  -atrium_fluid_lag_preconditioner 5 \
  -ksp_type gmres -ksp_rtol 1e-12 -ksp_max_it 500 \
  -ksp_norm_type unpreconditioned -pc_type lu \
  -pc_factor_mat_solver_type mumps -mat_mumps_icntl_20 0 \
  -ksp_diagonal_scale false -ksp_diagonal_scale_fix false
```

`-atrium_fsi_iterations`, `-atrium_fsi_robin`, `-atrium_fsi_rtol`, and
`-atrium_fsi_atol` configure the interface solve. The solid XDMF stores reference
geometry and `Displacement`; the fluid XDMF stores moving physical geometry,
`MeshDisplacement`, and `MeshVelocity`. Apply a displacement warp only to the
solid. The common diagnostics are written to `output/atrium-fsi/Solid.csv`.

```sh
python3 examples/Heart/Atrium/Tests/FSI.py --mpiexec mpiexec
```

The coupled checks use a conforming tetrahedron/wedge tube and examine rest,
contraction-driven flow, interface continuity, global flux balance, common
cavity volume, MPI consistency, and fluid feedback on the solid.

The fluid checks use a generated tetrahedral tube and require MPI and MUMPS:

```sh
python3 examples/Heart/Atrium/Tests/Fluid.py --mpiexec mpiexec
```

`--executable` selects a configured driver and `--output` selects a result
folder. The checks cover rest, pressure reversal, pressure-reference shifts,
MPI consistency, global flux balance, and outlet-flux self-convergence under a
smooth pressure source. They do not establish spatial convergence on the
anatomical mesh or periodic atrial hemodynamics.

`-atrium_mode fsi` selects both components on the same parent partition.
Use `--dt` for the solid and `-atrium_dt` for the fluid, with equal values in
seconds. Solid options use double dashes; fluid and PETSc solver options use
single dashes. The executables default to GMRES with overlapping additive
Schwarz, local ILU, and diagonal scaling with matrix restoration, configurable
through PETSc options. Output consists of a master XDMF file and per-rank HDF5
pieces.

`-atrium_fluid_lag_preconditioner N` refreshes the fluid preconditioner on the
first step and every $N$ timesteps thereafter. Its default is one. Larger
values require an iterative solver such as GMRES. FSI trials at the same
physical timestep retain the same scheduled factors; a lag of one refreshes
every trial. A failed solve with retained
factors is retried once with a fresh preconditioner before the timestep is
committed. Symbolic analysis remains reusable while sparsity is unchanged.
The policy takes precedence over `-ksp_reuse_preconditioner` for the fluid KSP.
The common CSV records Krylov iterations, preconditioner age, refresh count,
and whether a refresh or retry occurred.

Fluid arrays use PascalCase names: `Velocity`, `Pressure`, `WallShearStress`,
`TimeAveragedWallShearStress`, and `OscillatoryShearIndex`. The last two are
cycle statistics, updated after a completed cycle. `Thrombin`, `Fibrinogen`,
`Fibrin`, and `ThrombosisActivation` are exported only when transport kinetics
is enabled with `-atrium_kinetics true`; species evolution begins after the
configured flow warm-up cycles. `ThrombosisActivation` is a wall-shear-based
coagulation indicator, distinct from the solid's contractile activation.

The [visualization pipeline](Visualization/README.md) renders fluid speed as a
volume animation with a translucent cavity surface and physical timestamps.

A program linked to `AtriumModels` can configure the common source and advance
the selected components directly:

```cpp
using Atrium = Rodin::Examples::Heart::Atrium::Atrium;
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
