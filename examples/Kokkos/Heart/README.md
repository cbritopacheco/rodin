# Atrium on the device

`KokkosAtriumFlow` is the first step of the port of
`examples/Heart/Atrium/Viscoelastic` (branch `model/Atrium`) to the Kokkos
module. It advances a Newtonian fluid in the fixed atrial volume over cardiac
cycles, with the assembly and the linear solve on the device and one process.

The same source builds two executables that differ only in the linear algebra:

| Executable | Assembly | Linear solver |
|---|---|---|
| `KokkosAtriumFlow` | native Kokkos sparse matrix | cuDSS (`RODIN_USE_CUDSS`), else Kokkos-Kernels GMRES with Jacobi |
| `KokkosAtriumFlowPETSc` | COO assembly into a PETSc matrix of the Kokkos type | PETSc KSP on the device: BiCGStab with Jacobi by default, any PETSc option on the command line |

`KokkosAtriumFlowPETSc` is built when Rodin finds a PETSc configured with
Kokkos (`Rodin::Kokkos::PETSc`). Each solve starts from the previous step.

`KokkosAtriumRheologyPETSc` (`AtriumRheology.cpp`, PETSc path only) adds the
fluids of the CPU driver: see [Rheology](#rheology).

## Ported

- The domain: a tetrahedral MEDIT mesh of the fluid in millimetres
  (`resources/examples/Heart/Atrium/Coarse.mesh`), labels 1, 2, 3, 5 veins,
  6 mitral plane, 4 wall.
- The openings of Atrium/Fluid: venous traction $-p_{pv}(t)\,n$ with
  impedance, mitral traction $-(p_v+Z\,u\cdot n)\,n$ with the lumped ventricle
  `Rodin::Heart::CCMLC2014` and its valve law, no slip on the wall. The
  ventricle runs on the host.
- The time loop, the flow rates and mean pressures of the openings, CSV and
  XDMF output.

## Different from the CPU driver

- One linear Oseen step per time step (backward Euler, $w=u^n$), not Newton
  with BDF2 and step control.
- Skew-symmetric convection in place of the backflow term: energy neutral for
  either direction of the flow through an opening. The natural condition then
  carries $\tfrac12\rho(w\cdot n)u$, so the opening tractions are not
  identical to those of Atrium/Fluid.
- Streamline diffusion and a pressure Laplacian with the cellwise
  $\tau=[(2/\Delta t)^2+(2|w|/h)^2+(4\nu/h^2)^2]^{-1/2}$, evaluated on the device.
  They are not residual based: no orthogonal projection, no dynamic subscale,
  and a consistency error of order $\tau$.
- Venous impedance on the whole velocity; Atrium/Fluid uses equal normal and
  tangential coefficients, so the two coincide at its defaults.

- The ventricle reads $p_v+Q/K$ as its atrial pressure, the mean of the
  balance that the natural condition imposes on the mitral plane, and not the
  mean of $p$: with the skew-symmetric convection that balance holds for
  $p+\tfrac12\rho(u\cdot n)^2$, so the mean of $p$ falls below $p_v$ at the
  peak of the filling, which closed the valve in mid-diastole and left the
  ventricle with less than the volume that crossed the plane. On the test mesh
  $p_v$ at the end of the first cycle is now 220 Pa (CPU driver: 218 Pa; mean
  of $p$: 20 Pa).

## Rheology

`KokkosAtriumRheologyPETSc mesh.mesh fluid=N|CY|VE|TVE` runs the four fluids
of the CPU driver, which share the steady shear viscosity of Carreau-Yasuda:
Newtonian, Carreau-Yasuda, Owens with the aggregate size at equilibrium (VE)
and Owens with transported aggregate size (TVE). The file header states the
scheme and its references. For the elastic fluids, $\tau_p=G(c-I)$ and
$\psi=\log c$ is stored at the vertices; per time step:

1. velocity gradient $L^n$ of $u^n$ at the vertices: mean of the cell
   gradients, trace free, and $(\partial_n u)\,n^T$ on the wall;
2. $\psi$ and $n$ at the feet of the characteristics, $x-\Delta t\,u^{n+1/2}$,
   located by a walk through the cells and interpolated linearly: monotone,
   no linear system; developed blood enters by the veins;
3. prediction at every vertex with $L^n$: $n$ by the exact solution of its
   Riccati equation, $c$ by the exact solution of its equation with the
   gradient frozen;
4. momentum and continuity (BDF2, convection by $u^*=2u^n-u^{n-1}$, algebraic
   subgrid scales with quasi-static subscales) with the predicted stress in
   weak form and two terms that vanish for a resolved steady flow:
   $\eta_a(\nabla u^{n+1}-L^n)$, $\eta_a=G\lambda$ (DEVSS-G), and
   $\eta_d(\nabla u^{n+1}-\nabla u^n)$,
   $\eta_d=2G\lambda(1-e^{-\Delta t/\lambda})\lambda_{\max}(c)$, the implicit
   response of the stress to the velocity of the step;
5. correction: step 3 again with the gradient of $(u^n+u^{n+1})/2$.

Steps 1-3 and 5 are Kokkos kernels of the example; step 4 is one PETSc
matrix of the Kokkos type. There is no transport matrix. N and CY assemble
only their own terms. The elastic part is first order in time and space.

**Numbering.** The data of a P1 grid function are indexed by degree of
freedom, which the space numbers by cell: on the test mesh 1099 of 1100
vertices have a degree of freedom that is not their index. The vertex
kernels up to commit 2702be8 indexed the velocity by vertex, so the
gradients of CY, VE and TVE were those of a permuted velocity: of the right
size on a mesh of 2 mm, where they went unnoticed, and of order $U/h$ on
small cells, which is the divergence of VE and TVE on `Coarse.mesh`. The
tables of the kernels are now built by degree of freedom, and the run stops
if the numbering is not the expected one. Results of CY, VE and TVE of
earlier commits are void; N has no vertex kernel and is not affected.

**Checks**, all on the 3983-cell test mesh (2.2 mm) or scaled copies of it,
Kokkos Serial:

| check | result |
|---|---|
| steady flow at Wi 0.2 (`scale=1e-4 lv=0 p_in=5`), VE against CY | flow rate -1.4 %; +20 % with `devss=0` |
| the same without the wall form of the gradient | +1.9 %, and $\mathrm{tr}\,c$ to 3e5 on the wall in the test below |
| small cells, `scale=1e-4 lv=0` (0.22 mm, Wi to 28), 200 steps | VE and TVE run, velocity within 20 % of N; $\mathrm{tr}\,c$ peaks at 5.4e6 (VE, one transient) and 1.8e4 (TVE); `implicit=0` diverges at step 194 |
| `scale=5e-5 lv=0` (0.11 mm, Wi to 33), 200 steps | VE and TVE run, $\mathrm{tr}\,c$ to 2.1e3; `implicit=0` diverges at step 145 |
| 1.1 mm refinement (31 864 cells), TVE to 0.25 s | runs, $\mathrm{tr}\,c$ to 39 |
| one cycle with the ventricle, VE and TVE | run; valve closes and reopens as with N |
| time step, TVE, `lv=0`, 0.1 s, 2, 1, 0.5 ms against 0.25 ms | flow rate 3.0e-3, 1.3e-3, 4.3e-4: first order |
| N and CY before and after the split of the forms | equal digits |

Against the CPU driver on the test mesh (BDF2 at 10 ms, orthogonal
subscales), at $t$ = 0.15 s: flow rate 1.01e-4 against 0.95e-4 m^3/s for the
three fluids; Wi 1.6 against 1.8 (VE) and 2.3 against 3.1 (TVE); stress
$G(\mathrm{tr}\,c-3)$ 3.1 against 3.7 Pa (VE) and 7.1 against 8.1 Pa (TVE).

Not run with CUDA, nor on `Coarse.mesh`.

**Checkpoints.** Every `checkpoint_every` steps (by default at the end of
every cycle) and at its last step, a run writes `<output>.restart`
($u^n$, $u^{n-1}$, $\psi$, $n$, the state of the ventricle and the step).
`restart=<file>` continues from it on the same mesh and time step: time,
cycle numbers and the ventricle go on, and `cycles`/`t_end` count from there.
$\psi$ and $n$ come from the file only if the same fluid wrote it, otherwise
they start from $c=I$ and `n_init`, so that the four fluids can start from
one developed Newtonian flow. A run split at a checkpoint reproduces the
continuous run to the last digit (test mesh, N and TVE, 100 steps).

**Stabilisation parameters.** $\tau$ has no $(2/\Delta t)^2$ term: with it
the stabilised problem changes with the time step and the time
discretisation does not converge on a mesh that needs the stabilisation
(Bochev, Gunzburger & Shadid 2004). Without it the time derivative cannot
stay in the residual of the subscales (John & Novo 2011), so they are
quasi-static: consistency error of order $\tau$, spatial. `stab=split` keeps
the two diagonal terms alone.

**Time-scheme dissipation.** The CSV column `timeDissipation` is the kinetic
energy removed by the time scheme per unit time,
$\rho/(4\Delta t)\,\lVert u^{n+1}-2u^n+u^{n-1}\rVert^2$ for BDF2.

**Unbounded stress.** Owens' dumbbells are Hookean and their stress is
unbounded where $\lambda$ times the extension rate exceeds 1/2 (`tr(c)_max`
in the log); `ptt=0.01`, say, bounds it with a linear Phan-Thien-Tanner
factor that is not part of Owens' model.

**Linear solver.** BiCGStab(4) with Jacobi by default. On the 1.1 mm mesh
(Serial), a Schur field split with an algebraic multigrid on the pressure,
set up every ten steps and when the valve switches,

```
pc_lag=10 -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type schur
  -pc_fieldsplit_schur_fact_type lower -pc_fieldsplit_schur_precondition selfp
  -fieldsplit_u_ksp_type preonly -fieldsplit_u_pc_type jacobi
  -fieldsplit_p_ksp_type preonly -fieldsplit_p_pc_type gamg
```

took 17-38 iterations against 150-160 and half the time. Not tried with
CUDA.

## Not ported

gsPTT, the Eulerian age, wall shear indices, probes, the appendage metrics,
the periodic initialisation of the aggregate size (`resetCycle`), orthogonal
subscales with dynamic subscales, BDF2 with step control, and MPI.

## Memory

On an RTX 3060 (12 GB) the native path ran meshes of up to 113 104
unknowns (4.9 GB with GMRES, 5.6 GB with cuDSS) and ran out of memory on
`Coarse.mesh` (561 036 unknowns) in its symbolic assembly
(`"Rodin symbolic ranks"`), with either solver: about 43 kB per unknown. The
PETSc path assembles through the COO pattern of the Kokkos PETSc module,
which is sorted by ranges of rows to bound its memory: `KokkosAtriumFlowPETSc`
peaked at 6.3 GB on `Coarse.mesh` on that card (6.7 GB with ILU(0)).

## Status

`KokkosAtriumFlowPETSc` was compiled with Kokkos 5.2.2 (Serial) and PETSc
3.26.0 configured with Kokkos, and run on the same test mesh: with BiCGStab
and Jacobi (80-130 iterations per step at a relative tolerance of 1e-8) it
reproduces the native run to six digits at the end of the cycle.

`KokkosAtriumFlow` was compiled with Kokkos 5.2.2, Serial backend, and run on
a 3983-cell test mesh with the labels above for one cycle of 85 steps: the
flow rates through the veins and the mitral plane agree to within 0.5 % with
GMRES at every step, the valve closes and reopens with the ventricular
pressure, and areas and volume equal those of the CPU driver. With CUDA, on
an RTX 3060, it ran test meshes of 8 908 to 113 104 unknowns with cuDSS (flow
rates equal to all digits) and with GMRES and Jacobi; see Memory for
`Coarse.mesh`.

`KokkosAtriumFlowPETSc` with CUDA (PETSc 3.26.0, RTX 3060) on `Coarse.mesh`:
BiCGStab with Jacobi broke down at the first step; `-ksp_type bcgsl
-ksp_bcgsl_ell 4 -pc_type jacobi` takes about 420 iterations and 1.8 s per
step, plus 0.44 s of assembly; GMRES(100) with ILU(0)
(`-pc_factor_mat_solver_type kokkos`) about 340 iterations and 3 s. These
were measured before the change of the atrial pressure read by the ventricle.

## Run

```sh
cmake --build build-cuda --target KokkosAtriumFlow
build-cuda/examples/Kokkos/Heart/KokkosAtriumFlow /path/to/Coarse.mesh \
  cycles=2 steps=850 output=atrium
```

Keys: `cycles`, `steps` (per cycle), `period`, `rho`, `mu`, `p_in`,
`p_in_amplitude`, `ramp`, `z_in`, `lv` (0: prescribed `p_out`), `output`,
`output_every`, `tolerance`, `restart`, `max_iterations`, `scale`.

```sh
build-cuda/examples/Kokkos/Heart/KokkosAtriumFlowPETSc /path/to/Coarse.mesh \
  cycles=2 steps=850 output=atrium -ksp_type bcgsl -ksp_bcgsl_ell 4
build-cuda/examples/Kokkos/Heart/KokkosAtriumRheologyPETSc /path/to/Coarse.mesh \
  fluid=TVE cycles=2 steps=850 output=atrium_TVE
```

Keys of `KokkosAtriumRheologyPETSc`, besides those above: `fluid`, `stab`,
`order`, `inlet`, `t_end`, `n_in`, `n_init`, `devss`, `implicit`, `ptt`,
`eta_plasma`, `lambda_h`, `pc_lag`, `restart`, `checkpoint_every`. Its flow solver is BiCGStab(4) with
Jacobi unless PETSc options say otherwise.
