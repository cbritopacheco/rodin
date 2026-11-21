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
scheme; in short, per time step, with the extrapolations
$u^{n+1/2}=(3u^n-u^{n-1})/2$ and $u^*=2u^n-u^{n-1}$:

1. velocity gradient of $u^{n+1/2}$ at the vertices (mean of the cell
   gradients, trace free, and $(\partial_n u)\,n^T$ on the wall);
2. half a step at every vertex: the aggregate size by the exact solution of
   its Riccati equation, the conformation $c$ by the exact solution of its
   equation with the gradient frozen; $\psi=\log c$ is the stored and
   transported unknown;
3. $\psi$ (six components) and $n$ advected over the step by one scalar
   Crank-Nicolson SUPG matrix, with the inflow values imposed weakly on the
   faces where the flow enters;
4. the second half step of 2 (Strang splitting);
5. momentum and continuity by BDF2, with convection by $u^*$, $\tau_p=G(c-I)$
   of $t^{n+1}$ explicit and balanced by the viscosity $G\lambda$ on both sides,
   stabilised by the algebraic subgrid scale method with quasi-static
   subscales (`stab=asgs`: SUPG and PSPG with the residual of the momentum
   equation without its time derivative) or by the two terms of
   `KokkosAtriumFlow` (`stab=split`).

Steps 1, 2 and 4 are Kokkos kernels of the example; 3 and 5 are PETSc
matrices of the Kokkos type. `order=1` gives the first order version
(backward Euler throughout, one local step).

**Stabilisation parameters.** Neither $\tau$ has the $(2/\Delta t)^2$ term of
`KokkosAtriumFlow`. With it $\tau\to\Delta t/2$, the stabilised problem
changes with the time step, and on the test mesh neither order converged in
time (Bochev, Gunzburger & Shadid 2004). Without it, the time derivative
cannot stay in the residual of the subscales, since its term $\tau/\Delta t$
is destabilising for $\Delta t<\tau$ (John & Novo 2011; the run blew up). So
the subscales are quasi-static: consistency error of order $\tau$, spatial.
`tau_dt=1` restores the former version. Dynamic subscales, as in the CPU
driver, would remove both compromises; they need a subscale stored by cell.

**Time convergence**, test mesh, `lv=0`, $t=0.1$ s, against $\Delta t=0.125$ ms,
relative $L^2$ errors at $\Delta t$ = 2, 1, 0.5, 0.25 ms:

| | velocity | rate | $\psi$ | rate |
|---|---|---|---|---|
| N, order 2 | 3.9e-3 ... 4.7e-5 | 2.0 | | |
| VE, order 2 | 1.4e-3 ... 8.3e-5 | 1.2-1.5 | 3.1e-3 ... 9.5e-5 | 1.6-1.8 |
| TVE, order 2 | 1.2e-3 ... 6.1e-5 | 1.3-1.6 | 2.3e-3 ... 1.5e-4 | 1.2-1.5 |
| TVE, order 1 | 1.0e-2 ... 1.3e-3 | 1.0 | 2.4e-2 ... 3.0e-3 | 1.0 |

The elastic fluids miss rate 2 for a reason not identified, possibly the
inflow condition, which switches with the direction of the flow on parts of
the openings. With the ventricle the atrial pressure it reads is that of
$t^n$ (first order at the coupling); extrapolating it (`lv_extrapolate=1`)
made the valve chatter and $p_v$ diverge.

**Time-scheme dissipation.** The CSV column `timeDissipation` is the kinetic
energy removed by the time scheme per unit time,
$\rho/(4\Delta t)\,\lVert u^{n+1}-2u^n+u^{n-1}\rVert^2$ for BDF2. Against the
solvent and polymer dissipation, on that test: order 1, 3.3, 1.6 and 0.8 at
2, 1 and 0.5 ms (the scheme removed more energy than the fluid); order 2,
3e-3, 4e-4 and 5e-5. Over a full cycle at 1 ms with the ventricle, order 2:
about 1e-5 of the solvent and polymer dissipation in the E-wave.

**Against the CPU driver**, TVE, test mesh, one cycle: $p_v$ at the end 202 Pa
(CPU, BDF2 at 10 ms with orthogonal subscales: 218 Pa); in the E-wave,
polymer power 0.32 mW and elastic energy 1.4 uJ (CPU: 1.1-1.3 mW and
7.2-8.0 uJ). The stress is four to six times smaller on this 2.2 mm mesh:
nodal P1 stress from a recovered gradient, and SUPG diffusion of $\psi$, where
the CPU driver has orthogonal subscales. Compare on `Coarse.mesh` before
using the stress magnitudes. Checked besides: the vertex update against the
exact Oldroyd-B stress in steady shear, extension and start-up of shear;
`order=1 tau_dt=1 stab=split fluid=N` against `KokkosAtriumFlowPETSc` (equal
digits). Not run with CUDA, nor on `Coarse.mesh`: `stab=asgs` adds velocity
pressure matrices to the assembly, so its memory there is not known;
`stab=split` has the footprint of `KokkosAtriumFlowPETSc`.

Owens' dumbbells are Hookean and their stress is unbounded at a stagnation
point with strong extension (`tr(c)_max` in the log); `ptt=0.01`, say, bounds
it with a linear Phan-Thien-Tanner factor that is not part of Owens' model.

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
reproduces the native run to six digits at the end of the cycle. A Schur
field split (FGMRES, `selfp`, GAMG on the velocity, Jacobi on the pressure)
did not converge on that mesh with the options tried; it needs tuning.

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
`order`, `theta`, `tau_dt`, `lv_extrapolate`, `t_end`, `n_in`, `n_init`,
`bsd`, `ptt`, `eta_plasma`, `lambda_h`. Its flow solver is
BiCGStab(4) with Jacobi unless PETSc options say otherwise; the transport
solver takes options with the prefix `transport_`.
