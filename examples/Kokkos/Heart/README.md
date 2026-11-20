# Atrium on the device

`KokkosAtriumFlow` is the first step of the port of
`examples/Heart/Atrium/Viscoelastic` (branch `model/Atrium`) to the Kokkos
module. It advances a Newtonian fluid in the fixed atrial volume over cardiac
cycles, with the assembly and the linear solve on the device and one process.

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

## Not ported

Carreau-Yasuda, the log-conformation fluids (Owens with equilibrium or
transported structure, gsPTT), the Eulerian age, wall shear indices, volume
diagnostics and probes, and MPI.

What each needs from the Kokkos module, in order of effort:

1. **Carreau-Yasuda, age.** Nothing new. A device callable can read a grid
   function view and its Jacobian, as `tau` does here; the age is one more
   scalar problem.
2. **Orthogonal subscales.** $L^2$ projections of $(\nabla u)u$, $\nabla p$
   and $\nabla\cdot u$ are device problems with a mass matrix; the dynamic
   subscale needs a cellwise field (no P0 space on the device yet) or a P1
   surrogate.
3. **Log-conformation.** $\exp$, its derivative by Daleckii-Krein and the
   3x3 symmetric eigenproblem inside kernels, two more vector unknowns for
   $\psi$, a Newton loop, and storage at the quadrature points for the
   subscale of $\psi$. None exists on the device.
4. **Several devices.** The PETSc-Kokkos path of `KokkosMPIOseenScaling`.

## Status

Compiled with Kokkos 5.2.2, Serial backend, and run on a 3983-cell test mesh
with the labels above for one cycle of 85 steps: the flow rates through the
veins and the mitral plane agree to within 0.5 % at every step, the valve
closes and reopens with the ventricular pressure, and areas and volume equal
those of the CPU driver. It was not compiled with CUDA and not run on a GPU
or on `Coarse.mesh`. The GMRES path has no preconditioner; on a mesh of that
size the sparse direct device solver (`RODIN_USE_CUDSS`) is the intended one.

## Run

```sh
cmake --build build-cuda --target KokkosAtriumFlow
build-cuda/examples/Kokkos/Heart/KokkosAtriumFlow /path/to/Coarse.mesh \
  cycles=2 steps=850 output=atrium
```

Keys: `cycles`, `steps` (per cycle), `period`, `rho`, `mu`, `p_in`,
`p_in_amplitude`, `ramp`, `z_in`, `lv` (0: prescribed `p_out`), `output`,
`output_every`, `tolerance`, `restart`, `scale`.
