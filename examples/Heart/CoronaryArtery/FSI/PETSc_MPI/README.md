# Monolithic coronary ALE FSI prototype

This prototype solves fluid, wall, and mesh-motion corrections in a common
PETSc/SNES system. A single conforming mesh distinguishes fluid cells (label 1)
and solid cells (label 2). The same P1 displacement field represents physical
wall motion in the solid and ALE motion in the lumen; fluid velocity/pressure
use continuous P2/P1 spaces.

## Coupled model and correction system

The [fluid and Neo-Hookean balances](../README.md#component-models) are assembled
with conservative ALE time terms and Newmark solid inertia. Unknowns in the
linearized system are fluid corrections $\delta u_f,\delta p$ and displacement-
increment correction $\delta\eta$, where $d_s^{n+1}=d_s^n+\eta$.
With Newmark derivative $c_v=\gamma/(\beta\Delta t)$, the kinematic correction is

$$
\delta u_f-c_v\delta\eta=v_s-u_f\quad\text{on }\Gamma.
$$

This couples the new fluid velocity to wall velocity inside the correction
system. Interface stress terms transmit loading through the assembled residual;
there is no staggered Robin–Robin sweep. Harmonic mesh-motion equations extend
displacement into the fluid.

`newton` flow mode includes current ALE transport/geometry derivatives and the
Carreau–Yasuda viscosity derivative. `oseen` freezes fluid transport and viscosity
in the tangent. SNES remains necessary in both modes because the wall law is
nonlinear. Fluid time integration is BDF1 and the default Newmark parameters
are $\beta=1/4$, $\gamma=1/2$.

## Boundary data and limitations

Boundary/interface labels are 2 for the shared wall, 3 for the inlet, 4–9 for
outlets, and 100 for the solid support ring. The 0D heart supplies pressure data.
An outlet compartment is advanced by the driver relation

$$
C\frac{p_c^{n+1}-p_c^n}{\Delta t}
=Q-\frac{p_c^{n+1}-p_v^{n+1}}{R_d},\qquad
p_{\mathrm{out}}=p_c+R_pQ,
$$

where $p_v$ is the ventricular pressure used as downstream reference in this
update. This closure is distinct from the R–$\mu$ model in rigid-wall CFD and
the staggered P1/P1 example. The outer 0D update is not part of the monolithic
3D Newton unknown vector.

The driver is a prototype: a common correction system does not establish spatial
convergence, robustness on arbitrary partitions, or a fully implicit global
circulation coupling. Mesh Jacobians and nonlinear convergence remain requirements.
Source coordinates use the default metre conversion $10^{-3}$; solved quantities
use SI units. XDMF contains fluid and wall states; CSV contains lumped states,
flows, and mechanical/solver observables.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery_FSI_PETSc_MPI -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/FSI/PETSc_MPI/CoronaryArtery_FSI_PETSc_MPI
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default conforming mesh is `resources/examples/Heart/CoronaryArtery/MergedFSI.mesh`. Use `-coronary_fsi_mesh` to override it. Cell labels distinguish fluid and solid.
