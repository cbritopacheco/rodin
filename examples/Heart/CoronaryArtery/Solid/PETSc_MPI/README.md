# Distributed coronary wall dynamics

This example advances coronary wall displacement with PETSc/SNES on a reconciled
MPI mesh. A lumped heart model supplies a normal boundary-motion datum. The wall
is a compressible displacement-only solid; fluid and solid-pressure unknowns
are absent.

## Constitutive model and supports

With $F=I+\nabla_Xd$, $J=\det F$, and reference density $\rho_0$, the equations are

$$
\rho_0\ddot d-\mathrm{Div}_XP=0,\qquad
P=\mu_sF+(\lambda_s\log J-\mu_s)F^{-T}.
$$

`NeoHookean` uses the [stored energy](../../FSI/README.md#component-models) with
$E=5\times10^5\,\mathrm{Pa}$, $\nu=0.3$,
$\mu_s=E/[2(1+\nu)]$, and
$\lambda_s=E\nu/[(1+\nu)(1-2\nu)]$.
Continuous H1 order-two Lagrange functions approximate displacement;
the model does not use the mixed P1/P1 formulation of the local driver.

A ring is clamped. Heart-contact surfaces receive the normal-motion penalty
$k\int(d\cdot N-g(t))(w\cdot N)\,dS$.
Outlet supports add spring and dashpot terms

$$
\int_{\Gamma_{\mathrm{out}}}(a\,d+b\,v)\cdot w\,dS.
$$

The default coefficients are $k=a=10^5\,\mathrm{Pa/m}$ and
$b=5\times10^4\,\mathrm{Pa\,s/m}$. Dissipation here comes from these boundary
supports and the selected time integrator, rather than an atrial active or
finite-strain tissue-viscosity law.

## Time integration and architecture

The Newmark parameters are $\beta=0.3025$ and $\gamma=0.6$, derived from the
configured value $\alpha=-0.1$. The assembled forces are evaluated at the new
state; an HHT weighted-force balance is not added. The time increment is
$10^{-3}\,\mathrm{s}$. The retained nonlinear problem combines inertia, internal
virtual work, and support loads, with a consistent Newton correction.

The source mesh is loaded on rank zero, partitioned, sharded, and reconciled.
Ghost-aware grid-function updates preserve shared values. The 0D driver is
advanced identically on each rank to supply the same motion datum. Coordinates
are converted to metres; mechanical quantities use SI units. MPI XDMF contains
displacement, velocity, and acceleration; rank zero writes the lumped-state CSV.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArterySolid_PETSc_MPI -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/Solid/PETSc_MPI/CoronaryArterySolid_PETSc_MPI
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/CoronaryArtery/Solid.mesh`.
