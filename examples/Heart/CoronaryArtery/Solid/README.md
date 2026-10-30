# Coronary wall mechanics with mixed displacement and pressure

The local coronary wall example isolates deformation driven by a lumped heart
motion signal. The 0D model advances a radial displacement $y(t)$; its converted
normal-motion datum $g(t)$ enters surface penalties on the 3D wall. No fluid
velocity or ALE mesh-motion problem is solved.

## Mixed finite-strain model

Let $F=I+\nabla_Xd$ and $J=\det F$. Continuous P1 displacement and P1 solid
pressure $p_s$ approximate an incompressible Neo-Hookean wall. The implemented
first Piola stress is

$$
P=\mu_s(F-F^{-T})+p_sJF^{-T},\qquad
\rho_0\ddot d-\mathrm{Div}_XP=0.
$$

Here $p_s$ is a solid constraint multiplier, not the ventricular cavity pressure.
The stabilized incompressibility row is

$$
\int_{\Omega_{s,0}}(J-1)q\,dX
-\tau\int_{\Omega_{s,0}}\nabla_Xp_s\cdot\nabla_Xq\,dX=0.
$$

The code uses the fixed pressure-gradient coefficient $\tau=10^{-6}$.
The material shear modulus is obtained from $E=5\times10^5\,\mathrm{Pa}$ and
$\nu=0.3$ by $\mu_s=E/[2(1+\nu)]$; the separately computed Lamé coefficient is
not used in this mixed law. The assembled bulk model contains elastic stress
and inertia, without a bulk viscous term.

Ring label 100 is clamped. On labels 101 and 105, the weak load adds
$k\int(d\cdot N-g)(w\cdot N)\,dS$, with $k=10^4$ in SI units. This is prescribed
normal-motion loading, rather than a fluid traction transmission condition.

## Numerical procedure and output

Newmark average acceleration uses $\beta=1/4$, $\gamma=1/2$, and
$\Delta t=10^{-3}\,\mathrm{s}$. A block Newton method solves displacement and
pressure together using `InternalVirtualWorkResidual`, its displacement/pressure
tangent blocks, and `SparseLU`. Mesh coordinates are converted from millimetres
to metres. XDMF writes displacement, velocity, acceleration, and solid pressure;
CSV records the 0D state. The auxiliary homogeneous Laplace field is separate
from the mechanical unknowns.

The [PETSc/MPI driver](PETSc_MPI/README.md) uses a compressible displacement-only
law and quadratic displacement, so it represents a different discrete model.

No PETSc or MPI dependency is required for this local example.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronarySolid -j4
cd build
./examples/Heart/CoronaryArtery/Solid/CoronarySolid
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/CoronaryArtery/Solid.mesh`.

The [PETSc/MPI formulation](PETSc_MPI/README.md) provides a distributed solid
solve with Newmark integration and lumped heart forcing.
