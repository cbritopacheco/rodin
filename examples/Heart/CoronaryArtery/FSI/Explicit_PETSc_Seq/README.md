# Coronary FSI with Taylor–Hood flow

This sequential staggered example uses continuous P2 velocity, P1 pressure,
and P1 wall displacement on separate conforming meshes. Its
[Robin–Robin transmission](../README.md#transmission-conditions)
places velocity impedance and lagged traction on both subproblems.

## Model and numerical procedure

Carreau–Yasuda fluid stress and homogeneous compressible Yeoh wall stress satisfy the
[parent component balances](../README.md#component-models). Static follower-pressure
increments initialize the wall at a configured fraction of the initial arterial
pressure (default $0.975$), followed by a short dynamic pressure ramp. The dynamic wall load
combines the same follower-pressure operator with a transferred traction remainder.
The Yeoh constants are $c_1=7\times10^4$, $c_2=2\times10^5$,
$c_3=1.9\times10^7$, and $\kappa=1.4\times10^6$, all in Pa.
An optional reference-gradient viscosity adds
$\eta_s\int\nabla_Xv_s:\nabla_Xw\,dX$; its default coefficient is zero.

The solid Robin terms are weighted by the current-to-reference interface area
ratio $J_\Gamma$. The fluid target is the BDF1 wall displacement rate by default;
`-coronary_mesh_consistent_interface_vel false` selects the Newmark endpoint rate.
The two rates differ at finite time step. A clamped band anchors the cut ends,
and the fluid Robin form includes the interface convective correction.

Conservative BDF1 ALE/Oseen integration advances the fluid, with lagged
Carreau–Yasuda viscosity, Temam correction, orthogonal-subscale VMS, and grad-div
terms. The P2/P1 pair supplies pressure stability without the additional
equal-order pressure-gradient term of the P1/P1 variant. The Newmark defaults
are $\beta=0.3025$ and $\gamma=0.6$. One coupling sweep is the default;
`-coronary_coupling_iterations` controls further unrelaxed corrections. A
single sweep is not an exactly converged interface solve.

The 0D heart supplies inlet pressure and normal-motion data on heart-contact
patches. Each intramyocardial outlet has a
compartment pressure $p_c$, external pressure $p_{\mathrm{im}}$, and flow $Q$:

$$
C(\dot p_c-\dot p_{\mathrm{im}})=Q-q_d,\qquad
q_d=\mathcal Q_d(\max(p_c-p_{\mathrm{im}},0)),\qquad
p_{\mathrm{out}}=p_c+\Delta p_p(Q).
$$

Non-Newtonian branch-flow laws define $\mathcal Q_d$ and the proximal pressure
drop $\Delta p_p$. The full outlet pressure is lagged. A normal increment-damping impedance adds

$$
Z_{\mathrm{out}}\int_{\Gamma_{\mathrm{out}}}
[(u^{n+1}-u^n)\cdot n](v\cdot n)\,dS.
$$

This term vanishes at a steady velocity. It is distinct from the direct
proximal resistance term used by the P1/P1 R–$\mu$ closure.

## Geometry, units, and output

Wall label 2 identifies the interface before the cut-end band is promoted to
label 99. Opening caps stay fixed during weighted harmonic ALE extension; the extension
is solved in the P2 velocity space and coordinates are updated at mesh vertices.
Solid cut ends
and the ring band are clamped. Mesh coordinates are converted to metres with
factor $10^{-3}$ by default; all solved quantities use SI units.

XDMF records moving fluid and reference-solid fields. CSV records heart and
outlet states, fluxes, wall motion, and interface/coupling observables. Spatial
interface transfer and a finite correction count remain numerical approximations.

Requires `RODIN_USE_PETSC=ON`; MPI is not required.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery_FSI_Explicit_PETSc_Seq -j4
cd build
./examples/Heart/CoronaryArtery/FSI/Explicit_PETSc_Seq/CoronaryArtery_FSI_Explicit_PETSc_Seq
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default meshes are `resources/examples/Heart/CoronaryArtery/FSIFluid.mesh` and `FSISolid.mesh` in the same directory.
