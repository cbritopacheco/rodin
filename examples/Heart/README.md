# Heart examples

Cardiac mechanics, blood flow, and circulation involve different unknowns and
geometries. Lumped models retain chamber and vascular states, solid models solve
wall displacement, and fluid models solve velocity and pressure. Fluid–structure
interaction adds compatibility of motion and balance of traction at their shared
interface. Transport models use the computed flow to advance concentrations.

The corresponding balances have the forms

$$
\mathcal R_{0D}(z,\dot z,t)=0,\qquad
\rho_s\ddot d-\mathrm{Div}_X P_s=0,\qquad
\rho_f(\partial_tu+u\cdot\nabla u)-\nabla\cdot\sigma_f=0,\quad
\nabla\cdot u=0.
$$

Here $z$ collects lumped states, $d$ is reference-wall displacement, and
$u,p$ are spatial fluid velocity and pressure. Constitutive laws, loading, and
space/time discretization are specified by each example; these balances do not
imply that the examples share the same material or coupling algorithm.

| Example | Unknowns and driving data | Numerical model |
|---|---|---|
| [Atrium](Atrium/README.md) | Active wall, cavity flow, or both; contraction-driven FSI has open zero-pressure boundaries and no 0D driver. | MPI/Scotch, linear wedges and stabilized linear tetrahedral flow; Robin–Robin FSI. |
| [Ventricle](Ventricle/CCMLC2014_0D/README.md) | Chamber radius, pressures, and active variables. | Lumped CCMLC2014 model and two-compliance arterial circulation. |
| [Coronary artery](CoronaryArtery/README.md) | Rigid-wall flow, wall mechanics, or ALE FSI, with lumped heart/circulation loading. | Distinct local and PETSc/MPI formulations; coupling is specified per example. |
| [Ectatic tube](Thrombosis/EctaticTube/README.md) | Flow, coagulation concentrations, and a washout tracer. | Fixed-wall stabilized flow and one-way reaction–transport. |

Each model README defines its domains, constitutive assumptions, boundary data,
numerical procedure, inputs, outputs, and execution commands. XDMF fields describe
spatial states; CSV files describe scalar histories. A displayed deformation or
velocity scale does not change the solved physical field.
