# Coronary artery models

Coronary perfusion depends on ventricular loading and on the response of the
arterial wall and downstream microvasculature. These effects are represented at
different levels: a rigid lumen isolates hemodynamics, a solid problem isolates
wall motion, and FSI transfers motion and traction between them.

For spatial velocity $u$, pressure $p$, and wall displacement $d$, the three-dimensional
balances are

$$
\rho_f(\partial_tu+u\cdot\nabla u)-\nabla\cdot\sigma_f=0,\quad
\nabla\cdot u=0,\quad
\sigma_f=-pI+2\mu D(u),\quad
\rho_s\ddot d-\mathrm{Div}_X P_s=0.
$$

The fluid strain rate is $D(u)=(\nabla u+\nabla u^T)/2$. Carreau–Yasuda viscosity
accounts for shear thinning. The wall laws are Neo-Hookean or Yeoh, according to the driver; their
compressibility, support conditions, and damping differ. The 0D ventricular
model supplies pressure or boundary-motion data, rather than a resolved ventricular
geometry or an active material field in the coronary wall.

| Model | Spatial formulation | Loading and coupling |
|---|---|---|
| [CFD](CFD/README.md) | Rigid-wall Taylor–Hood flow, Oseen or Newton. | 0D arterial inlet and intramyocardial outlet compartments. |
| [Solid](Solid/README.md) | Local mixed linear displacement/pressure. | 0D-driven normal-motion penalties and a clamped ring. |
| [Solid/PETSc_MPI](Solid/PETSc_MPI/README.md) | Distributed quadratic displacement, compressible elasticity. | Normal-motion penalties, outlet springs/dashpots, and a clamped ring. |
| [FSI](FSI/README.md) | Moving fluid geometry and deformable solid. | Staggered Robin–Robin variants or a monolithic correction system. |

The CFD class `CoupledLV0DCoronary3D` retains the flow, projection, and transport
problems together with ventricular and outlet states. The standalone solid and
FSI drivers construct their component problems directly. Shared constitutive and
flow-assembly terms support these drivers. The coagulation header is used by
the atrial and ectatic-tube transport models. The linked
READMEs distinguish the actual discretizations and boundary closures.
