# Heart examples

These examples model lumped ventricular dynamics, coronary hemodynamics and
wall mechanics, atrial contraction, and coagulation transport.

| Example | Model | Requirements |
|---|---|---|
| [Atrium](Atrium/README.md) | Solid contraction and fluid flow on a common distributed geometry; FSI mode currently advances the components independently. | PETSc + MPI |
| [Atrium solid](Atrium/Solid/README.md) | Transient active contraction with passive and viscous tissue response. | PETSc + MPI |
| [CCMLC2014_0D](Ventricle/CCMLC2014_0D/README.md) | Lumped left-ventricle and Windkessel dynamics. | Local |
| [Coronary artery](CoronaryArtery/README.md) | Rigid-wall CFD, solid mechanics, and ALE FSI models. | Depends on formulation |
| [Ectatic tube thrombosis](EctaticTubeThrombosis/README.md) | Coagulation transport and washout in an ectatic tube. | PETSc |
