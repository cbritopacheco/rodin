# Heart examples

These examples cover lumped ventricular dynamics, coronary and atrial flow,
arterial fluid–structure interaction, coagulation, and active wall contraction.

| Example | Model | Build requirements |
|---|---|---|
| [Atrium](Atrium/README.md) | Atrial flow and coagulation with a lumped ventricular outlet. | PETSc + MPI |
| [AtriumContraction](AtriumContraction/README.md) | Transient wedge-wall contraction driven by prescribed activation. | Local |
| [CCMLC2014_0D](CCMLC2014_0D/README.md) | Lumped left-ventricle and Windkessel dynamics. | Local |
| [CoronaryArtery](CoronaryArtery/README.md) | Coronary flow coupled to ventricular and intramyocardial models. | PETSc + MPI |
| [CoronaryArterySolid_PETSc_MPI](CoronaryArterySolid_PETSc_MPI/README.md) | Distributed coronary wall dynamics. | PETSc + MPI |
| [CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq](CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq/README.md) | Sequential staggered FSI with linear velocity and pressure. | PETSc |
| [CoronaryArtery_FSI_Explicit_PETSc_MPI](CoronaryArtery_FSI_Explicit_PETSc_MPI/README.md) | Staggered ALE FSI with a single-rank interface map. | PETSc + MPI |
| [CoronaryArtery_FSI_Explicit_PETSc_Seq](CoronaryArtery_FSI_Explicit_PETSc_Seq/README.md) | Sequential staggered FSI with quadratic velocity and linear pressure. | PETSc |
| [CoronaryArtery_FSI_PETSc_MPI](CoronaryArtery_FSI_PETSc_MPI/README.md) | Monolithic coronary ALE FSI prototype. | PETSc + MPI |
| [CoronarySolid](CoronarySolid/README.md) | Coronary wall dynamics and a residual/tangent diagnostic. | Local |
| [ProteinTest](ProteinTest/README.md) | Coagulation and washout in an ectatic tube. | PETSc |
