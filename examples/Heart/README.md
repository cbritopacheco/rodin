# Heart examples

Each example has its own source folder, CMake target, and README with build
and run instructions. Executables are generated under
`build/examples/Heart/<example>/<example>`. Target names are unchanged.

| Example | Build requirements |
|---|---|
| [Atrium](Atrium/README.md) | PETSc + MPI |
| [AtriumContraction](AtriumContraction/README.md) | Local |
| [CCMLC2014_0D](CCMLC2014_0D/README.md) | Local |
| [CoronaryArtery](CoronaryArtery/README.md) | PETSc + MPI |
| [CoronaryArterySolid_PETSc_MPI](CoronaryArterySolid_PETSc_MPI/README.md) | PETSc + MPI |
| [CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq](CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq/README.md) | PETSc |
| [CoronaryArtery_FSI_Explicit_PETSc_MPI](CoronaryArtery_FSI_Explicit_PETSc_MPI/README.md) | PETSc + MPI |
| [CoronaryArtery_FSI_Explicit_PETSc_Seq](CoronaryArtery_FSI_Explicit_PETSc_Seq/README.md) | PETSc |
| [CoronaryArtery_FSI_PETSc_MPI](CoronaryArtery_FSI_PETSc_MPI/README.md) | PETSc + MPI |
| [CoronarySolid](CoronarySolid/README.md) | Local |
| [ProteinTest](ProteinTest/README.md) | PETSc |

Shared coronary helpers remain in `CoronaryArtery/`; the atrium, protein,
and sequential FSI examples include them through their private include path.
Meshes remain under `resources/examples/Heart/` and are stored in Git LFS.
