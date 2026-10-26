# CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq

Sequential staggered coronary ALE fluid–structure interaction using equal-order linear velocity and pressure.

Requires `RODIN_USE_PETSC=ON`; MPI is not required.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq -j4
cd build
./examples/Heart/CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq/CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default meshes are `resources/examples/Heart/CoronaryArtery/Stenosis80.mesh` and `Prismatic80.mesh` in the same directory.
