# CoronaryArtery_FSI_Explicit_PETSc_Seq

Sequential staggered coronary ALE fluid–structure interaction using a quadratic-velocity/linear-pressure fluid discretization, a hyperelastic wall, and lumped inlet/outlet models.

Requires `RODIN_USE_PETSC=ON`; MPI is not required.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery_FSI_Explicit_PETSc_Seq -j4
cd build
./examples/Heart/CoronaryArtery_FSI_Explicit_PETSc_Seq/CoronaryArtery_FSI_Explicit_PETSc_Seq
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default meshes are `resources/examples/Heart/CoronaryArtery/FSIFluid.mesh` and `FSISolid.mesh` in the same directory. Results are written under the configured results directory.

The driver source [`CoronaryArtery_FSI_Explicit_PETSc_Seq.cpp`](CoronaryArtery_FSI_Explicit_PETSc_Seq.cpp) contains the model parameters,
solver settings, and output configuration.
