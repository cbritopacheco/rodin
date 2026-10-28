# Coronary FSI: Explicit_PETSc_MPI

Staggered coronary ALE fluid–structure interaction with Robin coupling, a hyperelastic wall, and lumped inlet/outlet models. The interface map is not distributed, so this example must run on one MPI rank.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery_FSI_Explicit_PETSc_MPI -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/FSI/Explicit_PETSc_MPI/CoronaryArtery_FSI_Explicit_PETSc_MPI
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default meshes are `resources/examples/Heart/CoronaryArtery/FSIFluid.mesh` and `FSISolid.mesh` in the same directory.
