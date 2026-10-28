# Coronary FSI: PETSc_MPI

Monolithic coronary ALE fluid–structure interaction prototype using PETSc/MPI, a coronary fluid model, and a Neo-Hookean/Newmark wall.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery_FSI_PETSc_MPI -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/FSI/PETSc_MPI/CoronaryArtery_FSI_PETSc_MPI
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default conforming mesh is `resources/examples/Heart/CoronaryArtery/MergedFSI.mesh`. Use `-coronary_fsi_mesh` to override it. Cell labels distinguish fluid and solid.
