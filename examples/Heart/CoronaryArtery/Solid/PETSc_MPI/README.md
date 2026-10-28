# CoronaryArterySolid_PETSc_MPI

Distributed PETSc/MPI variant of the coronary solid dynamics example, with Newmark integration and the lumped heart model.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArterySolid_PETSc_MPI -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/Solid/PETSc_MPI/CoronaryArterySolid_PETSc_MPI
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/CoronaryArtery/Solid.mesh`.
