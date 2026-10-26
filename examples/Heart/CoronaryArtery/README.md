# CoronaryArtery

Coronary flow coupled to the lumped left-ventricle model and intramyocardial outlet beds. The driver supports Oseen and Newton flow solves. Shared coronary model headers and postprocessing scripts also live in this folder.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/CoronaryArtery
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/CoronaryArtery/Output3D.mesh`. See `CoronaryArtery.cpp` for the PETSc options, including `-coronary_dt`, `-coronary_nsteps`, and `-coronary_output_prefix`.

The driver source [`CoronaryArtery.cpp`](CoronaryArtery.cpp) contains the model parameters,
solver settings, and output configuration.
