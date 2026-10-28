# Coronary rigid-wall CFD

Three-dimensional coronary flow is coupled to a lumped left-ventricle model
and intramyocardial outlet beds through pressure and flow boundary data.
The Navier–Stokes problem supports Oseen and Newton linearizations. The arterial
wall is fixed with a no-slip velocity condition; no structural displacement
or ALE mesh motion is solved.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronaryArtery -j4
cd build
mpirun -n 1 ./examples/Heart/CoronaryArtery/CFD/CoronaryArtery
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/CoronaryArtery/Output3D.mesh`. See `CoronaryArtery.cpp` for the PETSc options, including `-coronary_dt`, `-coronary_nsteps`, and `-coronary_output_prefix`.
