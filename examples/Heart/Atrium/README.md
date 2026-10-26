# Atrium

Three-dimensional atrium flow with coagulation kinetics, pressure inlets, and an outlet coupled to the lumped left-ventricle model through a mitral diode.

Requires `RODIN_USE_PETSC=ON` and `RODIN_USE_MPI=ON`.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target Atrium -j4
cd build
mpirun -n 1 ./examples/Heart/Atrium/Atrium
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/Atrium/Coarse.mesh`. Override it with `-atrium_mesh`; `-atrium_dt`, `-atrium_flow_cycles`, and `-atrium_species_cycles` control the run.
