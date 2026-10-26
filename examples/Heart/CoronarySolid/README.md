# CoronarySolid

Coronary wall dynamics with a quasi-compressible Neo-Hookean material, Newmark time integration, damping, and a residual/tangent finite-difference diagnostic.

No PETSc or MPI requirement.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CoronarySolid -j4
cd build
./examples/Heart/CoronarySolid/CoronarySolid
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/CoronaryArtery/Solid.mesh`.
