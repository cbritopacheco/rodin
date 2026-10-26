# CCMLC2014_0D

Standalone lumped left-ventricle model coupled to a four-element Windkessel. Writes the pressure, volume, and flow history without a spatial mesh.

No PETSc or MPI requirement.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target CCMLC2014_0D -j4
cd build
./examples/Heart/CCMLC2014_0D/CCMLC2014_0D
```

Output files are written relative to the working directory.

## Input and configuration

No input mesh is required.
