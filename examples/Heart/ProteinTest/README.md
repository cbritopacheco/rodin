# ProteinTest

Sequential coagulation and washout testbed in a coronary-like tube with an ectatic segment. The outlet is traction-free. The mesh generators are included in this folder.

Requires `RODIN_USE_PETSC=ON`; MPI is not required.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target ProteinTest -j4
cd build
./examples/Heart/ProteinTest/ProteinTest
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/ProteinTest.mesh`. Generate it from the repository root with:

```sh
python3 examples/Heart/ProteinTest/make_test_mesh.py resources/examples/Heart/ProteinTest.mesh
```

The driver source [`ProteinTest.cpp`](ProteinTest.cpp) contains the model parameters,
solver settings, and output configuration.
