# Reconstruction

These examples reconstruct a circle in two dimensions or a sphere in three
dimensions. Cells are classified by their centroids; separating facets are
marked as the interface, and the outer box boundary is held fixed.

| Example | Displacement | Usage |
|---------|--------------|-------|
| `ReconstructionP1` | P1 | `SWIFT::Adapt` owns the solve and applies valid vertex displacement. |
| `ReconstructionP2` | P2 | `SWIFT::Problem` writes displacement; the application updates quadratic transformations. |
| `ReconstructionP3` | P3 | `SWIFT::Problem` writes displacement; the application updates cubic transformations. |

Build and run, for example:

```sh
cmake --build build --target ReconstructionP1 ReconstructionP2 ReconstructionP3 -j1
build/examples/Adaptation/SWIFT/ReconstructionP2 16 2
build/examples/Adaptation/SWIFT/ReconstructionP3 8 3
```

The optional arguments are grid points per axis and spatial dimension;
defaults are 16 and 2. Output is written to `swift/ReconstructionP*.xdmf`
and the corresponding HDF5 file in the working directory.

The examples use the production defaults and set the fixed reference spacing
through `parameters.model.h`. Additional controls can be set directly on that
parameter object before `setParameters()`. The printed report distinguishes
the geometric target from the sampled quality budget: a valid best-effort
result is not reported as a target hit.

Campaign parsing, response extraction and trajectory drivers remain in
`experiments/swift_calibration`, outside the examples.
