# Finite element comparison

Build the `FiniteElementSpaces` example target and run it in an output directory.
It constructs an eight-by-eight vertex grid, prints the global DOF counts for
several spaces, and interpolates the known function $f(x,y)=x^2+y^2$. The program compares field representations on a fixed mesh.

`FiniteElementSpaces.csv` contains physical-point evaluations, including the
cell and reference-sampling indices. `ReferenceP1.csv` contains the triangle
basis evaluations used by the Ciarlet guide. Coefficient storage order is not
used as a proxy for physical samples.

The `Connectivity` geometry example writes the actual vertex and edge indices
of its two-triangle mesh. The documentation figures can be reproduced from
both output directories with NumPy and Matplotlib installed:

```sh
python3 examples/Variational/visualize_spaces.py \
  --spaces /path/to/FiniteElementSpaces-output \
  --connectivity /path/to/Connectivity-output \
  --output /path/to/documentation-assets
```

Each interpolation-error panel has its own labeled color scale. The quadratic
field reproduces the chosen polynomial up to roundoff; this comparison is not
a convergence study for arbitrary PDE solutions. CSV sampling refinement preserves the finite element mesh and DOF counts.
