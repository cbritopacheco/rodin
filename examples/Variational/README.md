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
both output directories with NumPy and Matplotlib installed.

```sh
python3 examples/Variational/visualize_spaces.py \
  --spaces /path/to/FiniteElementSpaces-output \
  --connectivity /path/to/Connectivity-output \
  --output /path/to/documentation-assets
```

The three `SpaceInterpolation-p0/p1/p2.png` panels show the interpolated
fields on a common color scale from 0 to 2. P0 uses one constant value per
original triangle through a flat `PolyCollection`; the renderer verifies that
all exported samples in each cell agree. P1 and quadratic H1 retain the
within-cell sample lattice. The CSV also supports sampled error comparisons;
the quadratic field reproduces the chosen polynomial up to roundoff.
CSV sampling refinement preserves the finite element mesh and DOF counts.
A convergence study requires a sequence of refined meshes.

The three `P1Basis-0/1/2.png` files are separate reference-basis panels. The
documentation places two panels in the first row and centers the third below.
