# Signed Poisson distance normalization

`SpaldingTucker.cpp` retains the original three-disk geometry and generates
its unit-square background. MMG cuts the nodal level set with an explicit
inside/outside material split. A new P1 space supplies the auxiliary signed
Poisson solve and normalized output.

Build `SpaldingTucker` and run it in a scratch directory. The executable
checks CG convergence, the assembled relative residual, finite nodal values,
nonempty phases and the zero interface values. `DistanceNormalization.csv`
contains the measured geometric/normalization errors; the XDMF/HDF5 output
contains `potential`, `normalized` and `reference` fields.

The geometric and algebraic errors have different meanings. In particular,
Spalding–Tucker normalization has a bias even for the exact radial Poisson
potential in a disk. The generated **Constructing a distance approximation
from a Poisson field** tutorial derives that value and reports the measured
errors for this geometry.

Reproduce the figures from computed output.

```sh
pvpython examples/Models/Distance/visualize_paraview.py \
  /path/to/results doc/Images
```

The figures use the Turbo colormap, transparent backgrounds, Manrope for
prose and TeX for mathematical labels. The signed fields share a symmetric
color range. Individual panels stack at phone widths in the documentation.
