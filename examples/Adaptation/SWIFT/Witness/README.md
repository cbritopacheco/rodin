# Reference Witness Sets

`SWIFT_Witness` generates unrestricted Euclidean covering sets on Rodin reference
geometries. C++ performs continuous Voronoi evaluation, enclosing-ball refinement
and JSON serialization. Python draws the saved results. Neither stage changes
SWIFT's production quality witnesses or a running calibration campaign.

## Build

The executable is disabled by default. CGAL is discovered only when it is enabled:

```sh
cmake -S . -B build -DRODIN_BUILD_EXAMPLES=ON -DRODIN_BUILD_SWIFT_WITNESS=ON
cmake --build build --target SWIFT_Witness -j1
```

An external CGAL installation and Boost.JSON are required. For a nonstandard CGAL
installation, add `-DCGAL_DIR=/path/to/lib/cmake/CGAL`, or add its installation
prefix to `CMAKE_PREFIX_PATH`. CMake links `CGAL::CGAL` privately to
`SWIFT_Witness`; Rodin's library and reconstruction executables do not acquire
this dependency. No enclosing-ball library sources are vendored into Rodin.

CGAL's Bounding Volumes package is GPL-licensed, with commercial licensing also
available. Building an optional executable does not remove its distribution
obligations. Distribution of the combined Witness executable must comply with
the applicable CGAL license. See [CGAL licensing](https://www.cgal.org/license.html)
and the [package reference](https://doc.cgal.org/latest/Bounding_volumes/group__PkgBoundingVolumesRef.html).

The visualizers require only NumPy:

```sh
python3 -m pip install -r examples/Adaptation/SWIFT/Witness/requirements.txt
```

## Generate

Run from the repository root, replacing `build` with the configured build path:

```sh
# Default triangle; all witness locations are free.
build/examples/Adaptation/SWIFT/Witness/SWIFT_Witness 16

# One geometry and an explicit output path.
build/examples/Adaptation/SWIFT/Witness/SWIFT_Witness 16 \
  --geometry tetrahedron --output /tmp/witness/tetrahedron-16.json

# All eight geometries, executed serially.
build/examples/Adaptation/SWIFT/Witness/SWIFT_Witness 16 \
  --geometry all --output /tmp/witness/all.json
```

No vertices, edges or faces are prescribed. Every positive count is accepted,
except that the point geometry admits only one distinct location. `--geometry all`
uses one witness for that geometry. Both `--name value` and `--name=value` are
supported.

| Geometry | Reference domain | Vertices |
|---|---|---:|
| `point` | \(\{0\}\) | 1 |
| `segment` | \([0,1]\) | 2 |
| `triangle` | \(x,y\ge0,\ x+y\le1\) | 3 |
| `quadrilateral` | \([0,1]^2\) | 4 |
| `tetrahedron` | \(x,y,z\ge0,\ x+y+z\le1\) | 4 |
| `pyramid` | \(0\le z\le1,\ 0\le x,y\le1-z\) | 5 |
| `hexahedron` | \([0,1]^3\) | 8 |
| `wedge` | \(x,y\ge0,\ x+y\le1,\ 0\le z\le1\) | 6 |

Reference coordinates and half-spaces come from `Geometry::Polytope::Traits`.

| Setting | Meaning | Default |
|---|---|---|
| Positional count or `--n` | Number of distinct witnesses | `6` |
| `--geometry` | Reference geometry or `all` | `triangle` |
| `--output` | JSON path | `witness-<geometry>-<n>.json` |
| `--iterations` | Voronoi refinement sweep cap | `100` |
| `--step-tolerance` | Maximum accepted site-motion stopping threshold | `1e-6` |
| `--max-evaluations` | Candidate coverage evaluation cap | `0`: no extra cap |
| `--resolution` | Independent check lattice subdivision; the check uses twice this value | `12` |
| `--seed` | Enclosing-ball input-order seed | `13` |
| `--evaluate` | Evaluate supplied JSON coordinates without optimization | None |

The obsolete differential-evolution flags `--starts`, `--population`,
`--polish-iterations` and `--objective` are not supported. The tool is serial.

## Numerical Method

For a reference polytope \(K\), the covering objective is

\[
\min_{X\subset K,\ |X|=n}\eta_K(X),
\qquad
\eta_K(X)=\max_{x\in K}\min_{\xi\in X}\|x-\xi\|_2.
\]

Initialization uses the smallest uniform reference lattice containing at least
the requested number of points. Simplices use barycentric lattices; boxes use
tensor-product lattices; wedges use triangle lattices times a segment; pyramids
use shrinking square layers. If the count does not match a complete lattice,
deterministic farthest-point selection selects the requested number, starting
nearest the vertex barycenter. A single site starts at that barycenter. The
construction is recorded in JSON and is independent of `--resolution`.

Each site's Voronoi cell is clipped against the reference polytope and the other
sites' bisector half-spaces. Distance to a site is convex on its cell, so its
maximum is attained at a cell vertex. Continuous coverage is therefore evaluated
from these vertices in floating-point arithmetic, not from a sampling grid.

Each sweep freezes those cells and replaces every site by its cell's smallest
enclosing-ball center:

\[
\xi_i^{\mathrm{new}}
=\operatorname*{arg\,min}_{\xi}
\max_{a\in\operatorname{vertices}(V_i)}\|a-\xi\|_2.
\]

The center lies in the convex hull of its support points, hence within its cell
and the reference domain. In exact arithmetic this update cannot increase the
global covering radius. The Voronoi diagram is rebuilt after every sweep;
nonfinite or radius-increasing candidates are rejected.

The adapter translates and normalizes the cell vertices before calling
`CGAL::Min_sphere_of_spheres_d`, representing them as zero-radius spheres.
Containment and the active-support convex-hull condition are checked independently.
If floating-point degeneracy prevents those checks, a finite support search over
at most dimension plus one vertices is attempted. That fallback is combinatorial,
not a linear-time guarantee. Failed checks terminate with an error rather than
accepting an unchecked center.

## Termination and Limits

Refinement stops at the site-motion tolerance, sweep cap, evaluation cap or a
rejected candidate. Initial and final diagnostic evaluations do not count against
the candidate evaluation budget. JSON records the backend, sweep history,
termination reason, fallback count, local-solve cost and coverage-evaluation cost.
The interval solution is analytic.

Apart from the interval and point cases, global optimality is not certified.
Motion convergence only indicates a local fixed point; floating-point Voronoi
evaluation is not an interval certificate. An independent lattice check supplies
a lower bound on the returned radius, not an upper bound on numerical error.

This is geometric coverage, not constraint observability or quadrature design.
It uses neither \(Q\) nor \(j\), and assigns no hinge weights. Reference-space
coverage need not become uniform physical-space coverage under a curved or
anisotropic mapping.

## Visualize

Use `plot.py` with a single-geometry JSON object:

```sh
python3 examples/Adaptation/SWIFT/Witness/plot.py \
  /tmp/witness/tetrahedron-16.json --output /tmp/witness/tetrahedron-16.svg
```

The plot compares initial and final sites side by side, with clipped cells and
largest uncovered disks in two dimensions. Three-dimensional geometries retain
an oblique XYZ view together with
XY, YZ and XZ projections. Dashed circles in projected views are projections of
uncovered three-dimensional balls, not empty disks of the projected point set.
The output is vector SVG. Python does not optimize or recompute coverage.

## Verification

An independent SciPy/Qhull oracle is used only by the tests:

```sh
python3 -m pip install -r examples/Adaptation/SWIFT/Witness/tests/requirements.txt
SWIFT_WITNESS_BINARY="$(pwd)/build/examples/Adaptation/SWIFT/Witness/SWIFT_Witness" \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 -m unittest discover -s examples/Adaptation/SWIFT/Witness/tests -v
```

Checks cover all eight geometries, continuous radii against the independent
oracle, analytic interval solutions, initialization, determinism, monotone
refinement, evaluation budgets, invalid inputs and figure generation.
[The historical comparison](COMPARISON.md) retains the original enclosing-ball,
smooth-Newton and differential-evolution measurements; it is not a CGAL timing
benchmark.

Supplied sets can be evaluated without optimization:

```sh
# Input: {"geometry": "triangle", "points": [[0,0], [1,0], [0,1]]}.
build/examples/Adaptation/SWIFT/Witness/SWIFT_Witness \
  --evaluate /tmp/witness/input.json --output /tmp/witness/evaluated.json
```
