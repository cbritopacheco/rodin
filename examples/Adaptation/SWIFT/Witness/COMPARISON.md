# Local Minimax Comparison

This report records the historical experiment before CGAL integration. Its
timings must not be attributed to the current CGAL-backed executable.

## Protocol

The experiment compares continuous reference-space coverage, not SWIFT fitting
or physical-space quality. Every method starts from the same lattice-derived
configuration. One serial run is performed per geometry and method, without
warmups or repeated timing runs. The point geometry cannot contain 16 distinct
sites and is reported with its unique site. All other geometries use 16 sites.

| Setting | Value |
|---|---|
| Local enclosing-ball/Newton refinement cap | 100 Voronoi sweeps |
| Motion tolerance | \(10^{-6}\) |
| Current search | 2 seeded starts, 100 generations/start, up to 30 polling sweeps/start |
| Smooth maximum final squared-distance bias bound | \(10^{-8}\) |
| Newton cap | 50 corrections per temperature per cell |
| Temperature schedule | \(10^{-2}\), reduced by \(0.1\) to the prescribed bias bound |
| Newton gradient/step tolerance | \(10^{-10}\) |
| CPU policy | One optimizer worker; lowered scheduling priority |

The budgets are not equivalent work budgets. The current search retains its
existing defaults; the local methods use the same sweep cap and stopping rule.
The temperature bound describes smoothing error, not certified local-solve
accuracy when Newton reaches its correction cap.

Artifacts are in `tmp/witness-comparison-20261010/`. `results.json` contains all
24 completed pairs, initial/final coordinates, continuous Voronoi cells,
independent lattice lower bounds, sweep histories and costs. Each method also
has its own SVG gallery and article plate. No production defaults were changed.

## Coverage and Cost

Each entry is **continuous covering radius / total seconds**. Smaller radii
are better. Other agents' compilation and calibration were active; these
single-run wall times are indicative, not uncontaminated performance estimates.
Only the optimizer itself was serial; no calibration process was stopped.

| Geometry | Enclosing ball | Smooth Newton | Current search |
|---|---:|---:|---:|
| Point, one site | 0 / <0.001 | 0 / <0.001 | 0 / <0.001 |
| Segment | 0.03125000 / <0.001 | 0.03125000 / <0.001 | 0.03125000 / <0.001 |
| Triangle | 0.13067204 / 0.01042 | 0.13067212 / 0.11345 | 0.14137212 / 1.79688 |
| Quadrilateral | 0.17677759 / 0.00183 | 0.17677759 / 0.00327 | 0.23570226 / 2.25323 |
| Tetrahedron | 0.21826467 / 0.04207 | 0.21826462 / 0.30253 | 0.23696923 / 15.55507 |
| Pyramid | 0.26511257 / 0.03685 | 0.26511243 / 0.27520 | 0.28574463 / 19.16971 |
| Hexahedron | 0.36355812 / 0.04790 | 0.36355808 / 0.37930 | 0.41195138 / 26.23444 |
| Wedge | 0.28923699 / 0.04752 | 0.28923708 / 0.38660 | 0.35306775 / 21.80543 |
| Total measured seconds | 0.18673 | 1.46046 | 86.81478 |

The largest radius difference between enclosing-ball and Newton refinement is
\(1.415\times10^{-7}\). Both improve on the current search for every tested
two- and three-dimensional geometry. Point and segment use the same analytic
shortcuts, so their timings do not compare iterative solvers.

| Geometry | Sweeps, both local methods | Termination | Newton corrections, all cells/temperatures |
|---|---:|---|---:|
| Triangle | 100 | Sweep cap | 64,246 |
| Quadrilateral | 17 | Motion tolerance | 1,118 |
| Tetrahedron | 100 | Sweep cap | 74,142 |
| Pyramid | 82 | Motion tolerance | 62,280 |
| Hexahedron | 100 | Sweep cap | 77,693 |
| Wedge | 100 | Sweep cap | 82,971 |

The cap-limited local runs still move: the final maximum witness motions for
triangle, tetrahedron, hexahedron and wedge are respectively
\(4.78\times10^{-6}\), \(1.95\times10^{-6}\), \(1.32\times10^{-5}\), and
\(1.11\times10^{-4}\). No global optimum is certified. The quadrilateral radius
is close to the radius \(\sqrt{2}/8\) of a regular four-by-four cell-center
covering; that agreement alone is not a proof of global optimality.

## Why the Cost Differs

The enclosing-ball updates solve the finite cell minimax problems directly.
Newton repeatedly evaluates a smooth maximum while reducing its temperature;
the small linear systems are cheap, but the accumulated evaluations and
corrections are not free. This continuation is deliberately accurate, rather
than a benchmark of an aggressively loose surrogate.

In the enclosing-ball runs, Voronoi reconstruction accounts for approximately
50--83 percent of total time on the two- and three-dimensional geometries.
Newton instead spends most of its time in the local smooth solves, except for
the very short quadrilateral case. Reusing a single inverse would not remove
the changing partition or smooth weights.

The current search makes 41,026--63,938 objective calls per geometry. Every
continuous objective call clips all 16 Voronoi cells. Its quadrilateral result
essentially retains the initial radius even after both seeded starts and local
polling. A local coordinate-poll step tolerance is therefore not evidence that
the configuration is a minimax optimum: simultaneous cell-center motion improves
this same configuration substantially.

## Numerical Checks and Degeneracies

All 24 final radii agree with an independent SciPy/Qhull calculation to nine
decimal places. All local histories are nonincreasing within the recorded
floating-point tolerance. The 13 Python integration/audit tests pass. C++ checks
also pass for smooth gradients/Hessians by finite differences, local minimax
accuracy, and the observed degenerate support regression.

An initial trial of the Apache-licensed `hbf/miniball` implementation failed
explicit containment checks on symmetric lattice cells, including a triangle
cell. Its outputs were not accepted into the comparison. The successful
experiment uses Bernd Gartner's `Miniball.hpp` (2021 revision) from
`weddige/miniball`, under GPL-3.0-or-later. The header remains external in `/tmp`;
it is not a dependency of the standard generator or SWIFT.

The Gartner implementation can report negative support coefficients for
near-affinely-dependent boundary points even when its containment residual is
at floating-point roundoff. Such a report is not simply ignored: a separate
convex-combination certificate is sought using subsets of at most four active
boundary points. Containment plus this certificate establishes the local ball's
minimax optimality to the stated floating-point tolerance. Tetrahedron required
nine independent support checks and two checked input-order retries; wedge
required three independent checks and no retries. Failed intermediate attempts
were resumed from saved completed pairs, not pooled as completed runs.

| Provenance | SHA-256 |
|---|---|
| External Gartner header | `31e02559baf81a57b0d4be942a30f76da30afa841e3e98cd617acbe79289b93f` |
| Completed 3D comparison executable | `9dbd2b1ea7ac6a4eb21b3676aeedc412e084522ed449686511f842a4f0396eb3` |

Completed 2D pairs precede addition of the independent support-check fallback;
their local support tests passed without that fallback. No completed pair was
repeated. Subsequent executable changes add `--check` verification only.

## Conclusion

Checked enclosing-ball refinement is the preferred candidate for replacing the
expensive default search: it gives practically the same coverage as accurate
smooth Newton at much smaller measured cost, and improves all tested 2D/3D
results over the current search. Smooth Newton does not demonstrate a coverage
advantage sufficient to justify its additional work in this experiment.

Before integration, resolve the external library's licensing and keep the
containment/support checks. More initial configurations or further refinement
may improve local solutions, but are separate experiments. None of the local
results should be advertised as globally optimal witness sets.

## Canonical CGAL Integration

The canonical generator now discovers external CGAL only when
`RODIN_BUILD_SWIFT_WITNESS=ON` and links `CGAL::CGAL` privately to the Witness
executable. No enclosing-ball implementation is vendored. The conditional CMake
configuration and the GPL distribution implications are described in the
[usage documentation](README.md).

CGAL 6.1.1 was tested with Apple Clang 21 against the existing Rodin libraries,
using an isolated target build. The configured Clang 19 toolchain failed first
in GMP's use of standard floating-point macros; changing include order then
exposed missing infinity/NaN macros in libc++. No toolchain-specific workaround
was retained in the source. The running campaign executables were not rebuilt.

All 15 integration/audit tests pass, including independent Qhull coverage and
configuration with CGAL explicitly unavailable when the optional tool is disabled.
For the same eight geometries, initialization, count and sweep cap, the largest
radius difference from the historical enclosing-ball results is less than
\(9\times10^{-17}\). Artifacts are in `tmp/witness-cgal-20261010/all.json`.

| Geometry | Radius | Sweeps | Independent support fallbacks |
|---|---:|---:|---:|
| Point, one site | 0 | 0 | 0 |
| Segment | 0.03125000 | Analytic | 0 |
| Triangle | 0.13067204 | 100 | 0 |
| Quadrilateral | 0.17677759 | 17 | 0 |
| Tetrahedron | 0.21826467 | 100 | 1 |
| Pyramid | 0.26511257 | 82 | 2 |
| Hexahedron | 0.36355812 | 100 | 0 |
| Wedge | 0.28923699 | 100 | 2 |

Containment and support optimality remain independently checked. Five local
CGAL results required the finite-support fallback; they were not accepted
unchecked. No new uncontaminated timing comparison is claimed.
