# Atrial wall contraction demo

`AtriumContraction.cpp` solves quasistatic, finite-strain contraction on the
material-4 wedge wall with Rodin's `NeoHookean`, static `ActiveFiberLaw`,
`ActiveContraction`, `InternalVirtualWork`, `NewtonSolver`, and XDMF writer.
Activation is prescribed; no electrophysiology equation or evolving
Hill–Maxwell internal state is solved.

## A plane wave starting at the anatomical base

In the original `Coarse.mesh` volume, label 6 identifies the basal cap. The
original tetrahedra have volume material 1. In the wedge-only wall, the cap
is absent and its boundary label is preserved on the adjoining rim faces.
The wave direction is determined by the basal rim.

The inner vertices of rim 6 define a least-squares plane through their mean
position \(b\). Its unit normal \(d\) is oriented toward the wall centroid.
The projected coordinate and activation are

\[
s(X)=d\cdot(X-b),\qquad
 a(X,t)=\exp\!\left[-\left(\frac{s(X)-vt}{w}\right)^2\right],
\qquad c(X,t)=-c_{\max}a(X,t).
\]

The wave crest starts on the fitted basal plane at \(t=0\) and moves in the
positive \(d\)-direction. The first stored state is an equilibrium solve with
that basal activation. The last crest position is three widths past the
maximum wall projection.

This is a strict planar pulse. A portion of the appendage lies behind the
basal plane and therefore receives only the initial Gaussian tail. The
appendage-apex landmark has projection about \(-8.18\ \mathrm{mm}\) on this
mesh. Complete base-to-appendage propagation would require a different
activation construction, such as a wall-distance wave.

The fitted origin, inward normal, projected coordinate range, and landmark
information are written to `<stem>.geometry.json`. `WaveCoordinate` is a
static nodal XDMF field for checking the plane orientation.

## Simple landmark-based fibers

The default left-appendage apex landmark is

\[
X_{\mathrm{LAA}}=(46.39638367,187.44307667,-205.03317767)\ \mathrm{mm}.
\]

Coordinates are used instead of vertex or triangle numbers because MMG
renumbers entities. An apex patch is selected using inner-triangle centroids
within 3 mm of this point, together with the corresponding outer triangles.
A finite patch avoids imposing a point-only Dirichlet condition in the wall
volume. Its radius and center can be changed from the command line.

A single harmonic coordinate is solved with Rodin:

\[
-\Delta\psi=0\quad\text{in }\Omega_{\mathrm{wall}},\qquad
\psi=0\quad\text{on rim 6},\qquad
\psi=1\quad\text{on the appendage patch},
\]

with zero normal flux on the remaining boundary. A temporary boundary label
8 is used for this solve and then restored, preserving all original wall and
rim labels.

For each wedge, the known extrusion vectors supply a unit transmural
normal \(m\). The harmonic gradient is projected into the wall tangent
plane and the fiber is constructed as

\[
\ell=\frac{(I-mm^T)\nabla\psi}{\|(I-mm^T)\nabla\psi\|},
\qquad f=\ell\times m.
\]

Thus the fibers are tangent to the wall and follow contours of the
base-to-appendage coordinate. They are stored as one constant unit direction
per wedge. There is no transmural rotation in this one-layer model.

This is a simplified use of the Laplace-gradient frame in Piersanti et al.,
*Modeling cardiac muscle fibers in ventricular and atrial electrophysiology
simulations*, [arXiv:2101.10960v1](https://arxiv.org/abs/2101.10960v1), equations
(3) and (5). It does not implement their regional atrial bundle selection,
pulmonary-vein rules, or a validated anatomical fiber reconstruction.

Critical points are handled explicitly. If the harmonic tangent vanishes,
the projected direction toward the apex is used; at a remaining pole, a
Cartesian tangent is chosen. `FiberFallback` flags these cells: 0 means the
harmonic construction, 1 the geometric-apex fallback, and 2 the pole fallback.
The example mesh has one fallback cell, inside the constant-potential apex
patch. `FiberCoordinate`, `LongitudinalDirection`, `WallNormal`, and
`FiberDirection` are exported for inspection.

## Mechanics, parameters, and units

The active branch uses the prescribed Green–Lagrange target strain \(c\).
Stress depends on the actual fiber strain, and Newton's method solves the
mechanical equilibrium. All rim faces with labels 1, 2, 3, 5, and 6 are fixed.
The inner and outer wall faces, labeled 4 and 7, are free. The input must
contain only six-node wedges with volume material 4.

Coordinates and displacement are in mm, time in seconds, and stresses in
kPa. Defaults are demonstration values:

| Parameter | Default |
|---|---:|
| Passive Young's modulus | \(50\ \mathrm{kPa}\) |
| Passive Poisson ratio | \(0.3\) |
| Active fiber stiffness | \(20\ \mathrm{kPa}\) |
| Maximum target Green–Lagrange strain magnitude | \(0.03\) |
| Wave speed | \(200\ \mathrm{mm/s}\) |
| Gaussian width | \(6\ \mathrm{mm}\) |
| Basal boundary label | \(6\) |
| Appendage patch radius | \(3\ \mathrm{mm}\) |
| Time intervals | \(24\) |

The wall is compressible, with quasistatic equilibrium and no chamber
pressure or fluid coupling. Mechanical quadrature order 4 is an example
policy, not an exactness guarantee for hyperelasticity. Quantitative use
requires spatial, quadrature, and time refinement and material calibration.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target AtriumContraction -j4
build/examples/Heart/AtriumContraction/AtriumContraction \
  --mesh resources/examples/Heart/Atrium/AtriumWall_Repaired_Wedges_1mm.mesh \
  --output output/atrium-base6-plane/Contraction \
  --base-label 6 --apex 46.39638367 187.44307667 -205.03317767 \
  --apex-radius 3 --steps 24 --contraction 0.03 --speed 200 --width 6
```

`--young`, `--poisson`, and `--active-stiffness` configure the material.
`--stop-after` truncates the positive-time steps, including the initial
activated equilibrium in the stored series. `--fibers-only` exports the
geometry and fiber fields without solving mechanics; its zero displacement
is an inspection placeholder. `--help` lists the options.

The preflight check finds the exact minimum of each six-node wedge's
Jacobian determinant over its reference cell. The determinant is affine in
the triangle coordinates and quadratic through the thickness, so evaluating
the triangle vertices and thickness-quadratic minima is sufficient. A
nonpositive reference determinant rejects the input. The same test is applied
to each converged deformed cell before saving a mechanics snapshot.

## Mesh repair used for the demonstration

The source label-4 surface contains a small topological handle near
\((45.15,191.5,-193.5)\ \mathrm{mm}\). A direct normal extrusion folds inside
that handle, which persists after MMGS remeshing.

The example wall is constructed by replacing the local handle patch with
a disk, followed by Rodin's MMGS adaptation and normal extrusion:

1. Remesh the label-4 surface with tagged inlet/outlet boundary edges.
2. Replace the handle patch selected within an 8 mm radius by a disk.
3. Remesh with target size 5 mm, minimum size 1.5 mm, maximum size 7.5 mm,
   Hausdorff parameter 1 mm, gradation 1.3, and automatic ridge detection off.
4. Extrude one outward layer by exactly 1 mm using Rodin. Retain the
   inlet/outlet tags on rim quadrilaterals.

The repaired surface has 6,051 vertices and 11,924 triangles. The wall has
12,102 vertices, 11,924 wedges, and 184 rim quadrilaterals. Its minimum
reference determinant divided by base-triangle double area is approximately
\(0.0442645\ \mathrm{mm}\), and all cells have positive determinants
throughout their reference domain. Nodal extrusion distances equal 1 mm to
floating-point precision. This local positivity check does not certify the
absence of intersections between distant cells.

This repair changes the anatomy and topology. Sampled original-vertex to
repaired-surface distances have mean 0.221 mm and maximum 4.692 mm; sampled
repaired-vertex to original-surface distances have mean 0.024 mm and maximum
2.027 mm. These are vertex-to-triangle samples, not certified Hausdorff
bounds. Inspect the repaired patch before using the mesh for anatomical
conclusions.

## Reference run

With the default basal plane and harmonic fibers, all 25 stored equilibria
converged over \(0.452633\ \mathrm{s}\). The maximum nodal displacement was
\(0.284546\ \mathrm{mm}\). Independent ParaView readback found 11,924
six-node wedges and 12,102 shared vertices, zero displacement on every rim
vertex, and positive Jacobians throughout every deformed reference cell.
The minimum deformed determinant divided by the deformed base-triangle
double area was \(0.0438848\ \mathrm{mm}\). Activation values agreed with
the prescribed plane to within \(3.0\times10^{-15}\).

The harmonic linear-system relative residual was
\(3.56\times10^{-15}\). Fiber vectors were unit length and orthogonal to
the exported wall normals and longitudinal directions. One apex-patch cell
used the explicitly flagged geometric fallback. These checks verify the
specified demonstration, not anatomical fiber accuracy or numerical
convergence under mesh refinement.

## ParaView

Open `<stem>.xdmf`, keeping its companion HDF5 files in the same directory.
Apply **Warp By Vector**, choose `Displacement`, and use scale factor 1 for
physical displacement in mm. Color by `Activation` and play the time series.

For fibers, apply **Cell Centers**, then line glyphs oriented by the cell
vector `FiberDirection`. Color the wall by `FiberCoordinate` and inspect
`FiberFallback` near the apex patch. `ActiveExtension` stores the prescribed
target strain. The CSV records time, maximum displacement, activation
maximum, Newton iterations, final residual, and minimum scaled Jacobian.

The geometry repair above and the prescribed fibers are demonstration
choices. Inspect the repaired surface and the fiber field before drawing
anatomical conclusions.
