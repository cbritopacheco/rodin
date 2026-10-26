# Atrial wall contraction demo

`AtriumContraction.cpp` solves quasistatic, finite-strain contraction on the
material-4 wedge wall. It uses Rodin's `NeoHookean`, `ActiveFiberLaw`,
`ActiveContraction`, `InternalVirtualWork`, `NewtonSolver`, and native XDMF
writer. The delivered run was built against develop revision `6d0d0fe78b`.

The activation is prescribed. No electrophysiology equation or evolving
Hill–Maxwell internal state is solved. At reference position \(z\),

\[
a(z,t)=\exp\!\left[-\left(\frac{z-z_0-vt}{w}\right)^2\right],
\qquad c(z,t)=-c_{\max}a(z,t).
\]

The active fiber branch uses the prescribed Green–Lagrange target strain
\(c\). Its stress is computed from the actual mechanical fiber strain, and
the equilibrium displacement is solved by Newton's method. The initial
snapshot is exactly unactivated. The wave starts three widths below the wall
and finishes three widths above it.

The reference fiber direction is circumferential around the positive
\(z\)-axis, projected into each wedge's wall tangent plane. At a pole, the
projected Cartesian \(x\)-direction supplies a tangent direction. This is an
illustrative fiber field, rather than an anatomical fiber reconstruction.

All inlet/outlet rim faces with labels 1, 2, 3, 5, and 6 are fixed. The inner
and outer wall faces, labeled 4 and 7, are mechanically free. The input must
contain only six-node wedges with volume material label 4.

## Parameters and units

Coordinates and displacement are in mm, time is in seconds, and stresses are
in kPa. Defaults are demonstration values:

| Parameter | Default |
|---|---:|
| Passive Young's modulus | \(50\ \mathrm{kPa}\) |
| Passive Poisson ratio | \(0.3\) |
| Active fiber stiffness | \(20\ \mathrm{kPa}\) |
| Maximum target Green–Lagrange strain magnitude | \(0.03\) |
| Wave speed | \(200\ \mathrm{mm/s}\) |
| Gaussian width | \(6\ \mathrm{mm}\) |
| Time intervals over one complete passage | \(24\) |

The model uses a compressible wall, prescribed activation, and quasistatic
equilibrium without chamber pressure or fluid coupling. Quadrature order 4
is a numerical policy for this example, not an exactness guarantee for the
non-polynomial hyperelastic response. Quantitative use requires spatial,
quadrature, and time refinement plus material/fiber calibration.

## Build and run

With a recent develop checkout containing the XDMF interpolation-order fix:

```sh
cmake --build build --target AtriumContraction -j4
build/examples/Heart/AtriumContraction \
  --mesh resources/examples/Heart/Atrium/AtriumWall_Repaired_Wedges_1mm.mesh \
  --output output/atrium-contraction/Contraction \
  --steps 24 --contraction 0.03 --speed 200 --width 6
```

`--young`, `--poisson`, and `--active-stiffness` configure the material.
`--stop-after` truncates a run for a short numerical check while preserving
the complete-passage time step. `--help` lists options.

The preflight check finds the exact minimum of each six-node wedge's
Jacobian determinant over its reference cell. The determinant is affine in
the triangle coordinates and quadratic through the thickness, so checking
the three triangle vertices and the thickness quadratic's minima is
sufficient. A nonpositive reference determinant rejects the mechanics input.

## Mesh repair used for the demonstration

The original label-4 surface had a small topological handle near
\((45.15,191.5,-193.5)\ \mathrm{mm}\). The normal extrusion crossed inside
that handle. MMG remeshing reduced the number of folded wedges substantially
but retained the handle.

The separately saved repair candidate replaces the local handle patch with
a disk, then uses Rodin's MMGS adaptation and normal extrusion:

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
conclusions. Earlier source and combined meshes remain separate artifacts.

## ParaView

Open `output/atrium-contraction/Contraction.xdmf`, keeping all companion HDF5
files in the same directory. Apply **Warp By Vector**, choose `Displacement`,
and use scale factor 1 for physical displacement in mm. Color by `Activation`
and play the time series. `FiberDirection` is a cell vector field and can be
viewed with cell-center glyphs. `ActiveExtension` stores the prescribed target
strain. The CSV beside the XDMF records time, maximum displacement,
activation maximum, Newton iterations, and final residual for each saved step.


## Reference numerical run

The final run contains 25 snapshots spanning approximately 0.49708 s, on
11,924 wedges with 12,102 shared vertices. All 24 nonlinear equilibrium
steps converged. Peak nodal displacement is approximately 0.429981 mm.
Independent ParaView/VTK checks confirm zero motion on every supported rim
node and positive deformed Jacobians throughout every cell at every stored
time. The smallest normalized deformed determinant is approximately
0.0438773 mm. Activation-weighted centroid fiber strains are negative during
the active passage, confirming mechanical shortening.

These checks were performed independently on the stored XDMF fields using
ParaView/VTK. They establish solver and export consistency for this reference
run, without establishing spatial or quadrature convergence. The smaller
preview run was a workflow check, not a spatial convergence study.
