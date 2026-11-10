# Relative Jacobian and distortion

This prescribed-deformation example uses nine vertices and eight triangles on
the unit square. No optimization or linear solve is performed. The same field

\[
u(x,y)=\bigl(\gamma y,\;b x(1-x)\bigr),\qquad
\gamma=0.5,\quad b=0.4,
\]

is interpolated in continuous degree-one and degree-two displacement spaces.
The corresponding geometry has the same degree as the displacement. The
degree-one cells remain straight-sided; degree two represents the quadratic
bow exactly. Both configurations share the same nine vertex positions.

The diagnostics are relative to the undeformed mesh:

\[
F=I+\nabla u_h,\qquad j=\det F,\qquad
Q=\frac{\|F\|_F^2}{2j}.
\]

Thus `j` measures local area change, while `Q` measures distortion independently
of positive uniform scaling. These are not the absolute quality measures of
the reference triangle.

## Run

Build the `Adaptation_Quality` target in a configured Rodin build, then run:

```sh
./build/examples/Adaptation/Quality/Adaptation_Quality /tmp/rodin-quality
python3 examples/Adaptation/Quality/plot.py /tmp/rodin-quality
```

The plotting script requires NumPy and Matplotlib. The output directory contains:

- `P1.xdmf` and `P2.xdmf`, with background and moved blocks and their displacement
  fields, plus accompanying HDF5 files for ParaView.
- `P1.csv` and `P2.csv`, with reference and moved coordinates and pointwise `j`
  and `Q`. Each original cell is sampled independently on a 32-subdivision
  triangular lattice, including its boundary. This preserves jumps in the
  degree-one gradient rather than averaging across cells.
- `quality.png`, `quality.svg` and `quality.pdf`, showing both configurations
  with shared color scales. The SVG and PDF are vector figures.

The C++ example checks the installed geometry Jacobian against the displacement
gradient at every sample. The plotting script independently verifies the
coordinates and both quality measures against the analytic deformation (or its
piecewise-linear interpolant) before saving the figures. Color interpolation
between samples is only a visualization, not a certificate between points.
