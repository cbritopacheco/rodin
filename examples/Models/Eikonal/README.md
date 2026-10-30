# Arrival times with fast marching

`FMM.cpp` solves a point-source problem on the unit square with constant
positive speeds 1 and 2. It uses 16, 32 and 64 intervals per axis, placing an
exact mesh vertex at the center. The equation and analytic reference are

\[
F|\nabla T|=1,\qquad T(x)=\frac{|x-(0.5,0.5)|}{F}.
\]

Build `Eikonal_FMM` and run it in a scratch directory. No input data are
required. `FMMErrors.csv` records maximum and mean nodal errors for every
case. `FMMArrival.xdmf` and its HDF5 companions contain the finest unit-speed
computed and analytic nodal fields. All cases run sequentially.

The executable checks finite nonnegative values, a unique source, decreasing
maximum error with refinement and the loose regression envelope described in
the source. The observed unit-speed maximum errors are 0.0474507, 0.0299216
and 0.0182986. Speed 2 halves these errors. These measurements apply to the
generated triangular grids and central point source.

Render the actual field and contours with ParaView:

```sh
pvpython examples/Models/Eikonal/visualize_paraview.py \
  /path/to/results doc/Images/FMMArrival.png
```

Titles use the bundled Manrope font, mathematical labels use TeX, and the
background is transparent. Contours clarify the front without dense mesh edges.
