# PDE tutorials

The standalone `Poisson`, `Elasticity`, `Helmholtz`, `Periodic`, `SurfaceMesh`,
`Stokes`, `NavierStokes`, and `MixedBoundary` targets generate their own meshes and write XDMF
files with accompanying HDF5 data in the working directory. Enable examples
with `-DRODIN_BUILD_EXAMPLES=ON`. Run numerical examples one at a time.

To regenerate the tutorial figures and cavity video, run each executable from
its own directory named after the target, under a common output directory.

```sh
mkdir -p output/Helmholtz
cd output/Helmholtz
/path/to/build/examples/PDEs/Helmholtz
```

Repeat for `Poisson`, `Elasticity`, `Periodic`, `SurfaceMesh`, `Stokes`, and
`NavierStokes`, and `MixedBoundary`. Then use ParaView's `pvpython` from this
repository's root.

```sh
pvpython examples/PDEs/visualize_paraview.py /path/to/output /path/to/figures \
  --frames /path/to/temporary-frames
```

It loads the actual XDMF/HDF5 fields, plots the solution surfaces and cavity
streamlines, and exports transparent PNGs. FFmpeg with a VP9 encoder creates
an alpha-capable WebM from the cavity frames. Frames are temporary build
outputs and should stay outside the repository. Manrope typography is loaded
from `brand/fonts/static/Manrope-Regular.ttf` for titles and scalar bars.

The alternative Matplotlib renderer `visualize.py` remains available for
comparison, using `examples/rodin.mplstyle`. Neither renderer interprets finite
element coefficient storage as point values. Both consume the exported fields.
Streamlines are integrated through each instantaneous velocity field.

The transient cavity uses 2,048 triangles, quadratic velocity and linear
pressure, and a solver step of 0.01. Its smooth compact lid patch travels
horizontally between $x=0.28$ and $x=0.72$ with period 8, at peak-speed Reynolds
number 200 (peak lid speed 2 and viscosity 0.01). The 2,000 Oseen steps
cover simulation time 0–20. Output is saved
every five solver steps; every second output frame at 10 fps preserves
simulation-time playback. The example follows a transient flow. A convergence study requires independent
  mesh and time-step refinements.

The boundary-conditions guide uses `MixedBoundary`, whose exact solution on the
unit square is $u=x$. After running that executable under
`output/MixedBoundary`, its figure can be regenerated independently.

```sh
pvpython examples/PDEs/visualize_paraview.py /path/to/output /path/to/figures --example MixedBoundary
```

## Surface refinement check

`SurfaceMesh` accepts an optional vertex count per background axis; the default
is 11. To retain the checked 6/11/21 sequence, run the cases sequentially in
separate `refinement-6`, `refinement-11` and `refinement-21` directories.
Each writes `SurfaceMesh.csv` alongside its XDMF/HDF5 field output. The CSV
contains the assembled relative residual, surface area and relative error
against the ambient polynomial `x*y` on the actual flat-facet mesh.

Reproduce the refinement chart with the following command.

```sh
python3 examples/PDEs/visualize_surface_convergence.py \
  /path/to/refinement-root doc/Images/SurfaceMeshConvergence.png
```
