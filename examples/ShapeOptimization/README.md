# Cantilever examples

`SimpleCantilever2D` displaces a material boundary and remeshes it. Its XDMF
history exports the elastic `displacement` and unscaled regularized `descent`
at each solved state. The update scale is the source's `hmin`, currently 0.01.

`LevelSetCantilever2D` exports the current distance `dist`, normalized descent
`dJ` and transported field `advect` on the background `domain` block. The
`state` block contains the material submesh and elastic displacement `u`.
The tutorial pages document their input resources, controls and measured runs.

Run each executable in its own scratch directory. With two or more recorded
states available, render a single update from the repository root.

```sh
pvpython examples/ShapeOptimization/visualize_stages.py \
  /path/to/simple/SimpleCantilever.xdmf doc/Examples/ShapeOptimization \
  --mode simple --iteration 0 --step 0.01

pvpython examples/ShapeOptimization/visualize_stages.py \
  /path/to/level-set/out/ShapeOptimization.xdmf doc/Examples/ShapeOptimization \
  --mode levelset --iteration 0
```

The renderer extracts the selected state and the next state into temporary
XDMF views referring to the original HDF5 arrays. ParaView plots the fields
and meshes with transparent backgrounds, Manrope text and TeX labels.

For the direct method, the pre-remeshing mesh is obtained by applying the
exported P1 descent to its reference vertices using `Warp By Vector` and the
same normalization as the C++ update. Set `--step` to the run's movement
scale when changing `hmin`. The next snapshot supplies the mesh after MMG.
The color scale for descent is normalized by its sampled maximum; arrows
have equal display length to emphasize direction. Glyphs sample every Nth
node in mesh order with 50% opacity. The default stride targets at most
100 glyphs; `--glyph-stride N` sets the spacing explicitly.

For the level-set method, distance and transported-field colors are clipped
to the interval from -0.05 to 0.05. Their zero contours retain the computed
values. The last panel shows the material mesh used by the next state, after
reconstruction and the next iteration's geometry preparation.

`--iteration` selects a zero-based snapshot index and requires a following
snapshot. `visualize_paraview.py` renders the first/last elastic-state
comparison from a complete history. Generated fields and temporary data
remain outside the repository.
