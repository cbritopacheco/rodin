# Input and output examples

Run each executable in its own scratch directory. The programs write into the
current working directory.

| Target | Operation | Output |
| :--- | :--- | :--- |
| `Rodin_IO_MEDIT` | Save MEDIT, reload, convert to MFEM and check mesh counts | `Square.medit.mesh`, `Square.mfem.mesh` |
| `RodinIOMFEM` | Save and reload an affine scalar field; check cell-centroid values | `Square.mfem.mesh`, `Square.mfem.gf`, `Square.xdmf` and HDF5 data |
| `RodinIOXDMF` | Export a prescribed vector field at 41 times from 0 to 2 | `TimeEvolution.xdmf` and HDF5 data |

Keep each XDMF file beside its referenced HDF5 files. The time-series program
interpolates an analytic field on 1,922 triangles with 1,024 vertices.

To regenerate the tutorial figures, place the MFEM output in `OUTPUT_ROOT/MFEM`
and the time-series output in `OUTPUT_ROOT/XDMF`, then run:

```sh
pvpython examples/IO/visualize_paraview.py OUTPUT_ROOT IMAGE_DIRECTORY
```

The renderer uses the shared PDE figure palette and Manrope font, transparent
backgrounds, sparse vector glyphs and a fixed magnitude scale across snapshots.
