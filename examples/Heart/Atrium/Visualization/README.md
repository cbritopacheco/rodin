# Atrial flow visualization

The fluid component writes velocity and pressure at the requested output times.
`Volume.py` renders the velocity magnitude through the cavity, with a translucent
surface overlay to locate the wall and openings. An orthographic isometric
camera is centered on the cavity bounds and fitted to the viewport. The layout
shows back and front views side by side, with the front camera rotated by
$180^\circ$ about the vertical axis. Both views share the same velocity field,
colour scale, opacity function, and animation time. The translucent cavity
surface provides a visible anatomical outline around the flow.
The displayed scalar is

$$
s(X,t)=c\,\lVert u_h(X,t)\rVert_2,
\qquad
c=\begin{cases}100,&\text{for cm/s},\\1,&\text{for m/s}.\end{cases}
$$

Fluid velocity is stored in metres per second. The conversion changes the
displayed scalar and legend together. A single colour and opacity range is
used throughout the animation; its maximum is taken over the selected states
unless `--vmax` specifies it explicitly in the selected units.
MPI pieces are merged with coincident vertices before surface extraction,
so partition boundaries are excluded from the cavity outline.

The input may come from [fixed-domain fluid or coupled FSI](../README.md).
FSI fluid XDMF already contains moving physical coordinates; an additional
displacement warp would apply the motion twice. Solid output instead stores
reference coordinates and its separate `Displacement` field. Rendering computes
observables from saved states and adds no fluid or structural evolution.

## Rendering

ParaView's `pvpython` is required; movie encoding also requires `ffmpeg`.
From the repository root:

```sh
pvpython --force-offscreen-rendering \
  examples/Heart/Atrium/Visualization/Volume.py \
  output/atrium-fluid/Fluid.xdmf \
  output/atrium-fluid/Volume.mp4 \
  --units cm/s --resolution 128 --fps 12
```

On macOS, `pvpython` can be invoked through
`/Applications/ParaView-6.2.0.app/Contents/bin/pvpython`.
To generate an interactive scene without rendering a movie or requiring
`ffmpeg`, use:

```sh
pvpython examples/Heart/Atrium/Visualization/Volume.py \
  output/atrium-fluid/Fluid.xdmf \
  output/atrium-fluid/Volume.pvsm --state-only
```

Open the resulting `.pvsm` through **File > Load State** in ParaView and press
**Play** to advance through the saved simulation times.
Movie rendering also saves a ParaView state (`.pvsm`), the final frame (`.png`),
and rendering parameters and physical timestamps (`.json`) beside the movie.
The state can be opened with **File > Load State** to adjust the camera,
transfer functions, and surface opacity interactively.

`--resolution` sets the samples along each axis of a Cartesian rendering grid.
This interpolation controls visualization quality; it does not refine the
finite-element mesh or increase the accuracy of the computed flow.
`--size 1920 1080` sets the image dimensions, and `--max-frames` selects a
uniformly spaced subset of the available states.

Each frame represents a solved state and displays its physical time. `--fps`
sets playback speed. A short startup calculation remains a startup calculation
when played slowly; a full-cycle animation requires solving and saving that
cycle. Movies with comparisons should use the same `--vmax`, camera, and
opacity settings.

## Interactive ParaView pipeline

1. Open the fluid XDMF series, apply **Merge Blocks** with **Merge Points**
   enabled, then apply **Calculator**, with
   `100*mag(Velocity)` for speed in centimetres per second.
2. Apply **Resample To Image**, then select **Volume** representation and
   colour by the calculated scalar. Set the colour range for the full series
   and edit the opacity function to make low speeds transparent.
3. Add **Extract Surface** from the merged mesh, select a neutral colour,
   and reduce its opacity to show the cavity outline.
4. Add **Annotate Time Filter**, choose the camera, and use
   **File > Save Animation** to export the selected time range.

## Solid displacement in FSI

When `Solid.xdmf` exists beside the fluid series, both views also render the
complete deformed wedge wall. `--solid /path/to/Solid.xdmf` selects it explicitly.
The pipeline applies `Displacement` with scale factor 1 in mm, then converts
geometry to metres to align with the ALE fluid. The wall is translucent and
colored by displacement magnitude in mm, with a separate legend from fluid
speed. Fluid geometry is displayed directly. Cameras fit the union of fluid
and deformed-solid bounds over the selected saved states.

Only times saved by both components are included when the scene is generated.
For a running calculation, reload both XDMF readers in ParaView to discover new
states. Existing scenes must be regenerated to add the solid pipeline. No
visual displacement amplification is applied.
