# Examples

Build standalone programs with `-DRODIN_BUILD_EXAMPLES=ON`. Run numerical
programs one at a time; each example's README describes its inputs and output.

## Scientific figures

Figure titles, captions, axis labels and legends use the site's bundled
Manrope font. Matplotlib scripts load it through `plot_style.py` and apply
`rodin.mplstyle`. Mathematical expressions use TeX markup with mathematical
Computer Modern font rendering; prose retains Manrope. Mathematical labels
and legends use TeX expressions throughout. The ParaView renderer uses the same
bundled Manrope file for titles and scalar bars.

Use the Turbo colormap for solution fields and numerical plots, including
vector magnitudes. Line plots sample their curve colors from Turbo through
`plot_style.turbo_colors()`; markers and line styles distinguish overlapping
series. `apply_style()` sets Turbo as the default colormap and color cycle.
Retain common field ranges in comparison panels. Use categorical
colors for material regions and mesh diagrams.

Keep figure and axes backgrounds transparent. Include mesh edges where they
clarify geometry without hiding the field, and reserve space for color bars
and labels. HTML captions inherit Manrope from the documentation theme.

Generate figures from computed output and retain the script that reproduces
them. An axis background is transparent; colors representing the solution itself
remain opaque. Use alpha-capable WebM/VP9 for transparent animations. Encode documentation
videos at **30 or 60 fps**, with **30 fps** as the default. Preserve simulation
time when resampling saved snapshots; display interpolation does not add
solver timesteps. `PDEs/README.md` and `Variational/README.md` give reproduction commands
for the tutorial figures. Existing optimization illustrations already contain
transparent backgrounds.

Export individual panels for documentation comparisons. Wrap related figures
in `rodin-mesh-stages`. Desktop pages show two panels per row, center the
last panel of an odd-sized group, and stack panels on narrow screens. A standalone
panel can use up to 36rem of the content width. Keep captions with their
panels. Each pair shares a frame containing both images and their captions;
a single panel has its own centered frame. Composite images remain useful for video posters and external exports.

Plot cellwise fields from their cell values. In particular, a P0 field assigns
one flat color per original cell. Preserve the exported point/cell association
in ParaView; the shared `color()` helper selects it from the available arrays.
Discontinuous higher-order fields retain their cell-local samples and vertices.

Render light-theme figures with `RODIN_PLOT_THEME=light` and dark-theme figures
with `RODIN_PLOT_THEME=dark`. Their numerical ranges and Turbo colors remain
identical. Dark variants use a `-dark` filename suffix and are selected by the
documentation theme. Plot lettering uses the body font and high-contrast ink.
