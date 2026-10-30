# Examples

Build standalone programs with `-DRODIN_BUILD_EXAMPLES=ON`. Run numerical
programs one at a time; each example's README describes its inputs and output.

## Scientific figures

Figure titles, captions, axis labels and legends use the site's bundled
Manrope font. Matplotlib scripts load it through `plot_style.py` and apply
`rodin.mplstyle`. Mathematical expressions use TeX markup with mathematical
font rendering; prose retains Manrope. The ParaView renderer uses the same
bundled Manrope file for titles and scalar bars.

Keep figure and axes backgrounds transparent. Include mesh edges where they
clarify geometry without hiding the field, and reserve space for color bars
and labels. HTML captions inherit Manrope from the documentation theme.

Generate figures from computed output and retain the script that reproduces
them. An axis background is transparent; colors representing the solution itself
remain opaque. Use alpha-capable WebM/VP9 for transparent animations. `PDEs/README.md` and `Variational/README.md` give reproduction commands
for the tutorial figures. Existing optimization illustrations already contain
transparent backgrounds.

Export individual panels for documentation comparisons. Wrap related figures
in `rodin-mesh-stages`. Desktop pages show two panels per row, center the
last panel of an odd-sized group, and stack panels on narrow screens. A single
panel uses the same desktop width as one column. Keep captions with their
panels. Composite images remain useful for video posters and external exports.

Plot cellwise fields from their cell values. In particular, a P0 field assigns
one flat color per original cell. Preserve the exported point/cell association
in ParaView; the shared `color()` helper selects it from the available arrays.
Discontinuous higher-order fields retain their cell-local samples and vertices.
