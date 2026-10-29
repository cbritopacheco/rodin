# Rodin visual identity

Rodin's visual identity combines a fluid-line thinker, the Rodin wordmark,
and an R monogram. The identity is intended for scientific publications,
technical documentation, presentations, and software interfaces.

## Reference artwork

The kit contains SVG reference artwork and native Affinity Designer
files. Matching filenames identify the same design in each format.
Affinity documents retain editable Manrope text alongside vector curves,
fills, gradients, and layouts. SVG reference lettering is outlined from
Manrope so that placement does not depend on installed fonts.

| Folder | Contents |
| --- | --- |
| `logos/` | Horizontal and stacked signatures, standalone thinker, outlined wordmark, R monogram, and outline R |
| `icons/` | Light and dark app icons and avatars, micro R, monochrome mask icon, and maskable icon artwork |
| `designs/` | Visual identity board, repository social card, and presentation cover |
| `fonts/` | Bundled Manrope variable font, static desktop weights, and SIL Open Font License |
| `licenses/` | License for the GitHub icon used on the social card |

Start with [the visual identity board](designs/brand-board.svg).

## Marks

| Use | Reference SVG |
| --- | --- |
| Primary signature on a light background | [Horizontal color](logos/horizontal-color.svg) |
| Primary signature on a dark background | [Horizontal white](logos/horizontal-white.svg) |
| Single-color signature | [Horizontal navy](logos/horizontal-navy.svg) |
| Compact color mark | [R monogram](logos/monogram-color.svg) |
| Compact single-color mark | [R monogram in navy](logos/monogram-navy.svg) |
| Compact reversed mark | [R monogram in white](logos/monogram-white.svg) |
| Small interface mark | [Micro R](icons/favicon.svg) |

The color, navy, and white R monograms share the same geometry.
White variants have transparent backgrounds. The outline R is intended
for larger placements. Reference SVG lettering is stored as vector
outlines and does not require an installed font. Artwork is centered
within its placement area; horizontal signatures align the vertical
centers of the thinker and wordmark.

## Color and typography

| Color | Value | Role |
| --- | --- | --- |
| Navy | `#071F46` | Primary identity, body text, single-color artwork |
| Cyan | `#13B6D0` | Graphic accents and flowing contours |
| Blue | `#17679A` | Intermediate gradient color |
| Paper | `#F4F8FB` | Light background |
| White | `#FFFFFF` | Reversed identity |
| Muted | `#536477` | Secondary text |

Use navy for body copy and cyan for graphic accents. Avoid cyan for small
text on white. The identity uses Manrope Bold for the wordmark,
Manrope Semibold for headings and labels, and Manrope Regular for
supporting text. Install the three TTF files in `fonts/static/` before
editing live text in Affinity Designer. The variable TTF in `fonts/`
supports weights from 200 to 800. Avoid installing both the static and
variable versions simultaneously to prevent duplicate family entries.

Manrope is bundled under the SIL Open Font License 1.1, with its
copyright notice in [fonts/OFL.txt](fonts/OFL.txt). Static desktop
weights are instances of the bundled variable font. Retain the license
when redistributing the font files. The GitHub mark is supplied by
Octicons; its [MIT license](licenses/github-icon-MIT.txt) is included.

The caption is always:

> A modern C++ FEM library for shape and topology optimization.

Preserve its wording, capitalization, and final period, including when
it wraps onto two lines. The social card places the HTTPS website URL
above the GitHub icon and repository URL; its caption appears above
the link block.

## Placement and export

- Keep clear space of at least one quarter of the symbol's height around
  the mark. The SVG artboard inset does not provide this full clear space.
- As starting minimum heights, use 64 px for the primary signature, 48 px
  for the standalone thinker, and 64 px for the detailed R monogram.
  Check fine contours at the intended display size.
- Use the micro R at 16–32 px and the outline R at 128 px or larger.
- Use white artwork on dark, quiet backgrounds and navy artwork for
  single-color printing.
- Preserve proportions and contours. Avoid added shadows, rotation,
  stretching, or busy backgrounds.

Export delivery formats and pixel sizes from Affinity Designer as needed.
The social card uses a 1280 by 640 px canvas. App icons use a 512 px canvas;
the maskable artwork has a full-bleed background and keeps the micro R
inside the central safe area. Prepare the color profile and separations
required by the destination when exporting for print.
