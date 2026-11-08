# Scientific writing norms

Scientific prose should use a formal, impersonal, explanatory applied-
mathematics style. Geometric precision, visible logical progression, and
attention to numerical realization should be combined. The intended reader is
mathematically trained but may not know the construction; necessary reasoning
should therefore be explicit without replacing technical content with analogy.

## Voice and claims

- Prefer impersonal narration and passive constructions for authorial actions,
  modelling choices, derivations, and computational procedures. Direct
  statements remain appropriate when mathematical objects are the subjects.
- Do not use first-person claims, direct address, empty phrases such as “it can
  be seen,” or promotional language.
- Preserve the distinction between definitions, assumptions, proved results,
  formal derivations, approximations, and numerical observations. State the
  scope and hypotheses where they affect the claim.
- Do not turn an assumption into an unconditional fact, or a numerical
  similarity into equivalence, validation, or physical agreement.
- Keep attribution visible. A citation must support the statement to which it
  is attached.

## Functional tone

State the mathematical object, operation, conditions and result directly.
Use concise, explanatory prose and sentence-case headings. Repeated contrasts
of the form “this is X, not Y”, conversational warnings and commentary about
what the page does not claim should be replaced by a positive definition,
an explicit scope or an actionable requirement. Keep necessary mathematical
negations, assumptions, failure conditions and supported-backend limits.
State a qualification once, where it affects the construction or result.

Prefer connected prose and reserve colons for lists. Write clarifications,
definitions and consequences as complete sentences. Introduce code and equations
with grammatical prose rather than a label followed by a colon. Use bullet
points for genuinely parallel items. Keep colons required by mathematical
notation, qualified identifiers, URLs, citations and literal program output.
This rule applies to guides, tutorials, captions, headings and API descriptions.


## Exposition

Organize arguments around the mathematical need:

    setting → difficulty → consequence → construction → interpretation

Introduce the relevant objects and operations before explaining the obstruction
to a direct treatment. Introduce a construction as a response to that
obstruction, then state what it permits and what remains unresolved. Choose the
level of abstraction that exposes the mechanism with the least unnecessary
machinery. In extensions, distinguish unchanged structure from genuinely new
ingredients.

Introduce every object with its role and type. State domains, codomains,
spaces, geometry, boundary partitions, configurations, parameters, and
dependencies when they matter. Distinguish geometry from parametrization,
reference from current configurations, continuous fields from discrete
approximations, and local from global finite-element spaces. Keep notation and
terminology stable.

Displayed equations should be part of complete sentences. State their purpose
and explain the non-obvious consequence rather than repeating their symbols.
Explain transitions that use assumptions, cancellations, changes of variables,
derivative transfers, or elimination. Make the differentiation variable and
quantities held fixed explicit; absent explicit dependence is not necessarily
independence through a coupled state.

## Numerical realization and results

Connect continuous constructions to what is actually evaluated. Identify the
discrete spaces, representations, principal algorithmic steps, parameters,
tolerances, residuals, normalizations, units, and stopping conditions when
relevant. Distinguish continuous properties from behavior after
discretization, interpolation, quadrature, or remeshing.

Each experiment should answer an identifiable question. State what changes,
what is held fixed, the measured quantity, reference, sign convention,
normalization, and remaining confounders. Define norms on their actual spaces
and domains, and distinguish a norm of a difference from a difference of norms.
Separate observation from interpretation; explanations should be identified as
such unless the evidence establishes causation.

## Doxygen layout

Documentation comments place the brief on its own line, separate the longer
description with a blank line, and use Markdown section headings separated
from the surrounding text by blank lines. Use this layout for new and revised
documentation throughout the project:

```cpp
/**
 * @brief A concise description.
 *
 * The full description introduces the mathematical object and its role.
 *
 * ## Architecture
 *
 * The implementation stages and their invariants are described here.
 */
```

Keep focused code snippets beside their explanations. Link complete programs
to their generated example source pages with `@ref`; avoid repeating the full
program at the end of a tutorial.

Retain Doxygen mathematical delimiters for rendered expressions and explicit
links to public concepts. A layout change must not discard assumptions,
contracts, mathematical content, or existing cross-references.

## Figure typography

Use the site's bundled Manrope font for figure titles, captions, axes and
legends. Write mathematical expressions in TeX, including labels inside
figures. Mathematical glyphs use the Computer Modern LaTeX font; prose retains
Manrope. Use high-contrast lettering adapted to the light and dark themes. Use `examples/plot_style.py` for Matplotlib figures and the shared
ParaView renderer's Manrope font file for scientific field plots. Keep
backgrounds transparent, include mesh edges when they clarify the field,
and reserve space for labels and color bars. Use Turbo for solution fields
and numerical plots, including curve colors sampled from Turbo. Retain common
field ranges in comparisons. Categorical regions and mesh outlines keep their
consistent categorical colors.

Export comparisons as individual panels and group related figures with
`rodin-mesh-stages`. Use two panels per desktop row, center the final panel
when the count is odd, and stack panels on narrow screens. Frame each pair
together, including its captions. Standalone figures can use up to 36rem of
the content width. Trim empty figure margins so that plots and meshes use
the available space. Keep each caption with its panel.
Plot P0 fields using one value and flat color per original cell; preserve the
point/cell data association of all exported fields and the cell-local samples
of discontinuous higher-order fields.

Documentation videos must be encoded at **30 or 60 fps**, with **30 fps** as
the default. Preserve simulation-time playback when interpolating or resampling
recorded snapshots. Keep transparent backgrounds with alpha-capable WebM/VP9.

## Editing checklist

Before completion, verify that the main question and construction are visible,
the need for each substantial construction is explained, non-obvious steps are
justified, qualifications survive summaries and conclusions, comparisons have
clear measures and references, notation is consistent, and no unsupported
content has been introduced. Preserve mathematical meaning, assumptions,
attribution, and certainty during revision; flag substantive omissions instead
of inventing missing information.
