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

## Editing checklist

Before completion, verify that the main question and construction are visible,
the need for each substantial construction is explained, non-obvious steps are
justified, qualifications survive summaries and conclusions, comparisons have
clear measures and references, notation is consistent, and no unsupported
content has been introduced. Preserve mathematical meaning, assumptions,
attribution, and certainty during revision; flag substantive omissions instead
of inventing missing information.
