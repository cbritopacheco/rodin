# Rodin Verification Manual

## Scope and status

This is the source-controlled manual for Rodin's NAFEMS benchmark results.
It is currently a framework for future case reports: **no NAFEMS benchmark
result or certification is claimed here yet**. The independently specified
NAFEMS problems will test the complete modelling and solution workflow,
whereas `tests/convergence/` chiefly checks known finite-element rates for
manufactured fields.

The case inventory and implementation status live in [cases/README.md](cases/README.md).
Each completed result must link to its runnable test and documented input.

## Verification protocol

For each case, record:

1. The authoritative NAFEMS publication, case identifier, edition/revision,
   and any published corrections. Preserve provenance for every reference
   quantity and unit; do not reproduce copyrighted case sheets.
2. The continuum model and Rodin weak form, including dimensional reduction,
   sign conventions, constitutive law, boundary conditions and load
   application. State deliberate differences from the published idealization.
3. The geometry and mesh construction, element family and order, quadrature,
   solver and tolerances, backend, and Rodin commit/configuration needed to
   reproduce the result. Archive small, redistributable input fixtures in the
   case directory; provide a reproducible generator for anything larger.
4. The exact quantity of interest (QoI): component, coordinate, local frame,
   averaging or extrapolation rule, and units. A point stress and a nodal
   extrapolated stress are not interchangeable measurements.
5. At least a coarse, medium and fine discretization when refinement is
   meaningful. Report raw QoIs, relative deviation from the cited reference,
   and observed stabilization. Explain non-monotone results or dependence on
   element type; do not hide them behind a single passing tolerance.
6. A case-specific acceptance rule justified by the publication, reference
   uncertainty, singularity/regularity, and discretization error. Verify that
   a deliberately perturbed load, material parameter, or extraction point
   makes the automated check fail.

For a reference quantity $q_{\mathrm{ref}} \ne 0$, report the relative
deviation $|q_h-q_{\mathrm{ref}}|/|q_{\mathrm{ref}}|$. When the reference is zero
or near zero, use an explicitly stated absolute tolerance with the same units
as the QoI. Keep solver residual tolerance below the error budget being tested.

## Results register

| Case | Physics and QoI | Specification | Rodin implementation | Result |
| --- | --- | --- | --- | --- |
| None yet | — | — | — | Not evaluated |

Passing one benchmark supports only the documented model, element family,
geometry and backend. It is not a general certification of Rodin or of all
finite-element spaces. Claims should be expanded only as independent cases
and configurations are actually exercised.
