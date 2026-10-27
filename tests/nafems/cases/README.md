# Case catalog

The entries below are **candidates**, not implemented tests or verified
results. NAFEMS' [background to benchmarks](https://www.nafems.org/publications/resource_center/r0006/)
identifies the case names. The authoritative case specifications and their
revisions must be obtained and checked before implementing numerical targets.

| Candidate | Physical context | Intended first check | Status |
| --- | --- | --- | --- |
| LE1 — elliptic membrane | Linear elasticity | Displacement/stress extraction on curved geometry | Specification review pending |
| LE10 — thick plate under normal pressure | 3D linear elasticity | Through-thickness response and stress recovery | Specification review pending |
| T4 — two-dimensional heat transfer with convection | Steady thermal conduction | Conductivity plus convective boundary condition | Specification review pending |

This is a starting order, not a claim that all three cases are immediately
representable by existing Rodin APIs. In particular, confirm geometry,
boundary-condition and QoI-extraction support before choosing the first
executable. Name implemented directories by the published case identifier
(for example, `LE1/`) and keep each case's provenance and results beside its
code.
