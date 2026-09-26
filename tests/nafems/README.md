# NAFEMS verification benchmarks

This directory is the home for Rodin's implementations of published NAFEMS
verification problems and the accompanying [Rodin Verification Manual](VerificationManual.md).
These are **numerical accuracy benchmarks**, not the timing benchmarks in
`tests/benchmarks/`, and they complement rather than replace the
manufactured-solution and convergence suites.

The [case catalog](cases/README.md) records candidates and their status. No
NAFEMS case is implemented or certified by this scaffold. Cases should be
added only after the authoritative specification and permitted reference data
have been identified and checked. Do not infer a target value from a secondary
example or copy a publication's text, drawings, or meshes into this tree.

For each implemented case, keep a case-specific README with the source and
revision, idealization, units, material and load data, boundary conditions,
quantity-of-interest location and extraction rule, mesh/element family,
reference value, acceptance tolerance and its rationale, and known
limitations. Keep the executable, mesh-generation input, and any small
redistributable fixtures with that README. The manual should then link to the
case, record actual results, and distinguish a passing implementation from a
proposed one.

When tests are added, give them a separate `nafems` CTest label and an explicit
build option; do not silently add expensive benchmark solves to the regular
unit or convergence jobs.
