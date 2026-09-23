# How to build and run the tests

## Unit tests

## Numerical convergence tests

`convergence/` contains multi-resolution validation of theoretical finite-
element rates. Enable it with `RODIN_BUILD_CONVERGENCE_TESTS`; tests carry the
`convergence` CTest label. See `convergence/README.md` for the refinement axes
and acceptance rules.

## Benchmarks

`nafems/` is reserved for published numerical verification cases and the
Rodin Verification Manual. It is distinct from `benchmarks/`, which measures
performance. The NAFEMS case catalog currently lists candidates only; no
NAFEMS result is claimed yet.
