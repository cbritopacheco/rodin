# Numerical convergence tests

This suite validates asymptotic finite-element approximation properties rather
than accuracy on one fixed mesh. Each test solves the same problem on a nested
sequence of discretizations and measures the observed order

\[
  r_i = \frac{\log(e_{i-1}/e_i)}{\log(h_{i-1}/h_i)}.
\]

The suite is separate from `tests/manufactured`: manufactured tests protect a
particular formulation or regression at a useful fixed resolution, whereas
convergence tests require several resolutions and assert the mathematical rate.

The subdirectories split the four independent refinement mechanisms:

- `h`: decrease the mesh size while keeping the representing space degree fixed;
- `p`: increase the representing-space degree on a fixed mesh;
- `hp`: change mesh size and degree together;
- `isoparametric`: refine the discrete geometry map and verify that geometric
  approximation does not limit the field approximation.

Configure with `-DRODIN_BUILD_CONVERGENCE_TESTS=ON` and run with
`ctest --test-dir build/tests -L convergence --output-on-failure`.

`Convergence.h` contains the refinement-independent layer shared by these
modules: unit-box grid construction, direct L2/H1 error integration, error
histories, and algebraic or exponential rate calculation. Refinement-specific
directories add only the machinery unique to their refinement axis.
