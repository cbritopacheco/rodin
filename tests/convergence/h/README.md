# h-convergence

Tests in this directory keep the finite-element degree fixed and successively
halve the characteristic mesh size `h`. For a smooth solution represented by a
conforming degree-`K` space, Cea's lemma and interpolation estimates give

\[
  \lVert u-u_h\rVert_{H^1(\Omega)} = O(h^K), \qquad
  \lVert u-u_h\rVert_{L^2(\Omega)} = O(h^{K+1}).
\]

Every test records errors on at least three meshes, verifies monotone reduction,
and checks observed rates. Geometry-parameterized tests cover every cell type to
exercise reference-to-physical mappings and geometry-specific quadrature.

`Convergence.h` provides the h-specific mesh hierarchy:

- `UniformGridHierarchy` creates scaled unit-box mesh families and exposes the
  actual mesh size of every level.

The parent `tests/convergence/Convergence.h` supplies `UniformGrid`, `ErrorNorm`,
`ErrorHistory`, `ErrorNorms`, and `Rates` for both h- and p-convergence suites.

The `Poisson` suite verifies real P1 convergence, including exact affine
reproduction and homogeneous or nonhomogeneous Dirichlet data. The `Helmholtz`
suite verifies complex plane-wave convergence with both P1 and P2 elements.

Equation-specific suites should supply only the manufactured fields, assemble
and solve their weak formulation, then use these classes for measurement.
