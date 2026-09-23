# h-convergence

Tests in this directory keep the finite-element degree fixed and successively
halve the nominal coordinate spacing $h$. For a sufficiently regular solution
of a stable, consistent conforming degree $K$ problem, energy-norm estimates
and, for the L2 bound, the required dual regularity give

$$
  \lVert u-u_h\rVert_{H^1(\Omega)} = O(h^K), \qquad
  \lVert u-u_h\rVert_{L^2(\Omega)} = O(h^{K+1}).
$$

Every test records errors on at least three meshes, verifies monotone reduction,
and checks observed rates. Geometry-parameterized tests cover every cell type to
exercise reference-to-physical mappings and geometry-specific quadrature.

`Convergence.h` provides the h-specific mesh hierarchy:

- `UniformGridHierarchy` creates scaled unit-box mesh families and exposes the
  nominal coordinate spacing of every level.

The parent `tests/convergence/Convergence.h` supplies `UniformGrid`, `ErrorNorm`,
`ErrorHistory`, `ErrorNorms`, and `Rates` for both h- and p-convergence suites.

The `Poisson` suite verifies real P1, P2, and P3 convergence, including exact
affine reproduction, essential data, natural fluxes, Robin terms, and the
mean-constrained pure-Neumann problem. The `Helmholtz` suite verifies complex
plane-wave convergence with both P1 and P2 elements.

The `LinearElasticity` suite verifies vector P1 and P2 approximation, affine
patch exactness, full displacement and mixed traction boundaries, and a
divergence-free nearly incompressible regime.

The `Stokes` suite validates the mixed Taylor--Hood velocity-pressure pair,
including the pressure gauge, saddle-point coupling, and discrete
incompressibility.

`Conductivity` validates variable-coefficient diffusion with an exact affine
patch and P1/P2 smooth-field rates. `P0Projection` validates first-order L2
projection for real/complex scalar/vector piecewise constants; it also checks
exact reproduction of global constants by P0g.

`ReactionDiffusion` validates two coupled scalar fields, their separate P1/P2
rates, and a nonhomogeneous Dirichlet component on all cell geometries.
`NonlinearPoisson` validates a semilinear reaction law through a true Newton
solve and P1/P2 rates on the same geometry set.
When PETSc is enabled, `PETScPoisson` certifies local-context PETSc assembly
and CG through all-geometry P1 L2/H1 convergence rates.
When both MPI and PETSc are enabled, `PETScMPIPoisson` checks the same P1 rates
on distributed meshes, with globally reduced errors and one to four ranks.

Equation-specific suites should supply only the manufactured fields, assemble
and solve their weak formulation, then use these classes for measurement.
