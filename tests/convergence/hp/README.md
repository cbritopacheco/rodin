# hp-convergence

Tests in this directory refine the mesh and increase the representing-space
degree together. The expected rate depends on the refinement strategy and the
regularity distribution of the exact solution. Each suite documents its
strategy and separately measured norms. The first suite, `Poisson`, combines
three uniform-grid levels with degrees one through three for an analytic
manufactured solution on every cell geometry. `Helmholtz` follows the same
combined path for a complex analytic plane wave, checking both L2 and H1
errors on all seven geometries.
`Conductivity` uses the same path for a variable diffusion coefficient and
an analytic nonhomogeneous Dirichlet field.
