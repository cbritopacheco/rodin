# p-convergence

Tests in this directory keep the mesh fixed and increase the polynomial degree
`K` of the representing space. For a solution with finite Sobolev regularity,
the error decreases algebraically with `K`; for analytic solutions it decreases
exponentially until conditioning and floating-point error dominate. Each suite
states the solution regularity, norm, and expected degree dependence. The first
suite, `Poisson`, validates exact polynomial reproduction and successive
degree-one through degree-four decay for an analytic manufactured solution.
`Helmholtz` tests the same degree sequence for a complex plane wave. Both
rate studies cover all seven geometries.
`Conductivity` checks a variable diffusion coefficient through exact
degree-two polynomial reproduction and successive analytic-field degree
decay on the same geometry set.
`ReactionDiffusion` checks both fields in a symmetric coupled system on the
same fixed mesh, with separate L2 and H1-seminorm errors at degrees one
through four for all seven geometries.
