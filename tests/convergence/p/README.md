# p-convergence

Tests in this directory keep the mesh fixed and increase the polynomial degree
`K` of the representing space. For a solution with finite Sobolev regularity,
the error decreases algebraically with `K`; for analytic solutions it decreases
exponentially until conditioning and floating-point error dominate. Each suite
states the solution regularity, norm, and expected degree dependence. The first
suite, `Poisson`, validates exact polynomial reproduction and exponential decay
for an analytic manufactured solution.
