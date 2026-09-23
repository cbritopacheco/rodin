# Distributed PETSc Poisson h-convergence

This suite shards a unit-box UniformGrid, solves the P1 sine-product Poisson
problem with distributed PETSc CG, and measures L2 and H1-seminorm errors by
independent quadrature over owned cells. The local squared contributions are
summed across MPI ranks before taking square roots, so ghost cells cannot
double-count the norm. It checks three grid levels, expected rates two and
one, every seven cell geometries, and rank counts one through four.

The levels are 5, 9, and 17 points per axis. A two-cell-wide tetrahedral
coarsest grid (3 points per axis) was outside the stable asymptotic regime:
its first L2 ratio was 1.43. The certified hierarchy begins at four cells
per axis, as in the local P1 suite.

This tests distributed mesh partitioning, FE ownership, PETSc assembly and
solver behavior, and global error reduction together. It does not imply that
every other physics suite or space has distributed convergence coverage.
