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

P2 is not yet certified here. An exploratory three-level P2 run on
tetrahedra with $n=7,9,11$ yielded different errors at $n=9$ for one and
three ranks: respectively $E_{L^2}=0.000704082$ versus $0.000972264$ and
$E_{H^1}=0.0449821$ versus $0.0533736$. Errors at $n=7$ and $n=11$
agreed across those rank counts, and the owned-cell quadrature summed to
unit volume at every level. At $n=9$, interpolation of the exact field has
identical errors on one and three ranks, and the unconstrained stiffness and
load have identical Frobenius and Euclidean norms, respectively. Applying
the Dirichlet condition changes the system: 1,538 boundary DOFs are
constrained on one rank but only 1,534 DOFs are constrained on their owning
ranks in the three-rank run. The distributed boundary assembler visits
locally owned boundary faces, whereas PETSc eliminates only DOFs listed by
their owning ranks. The missing owned constraints are consistent with a
boundary-face/DOF ownership mismatch; the exact affected indices and a
corrective exchange have not yet been established. Distributed P2 rates
must not be certified until that mismatch is corrected and re-tested.
