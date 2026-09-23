# PETSc-backed Poisson h-convergence

This suite runs the same smooth sine-product Poisson problem as the local
Eigen-backed suite, but assembles the weak form into PETSc matrices and vectors
and solves with PETSc CG. It checks independently integrated L2 and H1 errors
and the P1 rates two and one on all seven UniformGrid cell geometries. Its
purpose is backend numerical parity, not another fixed-resolution residual
check. PETSc is initialized once per executable and the mesh context remains
local; distributed PETSc convergence requires a separate global-norm suite.
