# Distributed PETSc Poisson h-convergence

This suite shards a unit-box UniformGrid, solves the P1/P2 sine-product Poisson
problem with distributed PETSc CG, and measures L2 and H1-seminorm errors by
independent quadrature over owned cells. The local squared contributions are
summed across MPI ranks before taking square roots, so ghost cells cannot
double-count the norm. It checks three grid levels, expected L2/H1 rates
(2, 1) for P1 and (3, 2) for P2, all seven cell geometries, and rank counts
one through four. Each observed P2 L2 rate must lie in (2.4, 3.6), and each
H1-seminorm rate in (1.5, 2.5).

The boundary-identification rate regressions use the same differential
operator and independent error measurement. For P1, the exact field
$u_0(x)=\prod_{j=1}^{d}\sin(\pi x_j)$ satisfies the homogeneous
identification $u_0=-u_0$ on $\partial\Omega$. For P2, the shifted field
$u_1(x)=1+u_0(x)$ satisfies the affine identification $u_1=-u_1+2$
there. Both have source $f=d\pi^2u_0$; the shift changes the trace but not
the gradient or Laplacian. Three-level errors must decrease on every step,
with P1 L2/H1 rates respectively in (1.5, 2.5)/(0.7, 1.5) and P2 rates
in (2.4, 3.6)/(1.5, 2.5). These tests run on all seven cell geometries
and one through four ranks.

The smooth manufactured field is $u(x)=\prod_{j=1}^{d}\sin(\pi x_j)$ on
$\Omega=(0,1)^d$. It satisfies $-\Delta u=d\pi^2u$ and has zero trace on
$\partial\Omega$. With $h=1/(n-1)$, successive error ratios are compared
with $E_{L^2}=O(h^{p+1})$ and $|e|_{H^1}=O(h^p)$ for $p=1,2$. The norm
quadrature is separate from the assembled forms, and CG uses relative and
absolute residual tolerances $10^{-12}$ and $10^{-14}$, respectively.

The P1 levels are 5, 9, and 17 points per axis. A two-cell-wide tetrahedral
coarsest grid (3 points per axis) was outside the stable asymptotic regime:
its first L2 ratio was 1.43. The certified hierarchy begins at four cells
per axis, as in the local P1 suite. P2 uses 5, 7, and 9 points per axis,
except tetrahedra, which use 7, 9, and 11 to include the original failure.

A separate P2 patch test uses the exact quadratic field
$u(x)=1+\sum_{j=1}^{d}x_j^2$, source $f=-\Delta u=-2d$, and nonzero
Dirichlet trace $g=u|_{\partial\Omega}$. Its independently integrated L2
and H1-seminorm errors must each be below $10^{-8}$ for both
value-prescribed $u=g$ and affine-identified $u=-u+2g$ traces. The latter
tests spatially varying affine offsets through the P2 DOF functionals.
The same quadratic patch is solved after full-domain SubMesh extraction through
owned cells only, for both value and affine traces, on every geometry and
1–4 ranks. This fixed-mesh regression checks that finalization restores the
selected overlap and does not turn partition interfaces into Dirichlet boundaries;
it uses the same error norms and tolerances as the parent-mesh patch.
An ownership regression
checks that every boundary DOF found on any rank is also present in the
boundary-condition map on its owning rank. The tetrahedral case uses $n=9$,
where the unfixed implementation omitted four owned constraints on three
ranks; this regression fails with owned-face-only traversal.

This tests distributed mesh partitioning, FE ownership, PETSc assembly and
solver behavior, and global error reduction together. It does not imply that
every other physics suite or space has distributed convergence coverage.

The MPI assembler visits physical faces in the existing vertex-star
overlap and retains constraints for owned DOFs and DOFs used by owned cells.
This restriction excludes artificial outer-halo boundaries. Identification
rows and affine offsets select the same incident face by its smallest
distributed index. No P1/H1 constraint exchange is needed, and no coefficient
threshold or numerical duplicate comparison is used. The repair belongs to
MPI assembly independently of PETSc. Independent backend-free tests check
entity ownership, physical-boundary incidence, and the required constraint
sets across spaces and geometries. The PETSc MPI form
suite additionally solves the tetrahedral $n=9$ Poisson problem twice,
using $u=1$ and the mathematically equivalent $u=-u+2$ on the boundary,
and requires the solution-vector difference to be below $10^{-9}$ in the
Euclidean norm. Its mixed-field tests check nonzero affine offsets and
nontrivial master coefficients on two, three, and four ranks.

MPI affine evaluation uses the same shard-local functional selection as the
linear-row assembler; sequential evaluation uses local iterators. Global entity
counts must not be used as bounds for shard-local indices. The P2 affine solve
regression exercises this contract across partitions.
