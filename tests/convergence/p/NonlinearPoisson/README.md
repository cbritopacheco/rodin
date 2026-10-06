# Semilinear Poisson p-convergence

For $\Omega=(0,1)^d$, $d\in\lbrace 1,2,3\rbrace$, the problem is

$$
-\Delta u+u+u^3=f,\qquad u\rvert_{\partial\Omega}=0,
\qquad u=\prod_{j=0}^{d-1}\sin(\pi x_j).
$$

The source $f=(d\pi^2+1)u+u^3$ and the exact gradient are supplied by
`NonlinearPoissonData`, shared with the native h/hp and PETSc/SNES h suites.
The native algorithm is provided by `NonlinearPoissonProblem`. The reaction
derivative $1+3u^2\ge1$ yields a strongly monotone weak operator. The
analytic reference and regular unit-box meshes support degree-dependent
approximation improvement; a finite degree study does not establish an
asymptotic exponential theorem for nonlinear problems.

## Degree path and measurements

The physical mesh is fixed at `n=3` points per axis, $h=1/2$, while scalar
H1 degree increases through $p=1,2,3,4$. Segment, triangle, quadrilateral,
tetrahedron, pyramid, hexahedron, and wedge geometries are included.
Independent physical-cell quadrature measures

$$
E_0=\lVert u-u_p\rVert_{L^2(\Omega)},\qquad
E_1=\lVert\nabla u-\nabla u_p\rVert_{L^2(\Omega)}.
$$

Every adjacent interval requires finite positive errors, strict reduction,
and $\log(E_{\ell,p-1}/E_{\ell,p})>0.1$ for $\ell=0,1$. Four degrees
provide three measured intervals; no fixed-degree h order is asserted.

## Newton, tangent, and negative controls

Newton iteration begins at zero, uses SparseLU for each correction, and
stops at an absolute or relative assembled residual tolerance of $10^{-11}$,
with a limit of 20 iterations. Solver success, nonlinear convergence, and a
finite final residual are required. The tangent is

$$
J(u)[w;v]=(\nabla w,\nabla v)+((1+3u^2)w,v).
$$

At a nonzero interpolated state and direction, degrees 1–4 compare the
assembled tangent action with a central finite difference of the assembled
residual, using $\varepsilon=10^{-5}$. The relative action defect must be
below $10^{-6}$. A P4 negative control substitutes $1+u^2$ for $1+3u^2$
in the tangent and must produce a defect above $10^{-3}$.

A separate P4 control at `n=3` uses $u=4\phi$, where
$\phi=\prod_j\sin(\pi x_j)$, and derives its source with the same cubic
reaction. Omitting $u^3$ from both residual and tangent while retaining this
source gives a linear problem with solution $w$. Since
$(-\Delta+1)(w-u)=u^3$, its first sine-mode coefficient is

$$
c_d=\frac{4^3(3/4)^d}{d\pi^2+1}.
$$

Orthogonality gives $\lVert w-u\rVert_{L^2}\ge c_d2^{-d/2}$ and
$\lVert\nabla(w-u)\rVert_{L^2}\ge\pi\sqrt d\thinspace c_d2^{-d/2}$.
The smallest bounds over $d=1,2,3$ exceed $0.31$ and $1.69$, respectively.
The discrete incorrect solution must exceed $0.05$ in L2 and $0.2$ in H1
seminorm; the correct formulation on the same mesh must give both errors
below these bounds. This separates an analytically known missing-reaction
response from approximation error. Algebraic convergence alone is not the
acceptance oracle. The p-rate and tangent studies retain amplitude one.

Assembly and norm quadrature use order 16. The P4 sensitivity test raises
the order to 18 and tightens both Newton tolerances to $10^{-12}$; each
error must change by less than $10^{-6}$ relative to baseline.

## Execution scope

Native Eigen sequential/OpenMP assembly is selected by
`RODIN_MULTITHREADED` and tested in the local convergence CI matrix. Geometry
entries have `convergence;slow` labels and 600-second budgets. Run
`ctest --test-dir build/tests -R RodinConvergencePNonlinearPoisson --output-on-failure`.
PETSc/SNES, MPI, nonzero or mixed boundary data, curved maps, and arbitrary
degrees beyond four are not certified here. Point/0D has no spatial p-study.
