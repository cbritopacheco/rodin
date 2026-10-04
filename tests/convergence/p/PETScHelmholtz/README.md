# Complex PETSc Helmholtz degree refinement

On $\Omega=(0,1)^d$, the complex field satisfies

$$
-\Delta u-k^2u=f,\qquad k^2=\tfrac14,\qquad
u|_{\partial\Omega}=u_*,\qquad
u_*=e^{\mathrm{i}s},\quad s=\sum_jx_j,\quad f=(d-k^2)u_*.
$$

The complex weak form uses conjugate test functions. With full Dirichlet
traces, $k^2<\lambda_1(\Omega)=d\pi^2$ gives a coercive Hermitian form;
this is not a high-frequency or resonance study. The fixed $n=2$ grid and
degrees $p=1\to2\to3\to4$ match the native plane-wave degree suite.
For each independently integrated $L^2$ and $H^1$-seminorm error $E_p$,
every adjacent interval requires finite positive errors, strict reduction
and $\log(E_{p-1}/E_p)>0.25$. This is a finite-path decay floor, not a
universal exponential approximation constant.

## Architecture and controls

`HelmholtzData` supplies physical complex fields, gradients and sources.
`PETScHelmholtzProblem` constructs a fresh fixed-layout complex space and
system per measurement. `PETScHelmholtzRefinement` checks the degree/path
policy with the common `FieldConvergence` helper. The p/hp suites share
one entry point and the common geometry/rank registration. Norms use the
complex modulus, and MPI sums owned-cell squared errors once globally.

The quadratic patch is
$u_*=1+2\mathrm{i}+(1+\mathrm{i}/2)s^2$, with
$f=-2d(1+\mathrm{i}/2)-k^2u_*$. On $n=2$, degree one must have
$E_{L^2}>10^{-3}$ and $E_{H^1}>10^{-2}$, while degree two reproduces
the field and gradient within $10^{-9}$. Degree four reproduces this
patch separately on $n=3$. A resolved affine control
$u_*=1+2\mathrm{i}+(1+\mathrm{i}/2)s$ first satisfies the same exact-patch
budget; omitting the mass term while retaining the correct source and
trace must exceed both separation floors. Both real and imaginary parts
have nontrivial manufactured data.

On $n=3$ at degree four, assembly order $16\to18$, norm order $18\to20$,
and CG relative tolerance $10^{-13}\to10^{-14}$ vary separately.
Each error must change by less than $10^{-6}$ relatively. CG requires a
positive convergence reason, a finite reported residual below $10^{-9}$,
and an independent residual
$\|Ax-b\|_2/\max(1,\|b\|_2)<10^{-11}$.

All seven positive-dimensional cell families run with local meshes and
MPI ranks one through four, under sequential and OpenMP assembly. The
suite requires a native-complex PETSc installation; real PETSc does not
register it. A point has no nonconstant spatial degree rate. Tests are
slow-labelled, have 30-minute safety timeouts, and serialize pyramid
registrations through the shared degree-refinement resource lock.
