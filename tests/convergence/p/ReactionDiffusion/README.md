# Coupled reaction–diffusion p-convergence

On the fixed unit-box mesh, the two real fields solve

$$
-\Delta u+u+\alpha w=f_u,\qquad
-\Delta w+w+\alpha u=f_w,\qquad \alpha=0.2,
$$

with Dirichlet traces from the exact fields. In spatial dimension $d$, let
$s(x)=\sum_{j=1}^{d}x_j$, $u_*(x)=e^{s(x)}$ and $w_*(x)=e^{-s(x)}$. Then
$\nabla u_*=u_*(1,\ldots,1)$, $\nabla w_*=-w_*(1,\ldots,1)$ and the
manufactured loads are

$$
f_u=(1-d)u_*+\alpha w_*,\qquad
f_w=(1-d)w_*+\alpha u_*.
$$

Both fields use conforming `H1<K>` on the *same* mesh with two grid points
per axis. Degrees $K=1,2,3,4$ are solved, giving three successive p intervals
for each field. Independently integrated $L^2$ and $H^1$-seminorm errors must
decrease on every interval, with observed $\log(E_{K-1}/E_K)>0.25$. This
finite-range bound is evidence of degree decay for these analytic fields,
not a claim of an asymptotic constant. Assembly and error norms use
quadrature order 16; CG uses relative tolerance $10^{-13}$ and at most
20,000 iterations. Segment, triangle, quadrilateral, tetrahedron, pyramid,
hexahedron and wedge are all tested.

A second degree path retains $\alpha=0.2$ and instead uses the shared
`ReactionDiffusionData` variant
$\kappa=(1,2)$, $u_0=e^s$, $u_1=2e^{-s}$, and
$f_i=(1-d\kappa_i)u_i+\alpha u_{1-i}$. This matches the h/hp and
real-PETSc degree suites without replacing the original equal-diffusion
problem. It uses the same $n=2$, degrees and per-interval decay floors,
with assembly order 16 and independent norm order 18.
`FieldConvergence` checks both norms of both fields on every interval.

As a negative control, replacing the first load's coupling contribution
$\alpha w_*$ by $20\alpha w_*$ caused the adjacent-degree assertions to
fail on six of the seven geometries. The correct manufactured load was
restored before the reported passing run.
