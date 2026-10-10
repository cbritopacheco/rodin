# Coupled reaction–diffusion hp-convergence

This suite uses the same two-field system and `ReactionDiffusionData` as the
[PETSc h suite](../../h/PETScReactionDiffusion/README.md), with diffusion
coefficients $(1,2)$, reaction matrix

$$
R=\begin{pmatrix}1&0.2\cr 0.2&1\end{pmatrix},
$$

and full manufactured
Dirichlet traces. The analytic fields are $u_0=e^s$ and $u_1=2e^{-s}$,
where $s=\sum_jx_j$. Sources and gradients are derived independently of the
discrete fields; each component has its own L2 and H1 error history.

## Combined path and interpretation

The sequence is $(n,p)=(2,1)\to(3,2)\to(5,3)$ on all seven
positive-dimensional geometries, with $h=1/(n-1)$:

$$
(h,p)=(1,1)\to(1/2,2)\to(1/4,3).
$$

For each component and both adjacent intervals, errors must be finite,
positive, and strictly decreasing. Effective path rates are

$$
r_{\ell,i}=\frac{\log(E_{\ell,i}/E_{\ell+1,i})}
{\log(h_\ell/h_{\ell+1})},
$$

computed separately for the L2 norm and H1 seminorm. The accepted lower
bounds are $r_0>1.9$ and $r_1>0.9$, matching the existing analytic hp
studies. These are improvement criteria for this specified combined path,
not fixed-degree h-order estimates or a universal exponential hp theorem.
Analytic fields, conforming spaces, regular affine meshes, coercivity,
consistent quadrature, and sufficient solver accuracy are assumed.

## Higher-order controls and budget

At `n=3`, P3 must reproduce each quadratic field
$u_i=(i+1)(1+s^2)$ with L2 and H1-seminorm errors below $10^{-9}$.
Removing both off-diagonal reaction terms while preserving the correct
affine sources and traces must give errors above $10^{-3}$ in L2 and
$10^{-2}$ in H1 for both components. This checks that the high-order oracle
detects a missing physical coupling.

Assembly and independent error quadrature use order 16. Native CG uses
relative tolerance $10^{-13}$ and at most 50,000 iterations, and must
report success. At `n=3`, a P3 sensitivity study raises quadrature to 18
and tightens tolerance to $10^{-14}$; all four measured component errors
must change by less than $10^{-6}$ relative to baseline.

## Execution and exclusions

The suite uses local Eigen sequential/OpenMP assembly according to
`RODIN_MULTITHREADED`; both configurations are included in the existing
convergence CI matrix. All segment, triangle, quadrilateral, tetrahedron,
pyramid, hexahedron, and wedge cases have the `convergence;slow` labels
and a 600-second budget per test. Run
`ctest --test-dir build/tests -R RodinConvergenceHPReactionDiffusion --output-on-failure`.
PETSc/MPI hp, curved geometry, mixed boundaries, nonsymmetric reaction,
and complex fields are not certified by these entries. Point/0D has no
PDE hp-rate here.
