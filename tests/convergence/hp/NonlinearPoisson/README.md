# Semilinear Poisson hp-convergence

The continuous equation, analytic manufactured field, independent error
norms, and native Newton implementation are the same as the
[p suite](../../p/NonlinearPoisson/README.md). The fixed problem is
$-\Delta u+u+u^3=f$, with homogeneous Dirichlet trace and
$u=\prod_j\sin(\pi x_j)$ on the unit box.

## Combined refinement path

Both mesh size and H1 degree change along

$$
(n,p)=(2,1)\to(3,2)\to(5,3),\qquad
(h,p)=(1,1)\to(1/2,2)\to(1/4,3).
$$

Both adjacent intervals require finite positive errors, strict reduction,
and effective path orders $r_0>1.9$ and $r_1>0.9$, with

$$
r_{\ell,i}=\frac{\log(E_{\ell,i-1}/E_{\ell,i})}
{\log(h_{i-1}/h_i)},\qquad \ell\in\lbrace 0,1\rbrace .
$$

These are improvement criteria for this particular combined path, not
fixed-degree h orders or a universal hp decay theorem. All seven
positive-dimensional UniformGrid geometries are covered. On a coarsest P1
mesh with no free interior DOFs, the discrete homogeneous-trace field is
zero; this valid coarse approximation is not evidence of a nonlinear update.
The subsequent levels exercise nontrivial Newton iteration.

## Controls and numerical budget

Assembly and independent norm integration use order 16; Newton starts at
zero, uses SparseLU, and requires convergence within 20 iterations at
absolute or relative residual tolerance $10^{-11}$.
At `n=3`, P1–P3 tangent actions are checked against central differences with
$\varepsilon=10^{-5}$ and relative defect below $10^{-6}$.
A P3 tangent with the incorrect derivative $1+u^2$ must give defect above
$10^{-3}$. Omitting the cubic reaction while retaining the exact source must
give P3 errors above $0.05$ in L2 and $0.2$ in H1 seminorm at `n=5`, using
the amplitude-4 manufactured control and first-mode lower bounds derived
in the p specification.
The correct formulation on that same mesh must give both errors below these
bounds, isolating the missing physics from approximation error.
A P3 sensitivity check increases quadrature to 18 and tightens both Newton
tolerances to $10^{-12}$, requiring relative changes below $10^{-6}$ in
both errors.

Native Eigen sequential/OpenMP configurations are tested by the local
convergence CI matrix. Entries have `convergence;slow` labels and 600-second
budgets. Run
`ctest --test-dir build/tests -R RodinConvergenceHPNonlinearPoisson --output-on-failure`.
PETSc/SNES, MPI hp, curved geometry, nonzero/mixed boundaries, and degree
families beyond this path remain outside the claims. Point/0D is excluded.
