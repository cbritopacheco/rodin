# Semilinear Poisson h-convergence

The equation is $-\Delta u+u+u^3=f$, with
$u=\prod_i\sin(\pi x_i)$ and
$f=(d\pi^2+1)u+u^3$. Homogeneous Dirichlet data match the exact field.
Newton's tangent is $-\Delta+1+3u_k^2$; the solve starts at zero, so this
suite exercises nonlinear iteration rather than merely assembling a tangent
at an already converged state.

`NonlinearPoissonProblem` supplies the exact fields and native Newton solve
shared with the p/hp suites. This h suite uses quadrature order 12,
absolute and relative nonlinear tolerances $10^{-11}$, a limit of 20
iterations, and SparseLU corrections. The [p suite](../../p/NonlinearPoisson/README.md)
additionally checks residual/tangent consistency at degrees 1–4 and rejects
an incorrect derivative and an omitted cubic reaction.

The monotone cubic reaction and positive linear reaction give a unique
solution. For smooth data, conforming degree $K$ should attain L2 order
$K+1$ and H1-seminorm order $K$. P1 and P2 are checked on three grids
and on all seven UniformGrid cell geometries. Both nonlinear convergence and
the independent error norms must pass.
