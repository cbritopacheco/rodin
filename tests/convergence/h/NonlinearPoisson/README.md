# Semilinear Poisson h-convergence

The equation is $-\Delta u+u+u^3=f$, with
$u=\prod_i\sin(\pi x_i)$ and
$f=(d\pi^2+1)u+u^3$. Homogeneous Dirichlet data match the exact field.
Newton's tangent is $-\Delta+1+3u_k^2$; the solve starts at zero, so this
suite exercises nonlinear iteration rather than merely assembling a tangent
at an already converged state.

The monotone cubic reaction and positive linear reaction give a unique
solution. For smooth data, conforming degree $K$ should attain L2 order
$K+1$ and H1-seminorm order $K$. P1 and P2 are checked on three grids
and on all seven UniformGrid cell geometries. Both nonlinear convergence and
the independent error norms must pass.
