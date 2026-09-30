# Combined hp refinement for variable conductivity

The coefficient and manufactured field are

$$
\gamma(x)=1+\sum_{j=1}^{d}x_j,\qquad
u(x)=\exp\!\left(\sum_{j=1}^{d}x_j\right),\qquad
f=-\nabla\cdot(\gamma\nabla u)=-d(1+\gamma)u.
$$

The exact trace is imposed on the whole boundary of the unit box. The same
equation and data are solved on the combined path
$(n,p)=(2,1)\to(3,2)\to(5,3)$, where $n$ is the number of grid points
per coordinate axis and $h=1/(n-1)$. Thus $h$ takes the values
$1,1/2,1/4$ while the polynomial degree also increases.

The independently integrated errors
$E_{0,i}=\lVert u-u_i\rVert_{L^2(\Omega)}$ and
$E_{1,i}=|u-u_i|_{H^1(\Omega)}$ must strictly decrease. Each adjacent
ratio $\log(E_{j,i-1}/E_{j,i})/\log(h_{i-1}/h_i)$ must exceed $1.9$ in
$L^2$ and $0.9$ in the $H^1$ seminorm. These conservative bounds compare
the combined path with the minimum P1 h-only orders; the ratios are not
fixed-degree h convergence rates. The test covers all seven `UniformGrid`
geometries. Assembly and independent norm integration use order-twelve
quadrature; CG uses relative residual tolerance $10^{-13}$ and at most
20000 iterations.

As a negative control, replacing the spatially varying conductivity in the
assembled stiffness by $1$ while retaining the manufactured forcing caused
all seven geometry-parameterized assertions to fail. The correct
$\gamma(x)$ was restored after that check.
