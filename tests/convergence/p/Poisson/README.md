# Poisson p-convergence

This suite fixes a coarse mesh of the unit interval, square, or cube and raises
the degree `K` of the conforming `H1<K>` space. It validates two complementary
properties.

The polynomial manufactured solution

$$
  u(x)=\sum_{i=1}^d x_i^2, \qquad f=-2d
$$

is not represented exactly at degree one but belongs to every space from degree
two onward. The discrete degree-two solution must therefore reach roundoff-level
L2 and H1-seminorm error.

For the analytic manufactured solution

$$
  u(x)=\exp\!\left(\sum_{i=1}^d x_i\right),
  \qquad f=-\Delta u=-d u,
  \qquad \nabla u=(u,\ldots,u),
$$

the exact trace is imposed on the whole boundary. Analytic regularity permits
exponential best-approximation decay with degree, subject to the fixed mesh,
geometry, quadrature, and floating-point regime. With $u_K$ denoting the
discrete degree-$K$ solution on the fixed mesh, the test measures
$E_{0,K}=\lVert u-u_K\rVert_{L^2(\Omega)}$ and
$E_{1,K}=|u-u_K|_{H^1(\Omega)}$ independently at every degree
$K=1,2,3,4$. For $j\in\{0,1\}$ and each $K=2,3,4$, it computes

$$
  \alpha_{j,K}=\log\!\left(\frac{E_{j,K-1}}{E_{j,K}}\right),
$$

and requires strict error reduction and $\alpha_{j,K}>0.25$ in both norms.
The asymmetric exponential avoids the symmetry-induced plateau of the earlier
sine-product field between degrees two and three on some cell families; that
plateau made an adjacent-degree rate test inconclusive. The fixed mesh has two
grid points per axis. Assembly and error integration both use quadrature order
twelve, and CG uses relative tolerance $10^{-13}$ with at most 20000
iterations. All seven positive-dimensional `UniformGrid` cell geometries are
covered. These finite-range assertions do not prove exponential decay for
arbitrarily high degree.

As a negative control, setting $f=0$ while retaining the exact exponential
trace caused the adjacent-degree assertions to fail on all seven geometries.
The manufactured forcing $f=-d u$ was restored after that check.
