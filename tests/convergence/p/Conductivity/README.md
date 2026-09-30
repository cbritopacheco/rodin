# Variable-conductivity p-convergence

On the unit box $\Omega=(0,1)^d$, the diffusion coefficient and equation are

$$
\gamma(x)=1+\sum_{j=1}^{d}x_j,\qquad
-\nabla\cdot(\gamma\nabla u)=f\quad\text{in }\Omega,
$$

with the exact trace imposed on all of $\partial\Omega$. Since
$1\le\gamma\le 1+d$, the bilinear form is coercive. The mesh has two grid
points per coordinate axis and is held fixed while the conforming $H^1$
degree increases.

Two manufactured fields distinguish exact representation from finite-range
degree decay:

1. For $u_*(x)=\sum_j x_j^2$, the source is
   $f_*(x)=-2\{\sum_j x_j+d\gamma(x)\}$. The degree-one solution has
   nonzero $L^2$ and $H^1$-seminorm errors, while the degree-two solution
   reproduces $u_*$ to an absolute error below $10^{-10}$ in both norms.
2. For $u(x)=\exp(\sum_j x_j)$, one has $\nabla u=u(1,\ldots,1)$ and
   $f=-d(1+\gamma)u$. The errors
   $E_{0,K}=\lVert u-u_K\rVert_{L^2(\Omega)}$ and
   $E_{1,K}=|u-u_K|_{H^1(\Omega)}$ are integrated independently of the
   assembled residual. At every adjacent degree $K-1\to K$ for
   $K=2,3,4$, both errors must strictly decrease and satisfy
   $\log(E_{j,K-1}/E_{j,K})>0.25$ for $j=0,1$.

Assembly and error integration use quadrature order twelve. CG uses relative
residual tolerance $10^{-13}$ and at most 20000 iterations. The tests run on
all seven positive-dimensional `UniformGrid` cell geometries. The observed
finite-range degree decay does not establish an asymptotic exponential bound
for arbitrary degree.

As a negative control, replacing $\gamma(x)$ by the constant $1$ in the
assembled stiffness while retaining the manufactured source caused the
adjacent-degree assertions to fail on all seven geometries. The variable
coefficient was restored after this check.
