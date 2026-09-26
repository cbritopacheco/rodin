# Variable-conductivity h-convergence

The coefficient is $\gamma(x)=1+\sum_i x_i$, uniformly positive on the
unit box. We solve $-\nabla\cdot(\gamma\nabla u)=f$ with exact Dirichlet
data. An affine exact field tests consistency: its gradient is constant and
the source is $-d$. A smooth sine-product field tests approximation rates;
its manufactured source is

$$
f=d\pi^2\gamma u-\pi\sum_i\cos(\pi x_i)
       \prod_{j\ne i}\sin(\pi x_j).
$$

For conforming degree $K=1,2$, the expected H1-seminorm and L2 rates are
$K$ and $K+1$, respectively. The suite checks both successive rates on
three meshes for every seven UniformGrid cell geometries, using independent
high-order quadrature for the error norms.
