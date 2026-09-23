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

For the analytic homogeneous solution

$$
  u(x)=\prod_{i=1}^d\sin(\pi x_i), \qquad f=d\pi^2u,
$$

spectral approximation theory predicts exponential error decay
`e_K <= C exp(-alpha K)`. The test computes

$$
  \alpha_{1,4}=\frac{\log(e_1/e_4)}{4-1}
$$

in both L2 and the H1 seminorm, and requires positive decay. A degree-one to
degree-four interval is deliberately used: it is a non-degenerate fixed-mesh
enrichment for every supported cell family. Errors are integrated independently
at order twelve. All seven
positive-dimensional `UniformGrid` cell geometries are covered.
