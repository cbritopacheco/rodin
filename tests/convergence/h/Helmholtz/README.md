# Complex Helmholtz h-convergence

This suite solves the complex-valued Helmholtz equation

\[
  -\Delta u-\kappa^2u=f\quad\text{in }\Omega,
  \qquad u=g\quad\text{on }\partial\Omega
\]

on uniform meshes of the unit interval, square, and cube. The manufactured
plane wave

\[
  u(x)=\exp\!\left(i\sum_{j=1}^d x_j\right)
\]

satisfies

\[
  f=(d-\kappa^2)u,
  \qquad \nabla u=(iu,\ldots,iu).
\]

It is analytic and has nonzero real and imaginary parts on the boundary. For
conforming elements of degree `K`, standard approximation theory predicts

\[
  \lVert u-u_h\rVert_{H^1(\Omega)}=O(h^K),
  \qquad
  \lVert u-u_h\rVert_{L^2(\Omega)}=O(h^{K+1}).
\]

The tests verify rates `(2, 1)` for P1 and `(3, 2)` for P2 in `(L2, H1
seminorm)`. Errors are integrated independently with the complex modulus at
quadrature order twelve. The wavenumber is kept below the first Dirichlet
eigenfrequency so the homogeneous discrete operator remains coercive. All
seven positive-dimensional `UniformGrid` cell geometries are covered.
