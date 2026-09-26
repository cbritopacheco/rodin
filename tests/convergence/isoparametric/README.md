# Isoparametric convergence

Tests in this directory validate finite-element approximation when the physical
cell map is itself a finite-element field. Geometry and field spaces are
refined in a controlled way so that a geometry-map error cannot be mistaken for
a solution-space error.

`CurvedGeometry.cpp` uses the regular quadratic map

$$
  \Phi(\xi)_i=\xi_i\quad(i\lt d-1),\qquad
  \Phi(\xi)_{d-1}=\xi_{d-1}+0.1\xi_0^2.
$$

It installs P2 `ParametricTransformation`s on every positive-dimensional mesh
entity, including boundary faces. The map is exactly in the P2 geometry space;
the geometric patch test verifies that fact directly for every `UniformGrid`
cell type. A separate three-level study installs the P1 vertex interpolant
$\Phi_{1,h}$ and measures its map error on the *unwarped* unit box:

$$
G_h=\left(\int_{(0,1)^d}\lVert\Phi-\Phi_{1,h}\rVert^2\,\mathrm{d}\xi\right)^{1/2}.
$$

Both adjacent rates are required to lie between 1.7 and 2.3, consistent
with $G_h=O(h^2)$ for this regular quadratic map. The reference mesh is
retained separately from the curved mesh, so the norm is not contaminated
by a changing integration domain. P1 geometry uses $n=5,9,17$ points per
axis in 1D/2D and $n=3,5,9$ in 3D. A manufactured Poisson solution in
physical coordinates, $u(x)=\prod_i\sin(\pi x_i)$, then checks that P2
isoparametric fields retain
the mapped-element rates $O(h^3)$ in L2 and $O(h^2)$ in the H1 seminorm.

The distinction matters: an order-`K` isoparametric map controls the accuracy
of Jacobians, integration, traces, and curved-boundary placement. When the
geometry is represented at compatible order and remains regular, the usual
interpolation and Galerkin estimates apply on the mapped mesh. Later suites
will use non-polynomial boundary maps to measure a geometric-approximation
term that persists at every finite geometry order.

The same P2 map also supports a variable-conductivity study. In physical
coordinates, $\gamma(x)=1+\sum_{j=1}^{d}x_j$ and
$u_*(x)=\exp(\sum_{j=1}^{d}x_j)$ give the manufactured load

$$
-\nabla\cdot(\gamma\nabla u_*)=-d(1+\gamma)u_*.
$$

The exact trace is imposed on the mapped boundary. On all seven cell types,
both adjacent P2 field-error rates must match order three in $L^2$ and
order two in the $H^1$ seminorm. Poisson and conductivity use quadrature
order 16, an LU solve, and $n=5,9,17$ in 1D/2D or $n=3,4,5$ in 3D. The
quadratic map is exact in the P2 geometry space in these PDE studies; their
field-error rates must not be interpreted as convergence rates for an
under-resolved boundary.

As a negative control, replacing $\gamma\nabla u$ by $\nabla u$ in the
conductivity stiffness while retaining the manufactured load caused all
seven geometry rate tests to fail. The variable coefficient was restored
before the passing run.
