# Isoparametric convergence

Tests in this directory validate finite-element approximation when the physical
cell map is itself a finite-element field. Geometry and field spaces are
refined in a controlled way so that a geometry-map error cannot be mistaken for
a solution-space error.

`Poisson/` (implemented in `Poisson.cpp`) uses the regular quadratic map

\[
  \Phi(\xi)_i=\xi_i\quad(i<d-1),\qquad
  \Phi(\xi)_{d-1}=\xi_{d-1}+0.1\xi_0^2.
\]

It installs P2 `ParametricTransformation`s on every positive-dimensional mesh
entity, including boundary faces. The map is exactly in the P2 geometry space;
the geometric patch test verifies that fact directly for every `UniformGrid`
cell type. A manufactured Poisson solution in physical coordinates,
\(u(x)=\prod_i\sin(\pi x_i)\), then checks that P2 isoparametric fields retain
the mapped-element rates \(O(h^3)\) in L2 and \(O(h^2)\) in the H1 seminorm.

The distinction matters: an order-`K` isoparametric map controls the accuracy
of Jacobians, integration, traces, and curved-boundary placement. When the
geometry is represented at compatible order and remains regular, the usual
interpolation and Galerkin estimates apply on the mapped mesh. Later suites
will compare deliberately under-resolved geometry orders and use non-polynomial
boundary maps to measure the separate geometric-approximation term.
