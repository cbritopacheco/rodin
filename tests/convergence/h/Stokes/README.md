# Stokes h-convergence

This suite validates the Taylor--Hood discretization of the steady
incompressible Stokes equations,

\[
  -\Delta u + \nabla p = f, \qquad \nabla\cdot u = 0.
\]

The velocity is represented in vector-valued `H1<2>` and the pressure in
scalar `H1<1>`. A `P0g` Lagrange multiplier fixes the pressure mean, giving a
well-posed saddle-point system. For sufficiently smooth solutions and a
stable Taylor--Hood pair, the expected rates are

\[
  \lVert u-u_h\rVert_{H^1}=O(h^2),\quad
  \lVert u-u_h\rVert_{L^2}=O(h^3),\quad
  \lVert p-p_h\rVert_{H^1}=O(h),\quad
  \lVert p-p_h\rVert_{L^2}=O(h^2).
\]

An affine divergence-free velocity with affine zero-mean pressure is first
reproduced to roundoff. The rate study then uses
\(u=(x_1^3,0,\ldots,0)\) and \(p=x_0^2-\tfrac13\), with the manufactured
forcing \(f=(-6x_1+2x_0,0,\ldots,0)\). The velocity is exactly
divergence-free and the pressure has zero mean on the unit box.

The affine patch test additionally checks the L2 divergence. The rate study
intentionally does not require monotone divergence:
Taylor--Hood enforces incompressibility only weakly against the pressure
space, and the unresolved component of \(\nabla\cdot u_h\) has no monotonic
L2 convergence theorem.

Triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge grids
are covered. Segment is excluded because incompressible Stokes has no
non-degenerate velocity-pressure formulation in one spatial dimension.
