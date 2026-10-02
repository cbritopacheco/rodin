# Stokes h-convergence

This suite validates the Taylor--Hood discretization of the steady
incompressible Stokes equations,

$$
  -\Delta u + \nabla p = f, \qquad \nabla\cdot u = 0.
$$

The velocity is represented in vector-valued `H1<2>` and the pressure in
scalar `H1<1>`. A `P0g` Lagrange multiplier fixes the pressure mean, giving a
unique pressure gauge. Mixed well-posedness additionally requires discrete
inf-sup stability; fixing the mean does not eliminate other pressure null
modes. For sufficiently smooth solutions and a
stable Taylor--Hood pair, the expected rates are

$$
  \lVert u-u_h\rVert_{H^1}=O(h^2),\quad
  \lVert u-u_h\rVert_{L^2}=O(h^3),\quad
  \lVert p-p_h\rVert_{H^1}=O(h),\quad
  \lVert p-p_h\rVert_{L^2}=O(h^2).
$$

An affine divergence-free velocity with affine zero-mean pressure is first
reproduced to roundoff. The rate study then uses
$u=(x_1^3,0,\ldots,0)$ and $p=x_0^2-\tfrac13$, with the manufactured
forcing $f=(-6x_1+2x_0,0,\ldots,0)$. The velocity is exactly
divergence-free and the pressure has zero mean on the unit box.

The exact data are supplied by `StokesData`, shared with the
[PETSc local/MPI suite](../PETScStokes/README.md). `StokesProblem` supplies
the native mixed solve and independent field-error measurements, shared
with the [p](../../p/Stokes/README.md) and [hp](../../hp/Stokes/README.md)
studies. The native configuration uses SparseLU and quadrature order 12,
with patches at `n=3` and rates on `n=3→5→9`; $h=1/(n-1)$.
Solver success and normalized coefficient residual below $10^{-11}$ are
checked separately. Pressure mean is integrated at order 14 with absolute
bound $10^{-10}$. These finite-workload checks do not establish a uniform
inf-sup bound on every supported cell family.

The affine patch test additionally checks the L2 divergence. The rate study
intentionally does not require monotone divergence:
Taylor--Hood enforces incompressibility only weakly against the pressure
space, and the unresolved component of $\nabla\cdot u_h$ has no monotonic
L2 convergence theorem.

Triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge grids
are covered. Segment is excluded because incompressible Stokes has no
non-degenerate velocity-pressure formulation in one spatial dimension.
