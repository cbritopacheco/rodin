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

## Independent finite-mesh pressure stability

The separate `RodinConvergenceHStokesStability` target checks pressure modes
without solving a manufactured saddle-point system. Let $V_h^0$ be the
vector degree-two space with homogeneous essential trace and let $Q_h^0$
be the zero-mean subspace of the scalar degree-one pressure space. With
velocity energy $a(v,v)=\int_\Omega Dv:Dv\thinspace\mathrm{d}x$ and
pressure norm $\lVert p\rVert_{L^2(\Omega)}$, define

$$
\beta_h=\inf_{p\in Q_h^0\setminus\lbrace0\rbrace}
\sup_{v\in V_h^0\setminus\lbrace0\rbrace}
\frac{\int_\Omega p\,\operatorname{div}v\thinspace\mathrm{d}x}
{\lVert p\rVert_{L^2(\Omega)}\sqrt{a(v,v)}}.
$$

The real velocity energy matrix $A$, divergence matrix $B$, and pressure mass
matrix $M$ are assembled independently. Essential velocity DOFs are removed
using the actual logical boundary-index map, not an assumed nodal ordering.
With $n_v=\dim V_h^0$ and $n_p=\dim Q_h$, their dimensions are
$A\in\mathbb R^{n_v\times n_v}$, $B\in\mathbb R^{n_p\times n_v}$
and $M\in\mathbb R^{n_p\times n_p}$, where $Q_h$ is the full pressure space.
For coefficients $c$ of the independently interpolated pressure constant,
the mean functional is $m=Mc$. For a pressure coefficient vector
$\boldsymbol p\in\mathbb R^{n_p}$, a basis matrix
$T\in\mathbb R^{n_p\times(n_p-1)}$ for $m^\mathsf{T}\boldsymbol p=0$
uses columns $e_i-(m_i/m_r)e_r$, where $r$ selects a nonzero mean entry.
The pressure basis therefore removes only the mean constraint. It does not
discard modes through a numerical rank threshold.

The finite inf-sup quantity is obtained from

$$
S_0=T^\mathsf{T}BA^{-1}B^\mathsf{T}T,\qquad
M_0=T^\mathsf{T}MT,\qquad
S_0z=\lambda M_0z,\qquad \beta_h^2=\lambda_{\min}.
$$

The shared `MixedStability` measurement uses sparse LDLT for this Schur
route. Independently, factorizations $PAP^\mathsf{T}=LL^\mathsf{T}$ and
$M_0=CC^\mathsf{T}$ give the whitened divergence

$$
W=C^{-1}T^\mathsf{T}BP^\mathsf{T}L^{-\mathsf{T}},\qquad
\lambda_j(S_0,M_0)=\sigma_j(W)^2.
$$

All eigenvalues are compared with the squared singular-value spectrum,
including zeros required by rectangular dimensions. The normalized
eigenproblem, velocity-solve, mean-basis and spectral-comparison defects
must be below the dimensionless algebraic budget $10^{-10}$. Resolved
positivity requires $\lambda_{\min}>10^{-10}\max_j|\lambda_j|$;
this separates the positive finite spectrum from its algebraic consistency
scale, not from an independently proved mesh-uniform stability bound.
A synthetic unequal-mass example has the independent value
$\beta_h^2=4/3$ and checks logical constraint elimination and mean reduction.

All six applicable geometries use $n=2,3,5$ on affine and exact quadratic
maps. Assembly orders $16$ and $20$ are compared independently; each positive
minimum eigenvalue must change relatively by less than $10^{-8}$.
At $n=2$, every family except Pyramid has fewer free velocity DOFs than
zero-mean pressure DOFs. This exact dimension obstruction is checked
separately: $\operatorname{rank}B\le\dim V_h^0<\dim Q_h^0$ implies
additional pressure null modes even after fixing the mean. Pyramid at
$n=2$ and all families at $n=3,5$ require resolved positive spectra.
Omitting the divergence operator on the curved $n=3$ mesh must give an
exactly zero spectrum and fail the positivity predicate with unchanged
space dimensions.

Sequential and OpenMP native assembly are separate execution gates. This
target does not certify PETSc/MPI operator spectra or a mesh-uniform
inf-sup theorem. Geometry-specific registrations retain both tests, slow
labels and 1800-second watchdogs. Each dense velocity/pressure workspace
is bounded by 96 MiB before allocation; this is not a bound on all
simultaneously live matrices or factorization storage.
