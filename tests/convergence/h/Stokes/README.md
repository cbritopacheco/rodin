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
One same-precision residual correction reuses the direct LU factors;
the [algebraic-accuracy specification](../../isoparametric/Stokes/README.md#algebraic-accuracy-and-residual-correction)
states the construction and its limitations. The library default remains
zero correction steps; this setting is specific to the native test workload.
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
without solving a manufactured saddle-point system. For a velocity degree
$k\ge2$, let $V_h^0$ be the vector degree-$k$ space with homogeneous
essential trace and let $Q_h^0$ be the zero-mean subspace of the scalar
degree-$(k-1)$ pressure space. The pair-specific mesh sequences and
verification status are stated below. With
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

### Bounded-workspace spectral construction

The reduced divergence remains sparse. For a pressure-column block $T_J$,
the Schur route solves $AX_J=B^\mathsf{T}T_J$ and assigns
$S_{0,:,J}=T^\mathsf{T}(BX_J)$. The Frobenius residual is accumulated
over these disjoint column blocks; its normalization is unchanged.
The independent LLT route defines

$$
Z=L^{-1}PB^\mathsf{T}T.
$$

For a block $E_I$ of velocity coordinate columns, its rows are obtained
without constructing the full dense matrix $Z$:

$$
E_I^\mathsf{T}Z
=\bigl(T^\mathsf{T}BP^\mathsf{T}L^{-\mathsf{T}}E_I\bigr)^\mathsf{T}.
$$

Incremental Givens QR retains a pressure-sized upper-triangular factor $R$.
An orthogonal factor need not be stored because $Z=QR$ implies
$\sigma_j(ZC^{-\mathsf{T}})=\sigma_j(RC^{-\mathsf{T}})$.
Here the reduced factor has $\min(n_v,n_p-1)$ rows and $Q$ has orthonormal
columns. If $n_v<n_p-1$, the stored square factor contains additional zero
rows; the dimension-required zero singular values are padded explicitly.
Thus the independent SVD is retained rather than replaced by a
cross-product eigenproblem. Its pressure-sized factor is evaluated by
divide-and-conquer SVD; no singular vectors are requested. This changes the
spectral algorithm, not the compared quantity, padding rule or acceptance
budget. Algebraic known-spectrum controls and the independent generalized
Schur spectrum remain separate gates. Certification still requires the
complete geometry/backend matrix below, rather than agreement on small
algebraic cases alone. Identity rotations are skipped only for exact
zero entries. Rectangular zero padding depends on $n_v$ and $n_p-1$,
not on a numerical rank threshold.

Each rotation is constructed by Eigen's scaled Givens algorithm. For
subnormal entries $a,b$, forming $r=\operatorname{hypot}(a,b)$ first and
then taking $c=a/r$, $s=b/r$ can lose the identity $c^2+s^2=1$ through
rounding of $r$. The resulting rotation can corrupt ordinary-sized
trailing entries. The regression uses $A=I_2$, $M=I_3$, pressure-constant
coefficients $e_2$, and

$$
B=\begin{pmatrix}\eta&\eta\\1&0\\0&0\end{pmatrix},
$$

where $\eta$ takes two subnormal values and the smallest positive normal
value. At working precision the reduced squared spectrum is $(0,1)$.
The independent spectral comparison detects corruption of this value;
small pivots are retained without a rank threshold.

Blocks contain at most 64 columns, further limited by the 96 MiB allowance
for an individual blocked velocity workspace. Individual pressure-sized
matrices have a separate 128 MiB allowance. The finest $P_4/P_3$ pyramid
hierarchy has $n_p=3925$, so its full pressure matrix requires
$8n_p^2=123245000$ bytes in double precision. The 128 MiB ceiling is the
next binary allowance above that requirement; it does not change a
spectral or residual acceptance criterion. These policies do not bound
the total of all live workspaces,
sparse factors or backend storage. A separate slow algebraic regression
has $n_v=19999$ and $n_p=650$: the full dense coupling would exceed the
allowance, while $A=I$, $B=[I\;0]$ and $M=I$ give the independent spectrum
$\lambda_j=1$. Unequal mass, non-nodal constant coefficients, an exact
rectangular obstruction and omitted divergence are checked separately.
The algebraic regressions are locally verified; PDE-spectrum
recertification of this bounded-workspace implementation is pending.

For $k=2$, all six applicable geometries use $n=2,3,5$ on affine and exact quadratic
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

The $P_3/P_2$ extension was locally verified with the previous dense oracle
in sequential and OpenMP builds.
It retains all
six applicable geometries, affine and exact quadratic maps, and independent
assembly orders $16$ and $20$. Its three levels are $n=3,4,5$; every level
requires resolved positivity and both algebraic spectral routes to agree.
The smallest positive eigenvalue has the same $10^{-8}$ relative quadrature
contamination budget. A missing-divergence control on the curved $n=3$
mesh retains the dimensions and requires an exactly zero spectrum. The
pressure constant is obtained by interpolation in the actual degree-two
basis, not by assigning one to every coefficient. These tests do not assert
a power-law rate for $\beta_h$, a uniform lower bound, or behavior on the
unexamined single-interval $n=2$ higher-order mesh.

Sequential and OpenMP native assembly are separate execution gates. This
target does not certify PETSc/MPI operator spectra or a mesh-uniform
inf-sup theorem. Geometry-specific registrations retain all pair-specific tests, slow
labels and 1800-second watchdogs. The separate workspace regression is
also slow with a 1800-second watchdog. Each blocked velocity workspace
is bounded by 96 MiB and each pressure matrix by 128 MiB before allocation;
this is not a bound on all
simultaneously live matrices or factorization storage.

## Degree-four pressure-spectrum extension

The implemented $P_4/P_3$ extension applies the same finite-spectrum
criteria to all six geometries on affine and exact quadratic maps at
$n=3,4,5$. Each level compares orders $16$ and $20$, requires resolved
positivity, and checks both spectral routes and their residuals. The
curved $n=3$ control omits divergence while retaining all space dimensions
and requires an exactly zero spectrum. Pressure constants are interpolated
in the actual degree-three basis. This is implemented coverage with
numerical certification pending, not a mesh-uniform stability claim.

The complete hierarchy and its missing-divergence control have separate
slow registrations. The pyramid hierarchy has a 3600-second watchdog;
all other hierarchies and controls retain 1800 seconds. The existing
$P_2/P_1$ and $P_3/P_2$ geometry registrations exclude these new cases,
so every compiled case is selected once. No refinement level or algebraic
consistency budget is relaxed by the process partition.

CI runs this target separately from the baseline h/p/hp workload. Light
geometries (triangle, quadrilateral, hexahedron and wedge), tetrahedron,
and pyramid form three serial partitions per thread configuration.
Algebraic utility and workspace checks belong to the light partition.
Each job has a 180-minute build/execution budget; per-registration
watchdogs remain in force. The observed tetrahedral
$P_4/P_3$ hierarchy took approximately sixteen minutes locally. This
measurement motivates partitioning but is not a phase timing or a
performance regression threshold. Complete numerical certification of
the higher-order matrix remains pending. The native pyramid hierarchy
exceeded its earlier 1800-second limit without a reported numerical
assertion failure or memory stop. A later stack sample was inside the
independent pressure-sized singular-value calculation. The one-hour
allowance is an execution policy, not an established completed runtime
or a change to the spectral acceptance criteria.
