# Stokes combined mesh/degree verification

This suite uses the continuous problem, analytic data, pressure gauge, and
native mixed formulation specified in the
[degree-refinement suite](../../p/Stokes/README.md). On the unit box
$\Omega=(0,1)^d$, $d\in\lbrace 2,3\rbrace$, the viscosity is $\nu=1$ and

$$
u=\sin(\pi x_1)e_0,\qquad p=\cos(\pi x_0),\qquad
f=\bigl(\pi^2\sin(\pi x_1)-\pi\sin(\pi x_0)\bigr)e_0.
$$

The exact velocity trace is imposed on the full boundary. A global
$P_0^g$ multiplier enforces $\int_\Omega p_h\thinspace dx=0$. Both fields are
nonpolynomial, divergence vanishes, and the exact pressure has zero mean.
The forcing, viscosity, domain, and boundary data remain fixed while the
discrete spaces change.

## Refinement path and interpretation

With $n$ grid points per coordinate axis, nominal spacing $h=1/(n-1)$,
velocity degree $k$, and pressure degree $k-1$, the path is

| Level | $n$ | $h$ | Velocity/pressure degree |
| --- | --- | --- | --- |
| 0 | 3 | $1/2$ | $2/1$ |
| 1 | 4 | $1/3$ | $3/2$ |
| 2 | 5 | $1/4$ | $4/3$ |

The cell-appropriate scalar/vector H1 spaces use this degree offset on
every tested geometry. The coarsest mesh avoids the coarse-grid
pressure-rank obstruction documented and logically tested
in the p suite. Stability of a mixed pair remains a distinct hypothesis;
a degree offset alone does not prove a uniform inf-sup bound, particularly
on pyramid and wedge families.

For each field, independent physical-cell integration measures L2 and
H1-seminorm errors:

$$
E_{u,0,i}=\lVert u-u_i\rVert_{L^2},\quad
E_{u,1,i}=\lVert\nabla u-\nabla u_i\rVert_{L^2},\quad
E_{p,0,i}=\lVert p-p_i\rVert_{L^2},\quad
E_{p,1,i}=\lVert\nabla p-\nabla p_i\rVert_{L^2}.
$$

Every quantity must be finite and positive, and strictly decrease on
both adjacent intervals. The effective path rate is

$$
r_{w,j,i}=\frac{\log(E_{w,j,i-1}/E_{w,j,i})}
{\log(h_{i-1}/h_i)},\qquad w\in\lbrace u,p\rbrace,\quad j\in\lbrace 0,1\rbrace.
$$

The interval bounds are $r_{u,0,i}>2.5$, $r_{u,1,i}>1.5$,
$r_{p,0,i}>1.5$, and $r_{p,1,i}>0.5$. They are conservative finite-path
improvement bounds relative to the smooth lowest pair's nominal orders
$3,2,2,1$ under appropriate stability/regularity hypotheses. They are
not fixed-degree h orders or a claim about arbitrary combined refinement
paths. Separate field assertions prevent velocity accuracy from concealing
pressure error. Strong L2 divergence is not required to decrease
monotonically: the discrete divergence constraint is weak, not pointwise.

## Numerical budgets and controls

`StokesProblem` supplies the same fresh mixed system at every level and
is shared with native h and p tests. SparseLU factorization/solve success,
normalized coefficient residual below $10^{-11}$, and absolute pressure
mean below $10^{-10}$ are checked before field-error integration. Assembly
and error quadrature use order 16; the pressure-mean check uses order 18.
The native workload selects one same-precision residual correction with
the retained LU factors; the library default remains zero steps. The
[algebraic-accuracy methodology](../../isoparametric/Stokes/README.md#algebraic-accuracy-and-residual-correction)
states its construction and limitations without changing these budgets.
SparseLU is appropriate to the saddle-point operator; an SPD assumption is
not made.

Two independent controls on the coarsest `n=3` mesh exercise the highest
pair, $k=4$:

- The quartic patch $u=x_1^4e_0$, $p=x_0^3-1/4$, and
  $f=(-12x_1^2+3x_0^2)e_0$ must reproduce all four field quantities and
  strong L2 divergence with absolute errors below $10^{-9}$.
- With the quadratic reference $u=x_1^2e_0$, $p=x_0-1/2$, source
  $f=-e_0$, and unchanged trace, deliberately setting $\nu=2$ leaves
  the velocity exact but changes pressure to $3x_0-3/2$. The velocity
  errors remain below $10^{-9}$, while the pressure L2 and H1-seminorm
  errors exceed $0.1$ and $1$. Their analytic values are $1/\sqrt3$
  and $2$, respectively; the pressure oracle rejects the wrong physics.

A separate analytic $k=4$ sensitivity solve on `n=3` increases
assembly/error quadrature from 16 to 18. Each of the four errors must
change by less than $10^{-6}$ relative to baseline. This check does not
replace the three-level path or establish a quadrature bound at arbitrary
degree. SparseLU has no iterative tolerance to tighten.

## Execution scope

Triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge
families are tested with native Eigen sequential/OpenMP assembly.
Point and segment are excluded from this non-degenerate incompressible
workload. Curved geometry and PETSc/MPI hp solves have separate
[isoparametric](../../isoparametric/Stokes/README.md) and
[PETSc hp](../PETScStokes/README.md) specifications; their certificates
are not implied by this native study. Arbitrary higher degrees remain
outside these finite-workload tests.

Geometry entries are labeled `convergence;slow` with a 600-second limit.
Run `ctest --test-dir build/tests -R RodinConvergenceHPStokes --output-on-failure`.
