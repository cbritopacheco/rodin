# PETSc Stokes combined mesh/degree verification

The continuous problem and analytic data are those of the
[PETSc degree study](../../p/PETScStokes/README.md), with viscosity
$\nu=1$, full exact velocity trace, and zero-mean pressure on
$\Omega=(0,1)^d$, $d\in\lbrace 2,3\rbrace$:

$$
u=\sin(\pi x_1)e_0,\qquad p=\cos(\pi x_0),\qquad
f=\bigl(\pi^2\sin(\pi x_1)-\pi\sin(\pi x_0)\bigr)e_0.
$$

The geometry, coefficients, source, and boundary data remain fixed in
physical coordinates. Mesh spacing and both field degrees change along
the same path as the [native hp suite](../Stokes/README.md).

## Refinement and quantities

For $n$ grid points per coordinate axis, $h=1/(n-1)$, and velocity degree
$k$, scalar pressure has degree $k-1$ and a global $P_0^g$ multiplier
fixes its mean:

| Level | $n$ | $h$ | Velocity/pressure degree |
| --- | --- | --- | --- |
| 0 | 3 | $1/2$ | $2/1$ |
| 1 | 4 | $1/3$ | $3/2$ |
| 2 | 5 | $1/4$ | $4/3$ |

The coarsest mesh avoids the single-cell tensor-product pressure-rank
obstruction, not every possible mixed-space instability. Suitable inf-sup
stability and smooth-data regularity remain hypotheses; the finite path
does not prove uniform stability of the pyramid/wedge families.

Independent physical-cell integration measures velocity and pressure
L2 and H1-seminorm errors separately. For field $w\in\lbrace u,p\rbrace$ and
norm index $j\in\lbrace 0,1\rbrace$, the effective path rate is

$$
r_{w,j,i}=\frac{\log(E_{w,j,i-1}/E_{w,j,i})}
{\log(h_{i-1}/h_i)},\qquad i\in\lbrace 1,2\rbrace.
$$

Each error must be finite and positive and strictly decrease on both
intervals. The bounds are

$$
r_{u,0,i}>2.5,\quad r_{u,1,i}>1.5,\quad
r_{p,0,i}>1.5,\quad r_{p,1,i}>0.5.
$$

These conservative finite-path bounds are the native hp assertions,
relative to the smooth lowest pair's nominal orders $3,2,2,1$ under
the appropriate stability/regularity hypotheses. They are not fixed-degree
h orders, an asymptotic hp theorem, or evidence for arbitrary paths.
Strong divergence is not required to decrease monotonically, since the
discrete incompressibility constraint is weak rather than pointwise.

## Solver, integration, and controls

`PETScStokesProblem` is shared by h, p, and hp studies. Each solve creates
fresh spaces and a fresh PETSc layout; matrices are not resized between
levels. PREONLY/LU/MUMPS handles the saddle-point operator. Solver reason,
finite independently recomputed normalized coefficient residual below
$10^{-11}$, and absolute pressure integral below $10^{-10}$ are checked
before field-error measurement. Assembly/error quadrature use order 16;
the pressure integral uses order 18.

At `n=3`, the highest P4/P3 pair has separate controls:

- The quartic patch $u=x_1^4e_0$, $p=x_0^3-1/4$, and
  $f=(-12x_1^2+3x_0^2)e_0$ requires all four field errors and
  strong L2 divergence below $10^{-9}$.
- Retaining the quadratic reference source/trace while changing viscosity
  to $\nu=2$ leaves velocity exact but changes pressure from $x_0-1/2$
  to $3x_0-3/2$. Velocity errors stay below $10^{-9}$, while pressure
  L2/H1-seminorm errors exceed $0.1$ and $1$; their exact values are
  $1/\sqrt3$ and $2$.
- Raising assembly/error quadrature from 16 to 18 for the analytic P4/P3
  field changes each measured error by less than $10^{-6}$ relative.

These controls are not replacements for the three-level path, and the
coarse-mesh quadrature comparison is not an arbitrary-degree integration
bound. Direct factorization has no iterative tolerance parameter.

MPI norms sum owned-cell squared contributions globally before taking
square roots. Higher-order halo/empty-partition norm diagnostics are
registered in the shared-data PETSc p suite. The shared workload uses
`ICNTL(20)=0` for a centralized RHS and `ICNTL(14)=100` for a 100%
factor-workspace margin above the symbolic estimate. The latter accommodates
delayed pivots without changing the discrete problem or numerical acceptance;
the [p-suite specification](../../p/PETScStokes/README.md) records the
older-stack reproduction. MUMPS factorization and PETSc preconditioner
status are checked explicitly, with factor diagnostics reported on failure.
PETSc manages the distributed solution. Assembly/factorization remain
distributed, but iterative scalability is not claimed. Explicit PETSc options
retain precedence over the test defaults.

## Execution scope

Triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge
families are tested locally and at MPI ranks 1–4. Local PETSc
sequential/OpenMP assembly and MPI runs in both build configurations are
distinct evidence. Point and segment are excluded from this incompressible
workload. Real PETSc with MUMPS is required and checked by CMake; CI builds
the target explicitly. Complex fields, curved maps, arbitrary higher orders,
and scalable iterative preconditioning remain outside this suite's scope.

Geometry/rank entries have `convergence;petsc;slow` labels, a 600-second
limit, and declared MPI process counts. Run
`ctest --test-dir build/tests -R 'RodinConvergenceHPPETSc(MPI)?Stokes' --output-on-failure`.
