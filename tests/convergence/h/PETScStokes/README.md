# PETSc Stokes h-convergence

For $\Omega=(0,1)^d$, $d\in\lbrace 2,3\rbrace $, velocity
$u:\Omega\to\mathbb R^d$ and pressure $p:\Omega\to\mathbb R$ satisfy

$$
-\nu\Delta u+\nabla p=f,\qquad \nabla\cdot u=0,\qquad
u|_{\partial\Omega}=g,\qquad \int_\Omega p\,\mathrm dx=0,
$$

with $\nu=1$. The mixed spaces are vector H1 degree 2, scalar H1 degree 1,
and scalar P0g. The globally constant multiplier $\lambda$ fixes the
pressure gauge. For tests $(v,q,\mu)$, the implemented residual is

$$
\nu(\nabla u,\nabla v)-(p,\nabla\cdot v)
+(\nabla\cdot u,q)+(\lambda,q)+(p,\mu)-(f,v)=0.
$$

Only velocity has essential boundary constraints; pressure is not pinned
at an arbitrarily selected node. This is an indefinite mixed system, not
an SPD problem. PETSc PREONLY with LU/MUMPS is used in both mesh contexts.
The shared workload sets MUMPS `ICNTL(20)=0` for a centralized right-hand
side, avoiding a distributed-RHS scatter fault in the local MUMPS installation.
PETSc manages the distributed solution; `ICNTL(21)` is not a supported
runtime override in this interface. Matrix assembly and factorization remain
distributed; this is not a scalable iterative-solver benchmark.
The factor-workspace margin is `ICNTL(14)=100`, reserving 100% above the
symbolic estimate for delayed pivots in this indefinite system. This changes
storage provision, not the operator, pressure gauge, pivot threshold, or
acceptance tolerances. Explicit PETSc options retain precedence over these
test defaults. MUMPS `INFOG(1)` and PETSc preconditioner status are checked
alongside the convergence reason, with `INFOG(2)` reported on failure.
Real PETSc with MUMPS support is required; CMake reports exclusion when
that factorization is unavailable. CI builds the target explicitly so that
loss of the required support cannot silently remove this gate.

## Manufactured fields and controls

`StokesData` is shared with the native Eigen Stokes studies.
`PETScStokesProblem` supplies the PETSc mixed solve and diagnostics for
h, [p](../../p/PETScStokes/README.md), and
[hp](../../hp/PETScStokes/README.md) without degree-specific formulations.
Write
$x=(x_0,x_1,\ldots)$ and let $e_0$ denote the first coordinate vector.
The selected fields are

| Case | Velocity | Zero-mean pressure | Source |
| --- | --- | --- | --- |
| Affine patch | $x_1e_0$ | $x_0-1/2$ | $e_0$ |
| Quadratic patch | $x_1^2e_0$ | $x_0-1/2$ | $-e_0$ |
| Rate study | $x_1^3e_0$ | $x_0^2-1/3$ | $(-6x_1+2x_0)e_0$ |

Each velocity is divergence-free because its only nonzero component is
independent of $x_0$. Full nonzero traces are prescribed from these fields.
Both patches at `n=3` require velocity and pressure L2 and H1-seminorm
errors, and the L2 divergence, below $10^{-9}$.

A negative control uses $\nu=2$ with the quadratic source and trace still
derived for $\nu=1$. The exact velocity can remain unchanged: the incorrect
operator instead produces $p=3x_0-3/2$. Consequently, velocity-only acceptance
would miss this error. Pressure errors must exceed $0.1$ in L2 and $1$ in
H1 seminorm. The two fields are therefore measured separately.

## Refinement and acceptance

The rate path is `n=3→5→9`, with $h=1/(n-1)$, on triangle,
quadrilateral, tetrahedron, pyramid, hexahedron, and wedge grids. Under
smoothness, regular-mesh, stability, and dual-regularity assumptions, the
expected velocity L2/H1 orders are $3/2$ and pressure L2/H1 orders are $2/1$.
Every adjacent interval checks finite positive errors and strict reduction.
The accepted intervals are

$$
2.45<r_{u,0}<3.55,\quad 1.55<r_{u,1}<2.45,\quad
1.45<r_{p,0}<2.55,\quad 0.55<r_{p,1}<1.45.
$$

These are the existing native Stokes bounds. The finite mesh sequence does
not establish a uniform inf-sup theorem for arbitrary meshes or all element
families. Point and segment are excluded because the intended non-degenerate
incompressible velocity–pressure problem has dimension at least two.

Assembly and independent error integration use order 12. At `n=3`, order
14 must change each of the four errors by less than $10^{-8}$ relative.
The direct solve requires a positive PETSc convergence reason and an
independently recomputed residual

$$
\frac{\lVert A x-b\rVert_2}{\max(1,\lVert b\rVert_2)}<10^{-11}.
$$

The pressure integral is additionally checked below $10^{-10}$ using
Rodin's integral operator at order 14; this is a gauge diagnostic, whereas
the field errors are integrated independently. Weak incompressibility does
not imply monotone pointwise divergence errors, so no divergence rate is
asserted for the cubic case.

## Backend and execution scope

Local sequential/OpenMP PETSc assembly and MPI at ranks 1–4 are exercised.
Distributed norms integrate owned cells only before global reduction. An
additional MPI diagnostic interpolates $u=x_0e_0$ on `n=2` and checks
$\lVert\nabla\cdot u\rVert_{L^2(\Omega)}=1$; quadrilateral and hexahedron
cases include empty owned-cell partitions at multiple ranks. This checks
halo exclusion and the global divergence-norm reduction without a solve.

Each geometry entry has `convergence;petsc;slow` labels, a 600-second budget,
and the appropriate MPI processor count. Configure convergence and PETSc
support; enable MPI for distributed entries. Run
`ctest --test-dir build/tests -R 'RodinConvergenceHPETSc(MPI)?Stokes' --output-on-failure`.
These full-velocity-boundary studies do not cover complex fields, curved maps,
or scalable iterative preconditioning. Degree refinement is specified in the
linked p/hp suites; the distinct mixed-boundary formulation follows below.

## Physical traction and the pressure level

The `RodinConvergenceHPETScStokesBoundary` target uses
`PETScStokesTractionProblem`, with two fields rather than a pressure-mean
multiplier. Define $\Gamma_D=\lbrace x\in\partial\Omega:x_0=0\rbrace $ and
$\Gamma_N=\partial\Omega\setminus\Gamma_D$. With viscosity $\nu=1$,
the physical stress and boundary data are

$$
\sigma(u,p)=2\varepsilon(u)-pI,\qquad
u|_{\Gamma_D}=g,\qquad \sigma(u,p)n|_{\Gamma_N}=t.
$$

For velocity tests vanishing on $\Gamma_D$, the implemented residual is

$$
2(\varepsilon(u),\varepsilon(v))-(p,\nabla\cdot v)
+(\nabla\cdot u,q)-(f,v)-\langle t,v\rangle_{\Gamma_N}=0.
$$

This differs from the full-boundary vector-Laplacian formulation above.
For the manufactured divergence-free fields, their strong volume equations
coincide, since $\nabla\cdot(2\varepsilon(u))=\Delta u$.
Their natural boundary operators do not coincide. Prescribed traction fixes
the pressure constant: adding a mean-zero constraint would generally impose
an additional, incompatible condition.

All traction cases use $p=p_0+2$, where $p_0$ is the corresponding zero-mean
unit-box pressure. The nonzero pressure mean must be recovered, not removed.
Quadratic velocity/affine pressure patches use P2/P1; cubic velocity/quadratic
pressure patches use P3/P2. At `n=3`, velocity L2/H1 errors, pressure L2/H1
errors, and divergence L2 must be below $10^{-9}$; the pressure integral must
differ from $2$ by less than $10^{-9}$.

The pressure-sensitive control changes traction alone to $t+c n$, with $c=2$.
The exact solution is then $(u,p-c)$, so

$$
\Vert p_h-p\Vert_{L^2(\Omega)}=|c|=2,\qquad
|p_h-p|_{H^1(\Omega)}=0,\qquad \int_\Omega p_h\,\mathrm dx=0.
$$

The velocity remains exact. These quantities are checked within the patch
budget; velocity-only tests would miss an incorrect pressure level.

Both degree pairs use `n=3→5→9`, hence $h=1/(n-1)$. P2/P1 uses cubic
velocity and quadratic pressure; P3/P2 uses quartic velocity and cubic
pressure. Under the relevant approximation, mixed stability, and dual
regularity hypotheses, the expected L2/H1 orders are $(K+1,K)$ for velocity
and $(K,K-1)$ for pressure. Every adjacent interval requires finite positive
errors, strict reduction, and

$$
r_{u,0}>K+0.45,\quad r_{u,1}>K-0.45,\quad
r_{p,0}>K-0.55,\quad r_{p,1}>K-1.45.
$$

These are case-specific acceptance floors, not a uniform inf-sup or mixed
boundary regularity theorem for every element family. The measured pressure
integral also obeys the independent unit-volume bound
$|\int_\Omega p_h\,\mathrm dx-2|\leq\Vert p_h-p\Vert_{L^2(\Omega)}$,
up to the stated $10^{-9}$ integration/solve budget.

Assembly order 16 and independent norm order 18 are varied separately to
18 and 20 at `n=3`; each nonzero field error may change by at most $10^{-6}$
relative. The direct-solve residual budget is $10^{-11}$, as above.
Pointwise stress and source evaluation are local. Distributed assembly,
factorization, coefficient residuals, error norms, and pressure integrals have
global semantics and require participation by the mesh communicator.

The target registers six geometries, Local and MPI ranks 1–4, with sequential
and OpenMP assembly selected at configuration time. Point and segment are
excluded for the same incompressibility reason as the baseline suite.
Registrations carry `slow` labels and a 1800-second budget. Run
`ctest --test-dir build/tests -R '^RodinConvergenceHPETScStokesBoundary_' --output-on-failure -j 1`.
Registration and the mathematical specification do not by themselves certify
passing rates; backend/geometry validation is required separately.
