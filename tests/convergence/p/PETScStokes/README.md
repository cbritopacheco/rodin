# PETSc Stokes degree-refinement verification

This suite extends the [native degree study](../Stokes/README.md) to
real-scalar PETSc storage, local PETSc assembly, and distributed MPI meshes.
On $\Omega=(0,1)^d$, $d\in\{2,3\}$, the continuous problem is

$$
-\nu\Delta u+\nabla p=f,\qquad \nabla\cdot u=0,\qquad
u|_{\partial\Omega}=g,\qquad \int_\Omega p\,dx=0,
$$

with $\nu=1$, full manufactured velocity trace $g$, and analytic data

$$
u=\sin(\pi x_1)e_0,\qquad p=\cos(\pi x_0),\qquad
f=\bigl(\pi^2\sin(\pi x_1)-\pi\sin(\pi x_0)\bigr)e_0.
$$

The shared `StokesData` supplies these fields and their analytic derivatives
in physical coordinates. The source, trace, viscosity, domain, and mesh are
held fixed as the field degrees increase.

## Spaces, levels, and mathematical scope

On a fixed `n=3` mesh, $h=1/2$, the vector velocity and scalar pressure
use the cell-appropriate H1 families with degree pairs

$$
(k_u,k_p)=(2,1)\longrightarrow(3,2)\longrightarrow(4,3).
$$

A real $P_0^g$ multiplier fixes the global pressure mean. Pressure is not
pinned at a selected mesh node. `PETScStokesProblem` constructs a fresh
space tuple and PETSc system for each degree; assembled matrices are never
resized. The same mixed form is shared with the PETSc h and hp studies.

The fixed mesh avoids the single-cell tensor-product rank obstruction
proved by exact constrained-DOF counts in the native p suite. Under
appropriate mixed stability and regularity hypotheses, increasing degree
is expected to improve both fields. These finite-workload checks do not
establish a uniform discrete inf-sup bound, particularly for pyramid/wedge
families, or exponential convergence at arbitrary degree.

## Field errors and acceptance

Independent physical-cell quadrature measures

$$
E_{u,0}=\lVert u-u_h\rVert_{L^2(\Omega)},\quad
E_{u,1}=\lVert\nabla u-\nabla u_h\rVert_{L^2(\Omega)},\quad
E_{p,0}=\lVert p-p_h\rVert_{L^2(\Omega)},\quad
E_{p,1}=\lVert\nabla p-\nabla p_h\rVert_{L^2(\Omega)}.
$$

Each quantity must be finite and positive, strictly decrease on both
adjacent degree intervals, and have logarithmic decay

$$
\alpha_i=\frac{\log(E_{i-1}/E_i)}{k_i-k_{i-1}}>0.1.
$$

Velocity and pressure histories use their respective degrees; both
increments equal one. Each of the four errors is checked separately, so
accurate velocity cannot conceal inaccurate pressure. MPI norms integrate
only owned cells and globally sum squared contributions before taking
square roots; halo copies are not additional physical volume.

PREONLY/LU/MUMPS solves the saddle-point operator without an SPD assumption.
Positive PETSc convergence reason, finite coefficient residual
$\lVert Ax-b\rVert_2/\max(1,\lVert b\rVert_2)<10^{-11}$, and absolute
pressure integral below $10^{-10}$ are required. Assembly and field-error
quadrature use order 16; pressure mean uses order 18. The P4/P3 sensitivity
case raises assembly/error order to 18 and requires relative change below
$10^{-6}$ in every field quantity. The direct factorization has no
iterative tolerance to tighten.

MUMPS uses the shared h/p/hp defaults `ICNTL(20)=0` for a centralized
right-hand side and `ICNTL(14)=100` for a 100% factor-workspace margin above
the symbolic estimate. The latter accommodates delayed pivots without
altering the matrix, pressure gauge, pivot threshold, or error tolerances.
On the four-rank `n=3` pyramid P4/P3 system, PETSc 3.19.6/MUMPS 5.6.2
exhausts factor workspace at the solver's default 20% margin; the same
exported operator and right-hand side solve with the 100% margin and
normalized residual approximately $7.7\times10^{-14}$. This is a
solver-resource control, not evidence of a singular finite element system.
MUMPS `INFOG(1)` and PETSc preconditioner status must indicate success;
`INFOG(2)` and the degree/quadrature order identify failures.

PETSc manages the distributed solution; `ICNTL(21)` is not exposed as a
runtime override by this interface. Matrix assembly/factorization remain
distributed; no iterative-solver scalability claim is made. Explicit
PETSc options retain precedence over the shared test defaults.

## Polynomial, physical, and distributed controls

Each pair reproduces its highest-degree polynomial patch on `n=3`:

$$
u=x_1^m e_0,\qquad p=x_0^{m-1}-\frac1m,\qquad
f=\bigl(-m(m-1)x_1^{m-2}+(m-1)x_0^{m-2}\bigr)e_0,
\qquad m=k_u\in\{2,3,4\}.
$$

All four field errors and strong L2 divergence must be below $10^{-9}$.
These patches are distinct from the nonpolynomial degree study.

For every pair, a negative control retains the quadratic source/trace for
$u=x_1^2e_0$, $p=x_0-1/2$, but changes viscosity to $\nu=2$. Its velocity
remains exact while pressure becomes $3x_0-3/2$. The pressure-error
references are $1/\sqrt3$ in L2 and $2$ in H1 seminorm. Velocity errors
must remain below $10^{-9}$, while pressure errors must exceed $0.1$ and
$1$. Thus the field oracle rejects incorrect physics even when the
assembled residual is small.

An additional MPI norm diagnostic uses a P4 velocity $u_h=x_0e_0$,
reference velocity $u=2x_0e_0$, P3 pressure $p_h=0$, and reference pressure
$p=1$ on `n=2`. Without solving the singular coarse mixed problem, it checks

$$
E_{u,0}=1/\sqrt3,\qquad E_{u,1}=1,\qquad
E_{p,0}=1,\qquad E_{p,1}=0,\qquad
\lVert\nabla\cdot u_h\rVert_{L^2}=1.
$$

The absolute budget is $10^{-12}$. Quadrilateral/hexahedron cases at
multiple ranks include partitions with no owned cells. This verifies
higher-order global norm reduction and halo exclusion independently of
the PDE solve; it is not a tolerance-based mesh-ownership test.

## Execution scope

Triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge
families are registered locally and at MPI ranks 1, 2, 3, and 4. Local
sequential/OpenMP assembly and MPI runs in both build configurations are
distinct checks. Point and segment are excluded from the non-degenerate
incompressible workload. Complex PETSc, curved maps, arbitrary higher
degrees, and scalable iterative preconditioning are outside this claim.

Real PETSc with MUMPS is required. The shared convergence CMake capability
check excludes unsupported scalar/factorization configurations explicitly.
CI builds this target explicitly, preventing loss of required support
from silently removing the gate. Geometry/rank entries carry
`convergence;petsc;slow` labels, a 600-second limit, and MPI processor counts.
Run `ctest --test-dir build/tests -R 'RodinConvergencePPETSc(MPI)?Stokes' --output-on-failure`.
