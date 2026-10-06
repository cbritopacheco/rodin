# PETSc complex Helmholtz h-convergence

For $\Omega=(0,1)^d$, $d\in\lbrace 1,2,3\rbrace $, the complex field
$u:\Omega\to\mathbb C$ satisfies

$$
-\Delta u-k^2u=f,\qquad k^2=\tfrac14,
$$

with its exact trace prescribed on the full boundary. The homogeneous test
space is $H_0^1(\Omega;\mathbb C)$ and the trial-first form is

$$
a(u,v)=\int_\Omega\nabla u\cdot\overline{\nabla v}
-k^2u\overline v\,\mathrm{d}x
=\int_\Omega f\overline v\,\mathrm{d}x.
$$

On the unit box the first Dirichlet eigenvalue is $d\pi^2$. Since
$k^2<d\pi^2$, the form is Hermitian coercive and CG is applicable. This
is a low-wavenumber test, not certification of indefinite Helmholtz solvers
or high-frequency pollution behavior.

## Manufactured fields

Let $s=\sum_jx_j$ and $c=1+i/2$. `HelmholtzData` supplies:

| Case | Field | Gradient entries | Source |
| --- | --- | --- | --- |
| Affine P1 patch | $u=1+2i+cs$ | $\partial_j u=c$ | $f=-k^2u$ |
| Quadratic P2 patch | $u=1+2i+cs^2$ | $\partial_j u=2cs$ | $f=-2dc-k^2u$ |
| Smooth rates | $u=e^{is}$ | $\partial_j u=iu$ | $f=(d-k^2)u$ |

Both real and imaginary channels are nonzero. The plane wave has a varying
phase, so it is not merely a real field multiplied by a fixed complex scalar.
Nonzero manufactured traces exercise complex boundary values as well as
complex loads and the trial/test conjugation convention.

Independent quadrature measures

$$
E_0^2=\int_\Omega|u-u_h|^2\,\mathrm{d}x,\qquad
E_1^2=\int_\Omega\sum_{j=1}^{d}|\partial_ju-\partial_ju_h|^2\,\mathrm{d}x.
$$

Patches use `n=5` grid points per axis and require both errors below
$10^{-9}$. A negative control omits the mass term while keeping the correct
affine source and trace; it requires $E_0>10^{-3}$ and $E_1>10^{-2}$.
This demonstrates rejection of an incorrect Helmholtz operator.
An additional solver-independent segment regression interpolates the affine
field before assembly and requires $\lVert A x_{\mathrm{exact}}-b\rVert_2<10^{-11}$
after boundary elimination. In particular, the real mass coefficient must
remain present when the field and assembled operator are complex.

## Refinement and numerical budget

| Space | Levels | Expected L2/H1 orders | Accepted adjacent orders |
| --- | --- | --- | --- |
| Complex H1 P1 | `n=5→9→17` | $2/1$ | $1.65<r_0<2.35$, $0.75<r_1<1.25$ |
| Complex H1 P2 | `n=3→5→9` | $3/2$ | $2.45<r_0<3.55$, $1.55<r_1<2.45$ |

With $h=1/(n-1)$, every interval checks finite positive errors, strict
reduction, and both rate bounds. These match the existing Eigen Helmholtz
bounds and assume smooth fields, conforming spaces, regular affine meshes,
consistent quadrature, and the required dual regularity.

All forms and error integrals use order 12. CG uses relative tolerance
$10^{-13}$, absolute tolerance $10^{-14}$, divergence threshold $10^5$,
and at most 20,000 iterations. Each solve requires a positive PETSc
convergence reason and a finite reported residual below $10^{-9}$; the
reported residual is not an independently recomputed unpreconditioned one.
The P2 sensitivity check at `n=5` raises quadrature to 14 and tightens
relative tolerance to $10^{-14}$, requiring relative changes below
$10^{-6}$ in each error.

## Geometry and scalar/backend selection

Local PETSc and MPI PETSc entries cover segment, triangle, quadrilateral,
tetrahedron, pyramid, hexahedron, and wedge. MPI runs use ranks 1–4, with
empty owned-cell partitions on the coarsest P2 segment case. Squared norms
are integrated on owned cells only and globally reduced. A zero discrete
field compared with $1+i$ checks $E_0=\sqrt2$ and $E_1=0$ through both norm
interfaces, explicitly exercising complex magnitudes and halo exclusion.

Native complex PETSc is required: `PETSC_USE_COMPLEX` must be defined in
the selected headers and the library must be from the same configuration.
CMake probes the linked PETSc target and registers this suite only for a
complex build; existing real-scalar PETSc suites remain in the separate real
configuration. Use a fresh build directory when changing PETSc installations.
Configure with `RODIN_BUILD_CONVERGENCE_TESTS=ON`, `RODIN_USE_PETSC=ON`, and
`RODIN_USE_MPI=ON` for distributed entries. Pass `PETSc_DIR` or explicit
`PETSc_INCLUDE_DIR`/`PETSc_LIBRARY` paths to select the complex installation.
Run `ctest --test-dir build/tests -R 'RodinConvergenceHPETSc(MPI)?Helmholtz' --output-on-failure`.

The separate `PETScComplexConvergence` CI matrix uses complex PETSc 3.19
on Ubuntu 24.04 and sequential/OpenMP local assembly configurations, while
also running MPI at 1–4 ranks. Geometry suites have the `slow` label and
600-second budgets, with MPI processor counts declared for scheduling; the
assembly regression has a 60-second budget. Complex p/hp PETSc studies,
curved geometry, and indefinite/high-frequency
Helmholtz remain outside this suite's claims. Point/0D has no PDE h-rate here.

## Mixed Neumann and complex impedance target

`RodinConvergenceHPETScHelmholtzBoundary` uses the same physical-coordinate
`HelmholtzData` and shared PETSc problem as the degree-refinement studies.
On the unit box, the face $x_0=0$ carries essential data; the remaining
faces form $\Gamma_N$. With $k^2=1/4$, the natural condition is

$$
\partial_n u+\mathrm{i}\beta u=r,\qquad
r=\nabla u\cdot n+\mathrm{i}\beta u,\qquad \beta\in\lbrace 0,1\rbrace .
$$

The normal contraction uses no conjugation. The weak form is sesquilinear
in the field and test function:

$$
\int_\Omega\nabla u\cdot\overline{\nabla v}\,\mathrm{d}x
-k^2\int_\Omega u\overline v\,\mathrm{d}x
+\mathrm{i}\beta\int_{\Gamma_N}u\overline v\,\mathrm{d}s
=\int_\Omega f\overline v\,\mathrm{d}x
+\int_{\Gamma_N}r\overline v\,\mathrm{d}s.
$$

For functions vanishing on $x_0=0$, the unit-box Poincare inequality gives
$\Vert v\Vert _{L^2}^2\le4\Vert \nabla v\Vert _{L^2}^2/\pi^2$.
Consequently the real part of the form is coercive since
$k^2<\pi^2/4$. The mixed Neumann case uses CG; the non-Hermitian impedance
case uses GMRES with Jacobi preconditioning. This low-frequency coercive
configuration does not certify the indefinite or resonant regime.

The affine P1 and quadratic P2 patches retain their nonzero real and
imaginary parts and require L2/H1-seminorm errors below $10^{-9}$ at `n=3`.
The rate field is $u=\exp(\mathrm{i}\sum_jx_j)$, with
$f=(d-k^2)u$. Both boundary conditions use P1 on `5→9→17` and P2/P3 on
`3→5→9`. Every adjacent interval must reduce both errors and exceed the
L2/H1-seminorm floors $1.65/0.75$, $2.45/1.55$, and $3.45/2.55$.
The floors are case-specific numerical criteria; optimal L2 estimates
additionally require the corresponding adjoint regularity.

Two controls retain the exact forcing and traces while independently
removing the volume mass term or normal flux. Missing mass must violate
the affine patch with L2 error above $10^{-3}$ and H1-seminorm error above
$10^{-2}$; missing flux must increase each smooth-field error by a factor
greater than 5. The impedance datum and boundary mass remain present in
the missing-flux case.

Assembly order 16, independent norm order 18 and relative solver tolerance
$10^{-13}$ are varied independently to 18, 20 and $10^{-14}$, respectively;
each field error must change by less than $10^{-6}$ relatively. Every solve
checks positive solver status and an independently recomputed residual
$\Vert Au_h-b\Vert _2/\max(1,\Vert b\Vert _2)<10^{-11}$. Global MPI norms and residuals have
intentional collective semantics; boundary attributes are assigned before
partitioning, without coordinate-based entity reconciliation.

All seven positive-dimensional geometries have local and MPI rank 1–4
registrations, including an empty-shard case on the coarsest segment mesh.
Sequential and OpenMP boundary jobs have separate CI runtime budgets and
1800-second, slow-tagged geometry registrations. These registrations state
the verification scope; numerical certification requires passing results.
