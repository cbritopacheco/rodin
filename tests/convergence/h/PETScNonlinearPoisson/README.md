# PETSc/SNES semilinear Poisson h-convergence

On $\Omega=(0,1)^d$, $d\in\{1,2,3\}$, the homogeneous Dirichlet problem is

$$
-\Delta u+u+u^3=f,\qquad
u=A\prod_{j=0}^{d-1}\sin(\pi x_j),\qquad
f=(d\pi^2+1)u+u^3.
$$

`NonlinearPoissonData` supplies the continuous source and independent analytic
gradient shared with the native suites. Amplitude $A=1$ is used unless a
negative control states otherwise. The residual and tangent on
$V=H_0^1(\Omega)$ are

$$
F(u;v)=(\nabla u,\nabla v)+(u+u^3-f,v),\qquad
J(u)[w;v]=(\nabla w,\nabla v)+((1+3u^2)w,v).
$$

The reaction is strongly monotone because $1+3u^2\ge1$ for real $u$.
For this smooth reference on regular unit-box mesh families, the expected
finite-element errors satisfy $E_0=O(h^{p+1})$ and $E_1=O(h^p)$, where
$E_0=\lVert u-u_h\rVert_{L^2(\Omega)}$ and
$E_1=\lVert\nabla u-\nabla u_h\rVert_{L^2(\Omega)}$.
The finite studies check these orders; they do not prove an asymptotic theorem
for arbitrary nonlinear data or mesh families.

## Workload and acceptance

Scalar H1 degrees $p=1,2$ use `n=5→9→17` and `n=3→5→9`, respectively,
with $n$ points per axis and $h=1/(n-1)$. Every adjacent interval requires
finite positive errors, strict reduction, L2 rates in $(p+0.5,p+1.5)$,
and H1-seminorm rates in $(p-0.3,p+0.5)$. Assembly and independent
physical-cell error quadrature use order 12. MPI norms integrate owned cells
only and sum the squared errors globally.

`PETScNonlinearPoissonProblem` owns fixed spaces and PETSc layouts for one
mesh/degree. A separate state field receives the current SNES iterate,
including refreshed ghost values, before residual or tangent assembly.
Rodin assembles $J\delta u=-F$; the SNES residual callback returns the
negative assembled right-hand side. Zero correction trace is consistent with
the zero initial iterate and homogeneous reference trace.

SNES line-searched Newton starts at zero. Absolute and relative residual
tolerances are $10^{-11}$, step tolerance $10^{-14}$, with at most 20
iterations and 1000 function evaluations. The SPD tangent is solved by CG
with Jacobi preconditioning, relative tolerance $10^{-13}$, absolute
tolerance $10^{-14}$, and 50000-iteration limit. Positive SNES/KSP reasons
and a nonzero nonlinear iteration count are required. The final residual is
reassembled after explicit state synchronization, outside the SNES cache,
and must satisfy

$$
\frac{\lVert F_h(u_h)\rVert_2}
{\max(1,\lVert F_h(0)\rVert_2)}<10^{-10}.
$$

This is a coefficient-vector norm, distinct from the field-error norms.
The sensitivity case at $p=2$, `n=3` raises quadrature to 14 and tightens
SNES/CG tolerances by a factor ten. Each field error must change by less
than $10^{-6}$ relative to baseline.

## Independent and negative controls

At `n=3`, both degrees evaluate the SNES residual and Jacobian callbacks at
$u_h=I_h(u/2)$ in direction $w_h=I_h(u/4)$. With $\varepsilon=10^{-5}$,

$$
\frac{\left\lVert J_h(u_h)w_h-
\frac{F_h(u_h+\varepsilon w_h)-F_h(u_h-\varepsilon w_h)}{2\varepsilon}
\right\rVert_2}
{\left\lVert
\frac{F_h(u_h+\varepsilon w_h)-F_h(u_h-\varepsilon w_h)}{2\varepsilon}
\right\rVert_2}<10^{-6}.
$$

The perturbed states are distinct PETSc vectors. This also exercises iterate
identity in callback caching: a change counter belongs to one object and is
not a globally unique state key. A separate unit regression asserts exact
callback counts and exact residual equality for distinct vectors with equal
counters. A deliberately incorrect P2 tangent substitutes $1+u_h^2$ for
$1+3u_h^2$ and must give a relative defect above $10^{-3}$.

The physical negative control uses $A=4$, P2, and `n=9`. Omitting the cubic
reaction from both residual and tangent while retaining the correct source
gives a linear solution $z$ with $(-\Delta+1)(z-u)=u^3$. Its first sine-mode
coefficient is $c_d=64(3/4)^d/(d\pi^2+1)$. Thus

$$
\lVert z-u\rVert_{L^2}\ge c_d2^{-d/2},\qquad
\lVert\nabla(z-u)\rVert_{L^2}\ge\pi\sqrt d\,c_d2^{-d/2}.
$$

These lower bounds exceed $0.31$ and $1.69$ over the tested dimensions.
On the same mesh, the correct solution must have $E_0<0.05$, $E_1<0.2$;
the incorrect solution must exceed both budgets. A self-consistent but
physically incorrect residual–tangent pair therefore cannot pass solely
because SNES converges.

## Execution scope

Segment, triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge
meshes are covered locally and with MPI ranks 1–4. The coarse P2 segment
case includes ranks without owned cells. Sequential/OpenMP assembly is
selected by `RODIN_MULTITHREADED`; local and MPI contexts remain separate
tests. Each geometry/rank entry has labels `convergence;petsc;slow`
(plus `distributed` for MPI), a 600-second limit, and the MPI processor count.
Run `ctest --test-dir build/tests -R 'RodinConvergenceHPETSc(MPI)?NonlinearPoisson' --output-on-failure`.

The suite requires real-scalar PETSc. Complex cubic-reaction monotonicity,
0D spatial rates, nonzero/mixed traces, curved maps, PETSc p/hp paths,
nonlinear divergence behavior, and arbitrary high degrees are not certified
by this suite. CI uses PETSc 3.19's uncached callback path; local newer PETSc
also exercises identity/state caching. Passing one version does not certify
the other.
