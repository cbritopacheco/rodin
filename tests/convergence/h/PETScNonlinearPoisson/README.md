# PETSc/SNES semilinear Poisson h-convergence

On $\Omega=(0,1)^d$, $d\in\lbrace 1,2,3\rbrace$, the homogeneous Dirichlet problem is

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

The original target below retains its homogeneous Dirichlet configuration.
The separate natural-boundary target is specified in the next section.

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
and a nonzero nonlinear iteration count are required when the homogeneous
space has free DOFs. The final residual is
reassembled after explicit state synchronization, outside the SNES cache,
and must satisfy

$$
\frac{\lVert F_h(u_h)\rVert_2}
{\max(1,\lVert F_h(0)\rVert_2)}<10^{-10}.
$$

This is a coefficient-vector norm, distinct from the field-error norms.
The free dimension is determined from logical Dirichlet indices, not
residual magnitude. If $N$ is the global space size, $\mathcal I_r$ the
uniquely owned DOF range and $\mathcal C_r$ the constrained-index set
obtained from the assembled boundary map,

$$
N_{\mathrm{free}}=N-\sum_r\operatorname{card}(\mathcal C_r\cap\mathcal I_r).
$$

For $N_{\mathrm{free}}=0$, the homogeneous test space is $\lbrace 0\rbrace$;
zero initial residual and zero SNES iterations are required exactly.
The $n=2$, P1 regression checks this case on six families. Cube-centred
pyramid generation instead contributes one interior vertex and retains
the nontrivial solve contract. Refining to $n=3$ gives a positive free
dimension on every family. Owner-only integer counting also covers
MPI partitions without owned cells.

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
\lVert\nabla(z-u)\rVert_{L^2}\ge\pi\sqrt d\thinspace c_d2^{-d/2}.
$$

These lower bounds exceed $0.31$ and $1.69$ over the tested dimensions.
On the same mesh, the correct solution must have $E_0<0.05$, $E_1<0.2$;
the incorrect solution must exceed both budgets. A self-consistent but
physically incorrect residual–tangent pair therefore cannot pass solely
because SNES converges.

## Natural-boundary extension

`RodinConvergenceHPETScNonlinearPoissonBoundary` retains the same semilinear
operator. Mixed cases prescribe the manufactured trace on
$\Gamma_D=\lbrace x_0=0\rbrace$ and

$$
\partial_nu+\beta u=g,\qquad
g=\nabla u\cdot n+\beta u,\qquad \beta\in\lbrace 0,1\rbrace,
$$

on $\Gamma_N=\partial\Omega\setminus\Gamma_D$. Pure Neumann uses
$\Gamma_D=\varnothing$, $\Gamma_N=\partial\Omega$, and $\beta=0$.
For $V=\lbrace v\in H^1(\Omega):v\rvert_{\Gamma_D}=0\rbrace$, the physical residual
and its derivative are

$$
F_\beta(u;v)=\int_\Omega\bigl(\nabla u\cdot\nabla v
+(u+u^3-f)v\bigr)\thinspace \mathrm{d}x
+\int_{\Gamma_N}(\beta u-g)v\thinspace \mathrm{d}s,
$$

$$
J_\beta(u)[w;v]=\int_\Omega\bigl(\nabla w\cdot\nabla v
+(1+3u^2)wv\bigr)\thinspace \mathrm{d}x
+\beta\int_{\Gamma_N}wv\thinspace \mathrm{d}s.
$$

Since $(a+a^3-b-b^3)(a-b)\ge|a-b|^2$ for real $a,b$, the
reaction controls the constant mode even without essential data. No
compatibility projection, pressure-like multiplier, or mean constraint is
introduced. The global unconstrained dimension for pure Neumann is the
space's cached global size; a local size query does not introduce a reduction.

Nonzero traces use $u_h=I_hu+w_h$, with homogeneous correction on
$\Gamma_D$ and unrestricted correction on $\Gamma_N$. Every SNES callback
reconstructs this same physical state before residual or tangent assembly.
The manufactured flux is independent of the current iterate. Robin residual
and tangent both include their corresponding boundary term.

Sine fields with $A=1$ provide P1–P3 studies. P1 uses `n=5→9→17`, and
P2/P3 use `n=3→5→9`; every adjacent interval checks positive finite,
strictly decreasing errors. Pure-Neumann P1 tetrahedra instead use
`n=9→17→33`: the coarser `n=5→9` interval has an observed L2 rate
of approximately $1.602$ and does not meet the unchanged $1.65$ floor.
The finer hierarchy retains three levels and tests both adjacent intervals,
without claiming that the finite-resolution floor holds on the rejected
coarse hierarchy. Every accepted interval checks
strictly decreasing errors, with L2/H1-seminorm floors $1.65/0.75$,
$2.45/1.55$, and $3.45/2.55$. These policies retain the regularity
hypotheses stated above. At `n=3`, constant and affine P1 patches and the
quadratic P2 patch use $A=1/4$ and require both errors below $10^{-9}$:

$$
u=A,\qquad u=A\left(1+\sum_jx_j\right),\qquad
u=A\left(1+\sum_jx_j^2\right).
$$

For the last patch, $\nabla u=2Ax$, $\Delta u=2Ad$, and
$f=-2Ad+u+u^3$. Independent data tests check this formula on all seven
affine cell families. Removing the cubic term while retaining the affine
manufactured source must give $E_0,E_1>10^{-3}$. Independently removing
normal flux, but retaining the Robin exact-value load, is tested on the same
exactly representable affine P2 patch. The correct formulation must satisfy
$E_0,E_1<10^{-9}$, whereas the omitted-flux formulation must satisfy
$E_0,E_1>10^{-3}$ and exceed each corresponding correct error by a factor
greater than five. This control isolates the boundary term from coarse-mesh
approximation error; the sine field remains the rate-study solution.

Residual–Jacobian checks use the centered callback difference with
$\epsilon=10^{-5}$ described above. P1 and P2 require defect below
$10^{-6}$; the missing cubic derivative and missing Robin derivative
controls require defect above $10^{-3}$. The direction is proportional to
$x_0$: it vanishes on $\Gamma_D$ but has nonzero Robin trace. A direction
vanishing on the whole boundary could fail to detect the latter defect.

Independent sensitivity changes assembly order $16\to18$, norm order
$18\to20$, or SNES residual tolerance $10^{-12}\to10^{-13}$, one at a
time; each error must change by less than $10^{-6}$ relative to baseline.
The independently reassembled residual budget remains $10^{-10}$ relative
to $\max(1,\Vert F(u_h^0)\Vert_2)$. Inner KSP tolerances track SNES as above.
There are 37 cases per geometry/context, registered on all seven geometries
locally and at MPI ranks 1–4 in 35 `slow` CTest entries. The watchdog is
1800 seconds except for tetrahedron entries, which allow 3600 seconds:
the complete one-rank tetrahedron registration took approximately 1317
seconds with idle sleep prevented. The larger watchdog provides
platform/version headroom, rather than accommodating host suspension.
These are execution budgets, not numerical acceptance bounds. CI separates the
tetrahedron and pyramid registrations from the other geometries. Run
`ctest --test-dir build/tests -R '^RodinConvergenceHPETScNonlinearPoissonBoundary_' --output-on-failure`.

Space construction, distributed assembly, state synchronization, nonlinear
solving, residual norms, and error integration require communicator
participation. Manufactured field evaluation, normal contractions, and
cached layout queries are noncollective.

## Execution scope

Segment, triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge
meshes are covered locally and with MPI ranks 1–4. The coarse P2 segment
case includes ranks without owned cells. Sequential/OpenMP assembly is
selected by `RODIN_MULTITHREADED`; local and MPI contexts remain separate
tests. Each geometry/rank entry has labels `convergence;petsc;slow`
(plus `distributed` for MPI), a 600-second limit, and the MPI processor count.
Run `ctest --test-dir build/tests -R 'RodinConvergenceHPETSc(MPI)?NonlinearPoisson' --output-on-failure`.

Both targets require real-scalar PETSc. Complex cubic-reaction monotonicity,
0D spatial rates, curved maps, PETSc p/hp paths,
nonlinear divergence behavior, and arbitrary high degrees are not certified
by this suite. CI uses PETSc 3.19's uncached callback path; local newer PETSc
also exercises identity/state caching. Passing one version does not certify
the other.
