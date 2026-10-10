# Curved Taylor–Hood Stokes verification

This suite tests the velocity/pressure pair $P_2/P_1$, together with a
global $P_{0g}$ pressure-mean multiplier, on exact quadratic geometry and
quadratic approximations of a nonpolynomial geometry map. In the exact-map
study, the physical domain is fixed across all refinement levels:

$$
\Omega=\Phi((0,1)^d),\qquad
\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1},\qquad d\in\lbrace 2,3\rbrace.
$$

All six applicable cell families are covered: triangle, quadrilateral,
tetrahedron, pyramid, hexahedron and wedge. The incompressible spatial
velocity/pressure problem is not asserted in zero or one dimension.
Native Eigen and real-PETSc local/MPI workloads reuse the shared Stokes
problem classes; sequential and OpenMP assembly consume the same forms.
MPI registrations cover one through four ranks.

## Physical formulation and gauge

With unit viscosity, the manufactured problem is

$$
-\Delta u+\nabla p=f,\qquad \mathrm{div}u=0,
\qquad u\rvert_{\partial\Omega}=g,\qquad \int_\Omega p\thinspace \mathrm{d}x=0.
$$

Since $\det D\Phi=1$ and $x_0=\xi_0$, the pressure fields retain their
analytic mean under this map. In particular, the rate data are

$$
u(x)=x_1^3e_0,\qquad p(x)=x_0^2-\tfrac13,
\qquad f(x)=(-6x_1+2x_0)e_0.
$$

The pressure gauge follows by change of variables:

$$
\int_\Omega p\thinspace \mathrm{d}x
=\int_{(0,1)^d}(\xi_0^2-\tfrac13)\thinspace \mathrm{d}\xi=0.
$$

The exact volume and pressure mean are independently integrated on the
mapped cells. A fresh global mean-multiplier space enforces the computed
pressure gauge; no pressure point is pinned. Every solve checks the
algebraic residual and computed pressure mean before field-error integration.

## Refinement and acceptance

| Dimension | Grid points per axis | Nominal spacing | Fixed field/geometry degrees |
| --- | --- | --- | --- |
| 2 | $3\to5\to9$ | $h=1/(n-1)$ | Velocity 2, pressure 1, geometry 2 |
| 3 | $3\to4\to5$ | $h=1/(n-1)$ | Velocity 2, pressure 1, geometry 2 |

Each sequence has three levels, and every adjacent interval is checked.
The smaller three-dimensional grids bound the cost of mixed direct solves;
the rate denominator uses the actual spacing ratio, not an assumed factor two.

Subject to mapped-pair stability and the requisite solution/dual regularity,
the expected velocity L2/H1-seminorm orders are three/two, and the pressure
L2/H1-seminorm orders are two/one. Every field error must be finite,
positive and strictly decreasing. The L2 acceptance window is within 0.55
of its expected order, and the derivative window is within 0.45.
These are finite-resolution acceptance policies, not uniform inf-sup proofs
for arbitrary curved meshes or degrees.

The independently integrated divergence is controlled through

$$
\Vert\mathrm{div}u_h\Vert_{L^2(\Omega)}
\le\sqrt d\thinspace |u-u_h\rvert_{H^1(\Omega)}.
$$

This follows from $\mathrm{div}u=0$ and the Frobenius trace bound.
It ties the divergence error to the decreasing velocity derivative error.
A strict divergence rate is not imposed: cancellation can make this
quantity vanish or reach roundoff before the field errors do.
An MPI interpolation control with $u_h=x_0e_0$ independently checks
$\Vert\mathrm{div}u_h\Vert_{L^2(\Omega)}=1$, counting owned cells once.

## Algebraic accuracy and residual correction

Let $A\in\mathbb R^{N\times N}$ and $b\in\mathbb R^N$ denote the
assembled, constrained mixed system, and let $z_\star=A^{-1}b$ be its
exact-arithmetic solution. A small coefficient residual alone does not
bound the pressure derivative error independently of the inverse operator:

$$
z_\star-z_h=A^{-1}(b-Az_h).
$$

The native workload therefore applies one residual correction after the
direct solve. The original LU factors are retained, and the operator,
right-hand side, quadrature and field-error budgets remain unchanged:

$$
r_k=b-Az_k,\qquad A\delta z_k=r_k,\qquad
z_{k+1}=z_k+\delta z_k,\qquad k=0.
$$

Residual evaluation and correction use the same scalar precision as the
direct solve. `SparseLU::setRefinementSteps` is opt-in; its library default
is zero. No monotonic improvement theorem is asserted for arbitrary matrices.
This setting affects the native solver, not the PETSc solver or MPI ownership.

The fixed-mesh `NativeCubicWedgePressureForwardAccuracy` regression uses
cubic geometry, velocity degree three, pressure degree two and five grid
points per axis. Its represented and lifted pressure L2 norms and H1
seminorms must be below the dimensionless absolute budget $10^{-11}$.
The observed uncorrected pressure H1 error was $2.47\times10^{-10}$,
whereas one correction reduced it to $4.35\times10^{-13}$. This is an
accuracy regression, not a convergence-rate study or a performance benchmark.
On the nine-point endpoint, the same correction reduced the observed
pressure H1 error from $2.61\times10^{-9}$ to $1.71\times10^{-12}$ without
relaxing the separate three-level patch budget $10^{-9}$.

## Patch, negative and quadrature controls

The physical-affine patch uses $u=x_1e_0$ and $p=x_0-1/2$.
In two dimensions its velocity pulls back to
$(\xi_1+0.1\xi_0^2)e_0$, which belongs to the mapped P2 space.
The pressure pulls back to P1 in both dimensions. All velocity, pressure
and divergence errors must be below $10^{-9}$. A physical quadratic
velocity is not used as an exact P2 patch, because its pullback can have
degree four in two dimensions.

A wrong-viscosity control solves with viscosity two but retains the
unit-viscosity forcing and trace for $u=x_1^2e_0$, $p=x_0-1/2$.
Its pressure L2 error must exceed 0.1 and its pressure H1-seminorm error
must exceed one. In three dimensions the exact wrong pressure has gradient
$3e_0$ instead of $e_0$; in two dimensions an additional field-approximation
error is present. The control establishes rejection of an incorrect operator,
not a rate claim for that operator.

Assembly and field-norm integration use quadrature order 12. A separate
order-16 solve/integration must change each nonzero velocity and pressure
norm by less than $10^{-6}$ relatively. Mapped derivative integrands are
not assumed polynomial-exact. The coefficient residual obeys
$\Vert Az-b\Vert_2/\max(1,\Vert b\Vert_2)<10^{-11}$, where the coefficient vector
collects velocity, pressure and mean-multiplier unknowns. The computed pressure integral
has absolute magnitude below $10^{-10}$, in this dimensionless setting.
Independent analytic-volume/gauge checks use an absolute $10^{-12}$ bound.

Tests are labelled `convergence;slow`, with `petsc` and `distributed`
where applicable. Geometry-specific registrations have 30-minute safety
timeouts; pyramid registrations share a resource lock. Performance benchmarks
remain separate from these convergence assertions.

## Approximated domains and exact-domain lifts

The second study uses the exact reference map and its nodal quadratic
interpolant,

$$
\Phi(\xi)=\xi+0.1\sin(\pi\xi_0)e_{d-1},\qquad
\Phi_h=I_2\Phi,\qquad \Omega_h=\Phi_h((0,1)^d).
$$

Manufactured data are evaluated in physical coordinates on $\Omega_h$:

$$
u(x)=x_{d-1}^3e_0,\qquad p(x)=x_0^2-\tfrac13,\qquad
f(x)=(-6x_{d-1}+2x_0)e_0.
$$

Using the deformed coordinate in the shear profile gives a nonzero velocity
geometry defect in both dimensions. The reference domain has unit volume
and zero analytic pressure mean. Independent mapped-cell integrals check
these quantities on each represented domain too; no cellwise identity
$\det D\Phi_h=1$ is assumed. Map regularity is checked independently by
the shared geometry study and remains a hypothesis of the rate estimates.

Let $\ell_h=\Phi_h\circ\Phi^{-1}$ and $w_h^\ell=w_h\circ\ell_h$.
For either field $w\in\lbrace u,p\rbrace$, the measured defects on $\Omega$ are

$$
F_w=w\circ\ell_h-w_h^\ell,\qquad
G_w=w-w\circ\ell_h,\qquad T_w=w-w_h^\ell=F_w+G_w.
$$

Their L2 norms and H1 seminorms use the exact-domain Jacobian and the
chain rule $D(w_h^\ell)=(Dw_h\circ\ell_h)D\ell_h$.
Both triangle and reverse-triangle inequalities are checked numerically.
Velocity represented-domain and lifted field errors have orders three/two;
velocity geometry and total errors have orders three/two for this quadratic
geometry. Pressure represented-domain, lifted field and total errors have
orders two/one. Since $\ell_h$ preserves $x_0$, $G_p=0$ analytically;
its two measured norms must be below $10^{-10}$, rather than being assigned
a rate or used in a relative comparison.

| Cell family | Grid points per axis in the sine-map rate study |
| --- | --- |
| Triangle | $9\to17\to33$ |
| Quadrilateral | $3\to5\to9$ |
| Tetrahedron | $9\to11\to13$ |
| Wedge | $5\to7\to9$ |
| Pyramid, hexahedron | $3\to4\to5$ |

Every adjacent interval uses the actual spacing ratio. The finer triangle,
tetrahedron and wedge levels avoid coarse pressure transients; the same
fixed acceptance windows stated above apply without relaxation. These
finite-hierarchy observations do not certify a uniform inf-sup constant.

For each lifted velocity defect $E\in\lbrace F_u,G_u,T_u\rbrace$, independently
integrated divergence obeys

$$
\Vert\mathrm{tr}DE\Vert_{L^2(\Omega)}\le\sqrt d\thinspace |E\rvert_{H^1(\Omega)}.
$$

A lift of a physically divergence-free field on $\Omega_h$ need not remain
divergence-free on $\Omega$. The affine patch checks zero represented and
lifted field errors, but explicitly retains a nonzero velocity geometry
derivative and divergence defect. The wrong-viscosity control uses pressure
as its rejection oracle: for a shear flow, pressure can absorb an incorrect
viscosity without necessarily increasing velocity error.

Assembly order $12\to16$ and norm order $14\to18$ are varied separately
on the fixed $n=5$ mesh. Each nonzero field norm must change by less than
$10^{-6}$ relatively. Existing direct-solver residual and pressure-gauge
checks remain in force. Native and PETSc solvers expose const solution
views only during a solve-scoped observer callback; the common lifted-norm
helper evaluates owned quadrature points and reduces MPI contributions.
It neither duplicates the mixed solve nor retains solution references.

## Linear and cubic geometry with representable Taylor--Hood fields

The additional geometry-degree suites use $q\in\lbrace1,3\rbrace$ and
the same sine map, with velocity/pressure degrees $k/(k-1)$ for
$k=\max(2,q)$. The physical data on every represented domain are

$$
u_\ast(x)=x_{d-1}e_0,\qquad p_\ast(x)=x_0-\tfrac12,\qquad
f=e_0,\qquad \mathrm{div}u_\ast=0.
$$

The affine velocity pullback is represented at degree $k$; the pressure
pullback is represented at degree $k-1$. With the code's defect convention,

$$
e_{G,u}(\Phi(\xi))=
\left(\Phi_{h,d-1}(\xi)-\Phi_{d-1}(\xi)\right)e_0,
\qquad e_{G,p}=0,
\qquad e_{T,u}=e_{F,u}+e_{G,u}.
$$

The zero pressure geometry defect follows from preservation of $x_0$;
its derivative defect is zero as well. Represented and lifted-field
velocity errors, represented and lifted-field/total pressure errors,
represented divergence, and lifted-field divergence use the dimensionless
absolute budget $10^{-9}$. The pressure geometry norms retain the stricter
$10^{-10}$ budget. No pressure or field roundoff error is fitted to a rate.

| Geometry degree $q$ | Velocity/pressure degrees | Grid points per axis |
| --- | --- | --- |
| $1$ | $2/1$ | $3,5,9$ |
| $3$ | $3/2$ | $3,5,9$ |

Velocity geometry and total errors must decrease across both intervals,
with expected $L^2/H^1$ orders $(q+1,q)$ and the existing $0.55/0.45$
acceptance windows. These interpolation orders require a regular smooth
map and a nonvanishing leading defect. The
[Helmholtz specification](../Helmholtz/README.md) gives the map regularity
argument. The full mixed system is solved at every level: residual and
computed pressure-mean checks precede the solve-scoped lift observer.
An independent $n=3$ physical integration checks represented volume one
and analytic pressure mean zero for each geometry degree.

For each lifted defect, divergence retains the Frobenius trace bound.
Since the exact velocity is divergence-free, the field-divergence norm
also bounds the difference of total and geometry divergence norms:

$$
\left|\Vert\mathrm{div}e_{T,u}\Vert_{L^2}
-\Vert\mathrm{div}e_{G,u}\Vert_{L^2}\right|
\le\Vert\mathrm{div}e_{F,u}\Vert_{L^2}
\le\sqrt d\,|e_{F,u}|_{H^1}.
$$

At $n=5$, assembly order $12\to16$ and norm order $14\to18$ are
varied independently. Positive velocity geometry/total norms require
relative changes below $10^{-6}$; divergence comparisons use the absolute
$10^{-9}$ budget. Both backends use direct solvers, so no fictitious
iterative-tolerance variation is introduced. Every solve retains the
coefficient residual budget $10^{-11}$ and the pressure-mean checks.

The six applicable cell families are registered for native local and
real-PETSc local/MPI configurations with ranks one through four.
Sequential/OpenMP remains a build choice. These representable-field
certificates do not prove a uniform inf-sup constant or nonpolynomial
approximation rates for arbitrary mapped Taylor--Hood pairs.

Matched cubic-geometry native and PETSc pyramid registrations include the complete
rate hierarchy, independent quadrature study and pressure-gauge check.
They retain the `slow` label and pyramid resource lock, with a one-hour
CTest execution budget; the other matched-degree registrations retain
their 30-minute budgets. These are scheduling limits, not numerical error
budgets or performance regression thresholds. The cubic pyramid rate
case alone took about 25 minutes with MPI rank one, before its separate
quadrature study; the native rate hierarchy exceeded the 30-minute default.
These execution observations do not establish an end-to-end stage breakdown
or a controlled performance comparison.

## Physical traction on the exact quadratic domain

`RodinConvergenceIsoparametricPETScStokesBoundary` uses the same exact
quadratic map but a mixed physical-stress boundary condition. Set
$\Gamma_D=\Phi(\lbrace \xi_0=0\rbrace)$ and
$\Gamma_T=\partial\Omega\setminus\Gamma_D$. With unit viscosity,

$$
\sigma(u,p)=2\varepsilon(u)-pI,\qquad
\varepsilon(u)=\tfrac12(Du+Du^T),\qquad
-\mathrm{div}\sigma(u,p)=f,\qquad\mathrm{div}u=0,
$$

and prescribe $u=u_\star$ on $\Gamma_D$ and
$\sigma(u,p)n=t_\star$ on $\Gamma_T$. For zero-trace velocity tests $v$
and pressure tests $q$, the assembled equations are

$$
\int_\Omega 2\varepsilon(u):\varepsilon(v)-p\mathrm{div}v
+q\mathrm{div}u\thinspace \mathrm{d}x
=\int_\Omega f\cdot v\thinspace \mathrm{d}x
+\int_{\Gamma_T}t_\star\cdot v\thinspace \mathrm{d}s.
$$

Traction determines the pressure level, so this two-field system has no
pressure-mean multiplier. Here $\det D\Phi=1$ and $x_0=\xi_0$; the
manufactured pressure integral is exactly two. Adding $2n$ to the prescribed
traction changes the exact pressure to $p_\star-2$, with unchanged velocity.
The negative control checks zero pressure integral and pressure $L^2$ error
two, alongside velocity, pressure-gradient, and divergence patch budgets.

Rate studies use $(u_\star,p_\star)=(x_1^3e_0,x_0^2-1/3+2)$ for
P2/P1 and $(x_1^4e_0,x_0^3-1/4+2)$ for P3/P2. Both adjacent intervals
of $n=3,5,9$ are checked. Expected velocity $L^2/H^1$ orders are
$(K+1,K)$ and pressure orders $(K,K-1)$; the respective lower bounds are
$(K+0.45,K-0.45)$ and $(K-0.55,K-1.45)$. These finite-hierarchy
checks remain conditional on mixed-pair stability and solution regularity;
they do not certify a uniform inf-sup constant.

Exact physical quadratic and cubic shear patches use P4/P3 and P6/P5.
In two dimensions, $x_1=\xi_1+0.1\xi_0^2$, so their velocity pullbacks
can have degrees four and six. These choices also cover the undeformed
shear coordinate in three dimensions; all patch field/divergence errors
must be below $10^{-9}$. The pressure-level negative control uses the
quadratic P4/P3 patch. Independent assembly order $16\to18$ and norm
order $18\to20$ variations require changes below $10^{-6}$ relatively.
The existing coefficient residual budget $10^{-11}$ remains separate.

Six cases per geometry and context cover the six applicable geometries,
locally and at MPI ranks one through four in real-PETSc sequential/OpenMP
configurations. Reference boundary classification precedes rank-local map
installation. Pointwise stress and geometry evaluation remain noncollective;
mixed assembly, factorization, residual/norm reductions, and the pressure
integral are explicitly global operations.

CTest isolates each rate hierarchy in its own `P2P1Rates` or `P3P2Rates`
process and groups the four fixed-mesh cases in `Controls`, for each geometry
and rank count.
The flat fixture separately isolates its two rate hierarchies and groups
its four controls, as specified in the flat Stokes suite. Each group retains the slow label,
1800-second watchdog, MPI processor count and common pyramid resource lock.
No refinement level, finite-element degree, quadrature order, solver setting
or numerical assertion is changed. The distinction bounds process lifetime,
not the mathematical workload. An isolated OpenMP two-rank hexahedral P6/P5
patch passed with sampled aggregate participant RSS of approximately
1.01 GiB; the earlier combined process exceeded the 3 GiB validation limit
after preceding rate studies. This comparison motivates fresh-process
validation but does not identify the allocation or retention mechanism.
The current-source combined rate process additionally exceeded the same
memory guard at OpenMP MPI rank four on hexahedra. An isolated full P3/P2
hierarchy passed with sampled aggregate participant RSS of approximately
2.57 GiB, retaining all three levels and the same solver. This supports
separating the two rate processes without reducing their mathematical scope;
it does not identify a particular retained allocation. Earlier complete
matrix records are preserved, but current-source verification of the revised
process partition remains a separate gate and is not yet certified. Hosted
CI certification is also separate. Sampled RSS is a resource-guard observation,
not a continuous peak-memory measurement or a performance benchmark.

## Higher-order mixed fields on quadratic approximated geometry

The `ApproximatedP3Q2` extension uses velocity/pressure degrees $3/2$ on
the same quadratic-interpolated sine-map domain. The manufactured fields are

$$
u_\star(x)=x_{d-1}^4e_0,\qquad
p_\star(x)=x_0^3-\frac14,\qquad
f(x)=(-12x_{d-1}^2+3x_0^2)e_0,
\qquad d\in\lbrace2,3\rbrace.
$$

They satisfy $-\Delta u_\star+\nabla p_\star=f$ and
$\operatorname{div}u_\star=0$. The physical velocity trace is prescribed;
the existing global mean-pressure constraint is retained. The shear map
preserves $x_0$ and volume, so the physical cubic pressure has zero mean
on both domains and its geometry defect vanishes analytically. That defect
uses the existing absolute $10^{-10}$ budget, not a relative rate fit.

Three levels are $n=9,17,33$ for Triangle, $n=3,5,9$ for
Quadrilateral and $n=3,4,5$ in 3D, with
$h=(n-1)^{-1}$. The latter retains the module's three-dimensional mixed
hierarchy strategy; it does not change any older pair's refinement levels.
The original Triangle hierarchy $n=3,5,9$ gave pressure
$L^2$ slopes $4.05,3.42$ and $H^1$-seminorm slopes $3.14,2.55$.
These positive, strictly decreasing errors did not satisfy the stated
two-sided finite-resolution windows. An isolated finer native Triangle
study gave slopes $3.19,3.10$ and $2.31,2.19$, respectively, within
the original windows and approaching the expected orders $3$ and $2$.
The finer hierarchy therefore tests the resolved regime without widening
an acceptance margin. This isolated observation does not identify a separate
superconvergence theorem or certify another geometry or backend; the complete
registered matrix is verified separately below. The fixed-mesh
patch, quadrature and wrong-viscosity controls retain $n=3$ or $n=5$.
Under the stated regularity, map and mixed-stability hypotheses, represented
and lifted-field velocity errors have $L^2/H^1$ orders $4/3$ and geometry
errors have orders $3/2$. For norm index $j\in\lbrace0,1\rbrace$, the
shared mixed-order assertion checks both norm triangle inequalities and

$$
T_{j,f}<T_{j,c},\qquad
T_{j,f}\le F_{j,c}\rho^{4-j-\delta_j}
             +G_{j,c}\rho^{3-j-\delta_j}+10^{-11},
\qquad\rho=\frac{h_f}{h_c},\quad
\delta_0=0.55,\quad\delta_1=0.45.
$$

Every adjacent interval retains the independent velocity component windows.
Represented, lifted-field and total pressure errors retain strict decay and
the conditional $3/2$ rate floors with the same $0.55/0.45$ margins.
An upper slope window is not imposed on this coupled pressure error: its
velocity-induced contribution can dominate a finite hierarchy and decay
faster than the pressure approximation order. The independently separated
pressure-data study below retains two-sided pressure windows.
Divergence defects remain
bounded by $\sqrt d$ times the corresponding velocity derivative norm;
no exact divergence order is asserted. These are finite-hierarchy acceptance
policies, not a proof of uniform inf-sup stability or geometry dominance.

The uniform tetrahedral family also supplies an independent interpolation
reference for the pressure. On a coordinate slab $x_0\in[a,a+h]$, let
$t=(x_0-a)/h$. The quadratic nodal interpolant $I_2p_\star$ satisfies

$$
p_\star-I_2p_\star=h^3t\left(t-\frac12\right)(t-1),\qquad
\lVert p_\star-I_2p_\star\rVert_{L^2(\Omega_h)}
=\frac{h^3}{\sqrt{840}},\qquad
\lvert p_\star-I_2p_\star\rvert_{H^1(\Omega_h)}
=\frac{h^2}{\sqrt{20}}.
$$

These identities use the slab's vertex/edge-midpoint interpolation nodes,
preservation of $x_0$, and the unit-Jacobian shear map. They give a separate
approximation scale for the represented pressure space. Mixed-solve
pressure errors additionally depend on velocity approximation and finite
mixed stability; their observed rates cannot be inferred from interpolation
alone.

The additional `ApproximatedP3Q2SeparatedPressureRates` gate resolves this
coupling on the same mesh, spaces, quadrature and mean-pressure constraint.
At zero pressure offset, the original data are decomposed as

$$
(u_\star,p_\star)=(u_\star^{P},p_\star^{P})
                    +(u_\star^{V},p_\star^{V}),\qquad
(u_\star^{P},p_\star^{P})=\left(0,x_0^3-\frac14\right),\qquad
(u_\star^{V},p_\star^{V})=(x_{d-1}^4e_0,0).
$$

The corresponding body forces are $f^P=3x_0^2e_0$ and
$f^V=-12x_{d-1}^2e_0$. Linearity of the constrained discrete Stokes problem
then gives $p_h=p_h^P+p_h^V$ and $u_h=u_h^P+u_h^V$. For pressure defects
$e_h^P=p_h^P-p_\star$ and $e_h^V=p_h^V$, the actual combined error satisfies

$$
\lVert p_h-p_\star\rVert_j^2
=\lVert e_h^P\rVert_j^2+\lVert e_h^V\rVert_j^2
 +2\langle e_h^P,e_h^V\rangle_j,\qquad j\in\lbrace0,1\rbrace,
$$

where $\langle a,b\rangle_0=\int_{\Omega_h}ab$ and
$\langle a,b\rangle_1=\int_{\Omega_h}\nabla a\cdot\nabla b$.
The lifted variants use the exact-domain measure and lifted gradients instead.
Consequently, a finite-resolution slope of the combined pressure error need
not equal the approximation slope of its pressure-data component. A bound
$E(h)\le Ch^r$ does not supply the complementary lower bound needed to
exclude a faster observed slope.

The registered component gate retains the original coupled solve at every
level. Logical cell/DOF layouts are compared exactly before numerical
coefficient superposition is checked using
$\lVert c_h-c_h^P-c_h^V\rVert_2/\max(1,\lVert c_h\rVert_2)$ against the
existing dimensionless $10^{-9}$ roundoff budget, separately for velocity
and pressure. This coefficient check is rank-local, including
available halo coefficients, and introduces no collective. Represented,
lifted-field and total pressure norms are independently computed with the
existing owned-cell/global-norm contract. Both norm triangle inequalities
are checked for each decomposition. All three data sets retain strict error
reduction and rate floors $3-0.55$ and $2-0.45$; the pressure-data component
additionally retains the upper windows $3+0.55$ and $2+0.45$.
The same three-level geometry-specific hierarchy is used. The original coupled
study retains its velocity component windows, pressure decay and pressure rate
floors; the separated study additionally resolves the pressure-data order.
Neither registration alone certifies the complete numerical matrix. These
finite-hierarchy criteria do not establish an asymptotic lower bound for every
coupled pressure error or a uniform inf-sup constant.

At $n=5$, assembly order $12\to16$ and norm order $14\to18$ are varied
separately; positive error components retain the $10^{-6}$ relative
sensitivity budget. Direct factorization retains the existing independently
recomputed residual budget; no iterative tolerance variation is introduced.
At $n=3$, the physical affine shear patch must retain the existing $10^{-9}$
reproduction budgets. A separate viscosity control changes the operator
coefficient from one to two while retaining the original source and trace;
pressure errors must exceed the existing $0.1/1$ floors. The geometry defects
remain unchanged. Finally, physical volume and cubic-pressure mean are
integrated independently with the existing $10^{-12}$ gauge budget.

The original coupled hierarchy, the separated-pressure hierarchy and the four
fixed-mesh controls run in distinct `Rates`, `SeparatedRates` and `Controls`
process groups. Each rate group retains all three levels and both intervals;
process isolation changes neither the discrete problem nor its acceptance
criteria. These registrations cover six applicable geometries, native and
real-PETSc local execution, MPI ranks one through four and both thread
configurations. They retain slow labels, processor counts, pyramid locks and
1800-second watchdogs. The complete six-case finite matrix is locally verified:
216 process groups select 432 configurations with 864 successful rank-local
reports. Both thread configurations were compiled freshly for native and
real-PETSc execution. Source identity, dependency freshness, exact case
selection, runtime results and resource guards were audited independently.
The retained dependency libraries and finite mesh sequences do not constitute
a complete dependency-library rebuild, hosted-CI certification or a uniform
stability proof.

## Cubic-geometry pressure forward accuracy

Representable affine pressure on cubic geometry separates algebraic solve
error from approximation error. A small coefficient residual alone does not
bound the physical pressure-gradient error independently of conditioning.
The PETSc cubic-geometry fixture therefore requests two outer residual
corrections using the retained LU/MUMPS factors. For the assembled constrained
system $Ac=b$, unit-damped Richardson updates are

$$
c^{(0)}=c_{\mathrm{assembled}},\qquad
r^{(j)}=b-Ac^{(j)},\qquad
c^{(j+1)}=c^{(j)}+\widetilde A^{-1}r^{(j)},\qquad j=0,1,2,
$$

where $\widetilde A^{-1}$ denotes the computed factor solve and
$c_{\mathrm{assembled}}$ is the initial vector retained by the KSP wrapper.
The three updates comprise an initial factor-preconditioned update and two
corrections. This fixed count is a forward-accuracy policy, not a convergence
theorem. Solver type, norm policy and iteration count are checked exactly;
coefficient residual and physical pressure errors retain independent checks.
The MPI operation acts on the global distributed linear system. It adds no
collective to local mesh queries and does not rely on internal MUMPS refinement
for distributed solutions. Other geometry fixtures retain the default PREONLY
solve. Incompatible runtime overrides fail the explicit policy checks.

The `CubicPressureForwardAccuracy` regression uses $n=3$ on all six
applicable geometries, in Local and MPI rank-one-through-four contexts.
Represented, lifted-field and total pressure errors must each satisfy
$E_0<10^{-11}$ and $E_1<10^{-11}$, using dimensionless absolute budgets.
This fixed-mesh regression supplements the three-level $n=3,5,9$ geometry
study without relaxing its $10^{-9}$ reproduction budget. Fresh verification
of the latest-develop integration remains pending; prior detached-tree results
do not certify the newly combined implementation.
