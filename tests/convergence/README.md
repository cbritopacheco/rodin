# Numerical convergence tests

This suite examines discretization error as a finite-element space is refined.
Its reference is an exact manufactured solution or a field represented exactly
by the discrete space. It therefore provides **code-verification evidence** for
a specified formulation, space, geometry and backend; it does not compare a
model with physical observations. Fixed-mesh regressions belong in
`tests/manufactured/`, while independently specified NAFEMS benchmark
problems belong in `tests/nafems/`.

## Continuous and discrete problems

Let $\Omega\subset\mathbb{R}^d$, with
$d\in\lbrace 1,2,3\rbrace$, be the physical domain. Let
$u:\Omega\to\mathbb{F}^m$ be the exact field, where
$\mathbb{F}\in\lbrace\mathbb{R},\mathbb{C}\rbrace$ and $m$ is the number of
components. For a linear problem, a representative weak formulation is: find
$u\in g+V_0$ such that

$$
a(u,v)=\ell(v)\qquad\text{for every }v\in V_0.
$$

Here $V_0$ is the test space incorporating homogeneous essential-boundary
conditions, $g$ is a lifting into the domain of the prescribed boundary
trace, $g+V_0$ is the corresponding affine trial space, $a$ is the
continuous bilinear or sesquilinear form, and $\ell$ contains the forcing and
natural-boundary data. A nonlinear problem is stated instead as $F(u;v)=0$
for every $v\in V_0$; mixed problems specify all
fields, spaces, constraints and gauges. For each manufactured case, the exact
$u$ is selected first and the source and boundary data are derived from the
stated continuous equations. The domain, constitutive coefficients, forcing,
boundary partition and model parameters remain fixed across the refinement
sequence. An exactly representable polynomial or constant gives a separate
patch test; it does not replace a rate test for a non-polynomial field.

On refinement level $i$, let $\mathcal{T}_i$ be a mesh, $p_i$ the field degree,
and $V_i$ the corresponding discrete space. The computed field
$u_i\in g_i+V_i$ solves the discrete problem, including the documented
quadrature and algebraic-solver choices. The test changes only the declared
refinement variables. In particular, a p study keeps the physical mesh fixed,
while an h study keeps the element family and degree fixed.

## Error measurements

For conforming scalar or vector fields, the measured quantities are norms of
the *difference*, not differences of norms. On the physical domain
$\Omega_i$ represented by the mesh, they are defined by

$$
E_{0,i}=\left(\int_{\Omega_i}\lVert u-u_i\rVert^2\thinspace\mathrm{d}x\right)^{1/2},
\qquad
E_{1,i}=\left(\int_{\Omega_i}\lVert Du-Du_i\rVert_F^2\thinspace\mathrm{d}x\right)^{1/2}.
$$

Here $D$ is the spatial gradient for a scalar field and the spatial Jacobian
for a vector field. The value norm is the absolute value for a scalar and the
Euclidean norm for a vector; complex norms use conjugate magnitude. Thus
$E_{0,i}$ is an L2 error and $E_{1,i}$ is an H1-seminorm error for conforming
fields. Discontinuous P0 fields are assessed only in L2. A mixed formulation
reports errors for each field separately; pressure errors are measured after
the documented gauge has been imposed.

The norms are integrated independently of the assembled residual. In a local
calculation, each cell contribution is evaluated at physical quadrature
points $x_q=\Phi_K(\hat x_q)$ through its map
$\Phi_K:\hat K\to K\subset\Omega_i$, for example:

$$
E_{0,i}^2\approx\sum_{K\in\mathcal{T}_i}\sum_q
w_q\thinspace|\det D\Phi_K(\hat x_q)|\thinspace
\lVert u(x_q)-u_i(x_q)\rVert^2.
$$

The exact field is evaluated at those points; discrete coefficients are not
treated as point values. In MPI runs, squared contributions from owned cells
are summed globally before the square root. The norm-quadrature order, solve
tolerance and any geometry approximation are stated per suite. If the
represented domain $\Omega_i$ differs from the intended $\Omega$, the above
field error alone does not measure the domain error. In that case, the exact
field must be defined on $\Omega_i$ or compared through a stated pullback,
and a geometry-map check or an independently known map is also required.

## Refinement paths and rate interpretation

Each directory isolates a different question:

- `h`: $h_i$ decreases while the element family and $p_i=p$ remain fixed;
- `p`: $p_i$ increases on a fixed mesh and geometry;
- `hp`: both $h_i$ and $p_i$ change along a stated path;
- `isoparametric`: the geometry map $\Phi_i$ and the field approximation are
  checked together, with their errors distinguished where possible.

For h refinement, an error quantity $E_i\gt 0$ gives the adjacent-level order

$$
r_i^{(h)}=\frac{\log(E_{i-1}/E_i)}{\log(h_{i-1}/h_i)}.
$$

`UniformGridHierarchy` uses the nominal coordinate spacing
$h_i=1/(n_i-1)$ for $n_i$ points per axis. On a fixed cell family, the
geometry-dependent ratio between this spacing and element diameter cancels
between levels. For a sufficiently regular exact solution of a stable,
consistent conforming degree $p$ elliptic problem on regular meshes, the
standard expectations are $E_{1,i}=O(h_i^p)$ and, when the required dual
regularity holds, $E_{0,i}=O(h_i^{p+1})$. These orders are not asserted
unconditionally for singular solutions, mixed boundary corners, or
under-resolved geometric maps. Smooth-field P0 projection instead has an
L2 error of order one; P0g reproduces constants exactly and has no
nonconstant h rate.

For fixed-mesh p refinement, the observed adjacent-degree decay is

$$
\alpha_i=\frac{\log(E_{i-1}/E_i)}{p_i-p_{i-1}}.
$$

An analytic exact field may admit $E_i\le C\exp(-b p_i)$ over a range in
which geometry, quadrature, conditioning and floating-point errors do not
dominate. Positive measured $\alpha_i$ over a finite degree range supports
that case-specific observation; it is not a proof of exponential decay at
arbitrary degree. The Poisson, Helmholtz, conductivity, and coupled
reaction–diffusion p suites check every adjacent interval from degree 1
through degree 4. Native Stokes checks velocity/pressure degree pairs
$(2,1)\to(3,2)\to(4,3)$, measuring both fields separately.

For hp refinement, a log ratio against $h_i$ describes only the stated
combined path $(h_i,p_i)$; it is not a fixed-degree h order. In an
isoparametric study, the reference-to-physical map $\Phi_i$ has its own
approximation order and regularity. An exact-map patch test or a separately
measured geometry error is needed before a field-error rate can be attributed
to the field space.

## Acceptance and reproducibility

At least three discretizations are required for a new rate claim so that two
successive intervals are tested. The continuous data and measured quantity
are held fixed. Each interval checks finite positive errors, strict error
reduction, and a problem-specific rate interval or lower bound justified by
the expected regularity and numerical resolution. Exact-reproduction tests
instead use an absolute error bound compatible with quadrature and
floating-point roundoff. Solver convergence and sufficiently small algebraic
error are checked separately. When a rate approaches an assertion bound,
quadrature order and solver tolerance are varied to identify numerical
contamination; unexplained non-monotonicity or superconvergence is not
absorbed by widening the bound.

The present rate studies meet this minimum: h, hp and isoparametric paths
contain three meshes, while analytic p studies contain four degrees,
except Stokes, which contains three velocity/pressure pairs.
One-mesh patch tests and P0g exact-reproduction tests are not rate studies.

Each suite README states the equations, exact fields, derived data, discrete
spaces, mesh or degree sequence, quadrature and solver settings, norm domains,
expected rates and their hypotheses, assertion bounds, tested geometries and
backends, and known limitations. Failure traces report the individual errors
and rates. A negative control removes or perturbs a relevant term, parameter,
boundary condition or extraction rule and confirms that the acceptance test
fails; the correct formulation is then restored. The seven `UniformGrid`
cell types—segment, triangle, quadrilateral, tetrahedron, pyramid, hexahedron
and wedge—are tested where the formulation is meaningful, with any exclusion
recorded. Sequential, OpenMP, PETSc and MPI results are distinct evidence;
none implies the others. The PETSc conductivity suite additionally checks
nonzero-trace affine/P2 polynomial patches, P1/P2 rates on all seven geometries
at 1–4 ranks, a coefficient-omission negative control, and quadrature/solver
sensitivity. Its [suite specification](h/PETScConductivity/README.md) records
the exact data, levels, acceptance bounds, and backend scope.

## Coverage and execution

Current implemented coverage is summarized below; a dash means that no suite
exists yet.

| Context | h | p | hp | Isoparametric |
| --- | --- | --- | --- | --- |
| Poisson | P1–P3, boundary variants; PETSc local and MPI P1/P2 | P1/P2 patch; P1→P2→P3→P4 analytic | P1–P3 | P1/P2 on exact P2 maps; native local and real-PETSc local/MPI |
| Complex Helmholtz | P1/P2; native-complex PETSc local/MPI P1/P2 | P1–P4 | P1–P3 | P1/P2 fields on exact curved P2 maps; native and complex-PETSc local/MPI |
| Linear elasticity | Vector P1/P2, displacement and traction variants; PETSc local/MPI Dirichlet P1/P2 | Analytic vector P1→P2→P3→P4 | Analytic vector P1–P3 | P1/P2 displacement, strain and stress on exact P2 maps; native local and real-PETSc local/MPI |
| Stokes | Taylor–Hood P2/P1/P0g; native and PETSc local/MPI | Velocity/pressure pairs $2/1\to3/2\to4/3$; native and PETSc local/MPI | Analytic pairs $2/1\to3/2\to4/3$; native and PETSc local/MPI | P2/P1/P0g on exact curved P2 maps; native local and real-PETSc local/MPI |
| Variable conductivity | P1/P2; PETSc local/MPI P1/P2 | P1/P2 patch; P1→P2→P3→P4 analytic | P1–P3 | P1/P2 on exact P2 maps; native local and real-PETSc local/MPI |
| Coupled reaction–diffusion | P1/P2; PETSc local/MPI P1/P2 | P1→P2→P3→P4 analytic | Analytic two-field P1–P3 | P1/P2 on exact P2 maps; native local and real-PETSc local/MPI |
| Nonlinear Poisson | $P_1/P_2$; native and real-PETSc SNES local/MPI | Analytic P1→P2→P3→P4; tangent controls | Analytic P1–P3; tangent controls | P1/P2 on exact P2 maps; native Newton and real-PETSc SNES local/MPI |
| P0 projection | Real/complex scalar and vector, first-order L2 | Not applicable to fixed degree | Not applicable to fixed degree | Real/complex scalar/vector on exact P2 maps; native and PETSc local/MPI; cell-moment controls |
| P0g | Exact real/complex scalar and vector constants | Not applicable | Not applicable | Curved constant reproduction and analytic global means; no h-rate |
| Geometry approximation / Poisson patch | Map/derivative rates at fixed geometry degrees $q=1,2,3$ | Not yet implemented | Not yet implemented | Sine-map approximation and affine Poisson patches with $p=q$; native and real-PETSc local/MPI |

The main refinement sequences can be read with $n$ grid points per coordinate
axis, $h=1/(n-1)$, and field degree $p$:

- h: P1 commonly uses `n=5→9→17` ($h=1/4\to1/8\to1/16$); P2/P3 commonly
  use `n=3→5→9`. P0 projection and Taylor–Hood Stokes use `n=3→5→9`.
- p: Poisson, complex Helmholtz, variable conductivity, coupled
  reaction–diffusion, and linear elasticity check `p=1→2→3→4` on a fixed
  `n=2` mesh.
  Nonlinear Poisson instead uses a fixed `n=3` mesh. Stokes uses `n=3`
  with velocity degrees $2\to3\to4$ and pressure degrees $1\to2\to3$.
- hp: Poisson, Helmholtz, variable conductivity, linear elasticity,
  coupled reaction–diffusion, and nonlinear Poisson use
  `(n,p)=(2,1)→(3,2)→(5,3)`.
  Stokes instead uses $(n,k)=(3,2)\to(4,3)\to(5,4)$, with pressure
  degree $k-1$ and a global pressure-mean multiplier.
- isoparametric: curved P2 Poisson and conductivity use `n=5→9→17`
  in 1D/2D and `n=3→4→5` in 3D. A separate P1 geometry-map error
  study uses `n=5→9→17` in 1D/2D and `n=3→5→9` in 3D.
  The shared curved Poisson/conductivity backend suite uses `n=5→9→17`
  for P1 and `n=3→5→9` for P2 on all seven geometries; its
  [suite specification](isoparametric/Diffusion/README.md) states the
  physical fields, coefficient/source pairing, residuals and wrong-operator
  controls. Geometry remains exact and fixed in these field-rate studies.
  Curved complex Helmholtz instead uses `n=5→9→17` for P1 and
  `n=3→5→9` for P2 on every positive-dimensional geometry.
  Curved linear elasticity uses the same P1/P2 sequences and separately
  measures displacement, strain and stress; its
  [suite specification](isoparametric/LinearElasticity/README.md) states the
  physical-coordinate patch and omitted-volumetric-term controls.
  Curved Taylor–Hood Stokes uses `n=3→5→9` in 2D and `n=3→4→5`
  in 3D; its [suite specification](isoparametric/Stokes/README.md) states
  velocity/pressure rates, the physical pressure gauge, and divergence controls.
  Curved coupled reaction–diffusion uses `n=5→9→17` for P1 and
  `n=3→5→9` for P2 on all seven geometries; its
  [suite specification](isoparametric/ReactionDiffusion/README.md) states
  both component rates, representable patches and coupling controls.
  Curved nonlinear Poisson uses the same P1/P2 sequences; its
  [suite specification](isoparametric/NonlinearPoisson/README.md) describes
  nonzero-trace lifting, residual/tangent consistency and independent controls.
  Curved P0 projection uses `n=5→9→17` on all seven geometries; its
  [suite specification](isoparametric/P0Projection/README.md) describes
  physical cell moments, interpolation rejection and analytic P0g means.
  Nonpolynomial geometry approximation uses `n=3→5→9` separately for
  geometry degrees $q=1,2,3$, requiring both adjacent map-error orders
  $q+1$ and derivative-error orders $q$. Matched-degree affine Poisson
  patches measure field error on each approximated domain, not solution
  error on the exact domain. The
  [suite specification](isoparametric/GeometryApproximation/README.md)
  defines the independent sine-map oracle, logical chart correspondence,
  quadrature sensitivity, and separate wrong-map/wrong-trace controls.

These are rate-test levels, not one-mesh patch-test levels. The P1
mixed-traction elasticity test on tetrahedra uses `n=9→17→33` because its
coarsest mesh is pre-asymptotic. This is a coverage inventory, not a claim
that every supported space or physical context is certified; the suite
READMEs and test sources give all case-specific sequences and bounds.
Distributed P2 Poisson now checks three-level rates on all seven geometries
and one to four ranks. The ownership regression for the previously missing
tetrahedral boundary constraints and its mathematical patch test are recorded
in the `h/PETScMPIPoisson` README.
The distributed P2 identification regression checks every geometry; its
tetrahedral solve on two to four ranks also compares equivalent value and
affine-identification boundary conditions. Homogeneous P1 and affine P2
identification additionally have three-level rate checks on all seven
geometries and one to four ranks.

Configure with `-DRODIN_BUILD_CONVERGENCE_TESTS=ON` and run with
`ctest --test-dir build/tests -L convergence --output-on-failure`.

`Convergence.h` contains the refinement-independent layer shared by these
modules: unit-box grid construction and boundary partitioning, direct L2/H1
error integration, error histories, and algebraic or exponential rate
calculation. Refinement-specific directories add only the machinery unique to
their refinement axis. `StokesData` supplies common exact fields and sources;
`StokesProblem` supplies the native mixed solve, solver/residual and pressure
gauge checks, and separate field-error measurements for h, p, and hp.
`PETScStokesProblem` supplies the same workload with PETSc local/MPI
storage and direct factorization, sharing the continuous data and error oracles.
`CurvedGeometry` retains the unwarped vertices and installs the quadratic map
on cells, traces, and MPI halos. The
[curved Helmholtz specification](isoparametric/Helmholtz/README.md) distinguishes
superparametric P1 fields from strictly isoparametric P2 fields and records
the fixed-domain field rates and independent controls.

The CI convergence job runs the full local suite twice: once with sequential
assembly and once with OpenMP assembly (`RODIN_MULTITHREADED=OFF/ON`). A
separate real-PETSc sequential/OpenMP matrix checks local-context and
distributed P1/P2 Poisson,
conductivity, full-Dirichlet vector linear elasticity, and coupled
reaction–diffusion convergence with
PETSc assembly and CG, plus P2/P1 Stokes with a P0g pressure gauge and
LU/MUMPS factorization. Stokes additionally has PETSc local/MPI p and hp
paths through velocity/pressure degrees $4/3$ on all six applicable geometries.
Real-PETSc nonlinear Poisson additionally exercises SNES with CG/Jacobi
tangents on local and distributed $P_1/P_2$ spaces;
the [suite specification](h/PETScNonlinearPoisson/README.md) states the
residual/tangent callback checks, independently reassembled final residual,
and analytically bounded missing-reaction control.
The distributed suite uses mesh
families partitioned across one to four MPI ranks and globally reduced norms;
it is a separate check from the local-context suites.
A dedicated complex-PETSc job checks local and MPI P1/P2 Helmholtz with
sequential and OpenMP assembly configurations. It requires a native-complex
PETSc scalar build; a real PETSc installation does not register these tests.
The [Helmholtz specification](h/PETScHelmholtz/README.md) states the
coercive wave-number regime, complex manufactured data, patch and rate
oracles, and rank matrix. The
[reaction–diffusion hp specification](hp/ReactionDiffusion/README.md)
states the combined path and separate component errors; its path rates are
not fixed-degree h orders.
The [PETSc Stokes specification](h/PETScStokes/README.md) records the six
applicable geometries, pressure-sensitive negative control, global divergence
norm check, factorization settings, and separate velocity/pressure orders.
The nonlinear Poisson [p](p/NonlinearPoisson/README.md) and
[hp](hp/NonlinearPoisson/README.md) specifications record degree/path
improvement, higher-order residual/tangent checks, and numerical budgets.
The native Stokes [p](p/Stokes/README.md) and [hp](hp/Stokes/README.md)
specifications record mixed-field improvement, highest-degree polynomial
patches, pressure-sensitive viscosity controls, quadrature sensitivity, and
the exact coarse tensor-cell pressure-rank obstruction. They explicitly
distinguish finite-workload convergence from a uniform inf-sup theorem.
Their PETSc [p](p/PETScStokes/README.md) and
[hp](hp/PETScStokes/README.md) counterparts use the same levels and
acceptance bounds, globally reduced owned-cell norms, and ranks 1–4.
The p suite also checks known P4/P3 norm values on sparse/empty partitions.

## Remaining verification work

The coverage above includes the baseline merged in PR #333 and subsequent
suite additions. The following extensions
are planned; their presence in this workplan does not imply an implemented or
passing test. Completion is assessed per formulation, space, geometry,
refinement path, and backend, rather than by the presence of a directory.

| Priority | Extension | Required evidence |
| --- | --- | --- |
| 1 | PETSc local and MPI PDE coverage: remaining boundary/refinement variants of Poisson, Helmholtz, conductivity, linear elasticity, Stokes, coupled reaction–diffusion, and nonlinear Poisson | Independently integrated field errors and expected rates on each meaningful geometry; supported scalar/backend configurations stated explicitly; owned-cell global norms in MPI |
| 2 | Curved Poisson, conductivity, Helmholtz, linear-elasticity, Stokes, reaction–diffusion and nonlinear Poisson boundary/degree extensions | Physical-coordinate manufactured data, independent norm integration, regular maps, and case-specific field rates or exact reproduction |
| 3 | Smooth PDE field rates on approximated nonpolynomial geometry and exact-domain comparisons | Geometry degrees 1–3 now have map/derivative rate tests and matched-degree affine Poisson patches; remaining work must specify the solution comparison domain or pullback and separate field error from domain error |
| 4 | Remaining complex-vector and high-order structural combinations supported by the library | Exact index round trips, unique ownership, halo/incidence completeness, boundary and identification selection, and SubMesh restriction across geometries and rank counts |
| Last | Independent NAFEMS benchmarks, after the convergence/structural/backend batches | Authoritative specifications and usable reference data; independently defined quantities of interest, units, error budgets, and mesh studies in `tests/nafems` |
| Separate PR | Assembly performance across existing physical contexts, geometries, spaces, and backends ([PR #356](https://github.com/cbritopacheco/rodin/pull/356)) | Isolated stage timings, reproducible workload metadata, verified assembled operators, and controlled thread/rank scaling in `tests/benchmarks`; tracked independently from convergence certification |

### Coupled reaction–diffusion batch and continuation

The current implementation adds PETSc local/MPI P1/P2 h-convergence for
a two-component reaction–diffusion system on the unit box. For
$u=(u_1,u_2)^T$, use positive diffusion coefficients $\kappa_i$ and a
symmetric positive-definite reaction matrix $R$ with nonzero off-diagonal
entries. The [suite specification](h/PETScReactionDiffusion/README.md)
states the selected coefficients, fields, and acceptance bounds.
Manufactured sources are derived componentwise from

$$
f_i=-\kappa_i\Delta u_i+\sum_{j=1}^{2}R_{ij}u_j,
\qquad i\in\{1,2\},
$$

with full manufactured Dirichlet traces. Coupled-space and assembly/solver
support is exercised through two scalar H1 trial/test pairs in a coupled
problem. These equations do not imply that every backend combination is
supported.

- The suite includes exactly representable affine P1 and quadratic P2
  patches, followed by smooth fields with nonzero cross-coupling and
  independently derived gradients.
- L2 and H1-seminorm errors are measured separately for each component, so
  that agreement in one field cannot conceal failure in the other.
- Three levels are used: `n=5→9→17` for P1 and `n=3→5→9` for
  P2. Both intervals require error reduction and the expected L2/H1 orders
  $2/1$ and $3/2$, respectively, under the stated regularity assumptions.
- Entries cover all seven positive-dimensional geometries locally and with
  MPI ranks 1–4, including owned-cell global norms and empty-rank cases.
- A negative control omits cross-coupling while retaining the correct
  sources and traces; the oracle must reject that incorrect operator.
  Quadrature and solver-tolerance sensitivity checks have stated budgets.
- Equations, levels, bounds, and backend exclusions are documented in the
  suite specification; verification evidence belongs in the PR. Planned,
  implemented, locally verified, and CI-certified coverage remain distinct
  states.

After the reaction–diffusion h/hp, complex Helmholtz h, PETSc Stokes h,
and native nonlinear Poisson p/hp batches, real-PETSc/SNES nonlinear Poisson
adds the local/MPI h path. Native Stokes p/hp fills the remaining native
p/hp entries in the table. Priority 1 continues with missing
PETSc mixed-boundary/refinement variants; Stokes now has its local/MPI
p/hp counterparts. Priorities 2–4 address curved fields, approximated
nonpolynomial geometry, and exact-index MPI structural combinations.
The curved complex Helmholtz batch supplies P1/P2 field rates on exact
P2 maps, with native local and complex-PETSc local/MPI counterparts.
Independent NAFEMS implementation is the last phase, after these
verification batches. Assembly benchmarks advance separately in PR #356;
completion of that PR is not a prerequisite for implementing convergence
suites. Each convergence batch still requires its own passing CI evidence
before being described as CI-certified.

Backend extensions require an initial support check: a mathematically meaningful
formulation does not establish that every solver, scalar type, or assembly path
supports it. Unsupported combinations must be recorded explicitly. Structural
MPI regressions for a space do not establish PDE convergence for that space.
The existing structural rank matrix is 1, 2, 3, 4, and 8; distributed PDE
studies currently use 1–4 ranks. Extensions must state their own rank matrix,
including empty-rank or sparse-selection cases where relevant.
For mixed spaces, pressure-nullspace and inf-sup verification remain
distinct from convergence of selected manufactured fields. The Stokes
coarse-space rank regression establishes an obstruction, not stability of
every larger pair; uniform stability evidence for the pyramid/wedge
families remains unresolved.

Each new suite must document its continuous problem, derived data, discrete
spaces, geometry families, refinement levels, quadrature, algebraic error
budget, measured norms, expected rates, and exclusions. At least three levels
and two adjacent intervals are required for rate evidence. Negative controls
and quadrature/solver sensitivity checks must establish that the assertions
detect the targeted defect and that numerical integration or algebraic error
does not determine the observed rate. Expensive resolved hierarchies belong
in explicitly timed slow tests with the required backend coverage retained.

### Assembly performance workplan

Performance work is tracked in [PR #356](https://github.com/cbritopacheco/rodin/pull/356)
in the existing `tests/benchmarks` module, separately from numerical
convergence assertions. The following describes the intended final scope,
not a claim of complete implemented coverage. The physical contexts are
Poisson, variable conductivity, complex Helmholtz, vector linear elasticity,
coupled reaction–diffusion, Taylor–Hood Stokes, and nonlinear Poisson. Scalar
mass/projection forms additionally exercise supported real/complex scalar and
vector P0/P0g spaces; these are assembly workloads, not additional PDE rate
claims. H1 workloads begin with P1/P2 and extend through the orders already
covered by the convergence suites, using stable mixed pairs for Stokes.

Every meaningful formulation is to be exercised on segment, triangle,
quadrilateral, tetrahedron, pyramid, hexahedron, and wedge geometries.
Point/0D assembly is included only where a discrete form is meaningful;
unsupported or mathematically inapplicable combinations are recorded rather
than counted as measured coverage. Affine and supported curved geometry paths
are distinguished. Backend coverage includes Eigen sequential/OpenMP and
PETSc sequential/OpenMP/MPI assembly, with unsupported scalar/backend paths
identified explicitly. Thread and rank counts are varied independently.

For a fixed discrete form, distinguish setup (mesh, space, quadrature, sparsity
and allocation), element-kernel evaluation, global insertion/accumulation,
constraint application, and matrix/vector finalization. Measure both complete
assembly and isolated stages where instrumentation permits; nonlinear cases
measure residual and tangent assembly at a prescribed state. Solves and error
integration are timed separately and excluded from assembly throughput.
Cold construction and warmed repeated assembly are separate experiments;
repeated assembly must not silently accumulate previous contributions.

Report wall time $T_A$ in seconds, owned-cell throughput $N_K/T_A$, and time
per global DOF $T_A/N_D$, with global cell count $N_K$, global DOF count $N_D$,
matrix nonzeros, degree, quadrature order and point count, geometry, and field
components. MPI wall time is the maximum elapsed time over ranks, including
required communication/finalization. Strong scaling holds the global problem
fixed; weak scaling holds the owned workload approximately fixed per rank.
Both state the partition imbalance and actual thread/rank configuration.

At least three mesh sizes are required per size study, beginning with the
existing convergence hierarchies where feasible. Record build type, compiler,
dependencies, hardware, thread affinity, sanitizer status, repetitions, and
timing dispersion. Run isolated workloads without concurrent CTest jobs and
check registration-order/cache effects. Sanitized and optimized-build timings
are not interchangeable baselines. Numerical checks on operators, loads,
constraints, and resulting solutions accompany each benchmark; an optimization
requires baseline-equivalent numerical behavior, including solver iterations
and residuals. Performance regression thresholds are introduced only after
repeatability and variance have been established on a controlled runner.

Darcy remains deferred. Fixed-degree P0 has no p-refinement family; P0g has
no nonconstant approximation rate. These are mathematical exclusions, not
missing certification tasks. Additional physics can be added after the
existing formulations have their intended refinement and backend coverage.
