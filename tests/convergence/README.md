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
through degree 4.

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
contain three meshes, while analytic p studies contain four degrees.
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
none implies the others.

## Coverage and execution

Current implemented coverage is summarized below; a dash means that no suite
exists yet.

| Context | h | p | hp | Isoparametric |
| --- | --- | --- | --- | --- |
| Poisson | P1–P3, boundary variants; PETSc local and MPI P1/P2 | P1/P2 patch; P1→P2→P3→P4 analytic | P1–P3 | Curved P2 |
| Complex Helmholtz | P1/P2 | P1–P4 | P1–P3 | — |
| Linear elasticity, Stokes | Implemented | — | — | — |
| Variable conductivity | P1/P2 | P1/P2 patch; P1→P2→P3→P4 analytic | P1–P3 | Curved P2 |
| Coupled reaction–diffusion | P1/P2 | P1→P2→P3→P4 analytic | — | — |
| Nonlinear Poisson | P1/P2 | — | — | — |
| P0/P0g projection | Implemented | Not applicable to fixed degree | — | — |

The main refinement sequences can be read with $n$ grid points per coordinate
axis, $h=1/(n-1)$, and field degree $p$:

- h: P1 commonly uses `n=5→9→17` ($h=1/4\to1/8\to1/16$); P2/P3 commonly
  use `n=3→5→9`. P0 projection and Taylor–Hood Stokes use `n=3→5→9`.
- p: Poisson, complex Helmholtz, variable conductivity, and coupled
  reaction–diffusion check `p=1→2→3→4` on a fixed `n=2` mesh.
- hp: Poisson, Helmholtz, and variable conductivity use
  `(n,p)=(2,1)→(3,2)→(5,3)`.
- isoparametric: curved P2 Poisson and conductivity use `n=5→9→17`
  in 1D/2D and `n=3→4→5` in 3D. A separate P1 geometry-map error
  study uses `n=5→9→17` in 1D/2D and `n=3→5→9` in 3D.

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
their refinement axis.

The CI convergence job runs the full local suite twice: once with sequential
assembly and once with OpenMP assembly (`RODIN_MULTITHREADED=OFF/ON`). A
separate PETSc job checks local-context and distributed P1/P2 Poisson
convergence with PETSc assembly and CG. The distributed suite uses mesh
families partitioned across one to four MPI ranks and globally reduced norms;
it is a separate check from the local-context suites.
