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
$h_i=1/(n_i-1)$ for $n_i$ points per axis. On an affine uniform grid of a
fixed cell family, the geometry-dependent ratio between this spacing and
element diameter cancels between levels. For curved meshes, nominal spacing
remains the declared refinement parameter; its interpretation requires
uniformly regular maps and a shape-regular family with physical mesh diameter
$h_i^{\mathrm{phys}}=\max_{K\in\mathcal{T}_i}\operatorname{diam}(K)$ satisfying

$$
c h_i\le h_i^{\mathrm{phys}}\le C h_i,
\qquad 0<c\le C<\infty,
$$

with constants independent of the level. This comparability transfers
asymptotic orders, but does not make the physical and nominal adjacent-size
ratios identical on finite curved meshes. Each study must identify the size
parameter used in its reported rates. For a sufficiently regular exact
solution of a stable, consistent conforming degree $p$ elliptic problem on
regular meshes, the standard expectations are $E_{1,i}=O(h_i^p)$ and, when
the required dual regularity holds, $E_{0,i}=O(h_i^{p+1})$. These orders are not asserted
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

Shared test utilities keep error integration separate from acceptance policy.
`ErrorNorm` integrates physical scalar/vector errors; `LiftedErrorNorm`
integrates exact-domain field, geometry and total defects for real/complex
scalar/vector fields. `ErrorHistory` and `NormHistory` use the actual spacing
ratio on every adjacent interval. `LiftedConvergence` composes the scalar
decomposition, rate and independent numerical-sensitivity checks without
containing a PDE solve. Direct utility regressions include an invalid final
interval and both violated triangle inequalities. Physics workloads retain
their own manufactured data, solver policy and wrong-operator controls.

When field and geometry degrees differ, a single two-sided rate window for
the total error need not describe the resolved refinement range. For each
norm index $j\in\lbrace 0,1\rbrace$, let $F_j$ and $G_j$ be the
independently measured lifted field and geometry errors. With adjacent
spacing ratio $\rho=h_f/h_c<1$, the mixed-order acceptance rule uses

$$
T_{j,f}\le F_{j,c}\rho^{s_{F,j}-\delta_j}
             +G_{j,c}\rho^{s_{G,j}-\delta_j}+\varepsilon_j.
$$

The exponents $s_{F,j}$ and $s_{G,j}$ are specified by the workload's
field and geometry approximation estimates. The rate allowance $\delta_j$
and absolute floor $\varepsilon_j$ are acceptance policies, not proved
error constants. Independent component-rate windows and triangle-inequality
checks remain required. Strict decrease of $T_j$ is an additional numerical
policy; it does not follow from the triangle inequality when cancellation
changes between levels. This rule avoids presuming that the slower component
already dominates the total error on the tested meshes.

Representable-field geometry studies use a distinct acceptance path.
If $u_\ast\rvert_{\Omega_h}\in V_h$, the exact represented-domain
field error vanishes. Its measured counterpart is tested against a stated
absolute budget $\delta$, not fitted to a convergence rate. For
$j\in\lbrace 0,1\rbrace$, the protocol requires

$$
E_{\mathrm{represented},j}<\delta,\qquad E_{F,j}<\delta,\qquad
\left|E_{T,j}-E_{G,j}\right|\le\delta.
$$

Here $E_{F,j}$, $E_{G,j}$ and $E_{T,j}$ denote the lifted field,
geometry and total norms defined above. `LiftedConvergence` retains
only the geometry and total histories for this path; both must be finite,
positive and decreasing over at least three levels. For geometry degree
$q$, the expected orders are $q+1$ in L2 and $q$ in the H1 seminorm,
subject to smooth regular maps and a nonvanishing leading interpolation
defect. Every observed adjacent rate is checked directly. Quadrature and
solver sensitivity remain separate controls; relative comparisons apply
only to nonzero geometry and total norms. The budget $\delta$ is supplied
by the physical workload in the measured norm's units. These field-error
budgets do not identify mesh entities: chart correspondence and MPI
ownership checks remain based on logical indices.

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

Rate certification also requires a finite computed rate: finite positive
errors alone do not guarantee a finite intermediate quotient, and a lower
bound alone does not exclude an infinite result. Logarithmic reductions must
preserve resolved differences between nearby samples as well as handle the
full positive floating-point range. Consequently, the rate utilities require
independent algebraic regressions for extreme errors and scales, nearby
samples, and genuinely unrepresentable rates. These are arithmetic and
acceptance checks, not additional PDE convergence evidence.

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

## Qualification records

[Qualification records and work plans](qualification.md) are maintained
separately from this description of the convergence suite design.
