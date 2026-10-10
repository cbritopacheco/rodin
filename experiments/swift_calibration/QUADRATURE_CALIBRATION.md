# WNGIR quadrature calibration

This historical study predates the rename to SWIFT (Strain-distributed Welsch
Implicit-interface Fitting Technique). Recorded executable names, commands and
artifact paths below retain their original names for reproducibility. Current
sources and runners are in `experiments/swift_calibration` and use SWIFT names.

The screening study estimates a sufficient integration order for the canonical
robust WNGIR model. It does not establish polynomial exactness, certified
maximum geometric error, or a universal minimal order.

## Grid

| Setting | Values |
| --- | --- |
| Integration order | 2, 4, 6, 8, 12, 16, 24, 32 |
| FE-target cases | 2D, P1 displacement, n = 8 and 16, target P1/P2/P3 |
| Analytic-target controls | P1 2D n = 32; P1 3D n = 8; P2 isoparametric 2D n = 8 |
| Lobes | 4 |
| Radius / amplitude | 0.24 / 0.08 |
| Fitting / deviatoric / divergence weights | 1 / 0.0001 / 0.01 |
| Hinge weight | 10 |
| Outer / inner caps | 30 / 15 |
| Independent geometric validation order | 32, with the solver's vertex checks |
| Independent final quality sampling | 32 for 2D/P2; P1 affine 3D has cellwise constant deformation gradient |
| Backend / threads | MUMPS / 1 |
| Linear relative tolerance | 0.000001 |
| Attempts / timeout | 72 / 180 seconds per attempt |
| Warmups / repeats | None |

Order 32 and 24 references run first for every case. One process runs at a time;
other agents' workloads make wall-clock rankings provisional. The explicit
integration override changes both surface and volume sampling, so the study
calibrates that public setting, not the surface order in isolation.

## Finite-element targets

The experimental copy of `LevelSetWNGIRSweep.cpp` interpolates the analytic
target into H1 P1/P2/P3 and passes a scalar evaluation adapter and its native
`Grad` to WNGIR. The adapter is required by the current scalar-function API;
it evaluates the actual GridFunction, not the analytic target.
The gradient explicitly selects the classified interior trace for evaluations
on interface facets; target-cell gradients at relocated volume points retain
their native values. An unspecified facet trace is undefined in Rodin.
`WNGIR_TARGET_DEGREE` selects the degree only in the frozen experimental
binary. Production has no new target-selection flag. Classification remains
analytic and identical between integration orders. These runs measure fit to
the FE target, not exact distance to the original analytic target. In
particular, the target gradient can jump across cells; the composed integrand
need not be a single polynomial on a moved facet.

The source, compile/link commands and binary hashes are retained under
`tmp/wngir-native-quadrature-20261008`. The runner copies the supplied binaries
into its output directory before executing them. It refuses an existing
manifest rather than overwriting a previous campaign.

## Responses and selection

Retain every attempt, including failures and timeouts. Log the final solver
energy, sampled maximum geometric residual and normalized constant, target
hit, quality compliance, outer count, total/maximum inner corrections,
stationarity diagnostics, Jacobian, distortion, reason and total time. Raw
trace logs retain inner/outer evolution and assembly/solve timings. These are
the solver's final responses before the examples select their best-RMS output
mesh; exported example meshes must not be mistaken for those final states.

A case has a stable reference only if orders 24 and 32 return successfully,
respect sampled quality, agree on target hit, and their normalized geometric
constants differ by at most 5% of the order-32 value plus 0.001. Without that
agreement, no minimum is reported. Stability here is an empirical trajectory
test, not an integration-error certificate.

An observed acceptable order respects sampled quality, does not lose a
reference target hit, has a geometric constant no greater than 105% of the
reference plus 0.001, and takes at most two more outer iterations. Report all
acceptable orders, not just their minimum: nonlinear trajectories need not
vary monotonically with integration order. An improved geometric outcome is
allowed by this practical screening criterion; it does not imply agreement
of the assembled operators.

Energies are computed using different rules and possibly different sampled
normalizations. They are diagnostics, not common-state energy errors. A
production default reduction requires wider geometry/topology coverage and
fixed-state integral/force checks in addition to this screening study.

## Production order contract

`FittingTensor` reports zero order for a globally constant vector gradient.
The complete fitting integrand then propagates its order through Rodin's
expression algebra on affine entities. Explicit `quadrature.order` takes
precedence. Curved measures, arbitrary analytic targets and located FE
compositions remain unknown and use the established nonpolynomial policy.
Robust energy and force retain a common integration rule. Geometric validation
remains independently configured. General composed FE-order propagation needs
a certified single-target-cell restriction or integration-cell subdivision;
merely forwarding a target's element degree would be incorrect.

## Pilot results, 8 October 2026

The analytic controls are retained in
`tmp/wngir-quadrature-calibration-20261008`; the corrected FE targets are in
`tmp/wngir-quadrature-fe-calibration-20261008`. The first 48 FE setup attempts
were rejected for an unspecified gradient trace before returning solver
responses. They remain in the first directory and are not numerical results.
Only that subset was retried with an explicit interior trace. Thus the 72
logical configurations have 120 recorded attempts: 48 setup failures, 71
completed numerical attempts, and one timeout. No numerical warmups or repeats
were performed.

All corrected FE attempts returned successfully and passed independent
order-32 quality sampling. There were 35 target hits, 10 step-stagnation exits
and three outer-budget exits. The maximum inner count was three; the inner
cap was not reached.

| Target / displacement | Resolution | Smallest order passing the pilot criterion | Qualification |
| --- | --- | --- | --- |
| FE P1 / P1 | 16 | 8 | Target reached; stable 24/32 reference |
| FE P2 / P1 | 16 | 4 | Target reached; acceptable orders are not monotone |
| FE P3 / P1 | 16 | 4 | Reference itself misses the target; this reproduces its fit, not a successful target hit |
| FE P1/P2/P3 / P1 | 8 | Undetermined | 24/32 final geometric constants differ beyond the screening tolerance |
| Analytic 2D / P1 | 32 | 4 | Target reached in 11 outer iterations, as at orders 24 and 32 |
| Analytic 3D / P1 | 8 | Undetermined | Order 32 exceeded the 180-second timeout |
| Analytic 2D / curved P2 | 8 | Undetermined | Order 24 violates independently sampled quality |

For the FE targets, hit counts across the six geometries were:

| Integration order | Target hits |
| --- | --- |
| 2 | 2 / 6 |
| 4 | 4 / 6 |
| 6 | 5 / 6 |
| 8 | 5 / 6 |
| 12 | 4 / 6 |
| 16 | 5 / 6 |
| 24 | 5 / 6 |
| 32 | 5 / 6 |

These counts do not establish a universal minimal rule. The coarse references
stop at different fitting states, sometimes at different accepted iterations;
their disagreement is not a fixed-state integral error. Located gradients can
jump across target cells, and quadrature also changes sampled normalizations.
Consequently, increased order need not monotonically improve a trajectory.

The curved P2 quality discrepancy is measured on the same final displacement,
before the example selects its output mesh. At integration order 24,
independent order-32 sampling found maximum distortion 10.471 despite the
solver reporting less than 10. At integration order 8, the independent value
was 23.118. Order 32 passed its own sampled check, but this is not an
everywhere quality certificate. Increasing fitting quadrature alone does not
settle curved-mesh admissibility.

The 3D order-32 attempt accepted its first outer step after approximately
113 seconds, then timed out. The explicit global override also raises volume
quadrature; this is unnecessarily expensive for cellwise-constant P1
deformation gradients. Orders 6 and 24 completed in two outer iterations with
geometric constants approximately 0.659 and 0.662, respectively. That is a
promising lower-order comparison, not a passing 24/32 reference test.

The nonpolynomial production fallback remains unchanged. Orders 4, 6 and 8
are candidates for subsequent fixed-state and wider-geometry checks, not new
universal defaults. General FE-target order inference should resolve crossed
target-cell pieces before asserting polynomial order. Integration sampling
and admissibility validation should be distinguished before calibrating a
curved-mesh default. Concurrent external CPU workloads preclude reliable
wall-clock rankings in this pilot.

## Separated calibration, 8 October 2026

The next stage separates surface integration, volume integration, geometric
response sampling and actual-quality sampling. Model coefficients, robust
loss, directional Newton, Armijo and iteration budgets are unchanged.
The public common-order override remains supported; the more specific
surface/volume overrides take precedence. All new controls default to zero,
preserving the existing integration and acceptance policy. A positive
quality-validation order adds cell vertices to its Gaussian volume samples.

| Setting | Values |
| --- | --- |
| Frozen displacement degree | P1, P2, P3 |
| Dimension | 2, 3 |
| Frozen mesh | Affine simplices on a 3-point-per-axis unit grid |
| Frozen displacement states | Mild deformation; positive-Jacobian compression/shear |
| Target | Smooth analytic function; native H1 target of the displacement degree |
| Robust scale / gradient normalization | 0.2 / 1, frozen across rules |
| Surface orders | 2, 4, 6, 8, 12, 16, 24, 32, 64 |
| Volume orders | 2, 4, 6, 8, 12, 16, 24, 32 |
| Surface reference | 64, checked against 32 |
| Volume reference | 32, checked against 24 |
| Quality reference | Gaussian samples and vertices, all-cell simplex lattice of density 8, then critical-cell densities 16, 32, 64 |
| Geometric reference | Gaussian samples and vertices, all-interface lattices of density 8, 16, 32, 64 |
| Process / threads / repeats | Serial / 1 / none |

The volume comparison uses three native FE displacement probes and their
3-by-3 centered-distribution and affine-hinge tangent action matrices, not
the full assembled operators. The affine hinge is tested at a fixed
compressive increment with coefficient one; its activation is held to the
same model, not to the same sampled active set. Full surface matrices and
force vectors are assembled through Rodin. Surface integration passes if
relative energy error is at most 0.001 and force/matrix errors at most 0.01.
The 32/64 reference check is ten times stricter. Volume reference actions
must agree to 0.001; candidate actions to 0.01.

The 12 frozen states completed successfully, producing 24 surface studies.
Sources, hashes, commands, logs and `fixed-errors.csv` are preserved under
`tmp/wngir-fixed-order-calibration-20261008`. The initial driver was stopped
in `tmp/wngir-order-calibration-20261008`: copying the base reference returned
by `Grad(target).traceOf(...)` caused a `copy()` loop. Profiling identified
the problem; retaining the concrete Grad object and calling `traceOf`
separately fixed it. No completed numerical result from that setup is pooled.
Affine P1 volume rows all use evaluated order 2, explicitly logged; higher
order labels there do not constitute independent high-order checks. Exactness
is supported by the cellwise-constant integrand, not by those duplicate rows.

| Surface order | Worst relative energy error | Worst force error | Worst fitting-matrix error |
| --- | ---: | ---: | ---: |
| 4 | 0.00286 | 0.0674 | 0.373 |
| 6 | 0.000115 | 0.00107 | 0.00785 |
| 8 | 0.0000177 | 0.00108 | 0.00213 |
| 12 | 0.0000177 | 0.000241 | 0.000525 |
| 24 | 0.00000422 | 0.000155 | 0.000327 |
| 32 | 0.00000403 | 0.000102 | 0.000198 |

All surface references passed the stability criterion. Order 6 is the
smallest rule passing all these fixed-state surface tests. This is not a
universal order: the grid, target frequency, robust scale and deformations
are specified above. Order 4's worst error occurs for the stressed 3D P3
FE-target matrix. Located target gradients and tensor products need more
integration accuracy than the scalar energy alone suggests.

The smooth distribution actions converge rapidly: the worst relative error
at volume orders 4, 6 and 8 is respectively 0.00188, 0.0000810 and 0.00000414.
This does not settle active-hinge integration. In the stressed 2D P2 state,
the hinge tangent differs by 5.80% between orders 24 and 32, whereas its
energy differs by only 0.0446%. Its order-12 tangent error is 10.57%.
The stressed 2D P3 hinge reference also fails the stricter stability check.
The squared hinge is continuously differentiable but its tangent jumps at
activation. Different nonnested Gaussian rules sample different portions
of that region; increasing order does not monotonically reduce the tangent
error. A polynomial-degree rule cannot establish exactness for this model.

The boundary-inclusive reference found quality extrema more accurately than
Gaussian sampling alone. Across the manufactured P3 states, order 8 plus
vertices still differs by up to 1.81% in minimum Jacobian and 1.25% in maximum
distortion. Order 16 reduces these discrepancies below 0.35%. These are
empirical errors on this study, not safety margins valid for all curved
cells. Critical-cell refinement can miss extrema in cells not flagged by
the initial samples. Geometry order 8 plus vertices underestimates the
dense sampled maximum by up to 1.44%; order 32 reduces this to 0.213%.
No finite rule in this study certifies an everywhere supremum.

## Full-solve confirmation

Five geometries were tested with four settings: the old common order 12;
separate surface orders 8, 12 and 24 with volume order 2 for P1 or 8 for P2,
quality order 16 plus vertices, and geometric validation order 32. Eight
analytic attempts were rejected at setup because the runner set the
FE-target environment variable to zero instead of leaving it unset. They
are retained, and only those configurations were retried with unchanged
binary hashes. Thus 28 attempts cover 20 logical configurations: eight
setup failures and 20 completed numerical solves. Artifacts are in
`tmp/wngir-separated-quadrature-20261008` and
`tmp/wngir-separated-quadrature-analytic-20261008`.

| Displacement / target | Resolution | Common 12: final C / outer | Separate surface 8: C / outer | Separate surface 12: C / outer | Separate surface 24: C / outer |
| --- | ---: | ---: | ---: | ---: | ---: |
| P1 2D / analytic | 32 | 0.999816 / 11 | 0.999796 / 11 | 0.999816 / 11 | 0.999735 / 11 |
| P1 2D / FE P1 | 16 | 0.394101 / 2 | 0.413198 / 2 | 0.394101 / 2 | 0.405466 / 2 |
| P1 2D / FE P3 | 16 | 1.126514 / 18 | 1.127338 / 19 | 1.126514 / 18 | 1.128693 / 18 |
| P1 3D / analytic | 8 | 0.661598 / 2 | 0.666439 / 2 | 0.661598 / 2 | 0.661595 / 2 |
| Curved P2 2D / analytic | 8 | 2.740908 / 11 | 3.686147 / 13 | 2.772422 / 11 | 2.767848 / 12 |

Here C is the solver's sampled normalized geometric constant using background
grid spacing and exponent p+1: h squared for P1, h cubed for P2. The P1
analytic and FE-P1 cases hit their sampled target with every rule; the FE-P3
and curved-P2 cases do not. The latter stop by accepted-step stagnation.
Every P1 solve has zero inner corrections because its affine hinges remain
inactive. The curved-P2 common-order control uses eight total corrections,
at most one per outer iteration; the three independently validated variants
use zero. There are no inner-cap or outer-cap exits in these confirmations.

| Curved P2 setting | Energy | Reported maximum Q | Independent Gaussian order-32 maximum Q |
| --- | ---: | ---: | ---: |
| Common order 12 | 0.000000254 | 9.9999998 | 10.720763 |
| Surface 8 / volume 8 / quality 16 + vertices | 0.00000159 | 9.9999957 | 9.786731 |
| Surface 12 / volume 8 / quality 16 + vertices | 0.000000307 | 9.9999605 | 9.789304 |
| Surface 24 / volume 8 / quality 16 + vertices | 0.000000290 | 9.9999547 | 9.789729 |

The last column samples interior points only, and is not a certificate.
It catches the old control's violation but need not reproduce the vertex
maximum included by the new acceptance test. The independent-validation
settings avoid that detected violation; none justifies claiming target
attainment or everywhere quality compliance. The surface-8 P2 result is
materially worse than surface 12, so the P1 recommendation is not extended
to curved P2.

Eight additional P1 solves lower only independent quality order from 16 to
2, retaining vertices. Their energies, C, quality and iteration counts
reproduce the corresponding separate-rule results to printed precision.
These are in `tmp/wngir-p1-validation-quadrature-20261008`.
With surface order 12 unchanged, replacing common volume order 12 by 2
reproduces all four P1 controls to roundoff, including the stagnating case.
This establishes behavior preservation for those configurations. Reducing
surface order can change the trajectory, especially for located FE targets.

## Recommendation

| Use | Surface order | Volume order | Quality sampling | Geometric sampling | Status |
| --- | ---: | ---: | --- | --- | --- |
| Tested affine P1, analytic or located FE targets | 8 | 2 | 2 + vertices | 32 + facet vertices | Practical calibrated configuration; not universal |
| Preserve tested P1 common-order-12 results | 12 | 2 | 2 + vertices | 32 + facet vertices | Reproduces tested results to roundoff |
| Curved P2 | 12 | 8 | 16 + vertices, independently checked | 32 + facet vertices | Conservative screening candidate, not a calibrated universal rule |
| P3 displacement or active curved hinges | Retain nonpolynomial fallback | No universal minimum established | Boundary-inclusive validation required | Independently refine | Full-solve calibration remains open |

The generic production fallback is deliberately unchanged. Known globally
constant fitting tensors retain native order inference on affine entities;
unknown composed targets retain the explicit nonpolynomial policy. The new
settings permit inexpensive P1 volume/quality rules without reducing surface
or geometric-response sampling. Further curved-mesh work should refine
active-hinge integration and boundary-inclusive quality validation rather
than raising one common order indiscriminately. Other agents ran concurrent
MPI jobs, so wall-clock ratios are not controlled performance evidence.
The regression suite includes the override precedence and a P2 vertex
inversion missed by interior samples. Direct Geometry::Point evaluation is
used for vertices; an IntegrationPoint with a null quadrature formula must
not enter the native cached-Jacobian evaluation path.

Verification: 64 rebuilt tests pass (48 solver, 7 example/parser, 3 hinge,
2 distribution and 4 directional-Newton tests), Python runners compile,
and `git diff --check` passes. No broad parameter campaign is running.
No paper changes, numerical-default changes, commits or pushes were made.

## Adopted default policy

Following the calibration, the user requested adoption of the recommended
defaults. Automatic affine-P1 surface/volume/quality/geometric orders are
8/2/2/32; curved or P2/P3 orders are 12/8/16/32. Actual geometry order is
included in the selection. Cell and facet vertices are checked by default.
The low-order volume/quality branch is restricted to simplices: tensor-product
degree-one basis gradients are not generally cellwise constant.
Specific integration overrides remain available, while quality and geometric
sampling are independent of the common integration override. Native fitting
order inference is retained where its polynomial contract is valid.

These are practical screening policies, not uniform accuracy guarantees.
An existing coarse analytic-facet regression has relative energy error about
0.00280 and force/energy-vector error 0.00725 at order 8. The strict historical
0.0001 checks remain with an explicit order-12 override, and the automatic
screening rule is tested separately against a 1% allowance. The calibrated
default does not promise the stricter tolerance for arbitrary target data.
Higher-order defaults remain provisional, as established above.

Adoption verification: the 48 solver and 7 example/parser tests pass after
rebuilding with the new policy, including automatic vertex checks, curved
geometry and nonsimplicial policy selection, and independent overrides.
