# Complex Helmholtz on exact quadratic geometry

This suite isolates field approximation on curved cells from approximation
of the physical domain. Let $Q=(0,1)^d$, $d\in\{1,2,3\}$, and define

$$
\Phi(\xi)=\xi+a\xi_0^2e_{d-1},\qquad a=0.1,\qquad \Omega=\Phi(Q).
$$

The shared `CurvedGeometry` retains original grid vertices, maps geometry
control points once, and installs degree-two `ParametricTransformation`s
on every positive-dimensional entity, including boundary and halo entities.
The map belongs exactly to each cell's degree-two geometry family. Topology,
logical indices, and MPI ownership are unchanged; identical maps are evaluated
locally on shared entities after partitioning the unwarped grid.
Mesh entity iterators enumerate the local entities and halos used by these
setters. Global MPI entity counts are not local-index bounds.
The workload dimension is that of the prescribed grid (equivalently its
ambient dimension here), including on empty shards. Physical references
and the dimension passed to `getMeasure` are not inferred from the highest
dimension of an empty local incidence complex.

For $d\ge2$, $D\Phi=I+2a\xi_0e_{d-1}\otimes e_0$ and
$\det D\Phi=1$. In one dimension, $\Phi'=1+2a\xi_0\ge1$.
Thus the map is regular, and $|\Omega|=1$ for $d\ge2$, whereas
$\Omega=(0,1.1)$ and $|\Omega|=1.1$ for $d=1$. Geometry degree two
is held fixed at all refinement levels: there is no changing-domain error.
With field degree one this is a superparametric study; with field degree
two it is isoparametric in the strict sense.

## Continuous problem and manufactured data

The nondimensional complex Dirichlet problem is

$$
-\Delta u-\kappa^2u=f\quad\text{in }\Omega,\qquad
u=g\quad\text{on }\partial\Omega,\qquad \kappa^2=\frac14.
$$

The smooth physical field, source, gradient, and trace are

$$
s(x)=\sum_{j=0}^{d-1}x_j,\qquad u_*(x)=e^{is(x)},\qquad
f(x)=\left(d-\frac14\right)u_*(x),\qquad
\nabla_xu_*=iu_*(1,\ldots,1)^T,\qquad g=u_*|_{\partial\Omega}.
$$

`HelmholtzData` evaluates these quantities in physical coordinates, not
reference coordinates. In trial-first, test-second convention, the form is

$$
a(u,v)=\int_\Omega\nabla u\cdot\overline{\nabla v}\,dx
-\frac14\int_\Omega u\overline v\,dx,
\qquad \ell(v)=\int_\Omega f\overline v\,dx.
$$

The negative mass term does not make this selected workload indefinite.
Indeed, $\Omega\subset B=(0,1)^{d-1}\times(0,1.1)$ and zero extension
of $H_0^1(\Omega)$ functions into $B$ gives

$$
\lambda_1(\Omega)\ge\lambda_1(B)
=\pi^2\left(d-1+\frac1{1.1^2}\right)>\frac14.
$$

Consequently, $a(v,v)\ge(1-\kappa^2/\lambda_1(B))
\lVert\nabla v\rVert_{L^2(\Omega)}^2$ on homogeneous traces.
The Hermitian positive-definite constrained system is solved by CG;
this argument does not cover arbitrary wave numbers or resonance.

## Spaces, mesh levels, and field errors

The native degree-one path uses complex `P1`; degree two uses complex
`H1<2>`. PETSc local/MPI paths use complex `H1<1>` and `H1<2>`.
Each denotes its cell-appropriate scalar conforming family, rather than
total-degree simplex polynomials on every geometry. Fresh spaces and
linear systems are constructed per solve; assembled matrices are not resized.
`DirichletBC` applies the field-space DOF functionals to the analytic trace;
the discrete trace is $g_h=I_h^{\partial\Omega}g$. A nonpolynomial boundary
function is not claimed to belong exactly to the finite-element trace space.

With $n$ grid points per coordinate axis and nominal spacing $h=1/(n-1)$,
all seven positive-dimensional geometries use the same three-level paths:

| Field degree | $n$ | $h$ | Nominal $L^2/H^1$ orders |
| --- | --- | --- | --- |
| 1 | $5\to9\to17$ | $1/4\to1/8\to1/16$ | $2/1$ |
| 2 | $3\to5\to9$ | $1/2\to1/4\to1/8$ | $3/2$ |

The fixed regular map makes physical element diameters uniformly comparable
to this spacing. Independent physical-cell integration measures

$$
E_{0,h}=\lVert u_*-u_h\rVert_{L^2(\Omega;\mathbb C)},\qquad
E_{1,h}=\lVert\nabla_xu_*-\nabla_xu_h\rVert_{L^2(\Omega;\mathbb C^d)}.
$$

Both errors must be finite, positive, and strictly decrease on each of the
two adjacent intervals. The observed rates
$r_{j,i}=\log(E_{j,i-1}/E_{j,i})/\log(h_{i-1}/h_i)$ must satisfy

$$
\begin{array}{c|cc}
\text{degree}&r_{0,i}&r_{1,i}\\\hline
1&(1.65,2.35)&(0.75,1.25)\\
2&(2.45,3.55)&(1.55,2.45)
\end{array}
\qquad i\in\{1,2\}.
$$

Mapped-element approximation and coercivity provide the expected
$E_{1,h}=O(h^k)$; the usual $E_{0,h}=O(h^{k+1})$ estimate additionally
requires the relevant dual regularity. These finite smooth-data checks do
not prove domain regularity for arbitrary curved maps or certify all
Helmholtz regimes.

Assembly uses quadrature order 11, and field norms use order 13. The pyramid
dispatcher selects a positive conical product rule above order 10; order 11
uses $6^3=216$ points per cell, compared with $7^3=343$ at order 12.
Pyramid basis functions are rational, so polynomial exactness alone does
not establish mapped-integrand accuracy. The separate order-16 sensitivity
control below is required; an order-eight baseline does not satisfy its
$10^{-6}$ budget and is not used. Native CG
uses relative tolerance $10^{-13}$ and at most 20000 iterations. PETSc CG
uses relative tolerance $10^{-13}$, absolute tolerance $10^{-14}$,
divergence threshold $10^5$, and the same iteration limit. Solver success
or a positive PETSc convergence reason is required, together with an
independently recomputed finite coefficient residual
$\lVert Ax-b\rVert_2/\max(1,\lVert b\rVert_2)<10^{-11}$.

## Independent controls

- At $n=3$, degree one reproduces $u_*=1+2i$, with zero gradient and
  $f=-u_*/4$. Degree two reproduces the physical affine patch
  $u_*=1+2i+(1+i/2)s(x)$, whose quadratic pullback is exactly representable
  on the degree-two geometry. Both field errors must be below $10^{-9}$.
  A physical quadratic field is not claimed to lie in this curved P2 space.
- At $n=5$, the affine P2 source and trace are retained while the mass term
  is omitted. The physical-field oracle must reject the changed equation:
  $E_0>10^{-3}$ and $E_1>10^{-2}$, even though its algebraic solve succeeds.
- For each field degree at $n=5$, assembly/norm quadrature is raised
  independently from $11/13$ to $16/18$, and solver relative tolerance is
  separately tightened to $10^{-14}$ at the baseline quadrature. Each change
  must alter both measured errors by less than $10^{-6}$ relatively.
- At $n=2$, every positive-dimensional entity is checked against the
  independent analytic map at order-four reference quadrature points.
  Map discrepancies must be below $10^{-11}$; metric factors must be finite
  and positive. The known physical volume is checked within $10^{-12}$.
  This is a map/metric diagnostic, not an index-ownership comparison.
- The additional MPI P2 norm diagnostic sets $u_h=0$ and reference
  $u=1+i$, so $E_0=\sqrt{2|\Omega|}$ and $E_1=0$ on $n=2$ meshes.
  The absolute L2 budget is $10^{-12}$. Owned-cell squared contributions
  are globally reduced before the square root; halo copies do not add volume.
  Several geometry/rank pairs include ranks without owned cells.

These are nondimensional numerical budgets. Polynomial reproduction,
geometry checks, and sensitivity controls are distinct from the three-level
nonpolynomial field-rate studies.

## Execution and exclusions

Segment, triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and
wedge cases are registered for native local assembly and, with complex PETSc,
local and MPI ranks 1–4. Sequential and OpenMP configurations are distinct
verification paths. Real-scalar PETSc is explicitly excluded by the shared
CMake scalar-capability check. Point geometry has no positive-dimensional
Helmholtz gradient-refinement problem and is excluded.

Native and PETSc drivers retain their solver/storage policies; the
`HelmholtzTest` fixture shares geometry, rate, patch, negative-control, and
sensitivity assertions across them. CI builds the complex PETSc target
explicitly. Each geometry entry is labelled `convergence;slow`, with
additional PETSc/distributed labels and MPI processor counts where relevant,
and a 600-second limit, extended to 1800 seconds for curved pyramids. Run
`ctest --test-dir build/tests -R 'RodinConvergenceIsoparametric.*Helmholtz' --output-on-failure`.
Curved-pyramid entries share the CTest resource lock
`curved_helmholtz_pyramid`, preventing simultaneous memory-heavy mapped
quadrature workloads within one CTest run. This scheduling policy does not
change mesh levels, quadrature, solver settings, or numerical acceptance.
The larger pyramid limit accommodates repeated parametric-Jacobian
evaluation on the finest grid; it is not a relaxed field-error budget.

Mixed boundary conditions, arbitrary field/geometry degrees, nonpolynomial
domain approximation, complex coefficients, and resonant or high-frequency
Helmholtz workloads remain outside this suite's claim.
