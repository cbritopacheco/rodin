# Nonpolynomial geometry approximation

## Reference domain and map errors

The comparison domain is the fixed unit box $\widehat\Omega=(0,1)^d$,
$d\in\lbrace 1,2,3\rbrace$. The prescribed map is

$$
\Phi(\xi)=\xi+a\sin(\pi\xi_0)e_{d-1},\qquad a=0.1.
$$

Its derivative is $D\Phi=I+a\pi\cos(\pi\xi_0)e_{d-1}e_0^T$.
For $d=2,3$, its determinant is one; for $d=1$, it is bounded below
by $1-a\pi>0$. Finite geometry degrees do not represent this sine map
exactly. The represented physical domains are therefore
$\Omega_{h,q}=\Phi_{h,q}(\widehat\Omega)$, not an assumed fixed exact domain.
In one dimension, the endpoint images remain zero and one, so the physical
interval is unchanged: map errors measure the interior chart approximation,
not boundary displacement. In two and three dimensions, the curved boundary
is also approximated. A map error is not itself a PDE domain-error estimate.

For each original cell chart $T_K:\widehat K\to K\subset\widehat\Omega$,
geometry-node samples of $\Phi\circ T_K$ define a transformation of degree $q$
$X_{K,q}$. Geometry and physical field degrees are distinct variables.
On each original cell, the geometric errors are measured by

$$
G_0(h,q)^2=\sum_K\int_K
\Vert X_{K,q}\circ T_K^{-1}-\Phi\Vert^2\thinspace d\xi,
$$

$$
G_1(h,q)^2=\sum_K\int_K
\Vert D_rX_{K,q}(r)[D_rT_K(r)]^{-1}-D_\xi\Phi(\xi)\Vert_F^2\thinspace d\xi,
\qquad r=T_K^{-1}(\xi).
$$

Thus the derivative error is taken with respect to the global reference
coordinates, not cell coordinates. Quadrature weights include the original
chart distortion, not the approximated physical-domain distortion.
The exact sine map and its derivative are evaluated analytically in the oracle,
independently of the map-installation helper.

Each degree $q=1,2,3$ uses three meshes $n=3,5,9$, with $n$ grid points
per coordinate and $h=(n-1)^{-1}$. For smooth maps, shape-regular original
meshes and geometry spaces having the required approximation properties,
the targets are $G_0=O(h^{q+1})$ and $G_1=O(h^q)$.
Both adjacent intervals require finite positive decreasing errors and

$$
q+1-0.35<r_0<q+1+0.35,\qquad
q-0.35<r_1<q+0.35,\qquad
r_i=\frac{\log(G_i(h_{\ell-1},q)/G_i(h_\ell,q))}
{\log(h_{\ell-1}/h_\ell)}.
$$

These are finite-resolution policies, not proofs of asymptotics or global
invertibility. Positive finite determinant ratios are required at all sampled
cell quadrature points. Segment, Triangle, Quadrilateral, Tetrahedron, Pyramid,
Hexahedron and Wedge are tested. No spatial rate is defined for point geometry.

## Geometry degree and combined refinement

Two additional paths distinguish changing geometry degree from fixed-degree
mesh refinement. Degree refinement uses $n=3$ and $q=1,2,3$.
Combined refinement uses

$$
(n,q)=(2,1)\longrightarrow(3,2)\longrightarrow(5,3),
\qquad h=(n-1)^{-1}.
$$

Every path contains three discretizations and checks both adjacent intervals.
The same independent $G_0$ and $G_1$ integrals are evaluated at all levels.
All errors must be finite and positive. Each interval requires

$$
G_j(h_{\ell},q_{\ell})<\rho G_j(h_{\ell-1},q_{\ell-1}),
\qquad j\in\lbrace 0,1\rbrace,\quad\ell\in\lbrace 1,2\rbrace,
$$

with $\rho=1/2$ for degree refinement and $\rho=1/4$ for combined
refinement. These are finite-path improvement policies for the analytic sine
map, not asserted fixed-degree powers of $h$. Three degrees do not establish
an asymptotic exponential law in $q$. Geometry interpolation families,
reference domain, amplitude, map oracle and derivative convention remain
unchanged. The combined path changes both $h$ and $q$, so it does not isolate
either contribution.

At every level, integration orders fourteen and eighteen must agree within
$10^{-6}$ relatively for each map error. An independent physical affine
Poisson patch with $p=q$ must retain both $10^{-9}$ field-error bounds.
Thus improving geometry errors and exactly reproducing a physical field are
distinct assertions. On MPI partitions, the same reference-cell ownership
and globally reduced squared norms are used, including empty ranks on the
coarsest mesh. These paths inherit the existing native and real-PETSc
local/MPI rank and sequential/OpenMP matrix.

## Independent physical-field patch

To distinguish geometry error from field error, a separate Poisson problem
is solved on each represented domain:

$$
-\Delta u=0\quad\text{in }\Omega_{h,q},\qquad
u=1+\sum_jx_j\quad\text{on }\partial\Omega_{h,q}.
$$

Its exact physical solution is affine, with gradient $(1,\ldots,1)^T$.
Using field degree $p=q$, its pullback is represented by the same degree $q$
space that represents the coordinate map. Physical norms must satisfy

$$
E_0=\Vert u_h-u\Vert_{L^2(\Omega_{h,q})}<10^{-9},\qquad
E_1=\Vert\nabla u_h-\nabla u\Vert_{L^2(\Omega_{h,q})}<10^{-9},
$$

at $n=3$, while both geometry errors remain positive. This verifies field
evaluation, traces and assembly at degrees one through three without attributing
geometry approximation to the field solver. It is an exact patch test, not a
nonpolynomial PDE field-rate study or a domain-error estimate for Poisson on
the exact domain $\Phi(\widehat\Omega)$.

Two controls test independence. Omitting the sine deformation gives
$G_0=a/\sqrt2$ and $G_1=a\pi/\sqrt2$; the tests require
$G_0>0.05$ and $G_1>0.15$. The affine physical patch must still pass
on this wrong domain map: a field patch alone cannot validate the domain.
Conversely, replacing the correct trace by zero must give both physical
field errors above $0.1$. These absolute thresholds are dimensionless policies
on the prescribed unit-scale domains.

## Numerical realization and backend scope

Map integration order fourteen is checked against eighteen; each positive geometry
error must change relatively by less than $10^{-6}$ at $n=3$.
Patch assembly order twelve is separately increased to sixteen, and CG relative
tolerance $10^{-13}$ is separately reduced to $10^{-14}$; both field errors
must retain their absolute patch bounds. Physical error integration uses order
fourteen. CG permits at most $50000$ iterations. PETSc uses absolute tolerance
$10^{-14}$ and divergence tolerance $10^5$. Native solver success or a
positive PETSc convergence reason is required, together with

$$
\frac{\Vert A_hU_h-b_h\Vert_2}{\max(1,\Vert b_h\Vert_2)}<10^{-11},
$$

computed independently of the solver's convergence report.

Workload preserves the original charts by copying the mesh before installing
geometry on the second mesh. Logical cell indices and ordered vertex indices
are checked exactly; no coordinate-tolerance matching is used. Geometry::Point
supplies both chart Jacobians and metric factors. CurvedGeometry supplies
the existing node-installation machinery with an opt-in sine map; its quadratic
default remains unchanged. ErrorNorm and ErrorHistory supply the shared
scalar/vector magnitudes, physical norms and adjacent rates.

Native local and real-PETSc local/MPI configurations are registered, with
sequential/OpenMP assembly and one to four MPI ranks. The continuous dimension
comes from the declared grid family, including on empty shards. Geometry errors
sum owned original cells only, and physical errors sum owned represented cells
only; squared contributions are reduced globally before taking square roots.
Complex-PETSc builds do not register this real scalar suite. Registrations are
labelled slow, use 1800-second timeouts and share a pyramid resource lock.

Smooth represented-domain and lifted exact-domain PDE field studies are
specified separately in [diffusion](../Diffusion/README.md) and
[complex Helmholtz](../Helmholtz/README.md). Their field, geometry and total
errors are not interchangeable with the reference-domain map errors here.
