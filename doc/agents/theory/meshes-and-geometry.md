# Meshes and discrete geometry — theory ↔ code

The mathematics behind `src/Rodin/Geometry/`. Companion to
[../geometry.md](../geometry.md) (the API inventory).

## The mesh as a graded incidence complex

Abstractly, a mesh $\mathcal T_h$ is a cell complex: a family of polytopes
graded by dimension, closed under taking faces, glued along shared
sub-polytopes. What the data structure stores is the *combinatorics*: for
each pair $(d, d')$, the incidence relation "which $d'$-polytopes touch
each $d$-polytope". This is the FEniCS-style mesh model: all topology is
derivable from the cell→vertex relation by three primitive computations
(build, transpose, intersect), which is why
`Connectivity::compute(d, d')` exists as an explicit, on-demand operation
— each relation costs memory and time, so nothing is computed
speculatively, and code declares what it needs
(`RODIN_GEOMETRY_REQUIRE_*`).

Identity of a polytope is its vertex *set*, not an ordering —
`Polytope::Key` with symmetric hash/equality is the quotient by
permutations. Orientation (a sign on the ordering) matters separately:
facet normals (`FaceNormal`, `BoundaryNormal`) and jump conventions
($w^+ - w^-$) depend on a consistent orientation choice, which is fixed by
the stored vertex order and the incidence direction.

Attributes are labels $\mathcal T_h \to \mathbb N$ per polytope — the
discrete encoding of subdomains $\Omega_i$ and boundary parts $\Gamma_j$.
All region selection in the form language (`.over(attr)`, `.on(attr)`,
`traceOf`) is preimage-of-a-label; `trace({pair → attr})` labels the
interface $\partial\Omega_i \cap \partial\Omega_j$; `skin()` is the
boundary operator $\partial$ realized as a (d−1)-complex with parent maps.

## Geometry = charts attached to combinatorics

Each polytope carries a chart $x_\tau : K \to \mathbb R^{d_s}$
(`PolytopeTransformation`). The complex plus charts is a piecewise-smooth
manifold (possibly embedded: surface meshes have $d < d_s$, and
`isSurface()` is meaningful). All metric quantities derive from the chart:

- Jacobian $J(\hat x)$ → `Geometry::Point::getJacobian()`;
- volume measure $|\det J|\,d\hat x$; surface measure
  $\sqrt{\det(J^\top J)}\,d\hat x$ (Gram determinant — this is what makes
  integrals on embedded facets/surfaces correct);
- validity = orientation preservation: $\det J > 0$ everywhere. For affine
  simplices this is one number per cell; for curved (parametric) cells it
  must be *sampled* over a quadrature (the house check uses order 4) —
  a curved cell can be valid at its corners and folded in the interior.

Separating combinatorics from charts is what lets the same mesh be
"linear" or "curved P2" by swapping transformations
(`ParametricTransformation<FE>` over a `PointCloud` of control nodes,
`flush()` to demote to affine) with no topological change.

## Element quality and why it matters

Interpolation error constants degrade with element distortion: for the
standard $H^1$ estimate the constant scales with the aspect
ratio/condition of $J$ (shape-regularity). Practical quality functionals
for triangles: the normalized measure $q = 4\sqrt3\,|K| / \sum \ell_i^2$
(1 for equilateral, →0 for slivers) and Jacobian-based measures
$\mu(J W^{-1})$ against a target Jacobian $W$ (the target-matrix /
TMOP literature). Two regimes to keep apart conceptually:

- **Discrete (combinatorial) improvement** — split/collapse/swap change
  the complex. Correctness constraints are topological: the *link
  condition* guards collapses (merging two vertices must not pinch the
  manifold), and every rewritten cell is re-validated (orientation, area
  floor) before commit.
- **Continuous (geometric) improvement** — vertex relocation/smoothing and
  variational mesh optimization move charts, complex fixed. Correctness
  constraint is $\det J > 0$ (a barrier/line-search matter, not a
  topological one).

The optimizer's evaluate-then-commit structure (propose over all edges,
commit a vertex-disjoint independent set) is the standard way to make such
passes order-independent and parallelizable.

## Level sets and implicit domains

A domain can be represented implicitly: $\Omega = \{\phi < 0\}$ with
interface $\Gamma = \{\phi = 0\}$. The useful normalization is the
*signed distance function*, characterized by the eikonal equation
$|\nabla\phi| = 1$ (with $\phi = 0$ on $\Gamma$); then
$n = \nabla\phi$, curvature $= \Delta\phi$ on $\Gamma$. Code map:

- `Eikonal::FMM` — fast marching solves the eikonal equation in
  $O(N\log N)$ by upwind causality (Dijkstra-like); the standard
  redistancer.
- `Distance::{Eikonal, Poisson, SignedPoisson}` — PDE-based
  distance/approximate-distance solvers; `Rvachev`, `SpaldingTucker` —
  algebraic normalizations of a non-distance level set.
- Transport of $\phi$ under a velocity is the Hamilton–Jacobi/advection
  step (`Advection::Lagrangian`). Discrete transport destroys the distance
  property (gradients steepen/flatten), which is *why periodic
  redistancing is mandatory* — a modeling fact, not an implementation
  preference.
- Converting implicit → explicit geometry is discretization of the zero
  set: MMG's `LevelSetDiscretizer` produces a body-fitted complex whose
  interface facets carry an attribute (then everything above about
  attributes applies).

## Partitioning

Distributing a mesh = partitioning the cell-adjacency graph
(cells as nodes, shared facets as edges) to minimize cut (communication)
under balance (load). `MeshPartitioner` implementations: `Greedy`,
`BalancedCompact`, and Scotch's graph partitioner; `Sharder`/`Shard` carry
the resulting ownership + overlap (ghost) metadata that the MPI stack
consumes. Ghost layers exist because assembly of a facet/cell integral
needs both sides' DOFs — the overlap width is a discretization statement
(stencil support), not a networking detail.
