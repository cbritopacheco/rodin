# Geometry — the mesh model

The mesh is a **graded incidence complex with attached geometry**. Topology
(which polytopes exist, how they touch) and geometry (where they are, how
they curve) are separate layers throughout.

## Polytopes

- A polytope is identified by `(dimension d, Index i)`; the classes `Cell`
  (top dimension), `Face` (codimension 1), `Vertex` (0D) specialize
  `Polytope`. Each carries an `Attribute` label (subdomain/boundary tags).
- `Polytope::Key` = the vertex index multiset, stored in a fixed stack
  array of `RODIN_MAXIMUM_POLYTOPE_VERTICES` (8); `SymmetricHash` /
  `SymmetricEquality` make lookup independent of vertex ordering. Keys are
  the currency for "same polytope regardless of orientation".
- `Polytope::Type` enumerates reference geometries (Segment, Triangle,
  Quadrilateral, Tetrahedron, ...); reference nodes/centroids come from
  the type's traits.
  Their reference-coordinate dimension is exactly the topological dimension
  $d$, independently of physical embedding. In particular, the Point vertex
  and centroid belong to $\mathbb{R}^0$: both are empty coordinate tuples.
  `SpatialPoint{0}` instead constructs a one-dimensional coordinate and must
  not represent the zero-dimensional reference point.

## Connectivity — incidence d → d′

`Connectivity<Context::Local>` stores incidence relations between
dimensions. **Nothing is computed implicitly**: call
`mesh.getConnectivity().compute(d, dp)` before using an incidence, and
library code asserts the requirement via the `RODIN_GEOMETRY_REQUIRE_*`
macros in Mesh.h. `compute` has modes (e.g. `Discover` vs `Restrict`) for
whether unseen sub-polytopes get created. Typical prerequisites:
`compute(2,3)`/`compute(1,2)` for boundary/skin work, `compute(1,0)` for
edge vertices, `compute(d,d)` for neighbours.

Mutating topology means rebuilding: `Connectivity` has no element erase, so
algorithms that add/remove polytopes copy out, edit, and rebuild via the
Builder. Rebuilds re-index polytopes — never hold indices across a rebuild;
compare coordinate/key sets in tests.

## Geometry attachment — transformations

- Every polytope has a `PolytopeTransformation`: reference → physical map.
  Default is affine (from vertex coordinates); `IdentityTransformation`
  for trivial cases; `ParametricTransformation<FE>(PointCloud, fe)` for
  curved/isoparametric geometry (mesh takes ownership via
  `mesh.setPolytopeTransformation({d, i}, ...)`).
- `Geometry::Point(polytope, rc)` is the evaluation object: reference
  coordinates in, `getPhysicalCoordinates()` / `getJacobian()` (and
  inverse) out — including curvature if the transformation is parametric.
  **This is the only sanctioned way to get geometric Jacobians**; never
  reconstruct them from vertex coordinates in FES-generic code.
- `mesh.setVertexCoordinates(...)` flushes/invalidates transformations;
  `mesh.flush()` demotes everything back to affine.
- `PointCloud` — sdim × n coordinate container with Eigen views, the input
  for parametric transformations.
- `PolytopeQuadrature` caches quadrature data attached to mesh polytopes.
  Construction evaluates the cell's quadrature Jacobians in one bulk
  transformation call and stores them in the points' existing owned caches.
  Parametric P1 and H1 geometry borrow their reference tables only during this
  call; neither points nor transformations retain a table or formula pointer.
  Physical coordinates and derived metric quantities remain lazy. Mapped points
  own their reference coordinates and retain only the formula's logical identity
  and sample index until physical evaluation. A parametric map may immediately
  borrow an existing thread-local basis table by that identity; an absent or
  evicted table uses direct evaluation at the owned reference point. Formula
  destruction, assignment, cache eviction and transfer to another thread do not
  invalidate the owned sample. Rebinding clears provenance and geometric caches.
  Variational rules own only their bound polytope's mapped points. The explicit
  mesh cache is retained for callers requiring borrowed quadratures whose
  lifetime extends until geometry is flushed; variational quadrature rules
  do not populate it.

## Mesh algebra

- Constructors/factories: `Mesh::Builder` (reserve → initialize(sdim) →
  vertex(...)* → polytope(...)* → finalize; move-only), static
  `Mesh::UniformGrid(Polytope::Type, {nx, ny, ...})`, file `load`/`save`
  (MEDIT, MFEM formats; see integrations.md).
- Subset operations return `SubMesh<Context>` which remembers the parent
  and index maps: `skin()` (boundary as a (d−1)-mesh), `trim(attrs)`
  (drop cells with attributes), `keep(attrs)` (inverse). SubMesh supports
  round-tripping fields to/from the parent — prefer it over manual
  extraction.
- `trace({{attr1, attr2} → newAttr})` relabels interface facets between
  attribute pairs — the tool for tagging interfaces before FaceIntegral /
  DirichletBC on internal boundaries.
- Connected-component labeling (`ccl`) is built in (returns components as
  cell sets).

## Quadrature (QF/)

Formulas are defined on reference polytopes and selected by polytope type and
polynomial degree. `QF::PolytopeQuadratureFormula` dispatches to published
positive-interior tables where available and positive Gaussian product rules
otherwise. Variational integrators pick an order automatically but accept an
override; order 4 is the house choice for validity and diagnostic sampling of
curved cells. See quadrature.md for the exact ranges and fallback formulas.

## Partitioning and distribution

`MeshPartitioner` interface with `GreedyPartitioner` and
`BalancedCompactPartitioner`; `Sharder`/`Shard` carry partition, ownership
and overlap metadata for distributed meshes (used by the MPI stack;
Scotch provides a graph-partitioner implementation).

A local mesh or shard has the maximum dimension of its actual entities,
with dimension zero for an empty shard. A distributed mesh stores

$$
D = \max_r \dim S_r,
$$

where $S_r$ is the shard on rank $r$. Collective construction and loading
establish $D$; `getDimension()`, copy/move, and geometry-only `flush()` do not
communicate. Distributed submeshes use the same contract. Topology replacement
must go through collective construction or loading. Empty ranks participate
in reconciliation even when they have no local entities to reconcile.

## Classification and location utilities

- `MinSTCut<MeshType>` — mesh-aware binary Potts classifier, deduced from the
  constructor's mesh argument. Only the local specialization is implemented;
  the primary template is declaration-only, so MPI meshes have no serial
  fallback. Construct with a borrowed local mesh; explicitly compute
  face-to-cell incidence first.
  `classify(average)` evaluates a caller-supplied cell-average function once
  per cell and obtains volumes and interior-facet measures through
  `Polytope::getMeasure` and adjacency through `Connectivity`. `Parameters`
  are owned and configured with chainable `setParameters`; `fidelity` scales
  unary costs, and `smoothing(const Polytope&)` multiplies each interior
  facet's measure (default one, capturing lambdas supported). Results contain
  mesh-indexed `inside`, `outside`, and `cut`; applying attributes
  remains the caller's responsibility. The graph and Dinic solver are private.
  Cell values are averages, not integrals: volumes multiply them once.
  Boundary facets add no penalty; all labels are free. Nonmanifold facets
  are unsupported.
- `Location::AABB` — BVH point locator for point-in-mesh queries.

Mesh *optimization/remeshing* on this branch is MMG-based
(`MMG::MeshOptimizer`, `MMG::Adapt`); a native triangle optimizer exists
only on other branches — verify before referencing.
