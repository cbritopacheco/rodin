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

Formulas are defined on reference polytopes and selected by
(polytope type, order): `QF::GaussLegendre`, `QF::GaussLobatto`,
`QF::GrundmannMoller` (arbitrary-order simplex), `QF::Centroid` (1-point).
`QF::PolytopeQuadratureFormula` is the runtime dispatcher. Variational
integrators pick an order automatically but accept an override; order-4 is
the house choice for validity/diagnostic sampling of curved cells.

## Partitioning and distribution

`MeshPartitioner` interface with `GreedyPartitioner` and
`BalancedCompactPartitioner`; `Sharder`/`Shard` carry partition, ownership
and overlap metadata for distributed meshes (used by the MPI stack;
Scotch provides a graph-partitioner implementation).

## Classification and location utilities

- `MinSTCut` — serial binary Potts classifier via s-t min cut (converting
  a level set into discrete cell labels with perimeter regularization).
- `Location::AABB` — BVH point locator for point-in-mesh queries.

Mesh *optimization/remeshing* on this branch is MMG-based
(`MMG::MeshOptimizer`, `MMG::Adapt`); a native triangle optimizer exists
only on other branches — verify before referencing.
