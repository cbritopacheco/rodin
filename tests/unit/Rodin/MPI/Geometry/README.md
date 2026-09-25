# Distributed SubMesh verification

SubMesh finalization must preserve global selection, unique ownership, and the
parent's overlap needed for local incidence queries. For selected top-dimensional
entity IDs $S_r$ requested by rank $r$, and IDs $P_r$ present in its parent shard,
the completed local selection is

$$
S = \bigcup_r S_r, \qquad S_r^{\mathrm{local}} = S\cap P_r.
$$

Requests are sent to parent owners; the selected IDs are then sent to parent
holders. Existing builder operations include their subentities. Ownership is
resolved afterward, and Shared/Ghost membership is recomputed relative to owned
SubMesh cells. Parent IDs, coordinates, and restriction maps are retained.
This requires a parent with complete vertex-star overlap; it does not repair
an externally supplied parent with missing topology.

## Regression and coverage

`SelectionClosureAcrossGeometries` tests all seven positive-dimensional cell
geometries with 17 points for segments and 5 points per axis otherwise:

- Full-domain selection through owned cells only.
- Proper subregions selected by even distributed cell IDs, exposing new boundaries.
- Boundary skins, including 0D skins of segments.
- Single-cell selection from rank zero, including ranks without selected entities.
- A request submitted only by a nonowner.
- Nested extraction, preserving the completed local selection.

Global selection is gathered independently in test code. The tests compare exact
entity ID sets, require exactly one owner at every dimension, and classify owned
facets using independently gathered selected-cell incidence. On these manifold
meshes, a facet is a boundary precisely when it has one selected incident cell.
No coordinate or floating-point tolerance determines identity or boundaries.

The MPI Assembly suite additionally checks P1/P2 constraint index sets against
physical boundary IDs recorded before sharding, on all seven geometries. Its
9-by-9-point triangular case exercises boundary classification on three ranks.
The PETSc MPI Poisson suite also solves the exact quadratic patch on the extracted
SubMesh with value and affine traces on all seven geometries and 1–4 ranks.

Sparse nested extraction also tests two related dimension requirements: empty
shards need metadata capacity for the collective dimension, and MPI cell/face
iterators must use that dimension rather than an empty shard's local dimension.

Run `RodinMPIGeometrySubMeshTest_np*` on 1, 2, 3, 4, and 8 ranks. The existing
`RodinMPISubMeshGridFunctionTest_np*` suite supplies field-numbering and restriction
coverage on 1–4 ranks. PDE rates remain in the convergence module; these are
fixed-mesh structural regressions, not convergence-rate measurements.
