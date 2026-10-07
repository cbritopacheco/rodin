# Distributed constraint assembly

These tests exercise MPI assembly without a linear-algebra backend. The target
links Rodin::MPI, not Rodin::PETSc, and runs on 1, 2, 3, 4, and 8 ranks.

## Ownership and overlap contract

Each distributed entity (dimension, global index) has exactly one owner.
This does not require incident entities of different dimensions to have the
same owner. Constraint assembly must use the existing overlap rather than
restrict its traversal to owned boundary faces.

For rank $r$, let $O_r$ denote its owned DOFs and $C_r$ its owned cells.
The required constraint indices are

$$
R_r = O_r \cup \bigcup_{K\in C_r}\mathrm{DOFs}(K).
$$

For P1 and H1, the vertex-star overlap contains every physical boundary
face incident to these indices. Artificial faces on the outer edge of the
overlap cannot touch $R_r$: both cells adjacent to any such interior face
would already belong to that vertex-star overlap. Boundary functionals are
therefore selected from shard-local faces, including ghosts, restricted to
$R_r$. For each index, the eligible face with the smallest distributed index
is selected. Attribute-selected interior faces remain supported.

No MPI communication is performed during P1/H1 constraint assembly.
The same selection is used for prescribed values, identification rows, and
affine offsets. Ownership, connectivity, and global numbering are unchanged.
Selection uses integer indices, not coordinate matching or numerical tolerances.
The generic variational layer delegates distributed affine evaluation to the
MPI assembler and contains no shard-specific selection logic.

P0g is different: its globally supported DOFs may be owned on a rank with no
copy of the selected face. The smallest eligible owned face is selected
globally, and its evaluated payload is broadcast unchanged. This collective
operation is confined to P0g; it is not a general constraint exchange.

Source selection assumes compatible prescribed data. It does not establish
compatibility between separate boundary-condition objects or certify arbitrary
transitive identification graphs.

## Independent structural checks

The entity-owner test gathers identities independently and checks one owner
per entity, agreement of every holder's owner declaration, and exact equality
between the owner's halo and the gathered holder set.

The physical-boundary incidence test records physical boundary face indices
on the original mesh before sharding. For every index in $R_r$, its physical
boundary incidence set is compared with that oracle. Artificial shard-boundary
faces are checked to have no DOF in $R_r$. This test does not use the constraint
assembler to define its expected result.

The matrix covers segment, triangle, quadrilateral, tetrahedron, pyramid,
hexahedron, and wedge. UniformGrid sizes are 17 points for segments, 9 points
per axis for tetrahedra, and 5 points per axis otherwise. These are fixed-mesh
structural checks, not refinement-rate measurements.

| Spaces | Ranges | Contract |
| --- | --- | --- |
| P1, H1 orders 1–6, P0g | Real/complex scalar and real/complex vector with three components | Owned-value completeness and required affine-identification rows |
| P1, H1 orders 1–6, P0g | Real/complex 2-by-3 matrix | Exact value/identification index sets, matching master indices and empty P0g holders |
| P1, H1 order 2 | Real | Constraint assembly on one rank without collective communication |

Further regressions check remote P0g boundary selection, matching linear and
affine source selection, exact preservation of a nonzero $10^{-30}$ coefficient,
empty master rows, and distinct vector spaces with two and three components.
Slave indices must not be overwritten by a master lookup sharing a scratch
buffer. P1 local/global index round-trips are checked exactly.

Tagged interior-face selection and reselection to an absent attribute are
checked for P1, H1 order 2, and P0g. Expected DOF sets are obtained independently
from owned tagged faces and compared by integer indices. Reselection must clear
both the linear rows and the affine values, including on ranks with no selected
faces.

P0 exposes cell DOFs rather than face traces; its numbering is tested separately
in the MPI Variational suite. A 0D mesh has no codimension-one boundary and is
outside this boundary-condition test. These results concern the tested
vertex-overlap partitions, not every possible externally constructed shard.

## SubMesh extraction regression

Unique ownership does not imply complete overlap. Distributed SubMesh
finalization completes selection on existing parent holders before
resolving ownership. The selector still requires a parent with complete overlap;
it does not validate externally constructed incomplete parents at runtime.

The permanent P1/P2 regression compares exact constraint sets with the physical
boundary on all seven geometries after owned-cell-only extraction. Geometry-level
coverage includes proper subregions, skins, sparse/nonowner selections, and
nested extraction; see [SubMesh verification](../Geometry/README.md).

Global space sizes are established during collective construction. Size
queries, copies and moves are noncollective for the fixed-layout P0, P0g,
P1 and H1 spaces; a separate bounded test queries them on rank zero alone,
including empty shards and Point meshes. Index assertions are exact;
numerical PDE errors and refinement rates are tested separately in
[PETSc MPI Poisson](../../../../convergence/h/PETScMPIPoisson/README.md).

## Backend-independent build

Matrix-valued trace coverage uses real and complex 2-by-3 ranges in P0g,
P1 and H1 orders 1–6 on all seven positive-dimensional geometries. Value
constraints, affine-identification slave rows and their defect maps must
have exactly the required physical-boundary index set. Physical faces are
identified independently by their original parent IDs, retained before
sharding; owned DOFs and DOFs of owned cells define each rank's required
set for P1 and H1. Globally supported P0g instead retains all components
on every holder, including empty shards, under its collective payload
selection contract described above. The identification row must contain
the matching master index. No coordinate comparisons, floating-point rank
decisions or additional reconciliation exchanges participate in this oracle.
After collective field construction, P1/H1 matrix constraints are first
assembled on rank zero alone; a test-protocol barrier then precedes assembly
on all ranks. This checks noncollective evaluation and repeated assembly
without turning a local operation into a global one. P0g is excluded from
that single-rank phase. P0 has no conforming boundary trace and is not included;
Point metadata and restriction contracts
have separate tests.

The fixed-layout metadata test includes real/complex scalar, three-component
vector and 2-by-3 matrix ranges in P0, P0g, P1 and H1 orders 1–6. Each space
is constructed collectively. Its global size is then checked against an
independently reduced ownership count, while rank zero alone performs the
size query and copies and moves the space. Point meshes are included, with
the sole point placed on the last rank to exercise empty holders.

The matrix probes use three grid points per axis and release each H1 order
before constructing the next. Their geometry/rank CTest entries are named
`RodinMPIAssemblyMatrixTraceSlowTest_np*`, labelled `unit;distributed;slow`
with a 600-second protocol watchdog, on 1, 2, 3, 4 and 8 ranks. This limit
is not a numerical acceptance criterion. These probes supplement, rather
than replace, the larger tetrahedral ownership regression.

Configure with MPI enabled, PETSc disabled, and unit tests enabled. Build
RodinMPIAssemblyTest. The regular `RodinMPIAssemblyTest_np*` CTest entries
exclude the four 3D affine-identification cases, the two tetrahedral
boundary-incidence/constraint probes, the fixed-layout size-query tests and
the matrix-trace probes. The regular group has a 120-second limit.
The four 3D affine-identification cases run on each of 1, 2, 3, 4, and 8 ranks as
`RodinMPIAssemblyIdentificationSlowTest_np*` entries, labelled
`unit;distributed;slow` with a 600-second limit. CI's Release slow job runs
them separately. The two tetrahedral probes use
`RodinMPIAssemblyGeometrySlowTest_np*` with a 600-second limit; fixed-layout
size queries use `RodinMPISpaceSizeTest_np*` with a 30-second limit.
No geometry or rank is omitted. Run the complete assembly and size-query matrix with
`ctest --test-dir <build>/tests -R '^Rodin(MPIAssembly|MPISpaceSize)' --output-on-failure`.
