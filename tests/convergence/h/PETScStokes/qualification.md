# Local qualification records

## Retained historical observations

The finite local/MPI
matrix was locally verified with the previous dense oracle in
sequential/OpenMP configurations; this is not
hosted-CI certification or a mesh-uniform stability theorem.

The $P_3/P_2$ extension was locally verified with that oracle, independently
of the $P_2/P_1$ matrix.

Local and MPI ranks one
through four pass in sequential/OpenMP configurations; these are finite-mesh
checks, not hosted-CI certification or a mesh-uniform lower bound.

Certification of this extension and recertification of the existing
PDE spectra with the [bounded-workspace shared oracle](../Stokes/README.md#bounded-workspace-spectral-construction)
are pending. Its algebraic known-spectrum, non-nodal-constant and logical
constraint regressions are locally verified.

The finite flat-traction matrix is locally verified: 180 registrations,
360 configurations and 792 participant reports cover six geometries,
Local/MPI ranks one through four, and sequential/OpenMP assembly. Historical
full-matrix results are qualified by an exact formatting-only transformation
of the shared field-data header; a fresh Triangle Local/MPI-2 subset verifies
24 configurations and 36 reports against the formatted header. This is not
a second complete matrix run or hosted-CI certification.
