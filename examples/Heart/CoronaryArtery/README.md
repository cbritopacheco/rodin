# Coronary artery models

Coronary hemodynamics and arterial wall deformation involve different spatial
models and coupling assumptions. The examples provide rigid-wall flow, wall
mechanics, and fluid–structure interaction formulations.

| Model | Formulation |
|---|---|
| [CFD](CFD/README.md) | Three-dimensional rigid-wall flow coupled to lumped ventricular and intramyocardial circulation models. |
| [Solid](Solid/README.md) | Coronary wall dynamics driven by pressure and lumped heart motion. |
| [Solid/PETSc_MPI](Solid/PETSc_MPI/README.md) | Distributed coronary wall dynamics. |
| [FSI](FSI/README.md) | ALE fluid–structure formulations with staggered or monolithic coupling. |

The CFD and FSI models have distinct wall assumptions: CFD fixes the wall,
whereas FSI solves its displacement and transfers interface data between the
fluid and solid problems. Shared headers define constitutive, transport, and
flow-assembly terms used by the example drivers.
