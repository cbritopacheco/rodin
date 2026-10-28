# Coronary fluid–structure interaction

The fluid velocity and pressure are solved together with arterial wall motion.
An arbitrary Lagrangian–Eulerian description accounts for the changing fluid
geometry. Staggered variants advance the component problems separately;
the monolithic prototype assembles their coupled equations.

| Example | Formulation | Requirements |
|---|---|---|
| [Explicit_P1P1_PETSc_Seq](Explicit_P1P1_PETSc_Seq/README.md) | Sequential staggered coupling with stabilized linear velocity and pressure. | PETSc |
| [Explicit_PETSc_Seq](Explicit_PETSc_Seq/README.md) | Sequential staggered coupling with quadratic velocity and linear pressure. | PETSc |
| [Explicit_PETSc_MPI](Explicit_PETSc_MPI/README.md) | PETSc/MPI staggered ALE formulation; the interface map requires a single rank. | PETSc + MPI |
| [PETSc_MPI](PETSc_MPI/README.md) | Monolithic ALE prototype on a conforming fluid–solid mesh. | PETSc + MPI |

Each variant describes its mesh inputs and execution options in its own README.
