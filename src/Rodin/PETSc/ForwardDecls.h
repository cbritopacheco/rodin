#ifndef RODIN_PETSC_FORWARDDECLS_H
#define RODIN_PETSC_FORWARDDECLS_H

/**
 * @file
 * @brief Forward declarations for the @ref Rodin::PETSc module.
 */

namespace Rodin::PETSc
{
  /**
   * @namespace Rodin::PETSc
   * @brief PETSc integration module for Rodin.
   *
   * The PETSc module provides integration with the Portable, Extensible
   * Toolkit for Scientific Computation (PETSc), enabling the use of
   * distributed and parallel linear algebra objects, Krylov solvers, and
   * nonlinear solver frameworks within Rodin's finite element pipeline.
   *
   * ## Submodules
   *
   * - @ref Rodin::PETSc::Math "Math" — PETSc vector, matrix, and linear
   *   system wrappers
   * - @ref Rodin::PETSc::Solver "Solver" — KSP, SNES, CG, and GMRES
   *   solver aliases
   * - @ref Rodin::PETSc::Assembly "Assembly" — Assembly strategies for PETSc
   *   objects (sequential, MPI, OpenMP)
   * - @ref Rodin::PETSc::Variational "Variational" — Trial/test functions,
   *   forms, grid functions, and problems backed by PETSc
   *
   * @note Available only when Rodin is configured with PETSc support.
   *
   * @see Rodin::PETSc::Math, Rodin::PETSc::Solver
   */
}

#endif
