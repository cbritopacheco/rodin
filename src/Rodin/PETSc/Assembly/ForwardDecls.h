#ifndef RODIN_PETSC_ASSEMBLY_FORWARDDECLS_H
#define RODIN_PETSC_ASSEMBLY_FORWARDDECLS_H

/**
 * @file
 * @brief Forward declarations for PETSc assembly specializations.
 */

namespace Rodin::PETSc::Assembly
{
  /**
   * @namespace Rodin::PETSc::Assembly
   * @brief PETSc-specific assembly strategies.
   *
   * Provides sequential, MPI, and OpenMP assembly kernels that populate
   * PETSc vectors (`::Vec`) and matrices (`::Mat`) from Rodin variational
   * forms.
   *
   * @see Rodin::Assembly::Sequential, Rodin::Assembly::MPI
   */
}

#endif

