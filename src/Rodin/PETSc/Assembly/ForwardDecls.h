/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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
   * PETSc vectors (@c Vec) and matrices (@c Mat) from Rodin variational
   * forms.
   *
   * @see Rodin::Assembly::Sequential, Rodin::Assembly::MPI
   */
}

#endif

