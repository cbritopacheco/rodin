#ifndef RODIN_PETSC_ASSEMBLY_DEFAULT_H
#define RODIN_PETSC_ASSEMBLY_DEFAULT_H

/**
 * @file
 * @brief Default assembly selector for PETSc objects.
 */

#include "Rodin/Assembly/Default.h"

namespace Rodin::PETSc::Assembly
{
  /**
   * @brief Alias for @ref Rodin::Assembly::Default, selecting the
   * appropriate PETSc assembly strategy (sequential, MPI, or OpenMP)
   * based on the mesh context type.
   */
  template <class ... Ts>
  using Default = Rodin::Assembly::Default<Ts...>;
}

#endif
