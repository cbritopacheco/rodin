#ifndef RODIN_PETSC_ASSEMBLY_DEFAULT_H
#define RODIN_PETSC_ASSEMBLY_DEFAULT_H

/**
 * @file
 * @brief Default PETSc assembly selector.
 */

#include "Rodin/Assembly/Default.h"

namespace Rodin::PETSc::Assembly
{
  /**
   * @brief Alias to Rodin's generic assembly selector.
   *
   * This alias mirrors @ref Rodin::Assembly::Default and allows PETSc code
   * paths to select the matching assembly strategy (sequential, MPI, OpenMP).
   */
  template <class ... Ts>
  using Default = Rodin::Assembly::Default<Ts...>;
}

#endif
