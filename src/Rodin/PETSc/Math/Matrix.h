#ifndef RODIN_PETSC_MATH_MATRIX_H
#define RODIN_PETSC_MATH_MATRIX_H

/**
 * @file
 * @brief PETSc matrix aliases and traits integration.
 */

#include <boost/mpi/communicator.hpp>
#include <mpi.h>
#include <petsc.h>
#include <petscmat.h>
#include <petscsystypes.h>

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::PETSc::Math
{
  /**
   * @brief Alias to PETSc matrix handle type.
   */
  using Matrix = ::Mat;
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Form-language traits for PETSc matrices.
   */
  template <>
  struct Traits<::Mat>
  {
    using ScalarType = PetscScalar;
  };
}

#endif
