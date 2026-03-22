#ifndef RODIN_PETSC_MATH_MATRIX_H
#define RODIN_PETSC_MATH_MATRIX_H

/**
 * @file
 * @brief PETSc matrix type alias and form-language traits.
 *
 * Introduces @ref Rodin::PETSc::Math::Matrix as an alias for `::Mat` and
 * provides the @ref Rodin::FormLanguage::Traits specialization so that
 * Rodin's type-trait machinery recognises PETSc matrices.
 *
 * @see Rodin::PETSc::Math::Vector, Rodin::PETSc::Math::LinearSystem
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
   * @brief Alias for the PETSc sparse/dense matrix handle.
   *
   * Used as the operator type for @ref Rodin::Variational::BilinearForm
   * specializations and as the system matrix inside
   * @ref Rodin::PETSc::Math::LinearSystem.
   */
  using Matrix = ::Mat;
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Traits specialization for PETSc matrices.
   *
   * Allows the form language to deduce the scalar type of a `::Mat`.
   */
  template <>
  struct Traits<::Mat>
  {
    using ScalarType = PetscScalar;
  };
}

#endif
