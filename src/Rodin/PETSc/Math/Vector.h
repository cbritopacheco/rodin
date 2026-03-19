#ifndef RODIN_PETSC_MATH_VECTOR_H
#define RODIN_PETSC_MATH_VECTOR_H

/**
 * @file
 * @brief PETSc vector type alias and form-language traits.
 *
 * Introduces @ref Rodin::PETSc::Math::Vector as an alias for `::Vec` and
 * provides the @ref Rodin::FormLanguage::Traits specialization so that
 * Rodin's type-trait machinery recognises PETSc vectors.
 *
 * @see Rodin::PETSc::Math::Matrix, Rodin::PETSc::Math::LinearSystem
 */

#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <petsc.h>
#include <petscsys.h>
#include <petscsystypes.h>

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::PETSc::Math
{
  /**
   * @brief Alias for the PETSc distributed/sequential vector handle.
   *
   * Used as the data type for @ref Rodin::Variational::GridFunction and
   * @ref Rodin::Variational::LinearForm specializations throughout the
   * PETSc backend.
   */
  using Vector = ::Vec;
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Traits specialization for PETSc vectors.
   *
   * Allows the form language to deduce the scalar type of a `::Vec`.
   */
  template <>
  struct Traits<::Vec>
  {
    using ScalarType = PetscScalar;
  };
}

#endif
