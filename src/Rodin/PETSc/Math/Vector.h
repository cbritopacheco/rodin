#ifndef RODIN_PETSC_MATH_VECTOR_H
#define RODIN_PETSC_MATH_VECTOR_H

/**
 * @file
 * @brief PETSc vector aliases and traits integration.
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
   * @brief Alias to PETSc dense/distributed vector handle type.
   *
   * This alias is used across the PETSc backend so variational forms, assembly,
   * and solvers can target PETSc vectors through Rodin type traits.
   */
  using Vector = ::Vec;
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Form-language traits for PETSc vectors.
   */
  template <>
  struct Traits<::Vec>
  {
    using ScalarType = PetscScalar;
  };
}

#endif
