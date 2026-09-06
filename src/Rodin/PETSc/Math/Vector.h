/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_MATH_VECTOR_H
#define RODIN_PETSC_MATH_VECTOR_H

/**
 * @file
 * @brief PETSc vector type alias and form-language traits.
 *
 * Introduces @ref Rodin::PETSc::Math::Vector as an alias for @c Vec and
 * provides the @ref Rodin::FormLanguage::Traits specialization so that
 * Rodin's type-trait machinery recognises PETSc vectors.
 *
 * @see <a href="_p_e_t_sc_2_math_2_matrix_8h.html">Rodin::PETSc::Math::Matrix</a>
 * @see <a href="class_rodin_1_1_math_1_1_linear_system_3_1_1_mat_00_01_1_1_vec_01_4.html">Rodin::PETSc::Math::LinearSystem</a>
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
   * Allows the form language to deduce the scalar type of a @c Vec.
   */
  template <>
  struct Traits<::Vec>
  {
    /// @brief Scalar value type.
      using ScalarType = PetscScalar;
  };
}

#endif
