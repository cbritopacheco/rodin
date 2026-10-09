/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_MATH_MATRIX_H
#define RODIN_PETSC_MATH_MATRIX_H

/**
 * @file
 * @brief PETSc matrix type alias and form-language traits.
 *
 * Introduces @ref Rodin::PETSc::Math::Matrix as an alias for @c Mat and
 * provides the @ref Rodin::FormLanguage::Traits specialization so that
 * Rodin's type-trait machinery recognises PETSc matrices.
 *
 * @see <a href="_p_e_t_sc_2_math_2_vector_8h.html">Rodin::PETSc::Math::Vector</a>
 * @see <a href="class_rodin_1_1_math_1_1_linear_system_3_1_1_mat_00_01_1_1_vec_01_4.html">Rodin::PETSc::Math::LinearSystem</a>
 */

#include <boost/mpi/communicator.hpp>
#include <mpi.h>
#include <petsc.h>
#include <petscmat.h>
#include <petscsystypes.h>
#include <cassert>
#include <type_traits>
#include <utility>

#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Variational/NamedFormStorage.h"

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
   * Allows the form language to deduce the scalar type of a @c Mat.
   */
  template <>
  struct Traits<::Mat>
  {
      /// @brief Scalar value type.
      using ScalarType = PetscScalar;
  };
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Selects PETSc operators for PETSc trial coefficient storage.
   * @tparam Solution Trial solution type with PETSc vector storage.
   * @tparam Scalar Operator entry type.
   */
  template <class Solution, class Scalar>
    requires std::is_same_v<typename Traits<Solution>::DataType, ::Vec>
  struct NamedFormOperatorType<Solution, Scalar>
  {
      /// @brief PETSc matrix handle.
      using Type = ::Mat;
  };
}

namespace Rodin::Variational
{
  /// @brief Owns a named form's PETSc matrix, copying values and transferring moves.
  template <>
  class NamedFormStorage<::Mat>
  {
    public:
      /// @brief Constructs empty storage; the assembler creates the matrix.
      NamedFormStorage() = default;
      /**
       * @brief Deep-copies the matrix.
       * @param other Operator storage to copy.
       */
      NamedFormStorage(const NamedFormStorage& other)
      {
        if (other.m_operator)
        {
          const auto ierr = MatDuplicate(other.m_operator, MAT_COPY_VALUES, &m_operator);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
        }
      }
      /**
       * @brief Transfers matrix ownership.
       * @param other Operator storage to move.
       */
      NamedFormStorage(NamedFormStorage&& other) noexcept
        : m_operator(std::exchange(other.m_operator, nullptr))
      {}
      /// @brief Destroys the owned matrix.
      ~NamedFormStorage()
      {
        if (m_operator)
        {
          const auto ierr = MatDestroy(&m_operator);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
        }
      }
      /**
       * @brief Gets the writable matrix handle.
       * @returns Reference to the owned handle.
       */
      ::Mat& get()
      {
        return m_operator;
      }
      /**
       * @brief Inspects the matrix handle.
       * @returns Const reference to the owned handle.
       */
      const ::Mat& get() const
      {
        return m_operator;
      }

    private:
      ::Mat m_operator = nullptr;
  };
}

#endif
