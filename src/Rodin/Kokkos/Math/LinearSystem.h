/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MATH_LINEARSYSTEM_H
#define RODIN_KOKKOS_MATH_LINEARSYSTEM_H

/**
 * @file
 * @brief Linear system whose matrix stays in the memory of a Kokkos backend.
 */

#include "Rodin/Math/LinearSystem.h"

#include "SparseMatrix.h"
#include "Vector.h"

namespace Rodin::Math
{
  /**
   * @brief Linear system with a Kokkos sparse operator and selected vector storage.
   *
   * ## Architecture
   * Assembly fills the operator in the backend memory space. Backend vectors
   * keep the right-hand side and solution in that space as well; explicitly
   * selected host vectors retain compatibility with host coefficient clients.
   * Constraints are eliminated by Assembly, rather than by the algebra owner.
   *
   * A Problem assembles into this system when it is named in its template
   * arguments:
   * @code{.cpp}
   * using System = ::Rodin::Math::LinearSystem<
   *   Rodin::Kokkos::Math::SparseMatrix<Real, Backend>,
   *   Rodin::Kokkos::Math::Vector<Real, Backend>>;
   * Problem<System, decltype(u), decltype(v)> poisson(u, v);
   * @endcode
   *
   * @tparam MatrixScalar Coefficient type of the matrix.
   * @tparam Backend Kokkos backend of the matrix.
   * @tparam VectorType_ Coefficient vector container.
   */
  template <class MatrixScalar, class Backend, class VectorType_>
  class LinearSystem<::Rodin::Kokkos::Math::SparseMatrix<MatrixScalar, Backend>,
    VectorType_>
    : public LinearSystemBase<::Rodin::Kokkos::Math::SparseMatrix<MatrixScalar, Backend>,
        VectorType_,
        LinearSystem<::Rodin::Kokkos::Math::SparseMatrix<MatrixScalar, Backend>,
          VectorType_>>
  {
    public:
      /// @brief Device matrix type of the operator.
      using MatrixType = ::Rodin::Kokkos::Math::SparseMatrix<MatrixScalar, Backend>;

      /// @brief Storage type of the right-hand side and the solution.
      using VectorType = VectorType_;

      /// @brief Parent class type.
      using Parent =
        LinearSystemBase<MatrixType, VectorType, LinearSystem<MatrixType, VectorType>>;

      /// @brief Empty system.
      LinearSystem() = default;

      /**
       * @brief Copies the system, the matrix included.
       * @param other System to copy.
       */
      LinearSystem(const LinearSystem& other) = default;

      /**
       * @brief Moves the system.
       * @param other System to move.
       */
      LinearSystem(LinearSystem&& other) noexcept = default;

      /// @brief Destructor.
      virtual ~LinearSystem() = default;

      /**
       * @brief Copies the system, the matrix included.
       * @param other System to copy.
       * @return Reference to this system.
       */
      LinearSystem& operator=(const LinearSystem& other) = default;

      /**
       * @brief Moves the system.
       * @param other System to move.
       * @return Reference to this system.
       */
      LinearSystem& operator=(LinearSystem&& other) noexcept = default;

      /**
       * @brief Returns the operator @f$ A @f$.
       * @return Mutable reference to the device matrix.
       */
      MatrixType& getOperator()
      {
        return m_operator;
      }

      /**
       * @brief Returns the operator @f$ A @f$.
       * @return Read-only reference to the device matrix.
       */
      const MatrixType& getOperator() const
      {
        return m_operator;
      }

      /**
       * @brief Returns the right-hand side @f$ b @f$.
       * @return Mutable reference to the right-hand side.
       */
      VectorType& getVector()
      {
        return m_vector;
      }

      /**
       * @brief Returns the right-hand side @f$ b @f$.
       * @return Read-only reference to the right-hand side.
       */
      const VectorType& getVector() const
      {
        return m_vector;
      }

      /**
       * @brief Returns the solution @f$ x @f$.
       * @return Mutable reference to the solution.
       */
      VectorType& getSolution()
      {
        return m_solution;
      }

      /**
       * @brief Returns the solution @f$ x @f$.
       * @return Read-only reference to the solution.
       */
      const VectorType& getSolution() const
      {
        return m_solution;
      }

    private:
      MatrixType m_operator;
      VectorType m_vector;
      VectorType m_solution;
  };
}

#endif
