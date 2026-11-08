/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_LINEARSYSTEM_H
#define RODIN_KOKKOS_LINEARSYSTEM_H

/**
 * @file
 * @brief Linear system whose matrix stays in the memory of a Kokkos backend.
 */

#include "Rodin/Math/LinearSystem.h"

#include "SparseMatrix.h"

namespace Rodin::Math
{
  /**
   * @brief Linear system with a Kokkos::SparseMatrix operator and host
   * vectors.
   *
   * The device assembly of a Problem fills the matrix in the memory of the
   * backend, and a solver running on the device reads it there. The
   * right-hand side and the solution stay on the host, where the grid
   * functions keep their data. The constraints are eliminated during the
   * assembly, so the system offers no elimination of its own.
   *
   * A Problem assembles into this system when it is named in its template
   * arguments:
   * @code{.cpp}
   * using System = Math::LinearSystem<
   *   Rodin::Kokkos::SparseMatrix<Real, Backend>, Math::Vector<Real>>;
   * Problem<System, decltype(u), decltype(v)> poisson(u, v);
   * @endcode
   *
   * @tparam MatrixScalar Coefficient type of the matrix.
   * @tparam Backend Kokkos backend of the matrix.
   * @tparam VectorScalar Coefficient type of the vectors.
   */
  template <class MatrixScalar, class Backend, class VectorScalar>
  class LinearSystem<::Rodin::Kokkos::SparseMatrix<MatrixScalar, Backend>,
    Math::Vector<VectorScalar>>
    : public LinearSystemBase<::Rodin::Kokkos::SparseMatrix<MatrixScalar, Backend>,
        Math::Vector<VectorScalar>,
        LinearSystem<::Rodin::Kokkos::SparseMatrix<MatrixScalar, Backend>,
          Math::Vector<VectorScalar>>>
  {
    public:
      /// @brief Device matrix type of the operator.
      using MatrixType = ::Rodin::Kokkos::SparseMatrix<MatrixScalar, Backend>;

      /// @brief Host vector type of the right-hand side and the solution.
      using VectorType = Math::Vector<VectorScalar>;

      /// @brief Parent class type.
      using Parent =
        LinearSystemBase<MatrixType, VectorType, LinearSystem<MatrixType, VectorType>>;

      /// @brief Empty system.
      LinearSystem() = default;

      /// @brief Copies the system, the matrix included.
      LinearSystem(const LinearSystem&) = default;

      /// @brief Moves the system.
      LinearSystem(LinearSystem&&) noexcept = default;

      /// @brief Destructor.
      virtual ~LinearSystem() = default;

      /**
       * @brief Copies the system, the matrix included.
       * @return Reference to this system.
       */
      LinearSystem& operator=(const LinearSystem&) = default;

      /**
       * @brief Moves the system.
       * @return Reference to this system.
       */
      LinearSystem& operator=(LinearSystem&&) noexcept = default;

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
