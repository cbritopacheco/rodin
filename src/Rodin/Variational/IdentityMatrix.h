/**
 * @file
 * @brief Identity matrix function.
 */

#ifndef RODIN_VARIATIONAL_IDENTITYMATRIX_H
#define RODIN_VARIATIONAL_IDENTITYMATRIX_H

#include "MatrixFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Represents the identity matrix function @f$ I_n @f$.
   *
   * This class represents the matrix function that returns the @f$ n @f$ 
   * dimensional identity matrix at each point:
   * @f[
   *    F(x) = I_n = \begin{pmatrix}
   *      1 & 0 & \cdots & 0 \\
   *      0 & 1 & \cdots & 0 \\
   *      \vdots & \vdots & \ddots & \vdots \\
   *      0 & 0 & \cdots & 1
   *    \end{pmatrix}
   * @f]
   *
   * Commonly used in finite element formulations for identity operators
   * or as a reference configuration.
   *
   * @note The identity matrix is constant in space and symmetric.
   */
  class IdentityMatrix : public MatrixFunctionBase<Real, IdentityMatrix>
  {
    public:
      using ScalarType = Real;

      using Parent = MatrixFunctionBase<ScalarType, IdentityMatrix>;

      /**
       * @brief Constructs the identity matrix function.
       * @param n Dimension of identity matrix @f$ n \times n @f$
       */
      IdentityMatrix(size_t n)
        : m_n(n)
      {}

      /**
       * @brief Copy constructor.
       * @param other Identity matrix to copy
       */
      IdentityMatrix(const IdentityMatrix& other)
        : MatrixFunctionBase(other),
          m_n(other.m_n)
      {}

      /**
       * @brief Move constructor.
       * @param other Identity matrix to move
       */
      IdentityMatrix(IdentityMatrix&& other)
        : MatrixFunctionBase(std::move(other)),
          m_n(other.m_n)
      {}

      /**
       * @brief Gets number of rows.
       * @returns Dimension @f$ n @f$
       */
      size_t getRows() const
      {
        return m_n;
      }

      /**
       * @brief Gets number of columns.
       * @returns Dimension @f$ n @f$
       */
      size_t getColumns() const
      {
        return m_n;
      }

      /**
       * @brief Evaluates the identity matrix at a point.
       * @returns Identity matrix @f$ I_n @f$
       */
      auto getValue(const Geometry::Point&) const
      {
        return Math::Matrix<Real>::Identity(m_n, m_n);
      }

      /**
       * @brief Creates a polymorphic copy of the identity matrix.
       * @return Pointer to copy
       */
      IdentityMatrix* copy() const noexcept override
      {
        return new IdentityMatrix(*this);
      }

    private:
      const size_t m_n;
  };
}

#endif
