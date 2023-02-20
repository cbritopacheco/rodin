#ifndef RODIN_VARIATIONAL_IDENTITYMATRIX_H
#define RODIN_VARIATIONAL_IDENTITYMATRIX_H

#include "MatrixFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Represents the identity matrix function.
   *
   * This class represents the matrix function, which at each point, returns
   * the @f$ n @f$ dimensional identity matrix:
   *
   * @f$
   *   F(x) = I_n
   * @f$
   */
  class IdentityMatrix : public MatrixFunctionBase<IdentityMatrix>
  {
    public:
      /**
       * @brief Constructs the identity matrix function.
       * @param[in] n Dimension of identity matrix
       */
      constexpr
      IdentityMatrix(size_t n)
        : m_n(n)
      {}

      constexpr
      IdentityMatrix(const IdentityMatrix& other)
        : MatrixFunctionBase(other),
          m_n(other.m_n)
      {}

      constexpr
      IdentityMatrix(IdentityMatrix&& other)
        : MatrixFunctionBase(std::move(other)),
          m_n(other.m_n)
      {}

      inline
      constexpr
      size_t getRows() const
      {
        return m_n;
      }

      inline
      constexpr
      size_t getColumns() const
      {
        return m_n;
      }

      inline
      auto getValue(const Geometry::Point&) const
      {
        return Math::Matrix::Identity(m_n, m_n);
      }

    private:
      const size_t m_n;
  };
}

#endif
