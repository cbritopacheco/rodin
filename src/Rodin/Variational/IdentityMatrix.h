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
  class IdentityMatrix : public MatrixFunctionBase
  {
    public:
      /**
       * @brief Constructs the identity matrix function.
       * @param[in] n Dimension of identity matrix
       */
      IdentityMatrix(int n)
        : m_n(n)
      {
        assert(n > 0);
      }

      IdentityMatrix(const IdentityMatrix& other)
        :  MatrixFunctionBase(other),
          m_n(other.m_n)
      {}

      IdentityMatrix(IdentityMatrix&& other)
        :  MatrixFunctionBase(std::move(other)),
          m_n(other.m_n)
      {}

      int getRows() const override
      {
        return m_n;
      }

      int getColumns() const override
      {
        return m_n;
      }

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        Math::Matrix res = Math::Matrix::Identity(m_n, m_n);
        return res;
      }

      IdentityMatrix* copy() const noexcept override
      {
        return new IdentityMatrix(*this);
      }

    private:
      const int m_n;
  };
}

#endif
