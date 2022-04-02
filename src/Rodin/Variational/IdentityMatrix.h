#ifndef RODIN_VARIATIONAL_IDENTITYMATRIX_H
#define RODIN_VARIATIONAL_IDENTITYMATRIX_H

#include "MatrixFunction.h"

namespace Rodin::Variational
{
   class IdentityMatrix : public MatrixFunctionBase
   {
      public:
         IdentityMatrix(int n)
            : m_n(n)
         {}

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

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation&, const mfem::IntegrationPoint&) const override
         {
            value.SetSize(m_n);
            value = 0.0;
            for (int i = 0; i < m_n; i++)
               value(i, i) = 1.0;
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
