#ifndef RODIN_VARIATIONAL_IDENTITYMATRIX_H
#define RODIN_VARIATIONAL_IDENTITYMATRIX_H

#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   class IdentityMatrix : public MatrixCoefficientBase
   {
      public:
         IdentityMatrix(int n)
            : m_n(n)
         {}

         IdentityMatrix(const IdentityMatrix& other)
            :  MatrixCoefficientBase(other),
               m_n(other.m_n)
         {}

         IdentityMatrix(IdentityMatrix&& other)
            :  MatrixCoefficientBase(std::move(other)),
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
