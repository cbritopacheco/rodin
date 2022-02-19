#ifndef RODIN_VARIATIONAL_IDENTITYMATRIX_H
#define RODIN_VARIATIONAL_IDENTITYMATRIX_H

#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   class IdentityMatrix : public MatrixCoefficientBase
   {
      public:
         IdentityMatrix(int n)
            : m_n(n),
              m_mfemMatrixCoefficient(m_n)
         {}

         IdentityMatrix(const IdentityMatrix& other)
            : m_n(other.m_n),
              m_mfemMatrixCoefficient(other.m_n)
         {}

         int getRows() const override;

         int getColumns() const override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override
         {
            m_mfemMatrixCoefficient.Eval(value, trans, ip);
         }

         IdentityMatrix* copy() const noexcept override
         {
            return new IdentityMatrix(*this);
         }

      private:
         int m_n;
         mfem::IdentityMatrixCoefficient m_mfemMatrixCoefficient;
   };
}

#endif
