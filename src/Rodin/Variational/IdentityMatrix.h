#ifndef RODIN_VARIATIONAL_IDENTITYMATRIX_H
#define RODIN_VARIATIONAL_IDENTITYMATRIX_H

#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   class IdentityMatrix : public MatrixCoefficientBase
   {
      public:
         constexpr
         IdentityMatrix(int n)
            : m_n(n)
         {}

         constexpr
         IdentityMatrix(const IdentityMatrix& other)
            : m_n(other.m_n)
         {}

         int getRows() const override;

         int getColumns() const override;

         void build() override;

         mfem::MatrixCoefficient& get() override;

         IdentityMatrix* copy() const noexcept override
         {
            return new IdentityMatrix(*this);
         }

      private:
         int m_n;
         std::optional<mfem::IdentityMatrixCoefficient> m_mfemMatrixCoefficient;
   };
}

#endif
