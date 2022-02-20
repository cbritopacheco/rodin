#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SCALARMATRIXPRODUCT_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SCALARMATRIXPRODUCT_H

#include <memory>
#include <utility>
#include <optional>

#include "Rodin/Variational/MatrixCoefficient.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   class ScalarMatrixProduct : public MatrixCoefficientBase
   {
      public:
         ScalarMatrixProduct(
               const ScalarCoefficientBase& s, const MatrixCoefficientBase& m);

         ScalarMatrixProduct(const ScalarMatrixProduct& other);

         int getRows() const override;

         int getColumns() const override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override;

         MatrixCoefficientBase* copy() const noexcept override
         {
            return new ScalarMatrixProduct(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_scalar;
         std::unique_ptr<MatrixCoefficientBase> m_matrix;
   };

   ScalarMatrixProduct
   operator*(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs);

   ScalarMatrixProduct
   operator*(const MatrixCoefficientBase& lhs, const ScalarCoefficientBase& rhs);
}

#endif
