#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SCALARMATRIXPRODUCT_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SCALARMATRIXPRODUCT_H

#include <memory>
#include <utility>
#include <optional>

#include "Rodin/Variational/MatrixCoefficient.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::FormLanguage
{
   namespace Internal
   {
      class ScalarMatrixProductCoefficient : public mfem::MatrixCoefficient
      {
         public:
            ScalarMatrixProductCoefficient(
                  mfem::Coefficient& s, mfem::MatrixCoefficient& m)
               :  mfem::MatrixCoefficient(m.GetHeight(), m.GetWidth()),
                  m_scalar(s),
                  m_matrix(m)
            {}

            virtual void Eval(
                  mfem::DenseMatrix& K,
                  mfem::ElementTransformation& T,
                  const mfem::IntegrationPoint& ip) override
            {
               m_matrix.Eval(K, T, ip);
               K *= m_scalar.Eval(T, ip);
            }

         private:
            mfem::Coefficient& m_scalar;
            mfem::MatrixCoefficient& m_matrix;
      };
   }

   class ScalarMatrixProduct : public MatrixCoefficientBase
   {
      public:
         ScalarMatrixProduct(
               const ScalarCoefficientBase& s, const MatrixCoefficientBase& m);

         ScalarMatrixProduct(const ScalarMatrixProduct& other);

         int getRows() const override;

         int getColumns() const override;

         void buildMFEMMatrixCoefficient() override;

         mfem::MatrixCoefficient& getMFEMMatrixCoefficient() override;

         MatrixCoefficientBase* copy() const noexcept override
         {
            return new ScalarMatrixProduct(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_scalar;
         std::unique_ptr<MatrixCoefficientBase> m_matrix;
         std::optional<Internal::ScalarMatrixProductCoefficient> m_mfemMatrixCoefficient;
   };

   ScalarMatrixProduct
   operator*(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs);

   ScalarMatrixProduct
   operator*(const MatrixCoefficientBase& lhs, const ScalarCoefficientBase& rhs);
}

#endif
