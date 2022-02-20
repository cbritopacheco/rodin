#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SCALARPRODUCT_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SCALARPRODUCT_H

#include "Rodin/Variational/ScalarCoefficient.h"

namespace Rodin::Variational::FormLanguage
{
   class ScalarProduct : public ScalarCoefficientBase
   {
      public:
         ScalarProduct(const ScalarCoefficientBase& a, const ScalarCoefficientBase& b);

         ScalarProduct(const ScalarProduct& other);

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override;

         ScalarProduct* copy() const noexcept override
         {
            return new ScalarProduct(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_a;
         std::unique_ptr<ScalarCoefficientBase> m_b;
   };

   ScalarProduct
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);
}

#endif

