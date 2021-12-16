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

         void buildMFEMCoefficient() override;

         mfem::Coefficient& getMFEMCoefficient() override;

         ScalarProduct* copy() const noexcept override
         {
            return new ScalarProduct(*this);
         }

      private:
         std::unique_ptr<ScalarCoefficientBase> m_a;
         std::unique_ptr<ScalarCoefficientBase> m_b;
         std::optional<mfem::ProductCoefficient> m_mfemCoefficient;
   };

   ScalarProduct
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);
}

#endif

