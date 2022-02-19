#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <mfem.hpp>

#include "FormLanguage/ForwardDecls.h"

#include "ForwardDecls.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   template <class T>
   class Integral;

   template <>
   class Integral<FormLanguage::TrialTestProduct> : public BilinearFormIntegratorBase
   {
      public:
         Integral(const FormLanguage::TrialTestProduct& prod);

         Integral(const Integral& other);

         const std::set<int>& getAttributes() const override;

         IntegratorRegion getIntegratorRegion() const override;

         void build() override;

         mfem::BilinearFormIntegrator& get() override;

         mfem::BilinearFormIntegrator* release() override;

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
   };
}

#endif
