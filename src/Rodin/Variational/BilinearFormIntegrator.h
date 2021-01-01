#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H

#include <set>
#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class BilinearFormIntegratorBase : public FormLanguage::Base
   {
      public:
         virtual ~BilinearFormIntegratorBase() = default;

         virtual const std::set<int>& getAttributes() const = 0;

         virtual void buildMFEMBilinearFormIntegrator() = 0;

         virtual mfem::BilinearFormIntegrator& getMFEMBilinearFormIntegrator() = 0;

         virtual IntegratorRegion getIntegratorRegion() const = 0;

         /**
          * @internal
          * @brief Releases ownership of the mfem::BilinearFormIntegrator.
          *
          * @note After this call, calling getMFEMBilinearFormIntegrator()
          * will result in undefined behaviour.
          *
          * @warning The BilinearFormIntegratorBase instance must still be kept
          * in memory since it might contain objects which the
          * mfem::BilinearFormIntegrator instance refers to.
          */
         virtual mfem::BilinearFormIntegrator* releaseMFEMBilinearFormIntegrator() = 0;

         virtual BilinearFormIntegratorBase* copy() const noexcept override = 0;
   };

   class BilinearFormDomainIntegrator : public BilinearFormIntegratorBase
   {
      public:
         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Domain;
         }

         virtual BilinearFormDomainIntegrator& over(int attr) = 0;

         virtual BilinearFormDomainIntegrator& over(const std::set<int>& attrs) = 0;

         virtual BilinearFormDomainIntegrator* copy() const noexcept override = 0;
   };
}

#endif
