#ifndef RODIN_VARIATIONAL_BIlinearFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BIlinearFORMINTEGRATOR_H

#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class BilinearFormIntegratorBase : public FormLanguage::Base
   {
      public:
         BilinearFormIntegratorBase& over(int attr)
         {
            return over(std::vector{ attr });
         }

         BilinearFormIntegratorBase& over(const std::vector<int>& attrs)
         {
            m_attrs = attrs;
            return *this;
         }

         const std::vector<int>& getAttributes() const
         {
            return m_attrs;
         }

         virtual void buildMFEMBilinearFormIntegrator() = 0;

         virtual mfem::BilinearFormIntegrator& getMFEMBilinearFormIntegrator() = 0;

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

      private:
         std::vector<int> m_attrs;
   };

   class BilinearFormDomainIntegrator : public BilinearFormIntegratorBase
   {
      public:
         virtual BilinearFormDomainIntegrator* copy() const noexcept override = 0;
   };
}

#endif
