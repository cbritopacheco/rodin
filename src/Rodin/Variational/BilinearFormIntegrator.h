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
         virtual BilinearFormIntegratorBase& over(int attr)
         {
            return over(std::vector{ attr });
         }

         virtual BilinearFormIntegratorBase& over(const std::vector<int>& attrs)
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
         /**
          * @brief Specifies the attribute of the elements where the
          * integration should be done
          * @param[in] attrs Element attributes
          */
         BilinearFormDomainIntegrator& over(int attr) override
         {
            return BilinearFormDomainIntegrator::over(std::vector<int>{attr});
         }

         /**
          * @brief Specifies the attributes of the elements where the
          * integration should be done.
          * @param[in] attrs Element attributes
          */
         BilinearFormDomainIntegrator& over(const std::vector<int>& attrs) override
         {
            return static_cast<BilinearFormDomainIntegrator&>(
                  BilinearFormIntegratorBase::over(attrs));
         }

         virtual BilinearFormDomainIntegrator* copy() const noexcept override = 0;
   };
}

#endif
