#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class LinearFormIntegratorBase : public FormLanguage::Base
   {
      public:
         virtual ~LinearFormIntegratorBase() = default;

         virtual void buildMFEMLinearFormIntegrator() = 0;

         virtual mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() = 0;

         /**
          * @internal
          * @brief Releases ownership of the mfem::LinearFormIntegrator.
          *
          * @note After this call, calling getMFEMLinearFormIntegrator() will
          * result in undefined behaviour.
          *
          * @warning The LinearFormIntegratorBase instance must still be kept
          * in memory since it might contain objects which the
          * mfem::LinearFormIntegrator instance refers to.
          */
         virtual mfem::LinearFormIntegrator* releaseMFEMLinearFormIntegrator() = 0;

         virtual LinearFormIntegratorBase* copy() const noexcept override = 0;
   };

   class LinearFormDomainIntegrator : public LinearFormIntegratorBase
   {
      public:
         virtual LinearFormDomainIntegrator* copy() const noexcept override = 0;
   };

   class LinearFormBoundaryIntegrator : public LinearFormIntegratorBase
   {
      public:
         LinearFormBoundaryIntegrator(const std::vector<int>& bdrAttrs)
            : m_bdrAttrs(bdrAttrs)
         {}

         LinearFormBoundaryIntegrator(
               const LinearFormBoundaryIntegrator&) = default;

         LinearFormBoundaryIntegrator(LinearFormBoundaryIntegrator&&) = default;

         virtual ~LinearFormBoundaryIntegrator() = default;

         virtual LinearFormBoundaryIntegrator* copy() const noexcept override = 0;

         virtual const std::vector<int>& getBoundaryAttributes() const
         {
            return m_bdrAttrs;
         }

      private:
         std::vector<int> m_bdrAttrs;
   };

}

#endif
