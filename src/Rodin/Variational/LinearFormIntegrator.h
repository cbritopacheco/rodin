#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <set>
#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class LinearFormIntegratorBase : public FormLanguage::Base
   {
      public:
         virtual ~LinearFormIntegratorBase() = default;

         virtual const std::set<int>& getAttributes() const = 0;

         virtual void buildMFEMLinearFormIntegrator() = 0;

         virtual mfem::LinearFormIntegrator& getMFEMLinearFormIntegrator() = 0;

         virtual IntegratorRegion getIntegratorRegion() const = 0;

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
         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Domain;
         }

         virtual LinearFormDomainIntegrator& over(int attr) = 0;

         virtual LinearFormDomainIntegrator& over(const std::set<int>& attrs) = 0;

         virtual LinearFormDomainIntegrator* copy() const noexcept override = 0;
   };

   class LinearFormBoundaryIntegrator : public LinearFormIntegratorBase
   {
      public:
         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Boundary;
         }

         virtual LinearFormBoundaryIntegrator& over(int attr) = 0;

         virtual LinearFormBoundaryIntegrator& over(const std::set<int>& attrs) = 0;

         virtual LinearFormBoundaryIntegrator* copy() const noexcept override = 0;
   };

}

#endif
