#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <set>
#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class LinearFormIntegratorBase
      : public FormLanguage::Buildable<mfem::LinearFormIntegrator>
   {
      public:
         virtual ~LinearFormIntegratorBase() = default;

         /**
          * @brief Gets the attributes of the elements being integrated.
          */
         virtual const std::set<int>& getAttributes() const = 0;

         /**
          * @brief Gets the integration region.
          */
         virtual IntegratorRegion getIntegratorRegion() const = 0;

         /**
          * @internal
          */
         virtual void build() override = 0;

         /**
          * @internal
          */
         virtual mfem::LinearFormIntegrator& get() override = 0;

         /**
          * @internal
          * @brief Releases ownership of the mfem::LinearFormIntegrator.
          *
          * @note After this call, calling get() will result in undefined behaviour.
          *
          * @warning The LinearFormIntegratorBase instance must still be kept
          * in memory since it might contain objects which the
          * mfem::LinearFormIntegrator instance refers to.
          */
         virtual mfem::LinearFormIntegrator* release() override = 0;

         virtual LinearFormIntegratorBase* copy() const noexcept override = 0;
   };

   class LinearFormDomainIntegrator : public LinearFormIntegratorBase
   {
      public:
         /**
          * @brief Specifies the material reference over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material reference over which the integration should
          * take place.
          */
         virtual LinearFormDomainIntegrator& over(int attr) = 0;

         /**
          * @brief Specifies the material references over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material references over which the integration should
          * take place.
          */
         virtual LinearFormDomainIntegrator& over(const std::set<int>& attrs) = 0;

         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Domain;
         }

         virtual LinearFormDomainIntegrator* copy() const noexcept override = 0;
   };

   class LinearFormBoundaryIntegrator : public LinearFormIntegratorBase
   {
      public:
         /**
          * @brief Specifies the material reference over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material reference over which the integration should
          * take place.
          */
         virtual LinearFormBoundaryIntegrator& over(int attr) = 0;

         /**
          * @brief Specifies the material references over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material references over which the integration should
          * take place.
          */
         virtual LinearFormBoundaryIntegrator& over(const std::set<int>& attrs) = 0;

         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Boundary;
         }

         virtual LinearFormBoundaryIntegrator* copy() const noexcept override = 0;
   };

}

#endif
