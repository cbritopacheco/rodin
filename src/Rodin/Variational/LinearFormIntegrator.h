#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <set>
#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      class LinearFormIntegrator : public mfem::LinearFormIntegrator
      {
         public:
            LinearFormIntegrator(const LinearFormIntegratorBase& bfi);

            void AssembleRHSElementVect(
                  const mfem::FiniteElement& fe,
                  mfem::ElementTransformation& trans, mfem::Vector& vec) override;
         private:
            std::unique_ptr<LinearFormIntegratorBase> m_lfi;
      };
   }

   class LinearFormIntegratorBase
      : public FormLanguage::Buildable<Internal::LinearFormIntegrator>
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

         virtual void getElementVector(
               const mfem::FiniteElement& fe, mfem::ElementTransformation&
               trans, mfem::Vector& vec) = 0;

         /**
          * @internal
          */
         std::unique_ptr<Internal::LinearFormIntegrator> build() const override
         {
            return std::make_unique<Internal::LinearFormIntegrator>(*this);
         }

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
