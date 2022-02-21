#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <set>
#include <memory>
#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

namespace Rodin::Variational::Internal
{
   class LinearFormIntegrator : public mfem::LinearFormIntegrator
   {
      public:
         LinearFormIntegrator(const LinearFormIntegratorBase& bfi);

         LinearFormIntegrator(const LinearFormIntegrator& other);

         void AssembleRHSElementVect(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) override;
      private:
         std::unique_ptr<LinearFormIntegratorBase> m_lfi;
   };
}

namespace Rodin::Variational
{
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
         LinearFormDomainIntegrator() = default;
         LinearFormDomainIntegrator(const LinearFormDomainIntegrator&) = default;
         LinearFormDomainIntegrator(LinearFormDomainIntegrator&&) = default;

         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Domain;
         }

         const std::set<int>& getAttributes() const override
         {
            return m_attrs;
         }

         /**
          * @brief Specifies the material reference over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material reference over which the integration should
          * take place.
          */
         LinearFormDomainIntegrator& over(int attr)
         {
            return over(std::set<int>{attr});
         }

         /**
          * @brief Specifies the material references over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material references over which the integration should
          * take place.
          */
         LinearFormDomainIntegrator& over(const std::set<int>& attrs)
         {
            assert(attrs.size() > 0);
            m_attrs = attrs;
            return *this;
         }

         virtual LinearFormDomainIntegrator* copy() const noexcept override = 0;
      private:
         std::set<int> m_attrs;
   };

   class LinearFormBoundaryIntegrator : public LinearFormIntegratorBase
   {
      public:
         LinearFormBoundaryIntegrator() = default;
         LinearFormBoundaryIntegrator(const LinearFormBoundaryIntegrator&) = default;
         LinearFormBoundaryIntegrator(LinearFormBoundaryIntegrator&&) = default;

         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Boundary;
         }

         const std::set<int>& getAttributes() const override
         {
            return m_attrs;
         }

         /**
          * @brief Specifies the material reference over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material reference over which the integration should
          * take place.
          */
         LinearFormBoundaryIntegrator& over(int attr)
         {
            return over(std::set<int>{attr});
         }

         /**
          * @brief Specifies the material references over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material references over which the integration should
          * take place.
          */
         LinearFormBoundaryIntegrator& over(const std::set<int>& attrs)
         {
            assert(attrs.size() > 0);
            m_attrs = attrs;
            return *this;
         }


         virtual LinearFormBoundaryIntegrator* copy() const noexcept override = 0;
      private:
         std::set<int> m_attrs;
   };

}

#endif
