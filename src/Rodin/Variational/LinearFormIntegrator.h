/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <set>
#include <memory>
#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "TestFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Base class for linear form integrators.
    *
    * An instance of LinearFormIntegratorBase performs the assembly of the
    * element vector for each finite element.
    */
   class LinearFormIntegratorBase
      : public FormLanguage::Buildable<mfem::LinearFormIntegrator>
   {
      public:
         virtual ~LinearFormIntegratorBase() = default;

         virtual const ShapeFunctionBase<Test>& getTestFunction() const = 0;

         /**
          * @brief Gets the attributes of the elements being integrated.
          */
         virtual const std::set<int>& getAttributes() const = 0;

         /**
          * @brief Gets the integration region.
          */
         virtual IntegratorRegion getIntegratorRegion() const = 0;

         /**
          * @brief Performs the assembly of the element vector.
          * @param[in] fe Test space finite element
          * @param[in] trans Element transformation
          * @param[out] vec Element vector
          */
         virtual void getElementVector(
               const mfem::FiniteElement& fe, mfem::ElementTransformation&
               trans, mfem::Vector& vec) const = 0;

         std::unique_ptr<mfem::LinearFormIntegrator> build() const override;

         virtual LinearFormIntegratorBase* copy() const noexcept override = 0;
   };

   /**
    * @brief Represents a linear form integrator over the interior of the domain
    */
   class LinearFormDomainIntegrator : public LinearFormIntegratorBase
   {
      public:
         LinearFormDomainIntegrator(const ShapeFunctionBase<Test>& v)
            : m_v(v)
         {}

         LinearFormDomainIntegrator(const LinearFormDomainIntegrator& other)
            : LinearFormIntegratorBase(other),
              m_v(other.m_v),
              m_attrs(other.m_attrs)
         {}

         LinearFormDomainIntegrator(LinearFormDomainIntegrator&& other)
            : LinearFormIntegratorBase(std::move(other)),
              m_v(other.m_v),
              m_attrs(std::move(other.m_attrs))
         {}

         virtual ~LinearFormDomainIntegrator() = default;

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

         const ShapeFunctionBase<Test>& getTestFunction() const override
         {
            return m_v;
         }

         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Domain;
         }

         const std::set<int>& getAttributes() const override
         {
            return m_attrs;
         }

         virtual LinearFormDomainIntegrator* copy() const noexcept override = 0;
      private:
         const ShapeFunctionBase<Test>& m_v;
         std::set<int> m_attrs;
   };

   /**
    * @brief Represents a linear form integrator over the boundary of the domain
    */
   class LinearFormBoundaryIntegrator : public LinearFormIntegratorBase
   {
      public:
         LinearFormBoundaryIntegrator(const ShapeFunctionBase<Test>& v)
            : m_v(v)
         {}

         LinearFormBoundaryIntegrator(const LinearFormBoundaryIntegrator& other)
            : LinearFormIntegratorBase(other),
              m_v(other.m_v),
              m_attrs(other.m_attrs)
         {}

         LinearFormBoundaryIntegrator(LinearFormBoundaryIntegrator&& other)
            : LinearFormIntegratorBase(std::move(other)),
              m_v(other.m_v),
              m_attrs(std::move(other.m_attrs))
         {}

         virtual ~LinearFormBoundaryIntegrator() = default;

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

         const ShapeFunctionBase<Test>& getTestFunction() const override
         {
            return m_v;
         }

         IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Boundary;
         }

         const std::set<int>& getAttributes() const override
         {
            return m_attrs;
         }

         virtual LinearFormBoundaryIntegrator* copy() const noexcept override = 0;
      private:
         const ShapeFunctionBase<Test>& m_v;
         std::set<int> m_attrs;
   };
}

namespace Rodin::Variational::Internal
{
   class ProxyLinearFormIntegrator : public mfem::LinearFormIntegrator
   {
      public:
         ProxyLinearFormIntegrator(const LinearFormIntegratorBase& lfi)
            : m_lfi(lfi)
         {}

         ProxyLinearFormIntegrator(const ProxyLinearFormIntegrator& other)
            : mfem::LinearFormIntegrator(other),
              m_lfi(other.m_lfi)
         {}

         ProxyLinearFormIntegrator(ProxyLinearFormIntegrator&& other)
            : mfem::LinearFormIntegrator(std::move(other)),
              m_lfi(other.m_lfi)
         {}

         void AssembleRHSElementVect(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) override
         {
            m_lfi.getElementVector(fe, trans, vec);
         }
      private:
         const LinearFormIntegratorBase& m_lfi;
   };
}

#endif
