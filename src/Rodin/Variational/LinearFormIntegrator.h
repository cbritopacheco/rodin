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

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Assembly.h"
#include "TestFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Base class for linear form integrators.
    *
    * An instance of LinearFormIntegratorBase performs the assembly of the
    * element vector for each finite element.
    */
   class LinearFormIntegratorBase : public FormLanguage::Base
   {
      public:
         LinearFormIntegratorBase() = default;

         LinearFormIntegratorBase(const LinearFormIntegratorBase& other)
            : FormLanguage::Base(other)
         {}

         LinearFormIntegratorBase(LinearFormIntegratorBase&& other)
            : FormLanguage::Base(std::move(other))
         {}

         virtual ~LinearFormIntegratorBase() = default;

         virtual const ShapeFunctionBase<ShapeFunctionSpaceType::Test>& getTestFunction()
            const = 0;

         /**
          * @brief Gets the attributes of the elements being integrated.
          */
         virtual const std::set<int>& getAttributes() const = 0;

         /**
          * @brief Gets the integration region.
          */
         virtual IntegratorRegion getIntegratorRegion() const = 0;

         virtual void getElementVector(const Assembly::Common& as) const = 0;

         virtual void getElementVector(const Assembly::Device&) const
         {
            assert(false); // Unimplemented
         }

         virtual bool isSupported(Assembly::Type t) const
         {
            switch (t)
            {
               case Assembly::Type::Common:
                  return true;
               default:
                  return false;
            }
            return false;
         }

         std::unique_ptr<mfem::LinearFormIntegrator> build() const;

         virtual LinearFormIntegratorBase* copy() const noexcept override = 0;
   };

   /**
    * @brief Represents a linear form integrator over the interior of the domain
    */
   class LinearFormDomainIntegrator : public LinearFormIntegratorBase
   {
      public:
         LinearFormDomainIntegrator(const ShapeFunctionBase<ShapeFunctionSpaceType::Test>& v)
            : m_v(v.copy())
         {}

         LinearFormDomainIntegrator(const LinearFormDomainIntegrator& other)
            : LinearFormIntegratorBase(other),
              m_v(other.m_v->copy()),
              m_attrs(other.m_attrs)
         {}

         LinearFormDomainIntegrator(LinearFormDomainIntegrator&& other)
            : LinearFormIntegratorBase(std::move(other)),
              m_v(std::move(other.m_v)),
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


         const ShapeFunctionBase<ShapeFunctionSpaceType::Test>& getTestFunction() const override
         {
            return *m_v;
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
         std::unique_ptr<ShapeFunctionBase<ShapeFunctionSpaceType::Test>> m_v;
         std::set<int> m_attrs;
   };

   /**
    * @brief Represents a linear form integrator over the boundary of the domain
    */
   class LinearFormBoundaryIntegrator : public LinearFormIntegratorBase
   {
      public:
         LinearFormBoundaryIntegrator(const ShapeFunctionBase<ShapeFunctionSpaceType::Test>& v)
            : m_v(v.copy())
         {}

         LinearFormBoundaryIntegrator(const LinearFormBoundaryIntegrator& other)
            : LinearFormIntegratorBase(other),
              m_v(other.m_v->copy()),
              m_attrs(other.m_attrs)
         {}

         LinearFormBoundaryIntegrator(LinearFormBoundaryIntegrator&& other)
            : LinearFormIntegratorBase(std::move(other)),
              m_v(std::move(other.m_v)),
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

         const ShapeFunctionBase<ShapeFunctionSpaceType::Test>& getTestFunction() const override
         {
            return *m_v;
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
         std::unique_ptr<ShapeFunctionBase<ShapeFunctionSpaceType::Test>> m_v;
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
            m_lfi.getElementVector(Assembly::Common{fe, trans, vec});
         }

         void AssembleDevice(const mfem::FiniteElementSpace& fes,
               const mfem::Array<int>& markers, mfem::Vector& b) override
         {
            m_lfi.getElementVector(Assembly::Device{fes, markers, b});
         }

         bool SupportsDevice() override
         {
            return m_lfi.isSupported(Assembly::Type::Device);
         }

      private:
         const LinearFormIntegratorBase& m_lfi;
   };
}

#endif
