#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H

#include <set>
#include <memory>
#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "ShapeFunction.h"
#include "Assembly.h"

namespace Rodin::Variational
{
   class BilinearFormIntegratorBase : public FormLanguage::Base
   {
      public:
         BilinearFormIntegratorBase(
               const ShapeFunctionBase<TrialSpace>& u,
               const ShapeFunctionBase<TestSpace>& v)
            :  m_u(u.copy()), m_v(v.copy())
         {}

         BilinearFormIntegratorBase(const BilinearFormIntegratorBase& other)
            :  FormLanguage::Base(other),
               m_u(other.m_u->copy()), m_v(other.m_v->copy()),
               m_attrs(other.m_attrs)
         {}

         BilinearFormIntegratorBase(BilinearFormIntegratorBase&& other)
            :  FormLanguage::Base(std::move(other)),
               m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
               m_attrs(std::move(other.m_attrs))
         {}

         /**
          * @brief Gets the attributes of the elements being integrated.
          */
         const std::set<int>& getAttributes() const
         {
            return m_attrs;
         }

         /**
          * @brief Gets reference to trial function.
          */
         const ShapeFunctionBase<TrialSpace>& getTrialFunction() const
         {
            assert(m_u);
            return *m_u;
         }

         /**
          * @brief Gets reference to test function.
          */
         const ShapeFunctionBase<TestSpace>& getTestFunction() const
         {
            assert(m_v);
            return *m_v;
         }

         /**
          * @brief Specifies the material reference over which to integrate.
          * @returns Reference to self (for method chaining)
          *
          * Specifies the material reference over which the integration should
          * take place.
          */
         BilinearFormIntegratorBase& over(int attr)
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
         BilinearFormIntegratorBase& over(const std::set<int>& attrs)
         {
            assert(attrs.size() > 0);
            m_attrs = attrs;
            return *this;
         }

         std::unique_ptr<mfem::BilinearFormIntegrator> build() const;

         virtual ~BilinearFormIntegratorBase() = default;

         /**
          * @brief Gets the integration region.
          */
         virtual IntegratorRegion getIntegratorRegion() const = 0;

         virtual bool isSupported(Bilinear::Assembly::Type t) const
         {
            switch (t)
            {
               case Bilinear::Assembly::Type::Common:
                  return true;
               default:
                  return false;
            }
            return false;
         }

         virtual void getElementMatrix(const Bilinear::Assembly::Common& as) const = 0;

         virtual BilinearFormIntegratorBase* copy() const noexcept override = 0;

      private:
         std::unique_ptr<ShapeFunctionBase<TrialSpace>> m_u;
         std::unique_ptr<ShapeFunctionBase<TestSpace>>  m_v;
         std::set<int> m_attrs;
   };
}

namespace Rodin::Variational::Internal
{
   class ProxyBilinearFormIntegrator : public mfem::BilinearFormIntegrator
  {
      public:
         ProxyBilinearFormIntegrator(const BilinearFormIntegratorBase& bfi)
            : m_bfi(bfi)
         {}

         ProxyBilinearFormIntegrator(const ProxyBilinearFormIntegrator& other)
            : mfem::BilinearFormIntegrator(other),
              m_bfi(other.m_bfi)
         {}

         ProxyBilinearFormIntegrator(ProxyBilinearFormIntegrator&& other)
            : mfem::BilinearFormIntegrator(std::move(other)),
              m_bfi(other.m_bfi)
         {}

         void AssembleElementMatrix(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override
         {
            m_bfi.getElementMatrix(Bilinear::Assembly::Common{fe, fe, trans, mat});
         }

         void AssembleElementMatrix2(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
                  mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) override
         {
            m_bfi.getElementMatrix(Bilinear::Assembly::Common{trial, test, trans, mat});
         }
      private:
         const BilinearFormIntegratorBase& m_bfi;
   };
}

#endif
