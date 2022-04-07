/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_HPP
#define RODIN_VARIATIONAL_BILINEARFORM_HPP

#include <cassert>

#include "Rodin/Alert.h"
#include "FormLanguage/BilinearFormIntegratorSum.h"

#include "FiniteElementSpace.h"
#include "BilinearFormIntegrator.h"

#include "BilinearForm.h"

namespace Rodin::Variational
{
   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC>::BilinearForm(
         TrialFunction<TrialFEC>& u, TestFunction<TestFEC>& v)
      :  m_u(u), m_v(v),
         m_bf(new mfem::BilinearForm(&m_v.getFiniteElementSpace().getHandle()))
   {}

   template <class TrialFEC, class TestFEC>
   double BilinearForm<TrialFEC, TestFEC>::operator()(
         const GridFunction<TrialFEC>& u, const GridFunction<TestFEC>& v) const
   {
      return m_bf->InnerProduct(u.getHandle(), v.getHandle());
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC>&
   BilinearForm<TrialFEC, TestFEC>::operator=(const BilinearFormIntegratorBase& bfi)
   {
      from(bfi).assemble();
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC>&
   BilinearForm<TrialFEC, TestFEC>::operator=(const FormLanguage::BilinearFormIntegratorSum& bfi)
   {
      from(bfi).assemble();
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC>& BilinearForm<TrialFEC, TestFEC>::from(const BilinearFormIntegratorBase& bfi)
   {
      switch (bfi.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
         {
            m_bf.reset(new mfem::BilinearForm(&m_u.getFiniteElementSpace().getHandle()));
            m_domAttrMarkers.clear();
            add(bfi);
            break;
         }
         default:
            Alert::Exception() << "IntegratorRegion not supported." << Alert::Raise;
      }
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC>&
   BilinearForm<TrialFEC, TestFEC>::from(const FormLanguage::BilinearFormIntegratorSum& lsum)
   {
      m_bf.reset(new mfem::BilinearForm(&m_u.getFiniteElementSpace().getHandle()));
      m_bfiDomainList.clear();
      m_domAttrMarkers.clear();
      add(lsum);
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   void BilinearForm<TrialFEC, TestFEC>::assemble()
   {
      m_bf->Assemble();
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC>&
   BilinearForm<TrialFEC, TestFEC>::add(const FormLanguage::BilinearFormIntegratorSum& lsum)
   {
      for (const auto& p : lsum.getBilinearFormDomainIntegratorList())
         add(*p);
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC>& BilinearForm<TrialFEC, TestFEC>::add(
         const BilinearFormIntegratorBase& bfi)
   {
      assert(
         bfi.getTrialFunction().getRoot().getUUID()== getTrialFunction().getRoot().getUUID());
      assert(
         bfi.getTestFunction().getRoot().getUUID() == getTestFunction().getRoot().getUUID());

      switch (bfi.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
         {
            auto& l = m_bfiDomainList.emplace_back(bfi.copy());
            const auto& domAttrs = bfi.getAttributes();

            if (domAttrs.size() == 0)
            {
               m_bf->AddDomainIntegrator(l->build().release());
            }
            else
            {
               int size = m_u.getFiniteElementSpace().getMesh().getHandle().attributes.Max();
               auto data = std::make_unique<mfem::Array<int>>(size);
               *data = 0;
               for (const auto& b : domAttrs)
               {
                  // All domain attributes are one-indexed.
                  assert(b - 1 < size);
                  (*data)[b - 1] = 1;
               }
               m_bf->AddDomainIntegrator(
                     l->build().release(),
                     *m_domAttrMarkers.emplace_back(std::move(data)));
            }
            break;
         }
         default:
            Alert::Exception() << "IntegratorRegion not supported." << Alert::Raise;
      }
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   void BilinearForm<TrialFEC, TestFEC>::update()
   {
      m_bf->Update();
   }
}

#endif
