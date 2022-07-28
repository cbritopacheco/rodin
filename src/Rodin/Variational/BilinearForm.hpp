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

#include "FiniteElementSpace.h"
#include "BilinearFormIntegrator.h"
#include "BilinearFormIntegratorSum.h"

#include "BilinearForm.h"

namespace Rodin::Variational
{
   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>::BilinearForm(
         TrialFunction<TrialFEC, Traits::Serial>& u, TestFunction<TestFEC, Traits::Serial>& v)
      :  m_u(u), m_v(v),
         m_bf(new mfem::BilinearForm(&m_v.getFiniteElementSpace().getHandle()))
   {}

   template <class TrialFEC, class TestFEC>
   double BilinearForm<TrialFEC, TestFEC, Traits::Serial>::operator()(
         const GridFunction<TrialFEC, Traits::Serial>& u,
         const GridFunction<TestFEC, Traits::Serial>& v) const
   {
      return m_bf->InnerProduct(u.getHandle(), v.getHandle());
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>&
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>::operator=(const BilinearFormIntegratorBase& bfi)
   {
      from(bfi).assemble();
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>&
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>::operator=(const BilinearFormIntegratorSum& bfi)
   {
      from(bfi).assemble();
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>&
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>::from(const BilinearFormIntegratorBase& bfi)
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
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>&
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>::from(const BilinearFormIntegratorSum& lsum)
   {
      m_bf.reset(new mfem::BilinearForm(&m_u.getFiniteElementSpace().getHandle()));
      m_bfiDomainList.clear();
      m_domAttrMarkers.clear();
      add(lsum);
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>&
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>::add(const BilinearFormIntegratorSum& lsum)
   {
      for (const auto& p : lsum.getBilinearFormDomainIntegratorList())
         add(*p);
      return *this;
   }

   template <class TrialFEC, class TestFEC>
   BilinearForm<TrialFEC, TestFEC, Traits::Serial>& BilinearForm<TrialFEC, TestFEC, Traits::Serial>::add(
         const BilinearFormIntegratorBase& bfi)
   {
      assert(
         bfi.getTrialFunction().getLeaf().getUUID()== getTrialFunction().getLeaf().getUUID());
      assert(
         bfi.getTestFunction().getLeaf().getUUID() == getTestFunction().getLeaf().getUUID());

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
         case IntegratorRegion::Boundary:
         {
            auto& l = m_bfiBoundaryList.emplace_back(bfi.copy());
            const auto& bdrAttrs = bfi.getAttributes();

            if (bdrAttrs.size() == 0)
            {
               m_bf->AddBoundaryIntegrator(l->build().release());
            }
            else
            {
               int size = m_u.getFiniteElementSpace().getMesh().getHandle().bdr_attributes.Max();
               auto data = std::make_unique<mfem::Array<int>>(size);
               *data = 0;
               for (const auto& b : bdrAttrs)
               {
                  // All Boundary attributes are one-indexed.
                  assert(b - 1 < size);
                  (*data)[b - 1] = 1;
               }
               m_bf->AddBoundaryIntegrator(
                     l->build().release(),
                     *m_bdrAttrMarkers.emplace_back(std::move(data)));
            }
            break;
         }
         default:
            Alert::Exception() << "IntegratorRegion not supported." << Alert::Raise;
      }
      return *this;
   }
}

#endif
