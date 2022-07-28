/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_HPP
#define RODIN_VARIATIONAL_LINEARFORM_HPP

#include "Rodin/Alert.h"
#include "FiniteElementSpace.h"

#include "LinearForm.h"

namespace Rodin::Variational
{
   template <class FEC>
   LinearForm<FEC, Traits::Serial>::LinearForm(TestFunction<FEC, Traits::Serial>& v)
      :  m_v(v),
         m_lf(std::make_unique<mfem::LinearForm>(&v.getFiniteElementSpace().getHandle()))
   {}

   template <class FEC>
   LinearForm<FEC, Traits::Serial>&
   LinearForm<FEC, Traits::Serial>::operator=(const LinearFormIntegratorBase& lfi)
   {
      from(lfi).assemble();
      return *this;
   }

   template <class FEC>
   LinearForm<FEC, Traits::Serial>& LinearForm<FEC, Traits::Serial>
   ::operator=(const LinearFormIntegratorSum& lsum)
   {
      from(lsum).assemble();
      return *this;
   }

   template <class FEC>
   double LinearForm<FEC, Traits::Serial>::operator()(
         const GridFunction<FEC, Traits::Serial>& u) const
   {
      return *m_lf * u.getHandle();
   }

   template <class FEC>
   LinearForm<FEC, Traits::Serial>&
   LinearForm<FEC, Traits::Serial>::add(const LinearFormIntegratorBase& lfi)
   {
      assert(lfi.getTestFunction().getLeaf().getUUID() == getTestFunction().getLeaf().getUUID());
      switch (lfi.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
         {
            auto& l = m_lfiDomainList.emplace_back(lfi.copy());
            const auto& domAttrs = lfi.getAttributes();

            l->build();
            if (domAttrs.size() == 0)
            {
               m_lf->AddDomainIntegrator(l->build().release());
            }
            else
            {
               int size = m_v.getFiniteElementSpace().getMesh().getHandle().attributes.Max();
               auto data = std::make_unique<mfem::Array<int>>(size);
               *data = 0;
               for (const auto& b : domAttrs)
               {
                  // All domain attributes are one-indexed.
                  assert(b - 1 < size);
                  (*data)[b - 1] = 1;
               }
               m_lf->AddDomainIntegrator(
                     l->build().release(),
                     *m_domAttrMarkers.emplace_back(std::move(data)));
            }
            break;
         }
         case IntegratorRegion::Boundary:
         {
            auto& l = m_lfiBoundaryList.emplace_back(lfi.copy());
            const auto& bdrAttrs = lfi.getAttributes();

            if (bdrAttrs.size() == 0)
            {
               m_lf->AddBoundaryIntegrator(l->build().release());
            }
            else
            {
               int size = m_v.getFiniteElementSpace().getMesh().getHandle().bdr_attributes.Max();
               auto data = std::make_unique<mfem::Array<int>>(size);
               *data = 0;
               for (const auto& b : bdrAttrs)
               {
                  // All domain attributes are one-indexed.
                  assert(b - 1 < size);
                  (*data)[b - 1] = 1;
               }
               m_lf->AddBoundaryIntegrator(
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

   template <class FEC>
   LinearForm<FEC, Traits::Serial>&
   LinearForm<FEC, Traits::Serial>::add(const LinearFormIntegratorSum& lsum)
   {
      for (const auto& p : lsum.getLinearFormDomainIntegratorList())
         add(*p);
      for (const auto& p : lsum.getLinearFormBoundaryIntegratorList())
         add(*p);
      return *this;
   }

   template <class FEC>
   LinearForm<FEC, Traits::Serial>&
   LinearForm<FEC, Traits::Serial>::from(const LinearFormIntegratorSum& lsum)
   {
      m_lf.reset(new mfem::LinearForm(&m_v.getFiniteElementSpace().getHandle()));
      m_lfiDomainList.clear();
      m_lfiBoundaryList.clear();
      m_domAttrMarkers.clear();
      m_bdrAttrMarkers.clear();
      add(lsum);
      return *this;
   }

   template <class FEC>
   LinearForm<FEC, Traits::Serial>&
   LinearForm<FEC, Traits::Serial>::from(const LinearFormIntegratorBase& lfi)
   {
      assert(lfi.getTestFunction().getLeaf().getUUID() == getTestFunction().getLeaf().getUUID());
      switch (lfi.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
         {
            m_lf.reset(new mfem::LinearForm(&m_v.getFiniteElementSpace().getHandle()));
            m_lfiDomainList.clear();
            m_domAttrMarkers.clear();
            add(lfi);
            break;
         }
         case IntegratorRegion::Boundary:
         {
            m_lf.reset(new mfem::LinearForm(&m_v.getFiniteElementSpace().getHandle()));
            m_lfiBoundaryList.clear();
            m_bdrAttrMarkers.clear();
            add(lfi);
            break;
         }
      }
      return *this;
   }
}

#endif
