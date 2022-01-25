/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_HPP
#define RODIN_VARIATIONAL_BILINEARFORM_HPP

#include <cassert>

#include "FiniteElementSpace.h"
#include "BilinearFormIntegrator.h"

#include "BilinearForm.h"

namespace Rodin::Variational
{
   template <class FEC>
   BilinearForm<FEC>::BilinearForm(FiniteElementSpace<FEC>& fes)
      :  m_fes(fes),
         m_bf(new mfem::BilinearForm(&fes.getFES()))
   {}

   template <class FEC>
   double BilinearForm<FEC>::operator()(
         const GridFunction<FEC>& u, const GridFunction<FEC>& v) const
   {
      return m_bf->InnerProduct(u.getHandle(), v.getHandle());
   }

   template <class FEC>
   BilinearForm<FEC>& BilinearForm<FEC>::from(const BilinearFormDomainIntegrator& bfi)
   {
      m_bfiDomainList = FormLanguage::List<BilinearFormIntegratorBase>(bfi);
      m_bf.reset(new mfem::BilinearForm(&m_fes.getFES()));
      m_domAttrMarkers.clear();
      add(bfi);
      return *this;
   }

   template <class FEC>
   void BilinearForm<FEC>::assemble()
   {
      m_bf->Assemble();
   }

   template <class FEC>
   BilinearForm<FEC>& BilinearForm<FEC>::add(
         const BilinearFormDomainIntegrator& bfi)
   {
      auto& l = m_bfiDomainList.append(bfi);
      const auto& domAttrs = bfi.getAttributes();

      l.buildMFEMBilinearFormIntegrator();
      if (domAttrs.size() == 0)
      {
         m_bf->AddDomainIntegrator(l.releaseMFEMBilinearFormIntegrator());
      }
      else
      {
         int size = m_fes.getMesh().getHandle().attributes.Max();
         int* data = new int[size];
         std::fill(data, data + size, 0);
         for (size_t i = 0; i < domAttrs.size(); i++)
         {
            assert(domAttrs[i] < size);
            // All domain attributes are one-indexed.
            data[domAttrs[i] - 1] = 1;
         }

         auto& arr = m_domAttrMarkers.emplace_back(data, size);
         arr.MakeDataOwner();
         m_bf->AddDomainIntegrator(
               l.releaseMFEMBilinearFormIntegrator(),
               m_domAttrMarkers.back()
               );
      }
      return *this;
   }

   template <class FEC>
   void BilinearForm<FEC>::update()
   {
      m_bf->Update();
   }
}

#endif
