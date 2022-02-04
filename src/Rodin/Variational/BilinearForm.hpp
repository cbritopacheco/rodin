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
      FormLanguage::ProblemBody::BFIList newList;
      newList.emplace_back(bfi.copy());
      m_bfiDomainList = std::move(newList);
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
      auto& l = m_bfiDomainList.emplace_back(bfi.copy());
      const auto& domAttrs = bfi.getAttributes();

      l->buildMFEMBilinearFormIntegrator();
      if (domAttrs.size() == 0)
      {
         m_bf->AddDomainIntegrator(l->releaseMFEMBilinearFormIntegrator());
      }
      else
      {
         int size = m_fes.getMesh().getHandle().attributes.Max();
         auto data = std::make_unique<mfem::Array<int>>(size);
         *data = 0;
         for (const auto& b : domAttrs)
         {
            // All domain attributes are one-indexed.
            assert(b - 1 < size);
            (*data)[b - 1] = 1;
         }
         m_bf->AddDomainIntegrator(
               l->releaseMFEMBilinearFormIntegrator(),
               *m_domAttrMarkers.emplace_back(std::move(data)));
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
