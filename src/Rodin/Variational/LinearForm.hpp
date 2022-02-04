/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_HPP
#define RODIN_VARIATIONAL_LINEARFORM_HPP

#include "FiniteElementSpace.h"

#include "LinearForm.h"

namespace Rodin::Variational
{
   template <class FEC>
   LinearForm<FEC>::LinearForm(FiniteElementSpace<FEC>& fes)
      :  m_fes(fes),
         m_lf(std::make_unique<mfem::LinearForm>(&fes.getFES()))
   {}

   template <class FEC>
   double LinearForm<FEC>::operator()(const GridFunction<FEC>& u) const
   {
      return *m_lf * u.getHandle();
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::add(const LinearFormDomainIntegrator& lfi)
   {
      auto& l = m_lfiDomainList.emplace_back(lfi.copy());
      const auto& domAttrs = lfi.getAttributes();

      l->buildMFEMLinearFormIntegrator();
      if (domAttrs.size() == 0)
      {
         m_lf->AddDomainIntegrator(l->releaseMFEMLinearFormIntegrator());
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
         m_lf->AddDomainIntegrator(
               l->releaseMFEMLinearFormIntegrator(),
               *m_domAttrMarkers.emplace_back(std::move(data)));
      }
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::add(const LinearFormBoundaryIntegrator& lfi)
   {
      auto& l = m_lfiBoundaryList.emplace_back(lfi.copy());
      const auto& bdrAttrs = lfi.getAttributes();

      l->buildMFEMLinearFormIntegrator();
      if (bdrAttrs.size() == 0)
      {
         m_lf->AddBoundaryIntegrator(l->releaseMFEMLinearFormIntegrator());
      }
      else
      {
         int size = m_fes.getMesh().getHandle().bdr_attributes.Max();
         auto data = std::make_unique<mfem::Array<int>>(size);
         *data = 0;
         for (const auto& b : bdrAttrs)
         {
            // All domain attributes are one-indexed.
            assert(b - 1 < size);
            (*data)[b - 1] = 1;
         }
         m_lf->AddBoundaryIntegrator(
               l->releaseMFEMLinearFormIntegrator(),
               *m_bdrAttrMarkers.emplace_back(std::move(data)));
      }
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::from(const LinearFormDomainIntegrator& lfi)
   {
      FormLanguage::ProblemBody::LFIList newList;
      newList.emplace_back(lfi.copy());
      m_lfiDomainList = std::move(newList);
      m_lf.reset(new mfem::LinearForm(&m_fes.getFES()));
      m_domAttrMarkers.clear();
      add(lfi);
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::from(const LinearFormBoundaryIntegrator& lfi)
   {
      FormLanguage::ProblemBody::LFIList newList;
      newList.emplace_back(lfi.copy());
      m_lfiBoundaryList = std::move(newList);
      m_lf.reset(new mfem::LinearForm(&m_fes.getFES()));
      m_bdrAttrMarkers.clear();
      add(lfi);
      return *this;
   }

   template <class FEC>
   void LinearForm<FEC>::assemble()
   {
      m_lf->Assemble();
   }

   template <class FEC>
   void LinearForm<FEC>::update()
   {
      m_lf->Update();
   }
}

#endif
