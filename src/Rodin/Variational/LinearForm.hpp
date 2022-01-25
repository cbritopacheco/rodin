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
      auto& l = m_lfiBoundaryList.append(lfi);
      const auto& domAttrs = lfi.getAttributes();

      l.buildMFEMLinearFormIntegrator();
      if (domAttrs.size() == 0)
      {
         m_lf->AddDomainIntegrator(l.releaseMFEMLinearFormIntegrator());
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
         m_lf->AddDomainIntegrator(
               l.releaseMFEMLinearFormIntegrator(),
               m_domAttrMarkers.back()
               );
      }
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::add(const LinearFormBoundaryIntegrator& lfi)
   {
      auto& l = m_lfiBoundaryList.append(lfi);
      const auto& bdrAttrs = lfi.getAttributes();

      l.buildMFEMLinearFormIntegrator();
      if (bdrAttrs.empty())
      {
         m_lf->AddBoundaryIntegrator(l.releaseMFEMLinearFormIntegrator());
      }
      else
      {
         int size = m_fes.getMesh().getHandle().bdr_attributes.Max();
         int* data = new int[size];
         std::fill(data, data + size, 0);
         for (size_t i = 0; i < bdrAttrs.size(); i++)
         {
            assert(bdrAttrs[i] < size);
            // All boundary attributes are one-indexed.
            data[bdrAttrs[i] - 1] = 1;
         }

         auto& arr = m_bdrAttrMarkers.emplace_back(data, size);
         arr.MakeDataOwner();
         m_lf->AddBoundaryIntegrator(
               l.releaseMFEMLinearFormIntegrator(),
               m_bdrAttrMarkers.back()
               );
      }
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::from(const LinearFormDomainIntegrator& lfi)
   {
      m_lfiDomainList = FormLanguage::List<LinearFormIntegratorBase>(lfi);
      m_lf.reset(new mfem::LinearForm(&m_fes.getFES()));
      m_domAttrMarkers.clear();
      add(lfi);
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::from(const LinearFormBoundaryIntegrator& lfi)
   {
      m_lfiBoundaryList = FormLanguage::List<LinearFormIntegratorBase>(lfi);
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
