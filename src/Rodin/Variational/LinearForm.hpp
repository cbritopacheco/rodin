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
   LinearForm<FEC>& LinearForm<FEC>::operator=(const LinearFormIntegratorBase& lfi)
   {
      from(lfi).assemble();
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::operator=(
         const FormLanguage::LinearFormIntegratorSum& lsum)
   {
      from(lsum).assemble();
      return *this;
   }

   template <class FEC>
   double LinearForm<FEC>::operator()(const GridFunction<FEC>& u) const
   {
      return *m_lf * u.getHandle();
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::add(const LinearFormIntegratorBase& lfi)
   {
      switch (lfi.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
         {
            auto& l = m_lfiDomainList.emplace_back(lfi.copy());
            const auto& domAttrs = lfi.getAttributes();

            l->build();
            if (domAttrs.size() == 0)
            {
               m_lf->AddDomainIntegrator(l->release());
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
                     l->release(),
                     *m_domAttrMarkers.emplace_back(std::move(data)));
            }
            break;
         }
         case IntegratorRegion::Boundary:
         {
            auto& l = m_lfiBoundaryList.emplace_back(lfi.copy());
            const auto& bdrAttrs = lfi.getAttributes();

            l->build();
            if (bdrAttrs.size() == 0)
            {
               m_lf->AddBoundaryIntegrator(l->release());
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
                     l->release(),
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
   LinearForm<FEC>&
   LinearForm<FEC>::add(const FormLanguage::LinearFormIntegratorSum& lsum)
   {
      for (const auto& p : lsum.getLinearFormDomainIntegratorList())
         add(*p);
      for (const auto& p : lsum.getLinearFormBoundaryIntegratorList())
         add(*p);
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>&
   LinearForm<FEC>::from(const FormLanguage::LinearFormIntegratorSum& lsum)
   {
      m_lf.reset(new mfem::LinearForm(&m_fes.getFES()));
      m_lfiDomainList.clear();
      m_lfiBoundaryList.clear();
      m_domAttrMarkers.clear();
      m_bdrAttrMarkers.clear();
      add(lsum);
      return *this;
   }

   template <class FEC>
   LinearForm<FEC>& LinearForm<FEC>::from(const LinearFormIntegratorBase& lfi)
   {
      switch (lfi.getIntegratorRegion())
      {
         case IntegratorRegion::Domain:
         {
            m_lf.reset(new mfem::LinearForm(&m_fes.getFES()));
            m_lfiDomainList.clear();
            m_domAttrMarkers.clear();
            add(lfi);
            break;
         }
         case IntegratorRegion::Boundary:
         {
            m_lf.reset(new mfem::LinearForm(&m_fes.getFES()));
            m_lfiBoundaryList.clear();
            m_bdrAttrMarkers.clear();
            add(lfi);
            break;
         }
      }
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
