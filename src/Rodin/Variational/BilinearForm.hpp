/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_HPP
#define RODIN_VARIATIONAL_BILINEARFORM_HPP

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
      m_bfiDomainList = FormLanguage::List<BilinearFormDomainIntegrator>(bfi);
      (*m_bfiDomainList.begin()).buildMFEMBilinearFormIntegrator();
      m_bf.reset(new mfem::BilinearForm(&m_fes.getFES()));
      m_bf->AddDomainIntegrator(
            (*m_bfiDomainList.begin()).releaseMFEMBilinearFormIntegrator());
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
      m_bfiDomainList.append(bfi);
      m_bfiDomainList.back().buildMFEMBilinearFormIntegrator();
      m_bf->AddDomainIntegrator(
            m_bfiDomainList.back().releaseMFEMBilinearFormIntegrator());
      return *this;
   }

   template <class FEC>
   void BilinearForm<FEC>::update()
   {
      m_bf->Update();
   }
}

#endif
