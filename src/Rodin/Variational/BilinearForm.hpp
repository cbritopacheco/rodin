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
   BilinearForm<FEC>& BilinearForm<FEC>::operator=(
         const BilinearFormIntegratorBase& bfi)
   {
      from(bfi).assemble();
      return *this;
   }

   template <class FEC>
   BilinearForm<FEC>& BilinearForm<FEC>::from(const BilinearFormIntegratorBase& bfi)
   {
      m_bfi.reset(bfi.copy());
      m_bfi->buildMFEMBilinearFormIntegrator();
      m_bf.reset(new mfem::BilinearForm(&m_fes.getFES()));
      // TODO: Choose whether to add a Domain, Boundary, Face, etc. integrator
      m_bf->AddDomainIntegrator(m_bfi->releaseMFEMBilinearFormIntegrator());
      return *this;
   }

   template <class FEC>
   void BilinearForm<FEC>::assemble()
   {
      m_bf->Assemble();
   }
}

#endif
