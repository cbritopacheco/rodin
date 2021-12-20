/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_HPP
#define RODIN_VARIATIONAL_LINEARFORM_HPP

#include "FiniteElementSpace.h"
#include "LinearFormIntegrator.h"

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
      m_lfi.reset(lfi.copy());
      m_lfi->buildMFEMLinearFormIntegrator();
      m_lf.reset(new mfem::LinearForm(&m_fes.getFES()));
      // TODO: Choose whether to add a Domain, Boundary, Face, etc. integrator
      m_lf->AddDomainIntegrator(m_lfi->releaseMFEMLinearFormIntegrator());
      m_lf->Assemble();
      return *this;
   }

   template <class FEC>
   double LinearForm<FEC>::operator()(const GridFunction<FEC>& u) const
   {
      return *m_lf * u.getHandle();
   }
}

#endif
