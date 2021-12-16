/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FormLanguage/ScalarUnaryMinus.h"

#include "Problem.h"

#include "DomainLFIntegrator.h"

namespace Rodin::Variational
{
   DomainLFIntegrator::DomainLFIntegrator(const ScalarCoefficientBase& f)
      : m_f((-f).copy())
   {}

   DomainLFIntegrator::DomainLFIntegrator(const DomainLFIntegrator& other)
      : m_f(other.m_f->copy()), m_lf(other.m_lf)
   {}

   DomainLFIntegrator& DomainLFIntegrator::setLinearForm(LinearFormBase& lf)
   {
      m_lf.emplace(lf);
      return *this;
   }

   void DomainLFIntegrator::eval()
   {
      m_f->buildMFEMCoefficient();
      m_lf->get()
         .getHandle()
         .AddDomainIntegrator(
               new mfem::DomainLFIntegrator(
                  m_f->getMFEMCoefficient()));
   }

   DomainLFIntegrator& DomainLFIntegrator::toggleSign()
   {
      m_f.reset(new FormLanguage::ScalarCoefficientUnaryMinus(*m_f));
      return *this;
   }

   DomainLFIntegrator* DomainLFIntegrator::copy()
   const noexcept
   {
      return new DomainLFIntegrator(*this);
   }
}
