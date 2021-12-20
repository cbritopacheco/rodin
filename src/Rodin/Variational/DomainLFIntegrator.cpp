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
      : m_f(f.copy())
   {}

   DomainLFIntegrator::DomainLFIntegrator(const DomainLFIntegrator& other)
      : m_f(other.m_f->copy())
   {}

   void DomainLFIntegrator::buildMFEMLinearFormIntegrator()
   {
      m_f->buildMFEMCoefficient();
      m_mfemLFI = std::make_unique<mfem::DomainLFIntegrator>(
            m_f->getMFEMCoefficient());
   }

   mfem::LinearFormIntegrator& DomainLFIntegrator::getMFEMLinearFormIntegrator()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator* DomainLFIntegrator::releaseMFEMLinearFormIntegrator()
   {
      return m_mfemLFI.release();
   }

   DomainLFIntegrator* DomainLFIntegrator::copy() const noexcept
   {
      return new DomainLFIntegrator(*this);
   }
}
