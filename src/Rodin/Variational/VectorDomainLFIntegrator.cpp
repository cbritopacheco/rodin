/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FormLanguage/ScalarUnaryMinus.h"

#include "Problem.h"

#include "VectorDomainLFIntegrator.h"

namespace Rodin::Variational
{
   VectorDomainLFIntegrator::VectorDomainLFIntegrator(const VectorCoefficientBase& f)
      : m_f(f.copy())
   {}

   VectorDomainLFIntegrator::VectorDomainLFIntegrator(const VectorDomainLFIntegrator& other)
      :  m_attr(other.m_attr),
         m_f(other.m_f->copy())
   {}

   void VectorDomainLFIntegrator::build()
   {
      m_f->build();
      m_mfemLFI = std::make_unique<mfem::VectorDomainLFIntegrator>(m_f->get());
   }

   mfem::LinearFormIntegrator& VectorDomainLFIntegrator::get()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   mfem::LinearFormIntegrator* VectorDomainLFIntegrator::release()
   {
      return m_mfemLFI.release();
   }

   VectorDomainLFIntegrator* VectorDomainLFIntegrator::copy() const noexcept
   {
      return new VectorDomainLFIntegrator(*this);
   }
}

