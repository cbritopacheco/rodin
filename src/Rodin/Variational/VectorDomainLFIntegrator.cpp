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
      : m_f(f.copy()),
        m_mfemVector(f.build()),
        m_mfemLFI(*m_mfemVector)
   {}

   VectorDomainLFIntegrator::VectorDomainLFIntegrator(const VectorDomainLFIntegrator& other)
      :  m_attr(other.m_attr),
         m_f(other.m_f->copy()),
         m_mfemVector(other.m_f->build()),
         m_mfemLFI(*m_mfemVector)
   {}

   VectorDomainLFIntegrator* VectorDomainLFIntegrator::copy() const noexcept
   {
      return new VectorDomainLFIntegrator(*this);
   }
}

