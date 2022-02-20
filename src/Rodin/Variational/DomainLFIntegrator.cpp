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
      : m_f(f.copy()),
        m_mfemScalar(m_f->build()),
        m_mfemLFI(*m_mfemScalar)
   {}

   DomainLFIntegrator::DomainLFIntegrator(const DomainLFIntegrator& other)
      :  LinearFormDomainIntegrator(other),
         m_f(other.m_f->copy()),
         m_mfemScalar(other.m_f->build()),
         m_mfemLFI(*m_mfemScalar)
   {}
}
